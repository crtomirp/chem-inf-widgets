import json
import numpy as np
import torch
import torch.nn as nn
import torch.optim as optim
from sklearn.base import BaseEstimator, RegressorMixin
from AnyQt.QtWidgets import (QApplication, QTabWidget, QVBoxLayout, QWidget, 
                             QFileDialog, QTextBrowser, QLabel)
from AnyQt.QtCore import QThread, pyqtSignal
from Orange.widgets.widget import OWWidget, Input, Output
from Orange.widgets import gui
from Orange.widgets.settings import Setting
from Orange.data import Table, Domain, ContinuousVariable, DiscreteVariable
from sklearn.ensemble import RandomForestRegressor, GradientBoostingRegressor
from sklearn.svm import SVR
from sklearn.linear_model import ElasticNet, LogisticRegression, Lasso, Ridge
from sklearn.cross_decomposition import PLSRegression
from sklearn.tree import DecisionTreeRegressor
from sklearn.preprocessing import StandardScaler, MinMaxScaler
from sklearn.impute import SimpleImputer
from sklearn.model_selection import (cross_val_score, train_test_split, 
                                     GridSearchCV, RandomizedSearchCV)
from sklearn.pipeline import Pipeline
from sklearn.feature_selection import SelectKBest, f_regression, mutual_info_classif
from matplotlib.figure import Figure
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_pdf import PdfPages
# Additional metrics for regression evaluation
from sklearn.metrics import mean_squared_error, mean_absolute_error, median_absolute_error, explained_variance_score, r2_score

# ----------------------------------------------------------------------
# A simple PyTorch-based regressor, scikit-learn compatible.
class TorchRegressor(BaseEstimator, RegressorMixin):
    def __init__(self, hidden_layer_size=256, epochs=200, lr=0.01, batch_size=32, random_state=42):
        self.hidden_layer_size = hidden_layer_size
        self.epochs = epochs
        self.lr = lr
        self.batch_size = batch_size
        self.random_state = random_state
        self.model_ = None

    def fit(self, X, y):
        X_tensor = torch.tensor(X, dtype=torch.float32)
        y_tensor = torch.tensor(y, dtype=torch.float32).view(-1, 1)
        input_dim = X_tensor.shape[1]

        self.model_ = nn.Sequential(
            nn.Linear(input_dim, self.hidden_layer_size),
            nn.ReLU(),
            nn.Linear(self.hidden_layer_size, 1)
        )
        optimizer = optim.Adam(self.model_.parameters(), lr=self.lr)
        criterion = nn.MSELoss()

        dataset = torch.utils.data.TensorDataset(X_tensor, y_tensor)
        loader = torch.utils.data.DataLoader(dataset, batch_size=self.batch_size, shuffle=True)

        self.model_.train()
        for epoch in range(self.epochs):
            for batch_X, batch_y in loader:
                optimizer.zero_grad()
                outputs = self.model_(batch_X)
                loss = criterion(outputs, batch_y)
                loss.backward()
                optimizer.step()
        return self

    def predict(self, X):
        self.model_.eval()
        X_tensor = torch.tensor(X, dtype=torch.float32)
        with torch.no_grad():
            predictions = self.model_(X_tensor)
        return predictions.numpy().ravel()

# ----------------------------------------------------------------------
class OWQSARRegression(OWWidget):
    name = "QSAR Regression"
    description = ("Build QSAR regression models with flexible settings. "
                   "Splits the data into training and test sets, supports an optional external set, "
                   "advanced hyperparameter tuning (Grid Search and Randomized Search), and provides diagnostic plots "
                   "with a modern HTML5 report below.")
    icon = "icons/qsar_regression.png"  # Update this path as needed
    priority = 60
    keywords = ["QSAR", "Regression", "SMILES"]

    class Inputs:
        data = Input("Data", Table)
        external_data = Input("External Data", Table)

    class Outputs:
        model = Output("Model", object, auto_summary=False)
        train_results = Output("Train Results", Table)
        test_results = Output("Test Results", Table)
        external_results = Output("External Results", Table)

    want_main_area = True

    # Persistent settings
    selected_algorithm = Setting(0)
    normalization_method = Setting(0)
    imputation_method = Setting(1)
    cv_folds = Setting(5)
    test_size = Setting(0.3)
    tuning_method = Setting(0)
    n_iter = Setting(10)
    hyperparameters = Setting("")
    enable_feature_selection = Setting(False)
    num_features = Setting(10)

    # Updated list of available algorithms:
    algorithms = [
        ("Random Forest", RandomForestRegressor),
        ("Support Vector Regression", SVR),
        ("Gradient Boosting", GradientBoostingRegressor),
        ("PLS Regression", PLSRegression),
        ("Decision Tree Regression", DecisionTreeRegressor),
        ("Lasso Regression", Lasso),
        ("Ridge Regression", Ridge),
        ("Elastic Net", ElasticNet),
        ("Deep Learning Regression", TorchRegressor)
    ]
    normalization_options = [
        "None",
        "Standard Scaler",
        "MinMax Scaler"
    ]
    imputation_options = [
        "None",
        "Mean",
        "Median",
        "Most Frequent"
    ]
    tuning_options = [
        "None",
        "Grid Search",
        "Randomized Search"
    ]

    def __init__(self):
        super().__init__()
        self.data = None
        self.external_data = None
        self.model = None
        self.worker = None
        self.last_train_fig = None
        self.last_test_fig = None
        self.last_ext_fig = None
        self.last_model_name = ""

        # --- Control Panel Setup ---
        settings_box = gui.widgetBox(self.controlArea, "Model Settings")
        gui.comboBox(settings_box, self, "selected_algorithm",
                     label="Algorithm:",
                     items=[name for name, _ in self.algorithms],
                     callback=self.settings_changed)
        gui.comboBox(settings_box, self, "normalization_method",
                     label="Normalization:",
                     items=self.normalization_options,
                     callback=self.settings_changed)
        gui.comboBox(settings_box, self, "imputation_method",
                     label="Imputation:",
                     items=self.imputation_options,
                     callback=self.settings_changed)
        gui.spin(settings_box, self, "cv_folds", minv=2, maxv=20, step=1,
                 label="CV Folds:", callback=self.settings_changed)
        gui.doubleSpin(settings_box, self, "test_size", minv=0.1, maxv=0.9, step=0.05,
                       label="Test Set Fraction:", callback=self.settings_changed)
        gui.comboBox(settings_box, self, "tuning_method",
                     label="Hyperparameter Tuning:",
                     items=self.tuning_options,
                     callback=self.settings_changed)
        gui.spin(settings_box, self, "n_iter", minv=5, maxv=100, step=5,
                 label="Randomized Search Iterations:", callback=self.settings_changed)
        gui.lineEdit(settings_box, self, "hyperparameters",
                     label="Hyperparameters (JSON):", callback=self.settings_changed)
        gui.checkBox(settings_box, self, "enable_feature_selection",
                     "Enable Descriptor Subset Selection", callback=self.settings_changed)
        gui.spin(settings_box, self, "num_features", minv=1, maxv=100, step=1,
                 label="Number of Features:", callback=self.settings_changed)
        gui.button(settings_box, self, "Commit", callback=self.commit)
        gui.button(settings_box, self, "Export PDF", callback=self.export_pdf)

        self.info_box = gui.widgetBox(self.controlArea, "Model Performance")
        self.info_label = gui.label(self.info_box, self, "No model trained yet.")

        # --- Diagnostic Plots Setup in Main Area ---
        self.tabs = QTabWidget()
        self.mainArea.layout().addWidget(self.tabs)
        self.train_tab = QWidget()
        self.test_tab = QWidget()
        self.ext_tab = QWidget()
        self.tabs.addTab(self.train_tab, "Training Diagnostics")
        self.tabs.addTab(self.test_tab, "Test Diagnostics")
        self.tabs.addTab(self.ext_tab, "External Diagnostics")
        self.train_layout = QVBoxLayout(self.train_tab)
        self.test_layout = QVBoxLayout(self.test_tab)
        self.ext_layout = QVBoxLayout(self.ext_tab)

        # --- HTML5 Report Widget for the statistics report ---
        self.report_browser = QTextBrowser()
        self.report_browser.setStyleSheet(
            "background-color: #f9f9f9; padding: 10px; border: 1px solid #ccc; font-family: Arial, sans-serif;"
        )
        self.mainArea.layout().addWidget(self.report_browser)

    @Inputs.data
    def set_data(self, dataset):
        self.data = dataset
        self.commit()

    @Inputs.external_data
    def set_external_data(self, dataset):
        self.external_data = dataset
        self.commit()

    def settings_changed(self):
        pass

    def commit(self):
        if self.data is None:
            self.info_label.setText("No main data provided.")
            return

        if self.selected_algorithm < 0 or self.selected_algorithm >= len(self.algorithms):
            self.selected_algorithm = 0

        if self.worker is not None and self.worker.isRunning():
            self.worker.terminate()
            self.info_label.setText("Previous calculation terminated.")
        self.clear_layout(self.train_layout)
        self.clear_layout(self.test_layout)
        self.clear_layout(self.ext_layout)

        model_name = self.algorithms[self.selected_algorithm][0]
        self.last_model_name = model_name
        self.info_label.setText(f"Please wait calculation {model_name} is started")
        # Clear the HTML report widget and show a waiting message.
        self.report_browser.setHtml(
            '<div style="text-align: center; font-weight: bold; font-size: 14pt;">'
            'Please wait, calculation of the QSAR model in progress'
            '</div>'
        )

        config = {
            "selected_algorithm": self.selected_algorithm,
            "normalization_method": self.normalization_method,
            "imputation_method": self.imputation_method,
            "cv_folds": self.cv_folds,
            "test_size": self.test_size,
            "tuning_method": self.tuning_method,
            "n_iter": self.n_iter,
            "hyperparameters": self.hyperparameters,
            "enable_feature_selection": self.enable_feature_selection,
            "num_features": self.num_features,
            "algorithms": self.algorithms
        }

        self.worker = QSARWorker(self.data, self.external_data, config)
        self.worker.finished_signal.connect(self.handle_results)
        self.worker.error_signal.connect(self.handle_error)
        self.worker.start()

    def handle_results(self, result):
        self.model = result["model"]
        self.Outputs.model.send(self.model)
        self.Outputs.train_results.send(result["train_table"])
        self.Outputs.test_results.send(result["test_table"])
        if result.get("external_table") is not None:
            self.Outputs.external_results.send(result["external_table"])
        self.info_label.setText(f"Calculation {self.last_model_name} is completed.\n{result['performance_text']}")

        self.update_diagnostics("train", result["X_train"], result["y_train"], result["pipeline"],
                                  result["is_classification"])
        self.update_diagnostics("test", result["X_test"], result["y_test"], result["pipeline"],
                                  result["is_classification"])
        if self.external_data is not None:
            self.update_diagnostics("external", result["X_ext"], result["y_ext"], result["pipeline"],
                                    result["is_classification"])

        self.update_report_browser(result)

    def handle_error(self, error_msg):
        self.info_label.setText("Error: " + error_msg)

    def update_diagnostics(self, dataset_type, X, y, pipeline, is_classification=False):
        preds = pipeline.predict(X)
        fig = Figure(figsize=(10, 4))
        ax1 = fig.add_subplot(121)
        ax2 = fig.add_subplot(122)

        if not is_classification:
            residuals = y - preds
            threshold = 2 * np.std(residuals)
            inlier_mask = np.abs(residuals) <= threshold
            outlier_mask = np.abs(residuals) > threshold

            ax1.scatter(preds[inlier_mask], y[inlier_mask], alpha=0.7, c="blue", edgecolors="k", label="Inliers")
            if np.any(outlier_mask):
                ax1.scatter(preds[outlier_mask], y[outlier_mask], alpha=0.7, c="red", edgecolors="k", label="Outliers")
            ax1.plot([min(preds), max(preds)], [min(preds), max(preds)], 'r--', lw=2)
            ax1.set_title("Predicted vs Actual")
            ax1.set_xlabel("Predicted")
            ax1.set_ylabel("Actual")
            ax1.grid()
            ax1.legend()

            ax2.scatter(preds[inlier_mask], residuals[inlier_mask], alpha=0.7, c="blue", edgecolors="k", label="Inliers")
            if np.any(outlier_mask):
                ax2.scatter(preds[outlier_mask], residuals[outlier_mask], alpha=0.7, c="red", edgecolors="k", label="Outliers")
            ax2.axhline(0, color="r", linestyle="--", lw=2)
            ax2.set_title("Residuals vs Predicted")
            ax2.set_xlabel("Predicted")
            ax2.set_ylabel("Residuals")
            ax2.grid()
            ax2.legend()
        else:
            ax1.scatter(preds, y, alpha=0.7, c="green", edgecolors="k")
            ax1.plot([min(preds), max(preds)], [min(preds), max(preds)], 'r--', lw=2)
            ax1.set_title("Predicted vs Actual")
            ax1.set_xlabel("Predicted")
            ax1.set_ylabel("Actual")

            ax2.scatter(preds, (y != preds).astype(int), alpha=0.7, c="green", edgecolors="k")
            ax2.axhline(0, color="r", linestyle="--", lw=2)
            ax2.set_title("Misclassifications (1 if error)")
            ax2.set_xlabel("Predicted")
            ax2.set_ylabel("Error Indicator")

        canvas = FigureCanvas(fig)
        if dataset_type == "train":
            self.clear_layout(self.train_layout)
            self.train_layout.addWidget(canvas)
            self.last_train_fig = fig
        elif dataset_type == "test":
            self.clear_layout(self.test_layout)
            self.test_layout.addWidget(canvas)
            self.last_test_fig = fig
        elif dataset_type == "external":
            self.clear_layout(self.ext_layout)
            self.ext_layout.addWidget(canvas)
            self.last_ext_fig = fig

    def update_report_browser(self, result):
        total_desc = len(self.data.domain.attributes) if self.data is not None else 0
        used_desc = self.num_features if self.enable_feature_selection else total_desc
        cv_score = result.get("cv_score", None)
        cv_text = f"{cv_score:.3f}" if cv_score is not None else "N/A"
        train_metrics = result.get("train_metrics", {})
        test_metrics = result.get("test_metrics", {})
        external_metrics = result.get("external_metrics", {})

        metrics_list = ["R²", "RMSE", "MAE", "Median AE", "Explained Variance"]

        html = f"""
        <html>
          <head>
            <style>
              body {{ font-family: Arial, sans-serif; font-size: 12pt; color: #333; }}
              h2 {{ color: #444; }}
              table {{ width: 100%; border-collapse: collapse; margin-top: 10px; }}
              th, td {{ border: 1px solid #ccc; padding: 5px; text-align: center; }}
              th {{ background-color: #f0f0f0; }}
            </style>
          </head>
          <body>
            <h2>Model Report</h2>
            <p><b>Model:</b> {self.last_model_name}</p>
            <p><b>Total Descriptors:</b> {total_desc}</p>
            <p><b>Descriptors Used:</b> {used_desc}</p>
            <p><b>CV R²:</b> {cv_text}</p>
            <h3>Metrics</h3>
            <table>
              <tr>
                <th>Metric</th>
                <th>Training</th>
                <th>Test</th>
                <th>External</th>
              </tr>
        """
        for metric in metrics_list:
            train_val = f"{train_metrics[metric]:.3f}" if metric in train_metrics else "N/A"
            test_val = f"{test_metrics[metric]:.3f}" if metric in test_metrics else "N/A"
            ext_val = f"{external_metrics[metric]:.3f}" if external_metrics and metric in external_metrics else "N/A"
            html += f"""
              <tr>
                <td>{metric}</td>
                <td>{train_val}</td>
                <td>{test_val}</td>
                <td>{ext_val}</td>
              </tr>
            """
        html += """
            </table>
          </body>
        </html>
        """
        self.report_browser.setHtml(html)

    def clear_layout(self, layout):
        while layout.count():
            child = layout.takeAt(0)
            if child.widget():
                child.widget().deleteLater()

    def create_pdf_report_figure(self, result):
        fig = Figure(figsize=(8, 6))
        ax = fig.add_subplot(111)
        ax.axis("off")
        model_info = f"Model: {self.last_model_name}\n"
        total_desc = len(self.data.domain.attributes) if self.data is not None else 0
        used_desc = self.num_features if self.enable_feature_selection else total_desc
        descriptor_info = f"Total Descriptors: {total_desc}\nDescriptors Used: {used_desc}\n"
        cv_score = result.get("cv_score", None)
        cv_info = f"CV R²: {cv_score:.3f}\n\n" if cv_score is not None else "CV R²: N/A\n\n"
        report_text = model_info + descriptor_info + cv_info

        train_metrics = result.get("train_metrics", {})
        report_text += "Training Metrics:\n"
        for k, v in train_metrics.items():
            report_text += f"  {k}: {v:.3f}\n"
        report_text += "\nTest Metrics:\n"
        test_metrics = result.get("test_metrics", {})
        for k, v in test_metrics.items():
            report_text += f"  {k}: {v:.3f}\n"
        report_text += "\nExternal Metrics:\n"
        external_metrics = result.get("external_metrics", {})
        for k, v in external_metrics.items():
            report_text += f"  {k}: {v:.3f}\n"
        ax.text(0, 1, report_text, va="top", ha="left", fontsize=10, wrap=True)
        return fig

    def export_pdf(self):
        filename, _ = QFileDialog.getSaveFileName(self, "Export PDF", "", "PDF Files (*.pdf)")
        if filename:
            try:
                with PdfPages(filename) as pdf:
                    if self.worker is not None and self.worker.last_result is not None:
                        report_fig = self.create_pdf_report_figure(self.worker.last_result)
                        pdf.savefig(report_fig)
                    if self.last_train_fig is not None:
                        pdf.savefig(self.last_train_fig)
                    if self.last_test_fig is not None:
                        pdf.savefig(self.last_test_fig)
                    if self.last_ext_fig is not None:
                        pdf.savefig(self.last_ext_fig)
                self.info_label.setText("PDF Exported Successfully.")
            except Exception as e:
                self.info_label.setText("Error exporting PDF: " + str(e))

    def send_report(self):
        if self.model is not None:
            self.report_plot()
            self.report_caption("QSAR Regression Model\n" + self.info_label.text())

# ----------------------------------------------------------------------
class QSARWorker(QThread):
    finished_signal = pyqtSignal(dict)
    error_signal = pyqtSignal(str)

    def __init__(self, data, external_data, config, parent=None):
        super().__init__(parent)
        self.data = data
        self.external_data = external_data
        self.config = config
        self.last_result = None

    def run(self):
        try:
            sel_alg = self.config["selected_algorithm"]
            norm_method = self.config["normalization_method"]
            imp_method = self.config["imputation_method"]
            cv_folds = self.config["cv_folds"]
            test_size = self.config["test_size"]
            tuning_method = self.config["tuning_method"]
            n_iter = self.config["n_iter"]
            hyper_str = self.config["hyperparameters"]
            algorithms = self.config["algorithms"]
            enable_fs = self.config.get("enable_feature_selection", False)
            num_features = self.config.get("num_features", 10)

            X_all = np.array(self.data.X)
            y_all = np.array(self.data.Y).ravel()
            metas_all = np.array(self.data.metas)

            X_train, X_test, y_train, y_test, metas_train, metas_test = train_test_split(
                X_all, y_all, metas_all, test_size=test_size, random_state=42
            )

            algo_name, algo_class = algorithms[sel_alg]
            is_classification = False
            scoring = "r2"
            if algo_name == "Logistic Regression":
                if len(np.unique(y_train)) == 2:
                    is_classification = True
                    scoring = "accuracy"
                else:
                    scoring = "r2"

            steps = []
            if imp_method != 0:
                strat = {1: "mean", 2: "median", 3: "most_frequent"}.get(imp_method, "mean")
                steps.append(("imputer", SimpleImputer(strategy=strat)))
            if norm_method == 1:
                steps.append(("scaler", StandardScaler()))
            elif norm_method == 2:
                steps.append(("scaler", MinMaxScaler()))

            if enable_fs:
                score_func = mutual_info_classif if is_classification else f_regression
                steps.append(("feature_selection", SelectKBest(score_func=score_func, k=num_features)))

            model_instance = algo_class()
            steps.append(("regressor", model_instance))
            pipeline = Pipeline(steps)

            hp = {}
            if hyper_str.strip():
                try:
                    hp = json.loads(hyper_str)
                except Exception as e:
                    raise Exception("Error parsing hyperparameters: " + str(e))

            if tuning_method == 1 and hp:
                tuner = GridSearchCV(
                    pipeline, param_grid=hp, cv=cv_folds, scoring=scoring,
                    n_jobs=-1, error_score='raise'
                )
                tuner.fit(X_train, y_train)
                best_pipeline = tuner.best_estimator_
                cv_score = tuner.best_score_
                tuning_info = f"Grid Search best CV {scoring}: {cv_score:.3f}\n"
            elif tuning_method == 2 and hp:
                tuner = RandomizedSearchCV(
                    pipeline, param_distributions=hp, cv=cv_folds, scoring=scoring,
                    n_iter=n_iter, n_jobs=-1, random_state=42, error_score='raise'
                )
                tuner.fit(X_train, y_train)
                best_pipeline = tuner.best_estimator_
                cv_score = tuner.best_score_
                tuning_info = f"Randomized Search best CV {scoring}: {cv_score:.3f}\n"
            else:
                cv_scores = cross_val_score(pipeline, X_train, y_train, cv=cv_folds, scoring=scoring)
                cv_score = np.mean(cv_scores)
                best_pipeline = pipeline.fit(X_train, y_train)
                tuning_info = f"CV {scoring} (no tuning): {cv_score:.3f}\n"

            if not is_classification:
                train_preds = best_pipeline.predict(X_train)
                train_r2 = r2_score(y_train, train_preds)
                train_rmse = np.sqrt(mean_squared_error(y_train, train_preds))
                train_mae = mean_absolute_error(y_train, train_preds)
                train_median_ae = median_absolute_error(y_train, train_preds)
                train_explained_var = explained_variance_score(y_train, train_preds)

                test_preds = best_pipeline.predict(X_test)
                test_r2 = r2_score(y_test, test_preds)
                test_rmse = np.sqrt(mean_squared_error(y_test, test_preds))
                test_mae = mean_absolute_error(y_test, test_preds)
                test_median_ae = median_absolute_error(y_test, test_preds)
                test_explained_var = explained_variance_score(y_test, test_preds)

                performance_text = (
                    f"{algo_name}: \n \n {tuning_info}\n"
                    f"Train R²: {train_r2:.3f}, RMSE: {train_rmse:.3f}, MAE: {train_mae:.3f},"
                    f" MedAE: {train_median_ae:.3f}, Expl.Var: {train_explained_var:.3f}\n"
                    f"Test R²: {test_r2:.3f}, RMSE: {test_rmse:.3f}, MAE: {test_mae:.3f}, "
                    f"MedAE: {test_median_ae:.3f}, Expl.Var: {test_explained_var:.3f}\n"
                )
                train_metrics = {
                    "R²": train_r2,
                    "RMSE": train_rmse,
                    "MAE": train_mae,
                    "Median AE": train_median_ae,
                    "Explained Variance": train_explained_var
                }
                test_metrics = {
                    "R²": test_r2,
                    "RMSE": test_rmse,
                    "MAE": test_mae,
                    "Median AE": test_median_ae,
                    "Explained Variance": test_explained_var
                }
            else:
                test_score = best_pipeline.score(X_test, y_test)
                performance_text = f"{algo_name}: {tuning_info} | Test Accuracy: {test_score:.3f}"
                train_metrics = {}
                test_metrics = {"Accuracy": test_score}

            ext_table = None
            external_metrics = {}
            X_ext = y_ext = None
            if self.external_data is not None:
                X_ext = np.array(self.external_data.X)
                y_ext = np.array(self.external_data.Y).ravel()
                metas_ext = np.array(self.external_data.metas)
                ext_preds = best_pipeline.predict(X_ext).reshape(-1, 1)
                if not is_classification:
                    ext_domain = Domain(list(self.data.domain.attributes) + [ContinuousVariable("Predicted")],
                                        self.data.domain.class_vars,
                                        self.data.domain.metas)
                else:
                    ext_domain = Domain(list(self.data.domain.attributes) + [DiscreteVariable("Predicted")],
                                        self.data.domain.class_vars,
                                        self.data.domain.metas)
                ext_table = Table(ext_domain, np.hstack([X_ext, ext_preds]), y_ext.reshape(-1, 1), metas_ext)
                ext_preds_full = best_pipeline.predict(X_ext)
                ext_r2 = r2_score(y_ext, ext_preds_full)
                ext_rmse = np.sqrt(mean_squared_error(y_ext, ext_preds_full))
                ext_mae = mean_absolute_error(y_ext, ext_preds_full)
                ext_median_ae = median_absolute_error(y_ext, ext_preds_full)
                ext_explained_var = explained_variance_score(y_ext, ext_preds_full)
                external_metrics = {
                    "R²": ext_r2,
                    "RMSE": ext_rmse,
                    "MAE": ext_mae,
                    "Median AE": ext_median_ae,
                    "Explained Variance": ext_explained_var
                }

            if not is_classification:
                new_domain = Domain(list(self.data.domain.attributes) + [ContinuousVariable("Predicted")],
                                    self.data.domain.class_vars,
                                    self.data.domain.metas)
                train_table = Table(new_domain, np.hstack([X_train, train_preds.reshape(-1, 1)]),
                                    y_train.reshape(-1, 1), metas_train)
                test_table = Table(new_domain, np.hstack([X_test, test_preds.reshape(-1, 1)]),
                                   y_test.reshape(-1, 1), metas_test)
            else:
                new_domain = Domain(list(self.data.domain.attributes) + [DiscreteVariable("Predicted")],
                                    self.data.domain.class_vars,
                                    self.data.domain.metas)
                train_preds = best_pipeline.predict(X_train).reshape(-1, 1)
                test_preds = best_pipeline.predict(X_test).reshape(-1, 1)
                train_table = Table(new_domain, np.hstack([X_train, train_preds]),
                                    y_train.reshape(-1, 1), metas_train)
                test_table = Table(new_domain, np.hstack([X_test, test_preds]),
                                   y_test.reshape(-1, 1), metas_test)

            result = {
                "model": best_pipeline,
                "train_table": train_table,
                "test_table": test_table,
                "external_table": ext_table,
                "pipeline": best_pipeline,
                "X_train": X_train,
                "y_train": y_train,
                "X_test": X_test,
                "y_test": y_test,
                "is_classification": is_classification,
                "performance_text": performance_text,
                "train_metrics": train_metrics,
                "test_metrics": test_metrics,
                "external_metrics": external_metrics,
                "cv_score": cv_score
            }
            if self.external_data is not None:
                result["X_ext"] = X_ext
                result["y_ext"] = y_ext

            self.last_result = result
            self.finished_signal.emit(result)
        except Exception as ex:
            self.error_signal.emit(str(ex))

if __name__ == "__main__":
    app = QApplication([])
    ow = OWQSARRegression()
    ow.show()
    app.exec_()

