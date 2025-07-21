from Orange.widgets.widget import OWWidget, Input, Output
from Orange.widgets import gui
from Orange.data import Table, Domain, ContinuousVariable
from admet_ai import ADMETModel
from PyQt5.QtCore import QRunnable, QThreadPool, pyqtSignal, QObject

class WorkerSignals(QObject):
    finished = pyqtSignal(object)  # Will emit the predictions dataframe

class PredictionWorker(QRunnable):
    def __init__(self, model, smiles_list):
        super().__init__()
        self.model = model
        self.smiles_list = smiles_list
        self.signals = WorkerSignals()

    def run(self):
        # Run the prediction process in the background thread
        predictions = self.model.predict(smiles=self.smiles_list)
        predictions = predictions.round(4)
        self.signals.finished.emit(predictions)

class ADMETWidget(OWWidget):
    name = "ADMET Predictor"
    description = "Predicts ADMET properties from SMILES strings."
    icon = "icons/admet_predictor.png"  # Ensure you have an appropriate icon file
    priority = 81
    
    class Inputs:
        data = Input("Data", Table)

    class Outputs:
        data = Output("Data", Table)

    def __init__(self):
        super().__init__()
        self.data = None
        self.model = ADMETModel()
        self.infoLabel = gui.label(self.controlArea, self,
                                   "No data on input yet, waiting to get something.")
        # Use a global thread pool
        self.threadpool = QThreadPool.globalInstance()

    @Inputs.data
    def set_data(self, data):
        self.data = data
        if data is not None:
            self.infoLabel.setText(f"Received data with {len(data)} instances. Running predictions...")
            self.process_data()
        else:
            self.infoLabel.setText("No data on input yet, waiting to get something.")

    def process_data(self):
        # Find the 'SMILES' attribute in the data
        smiles_attr = next((var for var in self.data.domain if var.name.lower() == "smiles"), None)
        if smiles_attr is None:
            self.infoLabel.setText("No 'SMILES' column found in the input data.")
            return

        # Extract SMILES strings from the data
        smiles_list = [str(instance[smiles_attr]) for instance in self.data]

        # Create and start the worker for background prediction
        worker = PredictionWorker(self.model, smiles_list)
        worker.signals.finished.connect(self.handle_predictions)
        self.threadpool.start(worker)

    def handle_predictions(self, predictions):
        # Gather all original variable names from attributes, class_vars, and metas
        original_names = {var.name for var in self.data.domain.variables}

        # Create new ContinuousVariables for each ADMET property ensuring unique names
        new_attrs = []
        for col in predictions.columns:
            new_name = col
            suffix = " (predicted)"
            while new_name in original_names:
                new_name = col + suffix
                suffix += " (predicted)"
            new_var = ContinuousVariable(new_name)
            new_attrs.append(new_var)
            original_names.add(new_name)

        # Create a new domain by appending the new attributes to the original attributes.
        new_domain = Domain(
            self.data.domain.attributes + tuple(new_attrs),
            self.data.domain.class_vars,
            self.data.domain.metas
        )

        # Create a new table with the updated domain
        new_table = Table.from_table(new_domain, self.data)

        # Populate the new table with ADMET predictions
        for i, instance in enumerate(new_table):
            for new_var, col in zip(new_attrs, predictions.columns):
                instance[new_var.name] = predictions.iloc[i][col]

        # Update UI and send the new table as output
        self.infoLabel.setText("ADMET predictions completed.")
        self.Outputs.data.send(new_table)
