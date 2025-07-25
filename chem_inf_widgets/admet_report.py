import matplotlib.pyplot as plt
import numpy as np
import io
import base64
import json
import os

from Orange.widgets.widget import OWWidget, Input
from Orange.widgets import gui
from Orange.data import Table
from rdkit import Chem
from rdkit.Chem import Draw

from PyQt5.QtWidgets import QSpinBox, QFileDialog
from PyQt5.QtCore import Qt, QRunnable, QThreadPool, pyqtSignal, QObject
from PyQt5.QtWebEngineWidgets import QWebEngineView

# WorkerSignals: signals emitted by the worker thread.
class WorkerSignals(QObject):
    finished = pyqtSignal(str)  # Emits the generated HTML report when done.

# ReportWorker: runs the report generation function in a background thread.
class ReportWorker(QRunnable):
    def __init__(self, orange_data, generate_report_func, molecule_size):
        super().__init__()
        self.orange_data = orange_data
        self.generate_report_func = generate_report_func
        self.molecule_size = molecule_size
        self.signals = WorkerSignals()

    def run(self):
        # Execute the heavy computation.
        html = self.generate_report_func(self.orange_data, self.molecule_size)
        self.signals.finished.emit(html)

class AdmetReport(OWWidget):
    """
    The ADMET Report widget displays molecules with ADMET properties in a styled HTML report.
    For each molecule, it shows:
      - The 2D chemical structure.
      - A radial (spider) plot of key ADMET properties.
      - A table of ADMET properties grouped into two columns.
      
    Additionally, an option to export the generated report as a PDF (with one compound per A4 page)
    is provided. The PDF export includes a CSS rule that scales (shrinks) the entire report so that
    each compound’s content fits within an A4 page.
    """
    name = "ADMET Report"
    description = ("Displays molecules with ADMET properties in categorized sections "
                   "using an HTML5 report layout. Each molecule gets its own report card "
                   "showing its structure, a corresponding radial plot, and numerical values "
                   "organized in two property columns.")
    icon = "icons/admetreport.png"
    priority = 85

    class Inputs:
        orange_data = Input("Filtered Compounds", Table)

    def __init__(self):
        """Initialize the widget, UI elements, and load the ADMET configuration."""
        super().__init__()

        self.orange_data = None
        self.molecule_size = 350  # Default pixel size for molecule images

        # Load ADMET configuration from JSON file.
        json_path = os.path.join(os.path.dirname(__file__), "admet.json")
        self.admet_config = self.load_admet_json(json_path)

        # Main UI: QWebEngineView is used to display the HTML report.
        self.web_view = QWebEngineView()
        self.mainArea.layout().addWidget(self.web_view)

        # Create a section for Report Settings in the control area.
        gui.label(self.controlArea, self, "Report Settings")
        # A spin box to adjust the molecule image size.
        self.size_selector = QSpinBox()
        self.size_selector.setMinimum(100)
        self.size_selector.setMaximum(400)
        self.size_selector.setValue(self.molecule_size)
        self.size_selector.valueChanged.connect(self.update_molecule_size)
        gui.widgetBox(self.controlArea, orientation="horizontal").layout().addWidget(self.size_selector)

        # Button to export the report as a PDF.
        gui.button(self.controlArea, self, "Export as PDF", callback=self.export_pdf)

        # Info label to display status messages.
        self.info_label = gui.label(self.controlArea, self, "Awaiting input data...")

    def load_admet_json(self, filepath):
        """
        Load the ADMET JSON configuration file.
        This configuration defines which properties to display and their alert ranges.
        """
        if not os.path.exists(filepath):
            print(f"Warning: ADMET JSON file not found at {filepath}.")
            return {}
        try:
            with open(filepath, "r") as f:
                return json.load(f)
        except Exception as e:
            print(f"Error loading ADMET JSON: {e}")
            return {}

    @Inputs.orange_data
    def set_orange_data(self, data: Table):
        """
        Input handler that receives the filtered compounds data.
        Once data is available, update the info label and display the report.
        """
        self.orange_data = data
        if data is not None:
            self.info_label.setText(f"Received {len(data)} molecules.")
            self.display_report()
        else:
            self.info_label.setText("No data received.")

    def update_molecule_size(self, value):
        """
        Update the molecule image size (in pixels) when the user changes the spin box value,
        and refresh the report.
        """
        self.molecule_size = value
        self.display_report()

    def export_pdf(self):
        """
        Export the currently displayed report as a PDF file.
        Uses the QWebEngineView's built-in PDF export function.
        """
        file_name, _ = QFileDialog.getSaveFileName(self, "Save Report as PDF", "", "PDF Files (*.pdf)")
        if not file_name:
            return

        def pdf_print_finished(pdf_data):
            with open(file_name, 'wb') as f:
                f.write(pdf_data)
            self.info_label.setText("PDF Exported Successfully.")

        self.web_view.page().printToPdf(pdf_print_finished)

    def get_alert_icon(self, alert_color):
        """
        Returns an HTML <span> element styled as a colored circle based on the given alert color.
        """
        color = {
            "green": "green",
            "orange": "orange",
            "red": "red"
        }.get(alert_color, "lightgray")
        return f'<span style="display:inline-block; width:1em; height:1em; background-color:{color}; border-radius:50%;"></span>'

    def generate_radial_plot(self, instance, domain):
        """
        Generates a radial (spider) plot image based on five ADMET properties.
        Returns a Base64-encoded PNG image.
        """
        try:
            clin_var = next((v for v in (domain.variables + domain.metas)
                             if v.name == 'ClinTox_drugbank_approved_percentile'), None)
            hERG_var = next((v for v in (domain.variables + domain.metas)
                             if v.name == 'hERG_drugbank_approved_percentile'), None)
            solubility_var = next((v for v in (domain.variables + domain.metas)
                                   if v.name == 'Solubility_AqSolDB_drugbank_approved_percentile'), None)
            bioavail_var = next((v for v in (domain.variables + domain.metas)
                                 if v.name == 'Bioavailability_Ma_drugbank_approved_percentile'), None)
            bbb_var = next((v for v in (domain.variables + domain.metas)
                            if v.name == 'BBB_Martins_drugbank_approved_percentile'), None)

            if not all([clin_var, hERG_var, solubility_var, bioavail_var, bbb_var]):
                return ""
            clin = float(instance[clin_var])
            hERG = float(instance[hERG_var])
            solubility = float(instance[solubility_var])
            bioavail = float(instance[bioavail_var])
            bbb = float(instance[bbb_var])
        except Exception as e:
            print(f"Error generating radial plot data: {e}")
            return ""

        non_toxic = 100 - clin
        hERG_save = 100 - hERG
        souluble = solubility
        biavailble = bioavail
        bbb_save = 100 - bbb

        labels = ['Non-Toxic', 'hERG_Save', 'Soluble', 'Bioavailble', 'BBB_Save']
        values = [non_toxic, hERG_save, souluble, biavailble, bbb_save]
        values += values[:1]
        num_vars = len(labels)
        angles = [n / float(num_vars) * 2 * np.pi for n in range(num_vars)]
        angles += angles[:1]

        fig = plt.figure(figsize=(3, 4))
        ax = fig.add_subplot(111, polar=True)
        angles_deg = np.degrees(angles[:-1])
        ax.set_thetagrids(angles_deg, labels, fontsize=12)
        ax.tick_params(axis='x', pad=15)
        ax.plot(angles, values, linewidth=2, linestyle='solid')
        ax.fill(angles, values, color='skyblue', alpha=0.4)
        ax.set_rlabel_position(30)
        plt.yticks([20, 40, 60, 80, 100], ["20", "40", "60", "80", "100"], color="grey", size=8)
        plt.ylim(0, 100)

        buf = io.BytesIO()
        plt.savefig(buf, format='png', bbox_inches='tight')
        plt.close(fig)
        buf.seek(0)
        img_data = base64.b64encode(buf.getvalue()).decode("utf-8")
        return img_data

    def generate_report(self, data, molecule_size):
        """
        Build the complete HTML report using the input data.
        """
        domain = data.domain
        smiles_col = next((var for var in domain.variables + domain.metas
                           if var.name.lower() == "smiles"), None)
        if smiles_col is None:
            return "<html><body>No 'SMILES' column found in data.</body></html>"

        html = """
        <!DOCTYPE html>
        <html>
        <head>
            <meta charset="utf-8">
            <title>ADMET Properties Report</title>
            <style>
                body { 
                    font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif; 
                    margin: 20px; 
                    padding: 20px; 
                    background-color: #f0f2f5;
                }
                h1 { 
                    text-align: center; 
                    margin-bottom: 40px;
                    color: #333;
                }
                .molecule-card { 
                    margin-bottom: 40px; 
                    background-color: #fff; 
                    border-radius: 12px; 
                    box-shadow: 0 2px 5px rgba(0,0,0,0.1); 
                    overflow: hidden;
                }
                .card-header { 
                    background-color: #4a90e2; 
                    color: #fff; 
                    padding: 15px;
                    text-align: center;
                }
                .card-header h2 {
                    margin: 0;
                    font-size: 1.4em;
                }
                .card-header p {
                    margin: 5px 0 0;
                    font-size: 0.9em;
                }
                .card-body {
                    display: flex;
                    flex-wrap: wrap;
                    padding: 20px;
                    justify-content: space-around;
                    align-items: center;
                }
                .structure-image, .radial-plot {
                    flex: 1 1 300px;
                    text-align: center;
                    margin: 10px;
                }
                .structure-image img, .radial-plot img {
                    max-width: 100%;
                    border-radius: 8px;
                    border: 1px solid #ddd;
                }
                .properties-container {
                    display: flex;
                    gap: 20px;
                    padding: 20px;
                    background-color: #fafafa;
                }
                .properties-column {
                    flex: 1 1 calc(50% - 20px);
                    background-color: #fff;
                    border-radius: 10px;
                    padding: 15px;
                    box-shadow: 0 1px 3px rgba(0,0,0,0.1);
                }
                .properties-column h4 {
                    margin-top: 0;
                    border-bottom: 1px solid #ddd;
                    padding-bottom: 5px;
                    color: #4a90e2;
                }
                table {
                    width: 100%;
                    border-collapse: collapse;
                    margin-top: 10px;
                }
                table th, table td {
                    padding: 8px;
                    text-align: left;
                    border-bottom: 1px solid #ddd;
                }
                table th {
                    background-color: #f5f5f5;
                }
                @media print {
                    body {
                        zoom: 0.6;
                    }
                    .molecule-card { 
                        page-break-after: always; 
                    }
                }
                @page {
                    size: A4;
                    margin: 20mm;
                }
            </style>
        </head>
        <body>
            <h1>ADMET Properties Report</h1>
        """
        for idx, instance in enumerate(data):
            try:
                smiles = instance[smiles_col].value
                mol = Chem.MolFromSmiles(smiles)
                if mol:
                    img = Draw.MolToImage(mol, size=(molecule_size, molecule_size))
                    buffered = io.BytesIO()
                    img.save(buffered, format="PNG")
                    img_str = base64.b64encode(buffered.getvalue()).decode("utf-8")
                else:
                    img_str = ""

                radial_plot_img = self.generate_radial_plot(instance, domain)
                categorized_properties = self.get_categorized_properties(instance, domain)
                categories = list(categorized_properties.items())
                midpoint = 4  # Fixed midpoint for two columns.
                left_categories = categories[:midpoint]
                right_categories = categories[midpoint:]
                
                html += f'''
                <div class="molecule-card">
                    <div class="card-header">
                        <h2>Molecule {idx+1}</h2>
                        <p>SMILES: {smiles}</p>
                    </div>
                    <div class="card-body">
                        <div class="structure-image">
                            <img src="data:image/png;base64,{img_str}" width="{molecule_size}" height="{molecule_size}" />
                        </div>
                        <div class="radial-plot">
                            <img src="data:image/png;base64,{radial_plot_img}" height="350" />
                        </div>
                    </div>
                    <div class="properties-container">
                        <div class="properties-column">
                '''
                for category, props in left_categories:
                    html += f'<h4>{category}</h4>'
                    html += '<table>'
                    html += '<tr><th>Property</th><th>Value</th><th>Alert</th></tr>'
                    for prop in props:
                        alert_icon = self.get_alert_icon(prop["alert"])
                        html += f'<tr><td>{prop["key"]}</td><td>{prop["value_str"]}</td><td>{alert_icon}</td></tr>'
                    html += '</table>'
                html += '</div>'  # End left column.
                html += '<div class="properties-column">'
                for category, props in right_categories:
                    html += f'<h4>{category}</h4>'
                    html += '<table>'
                    html += '<tr><th>Property</th><th>Value</th><th>Alert</th></tr>'
                    for prop in props:
                        alert_icon = self.get_alert_icon(prop["alert"])
                        html += f'<tr><td>{prop["key"]}</td><td>{prop["value_str"]}</td><td>{alert_icon}</td></tr>'
                    html += '</table>'
                html += '</div>'  # End right column.
                html += '</div>'  # End properties-container.
                html += '</div>'  # End molecule-card.
            except Exception as e:
                print(f"Error processing molecule: {e}")
                continue

        html += "</body></html>"
        return html

    def get_categorized_properties(self, instance, domain):
        """
        Retrieve and group ADMET properties for the given molecule based on the
        ADMET JSON configuration.
        """
        categorized = {}
        for prop in self.admet_config:
            category = prop.get("Category", "Uncategorized")
            prop_key = prop.get("Property")
            if not prop_key:
                continue

            var = next((v for v in (domain.variables + domain.metas) if v.name == prop_key), None)
            if var is None:
                continue

            try:
                value = instance[var]
            except Exception as e:
                value = "N/A"

            try:
                num_value = float(value)
                value_str = str(round(num_value, 2))
            except (ValueError, TypeError):
                value_str = str(value)

            percentile_key = prop.get("percentile")
            if percentile_key:
                var_pct = next((v for v in (domain.variables + domain.metas) if v.name == percentile_key), None)
                if var_pct is not None:
                    try:
                        pct_val = float(instance[var_pct])
                        pct_val = round(pct_val, 2)
                        value_str = f"{value_str} ({pct_val}%)"
                    except Exception as e:
                        pass

            alert_color = self.get_alert_color(value, prop)
            if category not in categorized:
                categorized[category] = []
            categorized[category].append({
                "key": prop_key,
                "value_str": value_str,
                "alert": alert_color
            })
        return categorized

    def get_alert_color(self, value, prop_config):
        """
        Determine the alert color ("green", "orange", "red", or "white")
        based on the property value and its configured ranges.
        """
        try:
            num_value = float(value)
        except (ValueError, TypeError):
            return ""
        ranges = prop_config.get("ranges", {})
        if "range_1" in ranges:
            rng = ranges["range_1"]
            if rng.get("min") is not None and rng.get("max") is not None:
                if rng["min"] <= num_value <= rng["max"]:
                    return "green"
        if "range_2" in ranges:
            rng = ranges["range_2"]
            if rng.get("min") is not None and rng.get("max") is not None:
                if rng["min"] <= num_value <= rng["max"]:
                    return "orange"
        if "range_3" in ranges:
            rng = ranges["range_3"]
            if rng.get("min") is not None and rng.get("max") is not None:
                if rng["min"] <= num_value <= rng["max"]:
                    return "red"
        return "white"

    def display_report(self):
        """
        Launch the background worker to generate the HTML report.
        """
        if self.orange_data is None:
            self.info_label.setText("No input data available.")
            return

        worker = ReportWorker(self.orange_data, self.generate_report, self.molecule_size)
        worker.signals.finished.connect(self.update_web_view)
        QThreadPool.globalInstance().start(worker)
        self.info_label.setText("Generating report in background...")

    def update_web_view(self, html):
        """
        Update the QWebEngineView with the generated HTML report.
        """
        self.web_view.setHtml(html)
        self.info_label.setText("Report generated.")

if __name__ == "__main__":
    from Orange.widgets.utils.widgetpreview import WidgetPreview
    WidgetPreview(AdmetReport).run()
