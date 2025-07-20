from Orange.data import Table, Domain, StringVariable
from Orange.widgets.widget import OWWidget, Input, Output
from Orange.widgets import gui
from rdkit import Chem
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit.Chem import Draw, AllChem
from rdkit.Chem.Draw import rdMolDraw2D
from PyQt5.QtWidgets import (
    QCheckBox, QPushButton, QVBoxLayout, QFileDialog, QComboBox
)
from PyQt5.QtCore import Qt, QThread, pyqtSignal, QObject, QBuffer, QIODevice
from PyQt5.QtGui import QPixmap, QImage, QPainter, QFont
from PyQt5.QtPrintSupport import QPrinter
from PyQt5.QtWebEngineWidgets import QWebEngineView
import numpy as np
import base64
from io import BytesIO
from PIL import Image

# ---------------------------- Helper Functions ----------------------------

def pil_to_pixmap(pil_img):
    """Convert PIL Image to QPixmap"""
    if pil_img.mode == "RGBA":
        pass
    elif pil_img.mode == "RGB":
        pil_img = pil_img.convert("RGBA")
    data = pil_img.tobytes("raw", "RGBA")
    qimage = QImage(data, pil_img.size[0], pil_img.size[1], QImage.Format_RGBA8888)
    return QPixmap.fromImage(qimage)

def pixmap_to_base64(pixmap):
    """Convert QPixmap to base64 string"""
    buffer = QBuffer()
    buffer.open(QIODevice.ReadWrite)
    pixmap.save(buffer, "PNG")
    return base64.b64encode(buffer.data()).decode("utf-8")

def generate_placeholder(width=200, height=200):
    """Generate placeholder image for invalid structures"""
    pixmap = QPixmap(width, height)
    pixmap.fill(Qt.white)
    painter = QPainter(pixmap)
    painter.setPen(Qt.red)
    font = QFont("Arial", 20)
    painter.setFont(font)
    painter.drawText(pixmap.rect(), Qt.AlignCenter, "Invalid")
    painter.end()
    return pixmap

def break_line(text, max_len=30):
    """Break text every max_len characters."""
    return "\n".join(text[i:i+max_len] for i in range(0, len(text), max_len))

# ---------------------------- Standardization Worker ----------------------------

class StandardizationWorker(QObject):
    progress = pyqtSignal(int)
    finished = pyqtSignal(str, list, list, list)

    def __init__(self, smiles_list, operations, depiction_mode):
        super().__init__()
        self.smiles_list = smiles_list
        self.operations = operations
        self.depiction_mode = depiction_mode

    def process(self):
        orig_smiles = []
        std_smiles = []
        logs = []
        orig_images = []
        std_images = []

        for idx, smile in enumerate(self.smiles_list):
            orig_smiles.append(smile)
            log = []
            
            # Process original structure
            orig_mol = Chem.MolFromSmiles(smile, sanitize=False)
            if orig_mol:
                try:
                    # Partial sanitization excluding aromaticity and kekulization
                    #Chem.SanitizeMol(orig_mol, Chem.SanitizeFlags.SANITIZE_ALL ^ 
                    #Chem.SanitizeMol(orig_mol, 
                    #               (Chem.SanitizeFlags.SANITIZE_KEKULIZE |
                    #                Chem.SanitizeFlags.SANITIZE_SETAROMATICITY))
                    AllChem.Compute2DCoords(orig_mol)
                except:
                    orig_mol = None

            if orig_mol:
                drawer = rdMolDraw2D.MolDraw2DCairo(200, 200)
                
                # Original depiction preserves input bonds exactly
                drawer.drawOptions().useAromaticCircles = False
                drawer.drawOptions().prepareMolsBeforeDrawing = False
                
                try:
                    drawer.DrawMolecule(orig_mol)
                except:
                    orig_mol = None

            if orig_mol:
                drawer.FinishDrawing()
                orig_img = Image.open(BytesIO(drawer.GetDrawingText()))
                orig_pix = pil_to_pixmap(orig_img)
            else:
                orig_pix = generate_placeholder()
                log.append("Invalid original SMILES")

            # Standardization process
            std_mol = orig_mol
            if orig_mol:
                try:
                    # Full sanitization for standardization
                    Chem.SanitizeMol(orig_mol)
                    for op in self.operations:
                        initial_smiles = Chem.MolToSmiles(std_mol)
                        if op == "Cleanup":
                            std_mol = rdMolStandardize.Cleanup(std_mol)
                        elif op == "Normalize":
                            normalizer = rdMolStandardize.Normalizer()
                            std_mol = normalizer.normalize(std_mol)
                        elif op == "MetalDisconnector":
                            disconnector = rdMolStandardize.MetalDisconnector()
                            std_mol = disconnector.Disconnect(std_mol)
                        elif op == "LargestFragmentChooser":
                            chooser = rdMolStandardize.LargestFragmentChooser()
                            std_mol = chooser.choose(std_mol)
                        elif op == "Reionizer":
                            reionizer = rdMolStandardize.Reionizer()
                            std_mol = reionizer.reionize(std_mol)
                        elif op == "Uncharger":
                            uncharger = rdMolStandardize.Uncharger()
                            std_mol = uncharger.uncharge(std_mol)
                        
                        final_smiles = Chem.MolToSmiles(std_mol)
                        if initial_smiles != final_smiles:
                            log.append(f"{op} applied: {initial_smiles} → {final_smiles}")
                except Exception as e:
                    log.append(f"Standardization error: {str(e)}")
                    std_mol = None

            # Process standardized structure
            if std_mol:
                std_smi = Chem.MolToSmiles(std_mol)
                AllChem.Compute2DCoords(std_mol)
                drawer = rdMolDraw2D.MolDraw2DCairo(200, 200)

                # Handle user-selected depiction mode
                if self.depiction_mode == "As Is":
                    use_aromatic = any(c.islower() for c in std_smi)
                elif self.depiction_mode == "Aromatic":
                    use_aromatic = True
                else:  # Kekulized
                    use_aromatic = False

                if use_aromatic:
                    drawer.drawOptions().useAromaticCircles = True
                    drawer.drawOptions().prepareMolsBeforeDrawing = False
                else:
                    drawer.drawOptions().useAromaticCircles = False
                    drawer.drawOptions().prepareMolsBeforeDrawing = True

                try:
                    drawer.DrawMolecule(std_mol)
                    drawer.FinishDrawing()
                    std_img = Image.open(BytesIO(drawer.GetDrawingText()))
                    std_pix = pil_to_pixmap(std_img)
                    std_smiles.append(std_smi)
                except:
                    std_pix = generate_placeholder()
                    std_smiles.append("")
                    log.append("Rendering failed")
            else:
                std_pix = generate_placeholder()
                std_smiles.append("")
                log.append("Standardization failed")

            # Collect results
            orig_images.append(pixmap_to_base64(orig_pix))
            std_images.append(pixmap_to_base64(std_pix))
            logs.append("\n".join(log) if log else "No changes")
            self.progress.emit(idx + 1)

        # Generate HTML report
        html = self.generate_html(orig_smiles, std_smiles, logs, orig_images, std_images)
        self.finished.emit(html, orig_smiles, std_smiles, logs)

    def generate_html(self, orig_smiles, std_smiles, logs, orig_images, std_images):
        html = f"""
        <html>
        <head>
            <style>
                body {{ font-family: Arial, sans-serif; margin: 20px; }}
                table {{ width: 100%; border-collapse: collapse; }}
                th, td {{ border: 1px solid #ddd; padding: 10px; text-align: center; width: 250px; }}
                th {{ background-color: #4CAF50; color: white; }}
                img {{ max-width: 200px; max-height: 200px; }}
                .smiles {{ font-family: monospace; margin: 5px 0; }}
                .log {{ text-align: left; padding: 10px; }}
            </style>
        </head>
        <body>
            <h2>Molecular Standardization Report</h2>
            <table>
                <tr>
                    <th>Original Structure</th>
                    <th>Standardized Structure</th>
                    <th>Applied Changes</th>
                </tr>
        """
        
        for o_smi, s_smi, log, o_img, s_img in zip(orig_smiles, std_smiles, logs, orig_images, std_images):
            broken_log = break_line(log, 25)
            html += f"""
                <tr>
                    <td>
                        <img src="data:image/png;base64,{o_img}"><br>
                        <div class="smiles">{o_smi}</div>
                    </td>
                    <td>
                        <img src="data:image/png;base64,{s_img}"><br>
                        <div class="smiles">{s_smi if s_smi else 'N/A'}</div>
                    </td>
                    <td><pre>{broken_log}</pre></td>
                </tr>
            """
        
        html += """
            </table>
        </body>
        </html>
        """
        return html

# ---------------------------- Main Widget ----------------------------

class StandardizeMoleculesWidget(OWWidget):
    name = "Molecular Standardizer"
    description = "Standardizes molecules with exact input preservation and output options"
    icon = "icons/standardizer.png"
    priority = 5

    class Inputs:
        data = Input("Data", Table)

    class Outputs:
        data = Output("Standardized Data", Table)

    want_main_area = True
    operations = [
        "Cleanup",
        "Normalize",
        "MetalDisconnector",
        "LargestFragmentChooser",
        "Reionizer",
        "Uncharger"
    ]
    depiction_modes = ["As Is", "Aromatic", "Kekulized"]

    def __init__(self):
        super().__init__()
        self.smiles_data = []
        self.selected_ops = set()
        self.depiction_mode = "As Is"

        # GUI Setup
        self.controlArea.layout().addWidget(gui.label(self.controlArea, self, "Standardization Steps:"))
        self.checkboxes = {}
        for op in self.operations:
            cb = QCheckBox(op)
            cb.stateChanged.connect(self.update_selections)
            self.controlArea.layout().addWidget(cb)
            self.checkboxes[op] = cb

        self.depiction_combo = QComboBox()
        self.depiction_combo.addItems(self.depiction_modes)
        self.depiction_combo.currentTextChanged.connect(self.update_depiction_mode)
        self.controlArea.layout().addWidget(gui.label(self.controlArea, self, "Standardized Depiction:"))
        self.controlArea.layout().addWidget(self.depiction_combo)

        self.run_btn = QPushButton("Run Standardization")
        self.run_btn.clicked.connect(self.start_processing)
        self.controlArea.layout().addWidget(self.run_btn)

        self.export_btn = QPushButton("Export PDF Report")
        self.export_btn.clicked.connect(self.export_pdf)
        self.controlArea.layout().addWidget(self.export_btn)
        gui.rubber(self.controlArea)

        self.web_view = QWebEngineView()
        self.mainArea.layout().addWidget(self.web_view)

    @Inputs.data
    def set_data(self, data):
        self.smiles_data = []
        if data and "SMILES" in [var.name for var in data.domain.metas]:
            self.smiles_data = [str(row["SMILES"]) for row in data]

    def update_selections(self):
        self.selected_ops = {op for op, cb in self.checkboxes.items() if cb.isChecked()}

    def update_depiction_mode(self, mode):
        self.depiction_mode = mode

    def start_processing(self):
        if not self.smiles_data:
            self.web_view.setHtml("<h3 style='color: red;'>No input data!</h3>")
            return

        self.progressBarInit()
        self.thread = QThread()
        self.worker = StandardizationWorker(
            self.smiles_data, 
            list(self.selected_ops),
            self.depiction_mode
        )
        self.worker.moveToThread(self.thread)

        self.thread.started.connect(self.worker.process)
        self.worker.finished.connect(self.process_finished)
        self.worker.finished.connect(self.thread.quit)
        self.worker.progress.connect(self.update_progress)
        
        self.thread.start()
        self.run_btn.setEnabled(False)

    def process_finished(self, html, orig, std, logs):
        self.web_view.setHtml(html)
        self.run_btn.setEnabled(True)
        self.create_output_table(orig, std, logs)
        self.progressBarFinished()

    def create_output_table(self, orig, std, logs):
        domain = Domain([], metas=[
            StringVariable("Original SMILES"),
            StringVariable("Standardized SMILES"),
            StringVariable("Standardization Log")
        ])
        metas = np.array([[o, s, l] for o, s, l in zip(orig, std, logs)], dtype=object)
        output_table = Table.from_numpy(domain, X=np.empty((len(orig), 0)), metas=metas)
        self.Outputs.data.send(output_table)

    def update_progress(self, value):
        self.progressBarSet(100 * value / len(self.smiles_data))

    def export_pdf(self):
        filename, _ = QFileDialog.getSaveFileName(
            self, "Save Report", "", "PDF Files (*.pdf)")
        if filename:
            printer = QPrinter(QPrinter.HighResolution)
            printer.setOutputFormat(QPrinter.PdfFormat)
            printer.setOutputFileName(filename)
            self.web_view.page().print(printer, lambda ok: None)

if __name__ == "__main__":
    from Orange.widgets.utils.widgetpreview import WidgetPreview
    WidgetPreview(StandardizeMoleculesWidget).run()