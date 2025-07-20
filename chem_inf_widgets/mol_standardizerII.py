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
    """Convert PIL Image to QPixmap."""
    if pil_img.mode == "RGBA":
        pass
    elif pil_img.mode == "RGB":
        pil_img = pil_img.convert("RGBA")
    data = pil_img.tobytes("raw", "RGBA")
    qimage = QImage(data, pil_img.size[0], pil_img.size[1], QImage.Format_RGBA8888)
    return QPixmap.fromImage(qimage)

def pixmap_to_base64(pixmap):
    """Convert QPixmap to base64 string."""
    buffer = QBuffer()
    buffer.open(QIODevice.ReadWrite)
    pixmap.save(buffer, "PNG")
    return base64.b64encode(buffer.data()).decode("utf-8")

def generate_placeholder(width=200, height=200):
    """Generate placeholder image for invalid structures."""
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
        orig_smiles, std_smiles, logs = [], [], []
        orig_images, std_images = [], []

        for idx, smile in enumerate(self.smiles_list):
            orig_smiles.append(smile)
            log = []

            # --- Original molecule depiction ---
            orig_mol = Chem.MolFromSmiles(smile, sanitize=False)
            if orig_mol:
                try:
                    AllChem.Compute2DCoords(orig_mol)
                except:
                    orig_mol = None

            if orig_mol:
                drawer = rdMolDraw2D.MolDraw2DCairo(200, 200)
                # No aromatic circles for original
                drawer.drawOptions().circleAtoms = False
                drawer.drawOptions().prepareMolsBeforeDrawing = False
                try:
                    drawer.DrawMolecule(orig_mol)
                    drawer.FinishDrawing()
                    orig_img = Image.open(BytesIO(drawer.GetDrawingText()))
                    orig_pix = pil_to_pixmap(orig_img)
                except:
                    orig_pix = generate_placeholder()
                    log.append("Original rendering failed")
            else:
                orig_pix = generate_placeholder()
                log.append("Invalid original SMILES")

            # --- Standardization steps ---
            std_mol = orig_mol
            if std_mol:
                try:
                    Chem.SanitizeMol(std_mol)
                    for op in self.operations:
                        before = Chem.MolToSmiles(std_mol)
                        if op == "Cleanup":
                            std_mol = rdMolStandardize.Cleanup(std_mol)
                        elif op == "Normalize":
                            std_mol = rdMolStandardize.Normalizer().normalize(std_mol)
                        elif op == "MetalDisconnector":
                            std_mol = rdMolStandardize.MetalDisconnector().Disconnect(std_mol)
                        elif op == "LargestFragmentChooser":
                            std_mol = rdMolStandardize.LargestFragmentChooser().choose(std_mol)
                        elif op == "Reionizer":
                            std_mol = rdMolStandardize.Reionizer().reionize(std_mol)
                        elif op == "Uncharger":
                            std_mol = rdMolStandardize.Uncharger().uncharge(std_mol)
                        after = Chem.MolToSmiles(std_mol)
                        if before != after:
                            log.append(f"{op}: {before} → {after}")
                except Exception as e:
                    log.append(f"Standardization error: {e}")
                    std_mol = None

            # --- Standardized molecule depiction ---
            if std_mol:
                std_smi = Chem.MolToSmiles(std_mol)
                std_smiles.append(std_smi)
                AllChem.Compute2DCoords(std_mol)
                drawer = rdMolDraw2D.MolDraw2DCairo(200, 200)

                # Decide circleAtoms vs kekulize
                if self.depiction_mode == "As Is":
                    use_aromatic = any(c.islower() for c in std_smi)
                elif self.depiction_mode == "Aromatic":
                    use_aromatic = True
                else:  # Kekulized
                    use_aromatic = False

                drawer.drawOptions().circleAtoms = use_aromatic
                drawer.drawOptions().prepareMolsBeforeDrawing = not use_aromatic

                try:
                    drawer.DrawMolecule(std_mol)
                    drawer.FinishDrawing()
                    std_img = Image.open(BytesIO(drawer.GetDrawingText()))
                    std_pix = pil_to_pixmap(std_img)
                except:
                    std_pix = generate_placeholder()
                    log.append("Standardized rendering failed")
            else:
                std_smiles.append("")
                std_pix = generate_placeholder()
                log.append("Standardization failed")

            # Collect and emit progress
            orig_images.append(pixmap_to_base64(orig_pix))
            std_images.append(pixmap_to_base64(std_pix))
            logs.append("\n".join(log) or "No changes")
            self.progress.emit(idx + 1)

        html = self.generate_html(orig_smiles, std_smiles, logs, orig_images, std_images)
        self.finished.emit(html, orig_smiles, std_smiles, logs)

    def generate_html(self, orig_smiles, std_smiles, logs, orig_images, std_images):
        rows = ""
        for o, s, lg, oi, si in zip(orig_smiles, std_smiles, logs, orig_images, std_images):
            bl = break_line(lg, 25)
            rows += f"""
            <tr>
              <td><img src="data:image/png;base64,{oi}"><div class="smiles">{o}</div></td>
              <td><img src="data:image/png;base64,{si}"><div class="smiles">{s or 'N/A'}</div></td>
              <td><pre>{bl}</pre></td>
            </tr>"""
        return f"""
        <html><head><style>
        body{{font-family:Arial}}table{{width:100%;border-collapse:collapse}}
        th,td{{border:1px solid #ddd;padding:8px;text-align:center}}
        th{{background:#4CAF50;color:white}}.smiles{{font-family:monospace}}
        </style></head><body>
        <h2>Molecular Standardization Report</h2>
        <table>
          <tr><th>Original</th><th>Standardized</th><th>Log</th></tr>
          {rows}
        </table></body></html>
        """

# ---------------------------- Main Widget ----------------------------

class StandardizeMoleculesWidget(OWWidget):
    name = "Molecular StandardizerII"
    description = "Standardizes molecules with exact input preservation and output options"
    icon = "icons/standardizer.png"
    priority = 5

    class Inputs:
        data = Input("Data", Table)
    class Outputs:
        data = Output("Standardized Data", Table)

    operations = ["Cleanup", "Normalize", "MetalDisconnector",
                  "LargestFragmentChooser", "Reionizer", "Uncharger"]
    depiction_modes = ["As Is", "Aromatic", "Kekulized"]

    def __init__(self):
        super().__init__()
        self.smiles_data = []
        self.selected_ops = set()
        self.depiction_mode = "As Is"

        # Controls
        self.controlArea.layout().addWidget(
            gui.label(self.controlArea, self, "Steps:"))
        self.checkboxes = {}
        for op in self.operations:
            cb = QCheckBox(op)
            cb.stateChanged.connect(self.update_selections)
            self.controlArea.layout().addWidget(cb)
            self.checkboxes[op] = cb

        self.controlArea.layout().addWidget(
            gui.label(self.controlArea, self, "Depiction:"))
        self.depiction_combo = QComboBox()
        self.depiction_combo.addItems(self.depiction_modes)
        self.depiction_combo.currentTextChanged.connect(self.update_depiction_mode)
        self.controlArea.layout().addWidget(self.depiction_combo)

        self.run_btn = QPushButton("Run Standardization")
        self.run_btn.clicked.connect(self.start_processing)
        self.controlArea.layout().addWidget(self.run_btn)

        self.export_btn = QPushButton("Export PDF Report")
        self.export_btn.clicked.connect(self.export_pdf)
        self.controlArea.layout().addWidget(self.export_btn)
        gui.rubber(self.controlArea)

        # Web view
        self.web_view = QWebEngineView()
        self.mainArea.layout().addWidget(self.web_view)

    @Inputs.data
    def set_data(self, data):
        self.smiles_data = []
        if data and "SMILES" in [var.name for var in data.domain.metas]:
            self.smiles_data = [str(r["SMILES"]) for r in data]

    def update_selections(self):
        self.selected_ops = {op for op, cb in self.checkboxes.items() if cb.isChecked()}

    def update_depiction_mode(self, mode):
        self.depiction_mode = mode

    def start_processing(self):
        if not self.smiles_data:
            self.web_view.setHtml("<h3 style='color:red;'>No input data!</h3>")
            return
        self.progressBarInit()
        self.thread = QThread()
        self.worker = StandardizationWorker(self.smiles_data,
                                            list(self.selected_ops),
                                            self.depiction_mode)
        self.worker.moveToThread(self.thread)
        self.thread.started.connect(self.worker.process)
        self.worker.progress.connect(self.update_progress)
        self.worker.finished.connect(self.process_finished)
        self.worker.finished.connect(self.thread.quit)
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
            StringVariable("Log")
        ])
        metas = np.array([[o, s, l] for o, s, l in zip(orig, std, logs)], dtype=object)
        out = Table.from_numpy(domain, X=np.empty((len(orig), 0)), metas=metas)
        self.Outputs.data.send(out)

    def update_progress(self, val):
        self.progressBarSet(100 * val / len(self.smiles_data))

    def export_pdf(self):
        fname, _ = QFileDialog.getSaveFileName(self, "Save Report", "", "PDF Files (*.pdf)")
        if fname:
            printer = QPrinter(QPrinter.HighResolution)
            printer.setOutputFormat(QPrinter.PdfFormat)
            printer.setOutputFileName(fname)
            self.web_view.page().print(printer, lambda ok: None)
