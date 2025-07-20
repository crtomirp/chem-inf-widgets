from Orange.widgets import widget, gui, settings
from Orange.data import Table, Domain, ContinuousVariable, StringVariable
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors, AtomPairs, rdmolops, AllChem
from rdkit.Avalon import pyAvalonTools  # Updated import for Avalon
import numpy as np
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA

# Import PyQt5 threading support
from PyQt5.QtCore import QThread, pyqtSignal

# Define a QThread worker for fingerprint computation.
class FingerprintWorker(QThread):
    progress = pyqtSignal(int)      # Emits the current progress (molecule index)
    finished = pyqtSignal(list)     # Emits the list of computed fingerprint objects

    def __init__(self, molecules, fp_method, parent=None):
        super().__init__(parent)
        self.molecules = molecules
        self.fp_method = fp_method

    def run(self):
        results = []
        total = len(self.molecules)
        for i, mol in enumerate(self.molecules):
            # Compute fingerprint if molecule is not None.
            if mol is not None:
                fp = self.fp_method(mol)
            else:
                fp = None
            results.append(fp)
            self.progress.emit(i + 1)
        self.finished.emit(results)

class FingerprintWidget(widget.OWWidget):
    name = "Fingerprint Calculator_1"
    description = "Computes molecular fingerprints using RDKit and provides visualizations."
    icon = "icons/fingerprint.png"
    priority = 23

    inputs = [("Molecule Data", Table, "set_data")]
    outputs = [("Fingerprints", Table)]

    # User settings
    fp_type = settings.Setting(0)  # Default to Morgan Fingerprint
    bit_size = settings.Setting(1024)
    radius = settings.Setting(2)
    remove_low_variance = settings.Setting(False)  # Option to remove low variance descriptors
    variance_threshold = settings.Setting(0.01)     # Threshold for variance removal

    def __init__(self):
        super().__init__()
        self.mainArea.hide()  # Hides the main area
        self.data = None
        self.worker = None
        
        # User Interface
        box = gui.widgetBox(self.controlArea, "Fingerprint Settings")
        self.fp_radio = gui.radioButtons(box, self, "fp_type", btnLabels=[
            "Morgan Fingerprint", "RDKit Fingerprint", "MACCS Keys","Avalon Fingerprint"
        ])
        gui.spin(box, self, "bit_size", minv=128, maxv=4096, step=128, label="Bit Size")
        gui.spin(box, self, "radius", minv=1, maxv=5, step=1, label="Radius (Morgan)")
        
        # Controls for low variance removal.
        gui.checkBox(box, self, "remove_low_variance", label="Remove low variance descriptors")
        gui.doubleSpin(box, self, "variance_threshold", minv=0.0, maxv=0.25, step=0.005,
                        label="Variance Threshold", controlWidth=100)
        
        gui.button(self.controlArea, self, "Compute Fingerprint", callback=self.compute_fingerprints)
        gui.button(self.controlArea, self, "Show Histogram", callback=self.show_histogram)
        gui.button(self.controlArea, self, "Show PCA Projection", callback=self.show_pca_projection)
    
    def set_data(self, data):
        """Receives input data."""
        self.data = data
        if self.data:
            self.compute_fingerprints()
    
    def compute_fingerprints(self):
        """Computes the selected fingerprint type in a background thread."""
        if self.data is None:
            return
        
        # Assuming SMILES is the first meta column.
        smiles_col = self.data.domain.metas[0]
        molecules = [Chem.MolFromSmiles(str(row[smiles_col])) for row in self.data]
        
        # Define fingerprint methods.
        fp_methods = [
            (lambda mol: rdMolDescriptors.GetMorganFingerprintAsBitVect(mol, self.radius, self.bit_size), "Morgan"),
            (lambda mol: rdmolops.RDKFingerprint(mol, fpSize=self.bit_size), "RDKit"),
            (lambda mol: rdMolDescriptors.GetMACCSKeysFingerprint(mol), "MACCS"),
            (lambda mol: pyAvalonTools.GetAvalonFP(mol, nBits=self.bit_size)
                      if hasattr(pyAvalonTools, 'GetAvalonFP') else None, "Avalon")
        ]
        
        selected_fp_method, fp_name = fp_methods[self.fp_type]
        
        # Initialize the progress bar for fingerprint computation.
        self.progressBarInit(len(molecules))
        
        # Create and start the background worker.
        self.worker = FingerprintWorker(molecules, selected_fp_method)
        self.worker.progress.connect(self.progressBarSet)
        self.worker.finished.connect(lambda fp_values: self.handle_finished(fp_values, fp_name, smiles_col))
        self.worker.start()
    
    def handle_finished(self, fp_values, fp_name, smiles_col):
        """Processes computed fingerprints, drops molecules that failed computation,
        applies low variance filtering if enabled, and sends the output."""
        self.progressBarFinished()
        
        # Identify indices where fingerprint computation was successful.
        valid_indices = [i for i, fp in enumerate(fp_values) if fp is not None]
        if not valid_indices:
            print("Warning: No fingerprints were successfully computed.")
            return
        
        # Filter out invalid fingerprints and corresponding data.
        valid_fp_values = [fp_values[i] for i in valid_indices]
        valid_smiles = [str(row[smiles_col]) for i, row in enumerate(self.data) if i in valid_indices]
        
        # Convert fingerprint bit vectors into an array.
        array_fp = np.array([list(map(int, fp.ToBitString())) for fp in valid_fp_values])
        fingerprint_names = [f"{fp_name}_{i}" for i in range(array_fp.shape[1])]
        
        self.fingerprint_data = array_fp
        
        # If low variance removal is enabled, filter out descriptors below the threshold.
        if self.remove_low_variance:
            # Report progress for the filtering step.
            self.progressBarInit(100)
            variances = np.var(self.fingerprint_data, axis=0)
            self.progressBarSet(50)  # Midway progress update.
            keep = variances >= self.variance_threshold
            self.fingerprint_data = self.fingerprint_data[:, keep]
            fingerprint_names = [name for name, k in zip(fingerprint_names, keep) if k]
            self.progressBarSet(100)
            self.progressBarFinished()
        
        # Build the output domain and table.
        domain = Domain([ContinuousVariable(name) for name in fingerprint_names],
                        metas=[StringVariable("SMILES")])
        fingerprint_table = Table(
            domain,
            self.fingerprint_data,
            metas=np.array(valid_smiles, dtype=object).reshape(-1, 1)
        )
        self.send("Fingerprints", fingerprint_table)
    
    def show_histogram(self):
        """Displays a histogram of the top 20 most frequent fingerprint bits."""
        if hasattr(self, 'fingerprint_data'):
            bit_counts = np.sum(self.fingerprint_data, axis=0)
            top_indices = np.argsort(bit_counts)[-20:][::-1]
            plt.figure(figsize=(10, 5))
            plt.bar(range(20), bit_counts[top_indices],
                    tick_label=[f"{i}" for i in top_indices])
            plt.xlabel("Fingerprint Bit Index")
            plt.ylabel("Frequency")
            plt.title("Top 20 Most Frequent Fingerprint Bits")
            plt.show()
    
    def show_pca_projection(self):
        """Displays a PCA projection of the selected fingerprint data."""
        if hasattr(self, 'fingerprint_data'):
            pca = PCA(n_components=2)
            transformed = pca.fit_transform(self.fingerprint_data)
            plt.figure(figsize=(8, 6))
            plt.scatter(transformed[:, 0], transformed[:, 1], alpha=0.5)
            plt.xlabel("Principal Component 1")
            plt.ylabel("Principal Component 2")
            plt.title("PCA Projection of Fingerprint Data")
            plt.show()

