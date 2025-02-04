from Orange.widgets import widget, gui, settings
from Orange.data import Table, Domain, ContinuousVariable, StringVariable
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors, AtomPairs, rdmolops, AllChem
import numpy as np
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA

class FingerprintWidget(widget.OWWidget):
    name = "Molecular Fingerprint Calculator"
    description = "Computes molecular fingerprints using RDKit and provides visualizations."
    icon = "icons/fingerprint.svg"
    priority = 10

    inputs = [("Molecule Data", Table, "set_data")]
    outputs = [("Fingerprints", Table)]

    # User settings
    fp_type = settings.Setting(0)  # Default to Morgan
    bit_size = settings.Setting(1024)
    radius = settings.Setting(2)

    def __init__(self):
        super().__init__()
        self.mainArea.hide()  # Hides the main area
        self.data = None
        
        # User Interface
        box = gui.widgetBox(self.controlArea, "Fingerprint Settings")
        self.fp_radio = gui.radioButtons(box, self, "fp_type", btnLabels=[
            "Morgan Fingerprint", "RDKit Fingerprint", "MACCS Keys", 
            "Atom Pair Fingerprint", "Topological Torsion Fingerprint", "Avalon Fingerprint"
        ])
        gui.spin(box, self, "bit_size", minv=128, maxv=4096, step=128, label="Bit Size")
        gui.spin(box, self, "radius", minv=1, maxv=5, step=1, label="Radius (Morgan)")
        gui.button(self.controlArea, self, "Compute Fingerprint", callback=self.compute_fingerprints)
        gui.button(self.controlArea, self, "Show Histogram", callback=self.show_histogram)
        gui.button(self.controlArea, self, "Show PCA Projection", callback=self.show_pca_projection)
    
    def set_data(self, data):
        """Receives input data."""
        self.data = data
        if self.data:
            self.compute_fingerprints()
    
    def compute_fingerprints(self):
        """Computes selected fingerprint type for molecules in the dataset."""
        if self.data is None:
            return
        
        smiles_col = self.data.domain.metas[0]  # Assuming SMILES is the first meta column
        molecules = [Chem.MolFromSmiles(str(row[smiles_col])) for row in self.data]
        
        fingerprints = []
        fingerprint_names = []
        
        def process_fingerprint(fp_values, prefix):
            if fp_values:
                array_fp = np.array([list(map(int, fp.ToBitString())) for fp in fp_values])
                fingerprints.append(array_fp)
                fingerprint_names.extend([f"{prefix}_{i}" for i in range(array_fp.shape[1])])
        
        fp_methods = [
            (lambda mol: rdMolDescriptors.GetMorganFingerprintAsBitVect(mol, self.radius, self.bit_size), "Morgan"),
            (lambda mol: rdmolops.RDKFingerprint(mol, fpSize=self.bit_size), "RDKit"),
            (lambda mol: rdMolDescriptors.GetMACCSKeysFingerprint(mol), "MACCS"),
            (lambda mol: AtomPairs.GetAtomPairFingerprintAsBitVect(mol), "AtomPair"),
            (lambda mol: AtomPairs.GetTopologicalTorsionFingerprintAsBitVect(mol), "Torsion"),
            (lambda mol: AllChem.GetAvalonFP(mol, nBits=self.bit_size) if hasattr(AllChem, 'GetAvalonFP') else None, "Avalon")
        ]
        
        selected_fp_method, fp_name = fp_methods[self.fp_type]
        fp_values = [selected_fp_method(mol) for mol in molecules if mol is not None]
        
        if fp_values and all(fp is not None for fp in fp_values):
            process_fingerprint(fp_values, fp_name)
        else:
            print(f"Warning: {fp_name} fingerprint computation failed for some molecules.")
        
        if fingerprints:
            self.fingerprint_data = np.hstack(fingerprints)
            domain = Domain([ContinuousVariable(name) for name in fingerprint_names], metas=[StringVariable("SMILES")])
            fingerprint_table = Table(domain, self.fingerprint_data, metas=np.array([str(row[smiles_col]) for row in self.data], dtype=object).reshape(-1, 1))
            self.send("Fingerprints", fingerprint_table)
    
    def show_histogram(self):
        """Displays a histogram of the top 20 most frequent fingerprint bits."""
        if hasattr(self, 'fingerprint_data'):
            bit_counts = np.sum(self.fingerprint_data, axis=0)
            top_indices = np.argsort(bit_counts)[-20:][::-1]
            plt.figure(figsize=(10, 5))
            plt.bar(range(20), bit_counts[top_indices], tick_label=[f"{i}" for i in top_indices])
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
