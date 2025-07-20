import os
from Orange.widgets import gui, widget
from Orange.widgets.settings import Setting
from PyQt5.QtWebEngineWidgets import QWebEngineView, QWebEnginePage
from PyQt5.QtCore import QUrl, QTimer
from Orange.data import Table, Domain, StringVariable
from rdkit import Chem
from rdkit.Chem import AllChem

# Define search types for the widget
SEARCH_TYPES = ["substructure", "superstructure", "similarity", "exact"]

class CustomWebEnginePage(QWebEnginePage):
    """
    Custom QWebEnginePage to handle JavaScript console messages.
    """
    def javaScriptConsoleMessage(self, level, message, line_number, source_id):
        # Print JavaScript console messages for debugging
        print(f"JS Console: {message} (Line {line_number}) in {source_id}")

class OWSubstructureSearch(widget.OWWidget):
    """
    A widget for searching compounds based on substructure, superstructure, similarity, or exact match.
    """
    name = "Substructure Search"
    description = "Search compounds based on substructure, superstructure, similarity, or exact match."
    icon = "icons/substructuresearch.png"  # Ensure this icon exists in the 'icons' folder
    priority = 15

    # Define input and output channels
    inputs = [("Input Data", Table, "set_data")]
    outputs = [("Filtered Compounds", Table)]

    want_main_area = True  # Enable the main area for the widget
    resizing_enabled = True  # Allow resizing of the widget

    # Widget settings
    smiles_input = Setting("")  # Stores the input SMILES/SMARTS
    search_type = Setting(0)  # Default to "substructure"

    def __init__(self):
        """Initialize the widget and set up the GUI."""
        super().__init__()

        self.data = None  # Store input data

        # Set up the control panel
        self.setup_control_panel()

        # Set up the main area with JSME molecular editor
        self.setup_main_area()

        # Set up a timer to poll SMILES from the JSME editor
        self.polling_timer = QTimer()
        self.polling_timer.timeout.connect(self.poll_smiles)
        self.polling_timer.start(500)  # Poll every 500 ms

    def setup_control_panel(self):
        """Set up the control panel with search settings."""
        box = gui.widgetBox(self.controlArea, "Search Settings")
        gui.radioButtons(box, self, "search_type", SEARCH_TYPES, label="Search Type")
        gui.lineEdit(box, self, "smiles_input", label="Input SMILES/SMARTS:")
        gui.button(box, self, "Apply Search", callback=self.apply_search)

    def setup_main_area(self):
        """Set up the main area with the JSME molecular editor."""
        # Load the JSME HTML file
        html_path = os.path.abspath(os.path.join(os.path.dirname(__file__), "jsme/jsme_panel.html"))
        if not os.path.exists(html_path):
            raise FileNotFoundError(f"Cannot find the JSME HTML file at: {html_path}")

        # Create a QWebEngineView to display the JSME editor
        self.web_view = QWebEngineView(self.mainArea)
        self.web_view.setPage(CustomWebEnginePage(self.web_view))
        self.web_view.load(QUrl.fromLocalFile(html_path))
        self.mainArea.layout().addWidget(self.web_view)

    def poll_smiles(self):
        """Poll the SMILES string from the JSME editor."""
        self.web_view.page().runJavaScript("getSmiles();", self.update_smiles_field)

    def update_smiles_field(self, smiles):
        """Update the SMILES input field with the value from the JSME editor."""
        if smiles:
            self.smiles_input = smiles

    def set_data(self, data):
        """Handle input data."""
        self.data = data
        if data:
            print(f"Received input data with schema: {data.domain} and {len(data)} rows.")


    def apply_search(self):
        """Apply the search based on the selected search type and SMILES input."""
        if not self.data:
            self.error("No input data provided.")
            return

        # Detect the SMILES column in the input data
        smiles_col = next((var for var in self.data.domain.variables + self.data.domain.metas if var.name.lower() == "smiles"), None)
        if smiles_col is None:
            self.error("No SMILES column found in input data.")
            return

        # Parse the query SMILES/SMARTS
        smiles_query = self.smiles_input.strip()
        if not smiles_query:
            self.error("No query SMILES/SMARTS provided.")
            return

        try:
            search_type = SEARCH_TYPES[self.search_type]
            if search_type in ["substructure", "superstructure"]:
                query_mol = Chem.MolFromSmarts(smiles_query)
            else:
                query_mol = Chem.MolFromSmiles(smiles_query)

            if query_mol is None:
                raise ValueError("Invalid SMILES/SMARTS string for the selected search type.")
            
            query_mol.UpdatePropertyCache()
            _ = query_mol.GetRingInfo()  # Ensure ring info is updated

        except Exception as e:
            self.error(f"Error parsing query: {e}")
            return

        # Prepare filtered data and similarity scores
        filtered_data = []
        highlighted_atoms = []
        similarity_scores = []

        for row in self.data:
            smiles = row[smiles_col].value
            mol = Chem.MolFromSmiles(smiles)
            if not mol:
                continue

            match = self.match_molecule(mol, query_mol, search_type)
            similarity_score = None

            if match:
                if search_type == "substructure":
                    highlighted_atoms.append(",".join(map(str, mol.GetSubstructMatch(query_mol))))
                else:
                    highlighted_atoms.append("")

                if search_type == "similarity":
                    # Compute similarity index
                    fp1 = AllChem.GetMorganFingerprintAsBitVect(mol, 2)
                    fp2 = AllChem.GetMorganFingerprintAsBitVect(query_mol, 2)
                    similarity_score = AllChem.DataStructs.TanimotoSimilarity(fp1, fp2)

                similarity_scores.append(similarity_score if similarity_score is not None else "N/A")
                filtered_data.append(row)

        # Create and send the filtered data table
        if filtered_data:
            domain = self.data.domain
            highlighted_var = StringVariable("Highlighted Atoms")
            similarity_var = StringVariable("Similarity Index")

            new_domain = Domain(domain.attributes, domain.class_vars, list(domain.metas) + [highlighted_var, similarity_var])
            meta_data = [
                list(row.metas) + [highlighted_atoms[i], str(similarity_scores[i])]
                for i, row in enumerate(filtered_data)
            ]
            filtered_table = Table.from_numpy(new_domain, X=[row.x for row in filtered_data], metas=meta_data)
            self.send("Filtered Compounds", filtered_table)
        else:
            self.warning("No matching compounds found.")


    def match_molecule(self, mol, query_mol, search_type):
        """
        Match a molecule against the query molecule based on the search type.
        """
        if search_type == "substructure":
            return mol.HasSubstructMatch(query_mol)
        elif search_type == "superstructure":
            return query_mol.HasSubstructMatch(mol)
        elif search_type == "similarity":
            fp1 = AllChem.GetMorganFingerprintAsBitVect(mol, 2)
            fp2 = AllChem.GetMorganFingerprintAsBitVect(query_mol, 2)
            return AllChem.DataStructs.TanimotoSimilarity(fp1, fp2) > 0.1
        elif search_type == "exact":
            return Chem.MolToSmiles(mol) == Chem.MolToSmiles(query_mol)
        return False

if __name__ == "__main__":
    from Orange.widgets.utils.widgetpreview import WidgetPreview
    WidgetPreview(OWSubstructureSearch).run()

