from PyQt5.QtWidgets import QListWidget, QAbstractItemView, QFileDialog
from Orange.widgets.widget import OWWidget, Output
from Orange.widgets.settings import Setting
from Orange.widgets.gui import widgetBox, button, lineEdit, label
import os
from rdkit import Chem
from rdkit.Chem import Descriptors
from Orange.data import Table, Domain, ContinuousVariable, StringVariable
from rdkit import RDLogger

# Disable RDKit warnings to avoid cluttering the console
RDLogger.DisableLog('rdApp.*')

class OWSDFReader(OWWidget):
    """
    A widget for analyzing SDF files, selecting properties, and generating a data table.
    """
    name = "SDF Reader"
    description = "Analyzes an SDF file, allows selection of properties, and reads data based on the selection."
    icon = "icons/sdfreaderwidget.svg"  # Ensure this icon exists in the 'icons' folder
    priority = 2

    # Define outputs
    class Outputs:
        data = Output("Selected Data Table", Table)

    # Widget settings
    sdf_file_path = Setting("")  # Stores the path to the SDF file
    selected_properties = Setting([])  # Stores the list of selected properties

    def __init__(self):
        """Initialize the widget and set up the GUI."""
        super().__init__()

        # Create GUI elements
        self.setup_ui()
        self.mainArea.hide()  # Hides the main area

    def setup_ui(self):
        """Set up the user interface for the widget."""
        # File input section
        box = widgetBox(self.controlArea, "File Input")
        self.file_input = lineEdit(
            box, self, "sdf_file_path", label="SDF File Path:",
            orientation="horizontal", callback=self.analyze_file
        )
        button(box, self, "Browse...", callback=self.browse_file)

        # Property selection section
        self.property_selector = QListWidget(self.controlArea)
        self.property_selector.setSelectionMode(QAbstractItemView.MultiSelection)
        widgetBox(self.controlArea).layout().addWidget(self.property_selector)

        # Info label and load button
        self.info_label = label(self.controlArea, self, "Ready to analyze SDF file.")
        button(self.controlArea, self, "Read File", callback=self.load_file)

    def browse_file(self):
        """Open a file dialog to select an SDF file."""
        file_path, _ = QFileDialog.getOpenFileName(self, "Select SDF File", "", "SDF Files (*.sdf)")
        if file_path:
            self.sdf_file_path = file_path
            self.analyze_file()  # Analyze the file after selection

    def analyze_file(self):
        """Analyze the SDF file to extract available properties."""
        if not self.sdf_file_path or not os.path.isfile(self.sdf_file_path):
            self.info_label.setText("Invalid file path. Please select a valid file.")
            return

        try:
            suppl = Chem.SDMolSupplier(self.sdf_file_path)
            properties_set = set()

            # Collect all unique properties from the molecules in the SDF file
            for mol in suppl:
                if mol is None:
                    continue
                properties_set.update(mol.GetPropsAsDict().keys())

            # Update the property selector with the found properties
            self.property_selector.clear()
            for prop in sorted(properties_set):
                self.property_selector.addItem(prop)

            self.info_label.setText(f"Found {len(properties_set)} properties in the file.")
        except Exception as e:
            self.info_label.setText(f"Error analyzing file: {str(e)}")

    def load_file(self):
        """Load the SDF file and process molecules based on selected properties."""
        if not self.sdf_file_path or not os.path.isfile(self.sdf_file_path):
            self.info_label.setText("Invalid file path. Please select a valid file.")
            return

        # Get the selected properties from the list widget
        self.selected_properties = [item.text() for item in self.property_selector.selectedItems()]

        try:
            suppl = Chem.SDMolSupplier(self.sdf_file_path)
            data_list = []

            # Process each molecule in the SDF file
            for mol in suppl:
                if mol is None:
                    print("Invalid molecule detected, skipping...")
                    continue

                try:
                    # Extract basic molecular properties
                    atom_count = mol.GetNumAtoms()
                    mol_weight = Descriptors.MolWt(mol)
                    smiles = Chem.MolToSmiles(mol)

                    # Extract selected properties
                    properties = mol.GetPropsAsDict()
                    selected_properties = {prop: properties.get(prop, None) for prop in self.selected_properties}

                    # Append the data to the list
                    data_list.append([atom_count, mol_weight, smiles, selected_properties])
                except Exception as e:
                    print(f"Error processing molecule: {e}")
                    continue

            # Create an Orange data table from the processed data
            self.create_data_table(data_list)
            self.info_label.setText(f"Successfully processed {len(data_list)} molecules.")
        except Exception as e:
            self.info_label.setText(f"Error reading file: {str(e)}")

    def create_data_table(self, data_list):
        """Create an Orange data table based on the processed data."""
        # Define attribute and meta variables
        attribute_vars = [
            ContinuousVariable("Atom Count"), 
            ContinuousVariable("Molecular Weight")
        ]
        meta_vars = [
            StringVariable("SMILES")
        ]

        # Add selected properties as meta variables
        for prop in self.selected_properties:
            meta_vars.append(StringVariable(prop))

        # Create the domain for the Orange table
        domain = Domain(attribute_vars, metas=meta_vars)

        # Prepare data for the table
        attributes = []
        metas = []

        for row in data_list:
            atom_count, mol_weight, smiles, selected_properties = row
            attributes_row = [atom_count, mol_weight]
            meta_row = [smiles] + [selected_properties.get(var.name, None) for var in meta_vars if var.name != "SMILES"]

            attributes.append(attributes_row)
            metas.append(meta_row)

        # Create and send the Orange table
        table = Table.from_numpy(domain, X=attributes, metas=metas)
        self.Outputs.data.send(table)

if __name__ == "__main__":
    from Orange.canvas import main
    main()
