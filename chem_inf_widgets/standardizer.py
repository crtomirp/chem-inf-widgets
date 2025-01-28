from Orange.data import Table, Domain, StringVariable
from Orange.widgets.widget import OWWidget
from Orange.widgets import gui
from rdkit import Chem
from rdkit.Chem.MolStandardize import rdMolStandardize
import numpy as np
from PyQt5.QtWidgets import QCheckBox, QVBoxLayout, QPushButton

class StandardizeMoleculesWidget(OWWidget):
    name = "Molecular Standardization"
    description = "Reads SMILES from input table and outputs standardized SMILES."
    icon = "icons/standardizer.png"
    priority = 10

    inputs = [
        ("Input Table", Table, "set_data")
    ]

    outputs = [
        ("standardized_table", Table)  # Corrected output channel name
    ]

    want_main_area = False

    operations = [
        "ValidateSmiles",
        "Cleanup",
        "Normalize",
        "MetalDisconnector",
        "LargestFragmentChooser",
        "Reionizer",
        "Uncharger",
        "Tautomer"
    ]

    def __init__(self):
        super().__init__()
        self.data = None
        self.standardized_data = None
        self.selected_operations = set()

        # Add GUI elements
        self.info_label = gui.label(self.controlArea, self, "Awaiting SMILES data...")

        # Add checkboxes for operations
        self.group_box = gui.widgetBox(self.controlArea, "Select Operations")
        self.checkboxes = {}
        for operation in self.operations:
            checkbox = QCheckBox(operation)
            checkbox.setChecked(False)  # Default: unchecked
            checkbox.stateChanged.connect(self.update_selected_operations)
            self.group_box.layout().addWidget(checkbox)
            self.checkboxes[operation] = checkbox

        # Add a process button
        self.process_button = QPushButton("Apply Standardization")
        self.process_button.clicked.connect(self.standardize_smiles)
        self.controlArea.layout().addWidget(self.process_button)

    def update_selected_operations(self):
        """Update selected operations based on checkbox states."""
        self.selected_operations = {name for name, cb in self.checkboxes.items() if cb.isChecked()}
        self.info_label.setText(f"Selected operations: {', '.join(self.selected_operations) or 'None'}")

    def set_data(self, data):
        """Handle the input data."""
        self.data = data
        if self.data is not None:
            self.info_label.setText("Input data received. Ready for standardization.")
        else:
            self.info_label.setText("No data received.")

    def standardize_smiles(self):
        """Standardize SMILES strings from the input data."""
        smiles_column_name = "SMILES"

        if self.data is None:
            self.info_label.setText("No input data to standardize.")
            return

        # Find the index of the SMILES column in metas
        try:
            smiles_index = self.data.domain.metas.index(next(meta for meta in self.data.domain.metas if meta.name == smiles_column_name))
        except StopIteration:
            self.info_label.setText(f"Column '{smiles_column_name}' not found in input data metas. Please ensure the correct column name.")
            return

        smiles_column = self.data.metas[:, smiles_index]
        standardized_smiles = []
        change_logs = []

        for smile in smiles_column:
            try:
                mol = Chem.MolFromSmiles(smile)
                original_smile = smile
                log_entries = []
                if mol:
                    for operation in self.selected_operations:
                        initial_smile = Chem.MolToSmiles(mol) if mol else ""
                        if operation == "ValidateSmiles":
                            rdMolStandardize.ValidateSmiles(smile)
                            log_entries.append("Validation passed.")
                        elif operation == "Cleanup":
                            mol = rdMolStandardize.Cleanup(mol)
                            log_entries.append("Cleanup applied.")
                        elif operation == "Normalize":
                            normalizer = rdMolStandardize.Normalizer()
                            mol = normalizer.normalize(mol)
                            log_entries.append("Normalization applied.")
                        elif operation == "MetalDisconnector":
                            disconnector = rdMolStandardize.MetalDisconnector()
                            mol = disconnector.Disconnect(mol)
                            log_entries.append("Metal disconnection applied.")
                        elif operation == "LargestFragmentChooser":
                            chooser = rdMolStandardize.LargestFragmentChooser()
                            mol = chooser.choose(mol)
                            log_entries.append("Largest fragment chosen.")
                        elif operation == "Reionizer":
                            reionizer = rdMolStandardize.Reionizer()
                            mol = reionizer.reionize(mol)
                            log_entries.append("Reionization applied.")
                        elif operation == "Uncharger":
                            uncharger = rdMolStandardize.Uncharger()
                            mol = uncharger.uncharge(mol)
                            log_entries.append("Uncharging applied.")
                        elif operation == "Tautomer":
                            enumerator = rdMolStandardize.TautomerEnumerator()
                            mol = enumerator.Canonicalize(mol)
                            log_entries.append("Tautomer standardization applied.")

                        final_smile = Chem.MolToSmiles(mol) if mol else ""
                        if initial_smile != final_smile:
                            log_entries[-1] += f" (Modified: {initial_smile} -> {final_smile})"

                    standardized_smiles.append(Chem.MolToSmiles(mol))
                else:
                    standardized_smiles.append(None)
                    log_entries.append("Invalid SMILES.")
            except Exception as e:
                self.error(f"Error standardizing SMILES '{smile}': {e}")
                standardized_smiles.append(None)
                log_entries.append(f"Error: {e}")

            change_logs.append(" | ".join(log_entries))

        self.create_output_table(smiles_column, standardized_smiles, change_logs)

    def create_output_table(self, smiles_column, standardized_smiles, change_logs):
        """Create an Orange Table with the standardized SMILES."""
        domain = Domain([], metas=[
            StringVariable("Original SMILES"),
            StringVariable("SMILES"),
            StringVariable("Change Log")
        ])
        metas = [[original, standardized, log] for original, standardized, log in zip(smiles_column, standardized_smiles, change_logs)]

        # Ensure X is an empty array if no attributes exist
        X = np.empty((len(metas), 0))
        self.standardized_data = Table.from_numpy(domain, X, metas=np.array(metas, dtype=object))

        self.info_label.setText("Standardization complete. Output table ready.")
        self.send("standardized_table", self.standardized_data)

if __name__ == "__main__":
    from AnyQt.QtWidgets import QApplication
    import sys

    app = QApplication(sys.argv)
    widget = StandardizeMoleculesWidget()
    widget.show()
    sys.exit(app.exec_())

