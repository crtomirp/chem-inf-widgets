from PyQt5.QtWidgets import QListWidget, QAbstractItemView, QFileDialog, QComboBox, QCheckBox
from Orange.widgets.widget import OWWidget, Input, Output
from Orange.widgets.settings import Setting
from Orange.widgets.gui import widgetBox, button, lineEdit, label
from PyQt5.QtCore import Qt
import os
from rdkit import Chem
from rdkit.Chem import AllChem
from Orange.data import Table, Domain, ContinuousVariable, StringVariable
from rdkit import RDLogger

# Disable RDKit warnings
RDLogger.DisableLog('rdApp.*')


class OWSDFWriter(OWWidget):
    """
    An Orange widget that writes an SDF file from an input data table.
    The table must contain a SMILES column and optionally additional data.
    Coordinate generation options are provided: 2D, 3D, or 3D with minimisation.
    Optionally, hydrogens can be added.
    """
    name = "SDF Writer"
    description = "Creates an SDF file from a table with SMILES structures."
    icon = "icons/sdfwriter.png"  # Make sure this icon exists
    priority = 3

    class Inputs:
        data = Input("Data Table", Table)

    # Widget settings to store state between sessions
    sdf_file_path = Setting("")
    smiles_column = Setting("SMILES")
    coordinate_option = Setting("2D")
    add_hydrogens = Setting(False)  # New option to add hydrogens

    def __init__(self):
        super().__init__()
        self.data = None
        # Map column name to a tuple: (source, index) where source is "attribute" or "meta"
        self._column_mapping = {}

        # --- Output File Selection ---
        file_box = widgetBox(self.controlArea, "Output File")
        self.file_line = lineEdit(
            file_box, self, "sdf_file_path", label="SDF File Path:", orientation="horizontal"
        )
        button(file_box, self, "Browse...", callback=self.browse_file)

        # --- SMILES Column Selection ---
        smiles_box = widgetBox(self.controlArea, "SMILES Column")
        self.smiles_combo = QComboBox()
        smiles_box.layout().addWidget(self.smiles_combo)
        self.smiles_combo.currentTextChanged.connect(self.on_smiles_column_changed)

        # --- Coordinate Generation Options ---
        coord_box = widgetBox(self.controlArea, "Coordinate Generation")
        self.coord_option_combo = QComboBox()
        self.coord_option_combo.addItems(["2D", "3D", "3D with minimisation"])
        coord_box.layout().addWidget(self.coord_option_combo)

        # --- Molecule Options ---
        options_box = widgetBox(self.controlArea, "Molecule Options")
        self.add_hydrogens_cb = QCheckBox("Add Hydrogens")
        self.add_hydrogens_cb.setChecked(self.add_hydrogens)
        self.add_hydrogens_cb.stateChanged.connect(self.on_add_hydrogens_changed)
        options_box.layout().addWidget(self.add_hydrogens_cb)

        # --- Info and Action ---
        self.info_label = label(self.controlArea, self, "Waiting for data...")
        button(self.controlArea, self, "Save SDF", callback=self.write_sdf)

    @Inputs.data
    def set_data(self, data):
        """Receive the input data table and update the SMILES column options."""
        self.data = data
        if data is None:
            self.info_label.setText("No input data provided.")
            return

        # Clear previous mappings and combo box
        self.smiles_combo.clear()
        self._column_mapping = {}
        domain = data.domain

        # Process attributes: store with source "attribute"
        for i, var in enumerate(domain.attributes):
            if isinstance(var, StringVariable):
                self._column_mapping[var.name] = ("attribute", i)
                self.smiles_combo.addItem(var.name)
        # Process metas: store with source "meta"
        for j, var in enumerate(domain.metas):
            if isinstance(var, StringVariable):
                self._column_mapping[var.name] = ("meta", j)
                self.smiles_combo.addItem(var.name)

        if self.smiles_combo.count() == 0:
            self.info_label.setText("No string columns available for SMILES.")
        else:
            # Default to "SMILES" if present; otherwise, select the first available column.
            index = self.smiles_combo.findText("SMILES")
            if index >= 0:
                self.smiles_combo.setCurrentIndex(index)
                self.smiles_column = "SMILES"
            else:
                self.smiles_column = self.smiles_combo.itemText(0)
            self.info_label.setText("Data loaded. Select SMILES column, coordinate option, and molecule options.")

    def on_smiles_column_changed(self, text):
        """Update the SMILES column selection."""
        self.smiles_column = text

    def on_add_hydrogens_changed(self, state):
        """Update the add_hydrogens setting based on the check box."""
        self.add_hydrogens = (state == Qt.Checked)

    def browse_file(self):
        """Open a file dialog to select the output SDF file."""
        file_path, _ = QFileDialog.getSaveFileName(
            self, "Save SDF File", "", "SDF Files (*.sdf)"
        )
        if file_path:
            self.sdf_file_path = file_path

    def write_sdf(self):
        """Process the input table and write out an SDF file."""
        if self.data is None:
            self.info_label.setText("No input data to save.")
            return

        if not self.sdf_file_path:
            self.info_label.setText("Please specify a valid SDF file path.")
            return

        # Get the coordinate generation option
        coord_option = self.coord_option_combo.currentText()
        print(f"Selected coordinate generation option: {coord_option}")

        # Ensure the selected SMILES column exists in the mapping
        if self.smiles_column not in self._column_mapping:
            self.info_label.setText("Selected SMILES column not found in data.")
            return

        source, col_index = self._column_mapping[self.smiles_column]
        print(f"Using SMILES column '{self.smiles_column}' from {source} index {col_index}")

        writer = Chem.SDWriter(self.sdf_file_path)
        count = 0
        domain = self.data.domain

        for i, row in enumerate(self.data):
            try:
                # Access the SMILES value from the correct part of the row.
                if source == "attribute":
                    smiles = row[col_index]
                else:  # source == "meta"
                    smiles = row.metas[col_index]
                print(f"Row {i}: SMILES = {smiles}")
            except Exception as e:
                print(f"Row {i}: Error extracting SMILES: {e}")
                continue

            if not smiles:
                print(f"Row {i}: SMILES is empty, skipping.")
                continue

            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                print(f"Row {i}: Failed to create molecule from SMILES: {smiles}")
                continue

            # Add hydrogens if the option is enabled.
            if self.add_hydrogens:
                print(f"Row {i}: Adding hydrogens.")
                mol = Chem.AddHs(mol)

            # Generate coordinates according to the chosen option
            if coord_option == "2D":
                print(f"Row {i}: Generating 2D coordinates.")
                AllChem.Compute2DCoords(mol)
            elif coord_option == "3D":
                print(f"Row {i}: Generating 3D coordinates.")
                res = AllChem.EmbedMolecule(mol)
                print(f"Row {i}: 3D embedding result: {res}")
                if res != 0:
                    print(f"Row {i}: 3D embedding failed, skipping molecule.")
                    continue
            elif coord_option == "3D with minimisation":
                print(f"Row {i}: Generating 3D coordinates with minimisation.")
                res = AllChem.EmbedMolecule(mol)
                print(f"Row {i}: 3D embedding result: {res}")
                if res == 0:
                    AllChem.UFFOptimizeMolecule(mol)
                else:
                    print(f"Row {i}: 3D embedding failed, skipping molecule.")
                    continue

            # Attach all available properties from the table to the molecule
            # Process attributes
            for j, var in enumerate(domain.attributes):
                try:
                    value = row[j]
                    if value is not None:
                        mol.SetProp(var.name, str(value))
                except Exception as e:
                    print(f"Row {i}: Error setting attribute property '{var.name}': {e}")
            # Process metas
            for j, var in enumerate(domain.metas):
                try:
                    value = row.metas[j]
                    if value is not None:
                        mol.SetProp(var.name, str(value))
                except Exception as e:
                    print(f"Row {i}: Error setting meta property '{var.name}': {e}")

            try:
                writer.write(mol)
                print(f"Row {i}: Molecule written.")
                count += 1
            except Exception as e:
                print(f"Row {i}: Error writing molecule to SDF: {e}")

        writer.close()
        self.info_label.setText(f"SDF file saved successfully with {count} molecules.")
        print(f"SDF file saved successfully with {count} molecules at {self.sdf_file_path}.")


if __name__ == "__main__":
    from Orange.canvas import main
    main()

