from Orange.widgets.widget import OWWidget, Input
from Orange.widgets import gui
from Orange.data import Table
from rdkit import Chem
from rdkit.Chem import Draw, AllChem
from rdkit.Chem.rdCoordGen import AddCoords
import pandas as pd
from PyQt5.QtWidgets import QLabel, QVBoxLayout, QGridLayout, QWidget, QScrollArea, QSpinBox, QListWidget, QPushButton, QCheckBox
from PyQt5.QtGui import QPixmap, QImage
from PyQt5.QtCore import Qt
import io

class MolViewer(OWWidget):
    name = "Molecular Viewer"
    description = "Converts SMILES codes to molecule images and displays them in a grid."
    icon = "icons/molviewer.png"
    priority = 20

    class Inputs:
        orange_data = Input("Orange Data", Table)
        smiles_data = Input("SMILES Data", pd.DataFrame)

    def __init__(self):
        super().__init__()

        self.orange_data = None
        self.smiles_data = None

        # Default settings
        self.max_columns = 5
        self.selected_properties = []
        self.molecule_size = 300
        self.use_rdCoordGen = False

        # Create a scrollable area for the grid of molecule images
        self.scroll_area = QScrollArea()
        self.scroll_area.setWidgetResizable(True)
        self.scroll_widget = QWidget()
        self.grid_layout = QGridLayout()
        self.scroll_widget.setLayout(self.grid_layout)
        self.scroll_area.setWidget(self.scroll_widget)

        # Add scroll area to the main area
        self.mainArea.layout().addWidget(self.scroll_area)

        # Control area for settings
        gui.label(self.controlArea, self, "Grid Settings")

        # Molecule size selector
        size_label = gui.label(self.controlArea, self, "Molecule Image Size")
        self.size_selector = QSpinBox()
        self.size_selector.setMinimum(100)
        self.size_selector.setMaximum(500)
        self.size_selector.setValue(self.molecule_size)
        self.size_selector.valueChanged.connect(self.update_molecule_size)
        gui.widgetBox(self.controlArea, orientation="horizontal").layout().addWidget(self.size_selector)

        # Column selector
        column_label = gui.label(self.controlArea, self, "Grid Columns")
        self.column_selector = QSpinBox()
        self.column_selector.setMinimum(1)
        self.column_selector.setMaximum(10)
        self.column_selector.setValue(self.max_columns)
        self.column_selector.valueChanged.connect(self.update_columns)
        gui.widgetBox(self.controlArea, orientation="horizontal").layout().addWidget(self.column_selector)

        # Use rdCoordGen checkbox
        self.rdCoordGen_checkbox = QCheckBox("Use rdCoordGen for Better Structures")
        self.rdCoordGen_checkbox.stateChanged.connect(self.toggle_rdCoordGen)
        gui.widgetBox(self.controlArea, orientation="horizontal").layout().addWidget(self.rdCoordGen_checkbox)

        # Property selection
        gui.label(self.controlArea, self, "Select Properties to Display")
        self.property_selector = QListWidget()
        self.property_selector.setSelectionMode(QListWidget.MultiSelection)
        gui.widgetBox(self.controlArea, orientation="vertical").layout().addWidget(self.property_selector)

        # Apply button
        apply_button = QPushButton("Apply Selection")
        apply_button.clicked.connect(self.update_selected_properties)
        gui.widgetBox(self.controlArea, orientation="horizontal").layout().addWidget(apply_button)

        # Info label
        self.info_label = gui.label(self.controlArea, self, "Awaiting SMILES data...")

    @Inputs.orange_data
    def set_orange_data(self, data: Table):
        """Set the input data as Orange Table."""
        self.orange_data = data
        if data is not None:
            self.process_orange_data()

    @Inputs.smiles_data
    def set_smiles_data(self, data: pd.DataFrame):
        """Handles the input SMILES data as a pandas DataFrame."""
        self.smiles_data = data
        if data is not None:
            self.info_label.setText(f"Received {len(data)} SMILES entries.")
            self.update_property_selector()
            self.display_molecules()
        else:
            self.info_label.setText("No data received.")
            self.smiles_data = None

    def process_orange_data(self):
        """Process the Orange Table to extract SMILES data."""
        data = {attr.name: [str(instance[attr]) for instance in self.orange_data]
                for attr in self.orange_data.domain.attributes + self.orange_data.domain.metas}
        self.smiles_data = pd.DataFrame(data)
        self.info_label.setText(f"Extracted {len(self.smiles_data)} entries.")
        self.update_property_selector()
        self.display_molecules()

    def update_columns(self, value):
        """Update the number of columns in the grid."""
        self.max_columns = value
        self.display_molecules()

    def update_molecule_size(self, value):
        """Update the molecule image size."""
        self.molecule_size = value
        self.display_molecules()

    def toggle_rdCoordGen(self, state):
        """Toggle the use of rdCoordGen for generating 2D coordinates."""
        self.use_rdCoordGen = bool(state)
        self.display_molecules()

    def update_property_selector(self):
        """Populate the property selector with available columns."""
        self.property_selector.clear()
        if self.smiles_data is not None:
            for column in self.smiles_data.columns:
                self.property_selector.addItem(column)

    def update_selected_properties(self):
        """Update the selected properties list."""
        self.selected_properties = [
            item.text() for item in self.property_selector.selectedItems()
        ]
        self.display_molecules()

    def display_molecules(self):
        """Generates and displays images of molecules in a grid."""
        if self.smiles_data is None or "SMILES" not in self.smiles_data.columns:
            self.info_label.setText("No 'SMILES' column found in the data.")
            return

        smiles_column = self.smiles_data["SMILES"].dropna()

        # Clear the existing grid layout
        for i in reversed(range(self.grid_layout.count())):
            widget_to_remove = self.grid_layout.itemAt(i).widget()
            self.grid_layout.removeWidget(widget_to_remove)
            widget_to_remove.deleteLater()

        row, col = 0, 0

        for idx, smile in enumerate(smiles_column):
            try:
                mol = Chem.MolFromSmiles(smile)
                if mol:
                    # Generate 2D coordinates
                    if self.use_rdCoordGen:
                        AddCoords(mol)
                    else:
                        AllChem.Compute2DCoords(mol)

                    # Generate the image
                    img = Draw.MolToImage(mol, size=(self.molecule_size, self.molecule_size))
                    pixmap = self._convert_pil_to_pixmap(img)

                    # Create a QLabel for the image
                    image_label = QLabel()
                    image_label.setPixmap(pixmap)
                    image_label.setAlignment(Qt.AlignCenter)
                    image_label.setStyleSheet("border: 1px solid #ddd; padding: 5px; margin: 5px;")

                    # Create QLabels for the selected properties
                    property_labels = []
                    for prop in self.selected_properties:
                        prop_value = self.smiles_data.loc[idx, prop] if prop in self.smiles_data.columns else ""
                        text_label = QLabel(f"{prop}: {prop_value}")
                        text_label.setAlignment(Qt.AlignCenter)
                        text_label.setStyleSheet("font-size: 12px; color: #555;")
                        property_labels.append(text_label)

                    # Add all to a vertical layout
                    molecule_widget = QVBoxLayout()
                    molecule_widget.addWidget(image_label)
                    for label in property_labels:
                        molecule_widget.addWidget(label)

                    # Wrap in a QWidget for the grid layout
                    container = QWidget()
                    container.setLayout(molecule_widget)
                    container.setStyleSheet("background-color: #f9f9f9; border-radius: 5px;")

                    # Add to grid
                    self.grid_layout.addWidget(container, row, col)

                    col += 1
                    if col == self.max_columns:
                        col = 0
                        row += 1

            except Exception as e:
                print(f"Error processing SMILES {smile}: {e}")

        self.scroll_area.setWidget(self.scroll_widget)

    def _convert_pil_to_pixmap(self, pil_image):
        """Converts a PIL image to a QPixmap for display in QLabel."""
        buffer = io.BytesIO()
        pil_image.save(buffer, format="PNG")
        buffer.seek(0)
        qt_image = QImage.fromData(buffer.getvalue())
        return QPixmap.fromImage(qt_image)

if __name__ == "__main__":
    from Orange.widgets.utils.widgetpreview import WidgetPreview

    WidgetPreview(MolViewer).run()

