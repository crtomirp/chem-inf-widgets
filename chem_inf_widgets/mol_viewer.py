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
    description = "Displays molecules with optional substructure highlighting and property selection."
    icon = "icons/molviewer.png"
    priority = 10

    class Inputs:
        orange_data = Input("Filtered Compounds", Table)

    def __init__(self):
        super().__init__()

        self.orange_data = None
        self.max_columns = 5
        self.molecule_size = 300
        self.use_rdCoordGen = False
        self.highlight_enabled = True  # Default: Highlighting ON
        self.selected_properties = []  # Properties selected by the user

        # UI Components
        self.scroll_area = QScrollArea()
        self.scroll_area.setWidgetResizable(True)
        self.scroll_widget = QWidget()
        self.grid_layout = QGridLayout()
        self.scroll_widget.setLayout(self.grid_layout)
        self.scroll_area.setWidget(self.scroll_widget)
        self.mainArea.layout().addWidget(self.scroll_area)

        # Settings UI
        gui.label(self.controlArea, self, "Grid Settings")

        self.size_selector = QSpinBox()
        self.size_selector.setMinimum(100)
        self.size_selector.setMaximum(500)
        self.size_selector.setValue(self.molecule_size)
        self.size_selector.valueChanged.connect(self.update_molecule_size)
        gui.widgetBox(self.controlArea, orientation="horizontal").layout().addWidget(self.size_selector)

        self.column_selector = QSpinBox()
        self.column_selector.setMinimum(1)
        self.column_selector.setMaximum(10)
        self.column_selector.setValue(self.max_columns)
        self.column_selector.valueChanged.connect(self.update_columns)
        gui.widgetBox(self.controlArea, orientation="horizontal").layout().addWidget(self.column_selector)

        self.rdCoordGen_checkbox = QCheckBox("Use rdCoordGen for Better Structures")
        self.rdCoordGen_checkbox.stateChanged.connect(self.toggle_rdCoordGen)
        gui.widgetBox(self.controlArea, orientation="horizontal").layout().addWidget(self.rdCoordGen_checkbox)

        # New: Highlight Toggle Checkbox
        self.highlight_checkbox = QCheckBox("Highlight Substructures")
        self.highlight_checkbox.setChecked(True)  # Default: ON
        self.highlight_checkbox.stateChanged.connect(self.toggle_highlighting)
        gui.widgetBox(self.controlArea, orientation="horizontal").layout().addWidget(self.highlight_checkbox)

        # Property Selection UI
        gui.label(self.controlArea, self, "Select Properties to Display")
        self.property_selector = QListWidget()
        self.property_selector.setSelectionMode(QListWidget.MultiSelection)
        gui.widgetBox(self.controlArea, orientation="vertical").layout().addWidget(self.property_selector)

        # Apply Property Selection Button
        apply_button = QPushButton("Apply Selection")
        apply_button.clicked.connect(self.update_selected_properties)
        gui.widgetBox(self.controlArea, orientation="horizontal").layout().addWidget(apply_button)

        self.info_label = gui.label(self.controlArea, self, "Awaiting input data...")

    @Inputs.orange_data
    def set_orange_data(self, data: Table):
        """Receives the processed molecules with highlighted atoms."""
        self.orange_data = data
        if data is not None:
            self.info_label.setText(f"Received {len(data)} molecules.")
            self.update_property_selector()
            self.display_molecules()
        else:
            self.info_label.setText("No data received.")

    def update_columns(self, value):
        """Update number of columns in the grid."""
        self.max_columns = value
        self.display_molecules()

    def update_molecule_size(self, value):
        """Update molecule image size."""
        self.molecule_size = value
        self.display_molecules()

    def toggle_rdCoordGen(self, state):
        """Toggle rdCoordGen for 2D coordinates."""
        self.use_rdCoordGen = bool(state)
        self.display_molecules()

    def toggle_highlighting(self, state):
        """Enable or disable substructure highlighting."""
        self.highlight_enabled = bool(state)
        self.display_molecules()

    def update_property_selector(self):
        """Populate the property selector with available columns."""
        self.property_selector.clear()
        if self.orange_data is not None:
            for var in self.orange_data.domain.variables + self.orange_data.domain.metas:
                if var.name.lower() != "smiles" and var.name != "Highlighted Atoms":
                    self.property_selector.addItem(var.name)

    def update_selected_properties(self):
        """Update the selected properties list."""
        self.selected_properties = [
            item.text() for item in self.property_selector.selectedItems()
        ]
        self.display_molecules()

    def display_molecules(self):
        """Generate and display molecule images with optional highlighting and property selection."""
        if self.orange_data is None:
            self.info_label.setText("No input data available.")
            return

        domain = self.orange_data.domain
        smiles_col = next((var for var in domain.variables + domain.metas if var.name.lower() == "smiles"), None)
        highlight_col = next((var for var in domain.variables + domain.metas if var.name == "Highlighted Atoms"), None)

        if smiles_col is None:
            self.info_label.setText("No 'SMILES' column found in data.")
            return

        # Clear the existing grid
        for i in reversed(range(self.grid_layout.count())):
            widget_to_remove = self.grid_layout.itemAt(i).widget()
            self.grid_layout.removeWidget(widget_to_remove)
            widget_to_remove.deleteLater()

        row, col = 0, 0

        for instance in self.orange_data:
            smiles = instance[smiles_col].value
            highlighted_atoms = []
            if self.highlight_enabled and highlight_col:
                highlight_data = instance[highlight_col].value
                if highlight_data:
                    highlighted_atoms = list(map(int, highlight_data.split(",")))  # Convert to list of atom indices

            try:
                mol = Chem.MolFromSmiles(smiles)
                if mol:
                    if self.use_rdCoordGen:
                        AddCoords(mol)
                    else:
                        AllChem.Compute2DCoords(mol)

                    img = Draw.MolToImage(mol, size=(self.molecule_size, self.molecule_size),
                                          highlightAtoms=highlighted_atoms if self.highlight_enabled else [])

                    pixmap = self._convert_pil_to_pixmap(img)
                    image_label = QLabel()
                    image_label.setPixmap(pixmap)
                    image_label.setAlignment(Qt.AlignCenter)
                    image_label.setStyleSheet("border: 1px solid #ddd; padding: 5px; margin: 5px;")

                    # Property labels
                    property_labels = []
                    for prop in self.selected_properties:
                        prop_value = instance[domain[prop]].value if prop in domain else ""
                        text_label = QLabel(f"{prop}: {prop_value}")
                        text_label.setAlignment(Qt.AlignCenter)
                        text_label.setStyleSheet("font-size: 12px; color: #555;")
                        property_labels.append(text_label)

                    # Layout for molecule + properties
                    molecule_widget = QVBoxLayout()
                    molecule_widget.addWidget(image_label)
                    for label in property_labels:
                        molecule_widget.addWidget(label)

                    container = QWidget()
                    container.setLayout(molecule_widget)
                    container.setStyleSheet("background-color: #f9f9f9; border-radius: 5px;")

                    self.grid_layout.addWidget(container, row, col)
                    col += 1
                    if col == self.max_columns:
                        col = 0
                        row += 1

            except Exception as e:
                print(f"Error processing SMILES {smiles}: {e}")

        self.scroll_area.setWidget(self.scroll_widget)

    def _convert_pil_to_pixmap(self, pil_image):
        buffer = io.BytesIO()
        pil_image.save(buffer, format="PNG")
        buffer.seek(0)
        qt_image = QImage.fromData(buffer.getvalue())
        return QPixmap.fromImage(qt_image)

if __name__ == "__main__":
    from Orange.widgets.utils.widgetpreview import WidgetPreview
    WidgetPreview(MolViewer).run()
