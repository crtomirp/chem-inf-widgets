import os
import json
import numpy as np
from Orange.widgets import gui, widget
from Orange.widgets.settings import Setting
from PyQt5.QtWidgets import QFileDialog, QLabel, QMessageBox, QSizePolicy
from PyQt5.QtWebEngineWidgets import QWebEngineView, QWebEnginePage
from PyQt5.QtCore import QUrl, Qt
from PyQt5.QtGui import QColor
from Orange.data import Table, Domain, StringVariable, ContinuousVariable
import traceback

# RDKit configuration
try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors, Lipinski
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False

class CustomWebEnginePage(QWebEnginePage):
    def javaScriptConsoleMessage(self, level, message, line_number, source_id):
        print(f"JS Console: {message} (Line {line_number}) in {source_id}")

class OWJSMEMolecularSketcher(widget.OWWidget):
    name = "MolSketcher"
    description = "JSME-based molecular editor for DataBase creation"
    icon = "icons/molsketcher.png"
    priority = 11
    keywords = ["chemistry", "molecule", "sketcher"]

    outputs = [("Compounds", Table)]

    want_main_area = True
    resizing_enabled = True

    # Settings
    json_path = Setting("")
    data = Setting([])

    PROPERTY_MAP = {
        'mw': ('Molecular Weight', Descriptors.MolWt),
        'logp': ('LogP', Descriptors.MolLogP),
        'tpsa': ('TPSA', Descriptors.TPSA),
        'hbd': ('H-Bond Donors', Lipinski.NumHDonors),
        'hba': ('H-Bond Acceptors', Lipinski.NumHAcceptors),
        'rotatable_bonds': ('Rotatable Bonds', Lipinski.NumRotatableBonds),
        'inchi': ('InChI', Chem.MolToInchi),
        'inchikey': ('InChI Key', Chem.MolToInchiKey)
    }

    def __init__(self):
        super().__init__()
        
        self.fields = []           # List of compound properties (from JSON config)
        self.user_metadata = []    # List of metadata fields (from JSON config)
        self.current_smiles = ""
        self.dbkey_config = None   # Optional dbkey configuration (if provided)
        self._dbkey_counter = 1    # Counter for auto-incrementing DB key (modified with "initial" value if provided)
        
        self._setup_jsme()
        self._setup_ui()
        
        if self.json_path and os.path.exists(self.json_path):
            self._load_config()

    def _setup_jsme(self):
        html_path = os.path.join(os.path.dirname(__file__), "jsme/jsme_panel.html")
        if not os.path.exists(html_path):
            raise FileNotFoundError(f"JSME HTML not found: {html_path}")
            
        self.web_view = QWebEngineView(self.mainArea)
        self.web_view.setPage(CustomWebEnginePage(self.web_view))
        self.web_view.load(QUrl.fromLocalFile(html_path))
        self.web_view.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        self.mainArea.layout().addWidget(self.web_view)

    def _setup_ui(self):
        self.setStyleSheet("""
            QWidget { font-size: 12px; }
            QPushButton { 
                background-color: #4CAF50;
                color: white; 
                border: none;
                padding: 5px 10px;
                border-radius: 4px;
                min-width: 80px;
            }
            QPushButton:hover { background-color: #45a049; }
            QPushButton:disabled { background-color: #cccccc; }
            QLabel#status { 
                color: #666; 
                font-style: italic;
                margin-top: 10px;
            }
        """)

        # Configuration section
        config_box = gui.widgetBox(self.controlArea, "Configuration", orientation=Qt.Vertical)
        gui.button(
            config_box, self, "Load Config...",
            callback=self._load_config_dialog,
            tooltip="Load JSON configuration file"
        )
        self.config_label = QLabel("No configuration loaded")
        self.config_label.setWordWrap(True)
        self.config_label.setStyleSheet("color: #666; margin-top: 5px;")
        config_box.layout().addWidget(self.config_label)
        
        # Metadata inputs
        self.metadata_box = gui.widgetBox(self.controlArea, "Sample Metadata")
        
        # Status label
        self.info_label = QLabel()
        self.info_label.setObjectName("status")
        self.controlArea.layout().addWidget(self.info_label)

        # Action buttons
        btn_box = gui.widgetBox(self.controlArea, "Actions", orientation=Qt.Vertical)
        self.add_btn = gui.button(
            btn_box, self, "Add Compound",
            callback=self._add_compound,
            tooltip="Add current molecule to the table"
        )
        self.clear_btn = gui.button(
            btn_box, self, "Clear All",
            callback=self._clear_table,
            tooltip="Remove all compounds from the table"
        )
        self.clear_btn.setEnabled(False)

        if not RDKIT_AVAILABLE:
            self.warning("RDKit not installed. Chemical properties unavailable.")

    def _load_config_dialog(self):
        path, _ = QFileDialog.getOpenFileName(
            self, "Load Configuration", "", "JSON Files (*.json)"
        )
        if path:
            self.json_path = path
            self._load_config()

    def _load_config(self):
        try:
            with open(self.json_path, 'r') as f:
                config = json.load(f)
            
            # --- Handle optional iterative dbkey ---
            if 'dbkey' in config:
                self.dbkey_config = config['dbkey']
                # Initialize the counter with the provided "initial" value, or default to 1
                self._dbkey_counter = self.dbkey_config.get("initial", 1)
                # If the dbkey field is not already in the list of fields, insert it as the first field
                fields_from_config = config.get('fields', [])
                if not any(field.get('name') == self.dbkey_config.get('name') for field in fields_from_config):
                    self.fields = [self.dbkey_config] + fields_from_config
                else:
                    self.fields = fields_from_config
            else:
                self.dbkey_config = None
                self.fields = config.get('fields', [])
            
            self.user_metadata = config.get('user_metadata', [])
            
            # Clear existing metadata inputs
            while self.metadata_box.layout().count():
                item = self.metadata_box.layout().takeAt(0)
                if widget_inst := item.widget():
                    widget_inst.deleteLater()
            
            # Create new metadata inputs
            for meta in self.user_metadata:
                field_name = meta['name']
                field_label = meta.get('label', field_name)
                setattr(self, field_name, "")
                gui.lineEdit(
                    self.metadata_box, self, field_name,
                    label=field_label + ":",
                    tooltip=meta.get('description', ''),
                    callback=self._update_metadata
                )
            
            self.config_label.setText(f"Loaded: {os.path.basename(self.json_path)}")
            self.info_label.setText("Configuration loaded successfully")
            self.info_label.setStyleSheet("color: #4CAF50;")
            
        except Exception as e:
            QMessageBox.critical(self, "Configuration Error", 
                               f"Failed to load configuration:\n{str(e)}")
            self.config_label.setText("Invalid configuration file")
            self.config_label.setStyleSheet("color: #f44336;")
            traceback.print_exc()

    def _update_metadata(self):
        self.info_label.setText("Metadata fields updated")
        self.info_label.setStyleSheet("color: #2196F3;")

    def _add_compound(self):
        self.web_view.page().runJavaScript("getSmiles();", self._process_smiles)

    def _process_smiles(self, smiles):
        try:
            if not smiles:
                self.error("Please draw a molecule first")
                return

            compound = {'smiles': smiles}
            if RDKIT_AVAILABLE:
                mol = Chem.MolFromSmiles(smiles)
                if not mol:
                    raise ValueError("Invalid SMILES structure")
                
                for field in self.fields:
                    prop_name = field['name']
                    if prop_name in self.PROPERTY_MAP:
                        label, calculator = self.PROPERTY_MAP[prop_name]
                        try:
                            compound[label] = calculator(mol)
                        except Exception as e:
                            compound[label] = f"Error: {str(e)}"
                    elif prop_name == 'smiles':
                        compound[field['label']] = smiles

            # --- Add iterative DB key if configured ---
            if self.dbkey_config is not None:
                compound[self.dbkey_config['label']] = self._dbkey_counter
                self._dbkey_counter += 1

            for meta in self.user_metadata:
                compound[meta['label']] = getattr(self, meta['name'], "")
            
            self.data.append(compound)
            self._update_output()
            self.clear_btn.setEnabled(True)
            
            self.info_label.setText(
                f"✓ Added: {smiles[:20]}...\nTotal compounds: {len(self.data)}"
            )
            self.info_label.setStyleSheet("color: #4CAF50;")

        except Exception as e:
            self.error(f"Error processing compound: {str(e)}")
            self.info_label.setText(f"✗ Error: {str(e)}")
            self.info_label.setStyleSheet("color: #f44336;")
            traceback.print_exc()

    def _update_output(self):
        if not self.data:
            self.send("Compounds", None)
            return
        
        attributes = []
        metas = []
        
        # Create variables from fields.
        # Here we check if a field is numeric by looking at its type (for example, "float" or "int").
        for field in self.fields:
            label = field['label']
            if field.get('type') in ['float', 'int']:
                var = ContinuousVariable(label)
                attributes.append(var)
            else:
                var = StringVariable(label)
                metas.append(var)
        
        for meta in self.user_metadata:
            metas.append(StringVariable(meta['label']))
        
        domain = Domain(attributes, metas=metas)
        
        X, M = [], []
        for comp in self.data:
            x_row = [comp.get(var.name, np.nan) for var in attributes]
            meta_row = [str(comp.get(var.name, '')) for var in metas]
            X.append(x_row)
            M.append(meta_row)
        
        table = Table.from_numpy(
            domain, 
            X=np.array(X, dtype=float) if X else np.empty((len(X), 0)),
            metas=np.array(M, dtype=object)
        )
        table.name = "Compounds"
        self.send("Compounds", table)

    def _clear_table(self):
        self.data.clear()
        self._update_output()
        self.clear_btn.setEnabled(False)
        self.info_label.setText("✓ All compounds cleared")
        self.info_label.setStyleSheet("color: #4CAF50;")

if __name__ == "__main__":
    from Orange.widgets.utils.widgetpreview import WidgetPreview
    WidgetPreview(OWJSMEMolecularSketcher).run()
