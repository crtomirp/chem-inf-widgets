import requests
import pandas as pd
import numpy as np
import logging
from datetime import datetime

from Orange.widgets.widget import OWWidget, Input, Output
from Orange.widgets import gui
from Orange.data import Table, Domain, StringVariable, ContinuousVariable

# PyQt imports
from PyQt5.QtWidgets import QPlainTextEdit, QSizePolicy, QApplication
from PyQt5.QtCore import Qt, pyqtSlot, pyqtSignal
from Orange.widgets.utils.concurrent import ThreadExecutor, FutureWatcher, methodinvoke

class ChEMBLBioactivityWidget(OWWidget):
    """Orange widget to fetch ChEMBL bioactivity data with drug properties."""
    
    name = "ChEMBL Bioactivity Retriever"
    description = "Fetches bioactivity data with drug design properties"
    icon = "icons/chembl.png"
    priority = 4

    # Declare a signal to send log messages from any thread.
    logMessage = pyqtSignal(str)

    class Outputs:
        output_data = Output("Bioactivity Data", Table)

    def __init__(self):
        super().__init__()
        # Hide the main area so everything is in the control panel
        self.mainArea.hide()

        self.target_id = ""
        # Set up a logger for detailed debug/info messages.
        self.logger = logging.getLogger(__name__)
        self.logger.setLevel(logging.INFO)
        if not self.logger.handlers:
            handler = logging.StreamHandler()
            formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
            handler.setFormatter(formatter)
            self.logger.addHandler(handler)

        # Connect the logMessage signal to the log_status slot
        self.logMessage.connect(self.log_status)

        # Build the entire UI in the control area
        self._build_ui()
        # Set up an executor for background processing
        self.executor = ThreadExecutor(self)
        self._future = None

    def _build_ui(self):
        """Construct the user interface in the control area."""
        control_box = gui.widgetBox(self.controlArea, orientation="vertical", spacing=6)

        input_box = gui.widgetBox(control_box, "Retrieve Bioactivity Data", orientation="vertical")
        gui.label(input_box, self, "Enter ChEMBL Target ID (e.g., CHEMBL2095150):")
        gui.lineEdit(input_box, self, "target_id", placeholderText="CHEMBLxxxxxx")
        self.fetch_button = gui.button(input_box, self, "Fetch Data", callback=self.fetch_bioactivity_data)

        self.status_label = gui.label(input_box, self, "Status: Awaiting input.")
        control_box.layout().addStretch(1)

        log_box = gui.widgetBox(control_box, "Status Log", orientation="vertical")
        self.log_text = QPlainTextEdit()
        self.log_text.setReadOnly(True)
        self.log_text.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        log_box.layout().addWidget(self.log_text)

    @pyqtSlot(str)
    def log_status(self, message: str):
        """Update the status label and log widget with a timestamped message."""
        timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S,%f")[:-3]
        formatted_message = f"{timestamp} - INFO - {message}"
        self.status_label.setText(message)
        self.logger.info(message)
        current_text = self.log_text.toPlainText()
        new_text = f"{current_text}\n{formatted_message}" if current_text else formatted_message
        self.log_text.setPlainText(new_text)
        self.log_text.verticalScrollBar().setValue(self.log_text.verticalScrollBar().maximum())
        QApplication.processEvents()

    @pyqtSlot(bool)
    def set_fetch_button_enabled(self, enabled: bool):
        """Helper to enable or disable the fetch button (called from main thread)."""
        self.fetch_button.setEnabled(enabled)

    def fetch_bioactivity_data(self):
        """Start background data retrieval and processing."""
        target_id = self.target_id.strip()
        if not target_id:
            self.logMessage.emit("Status: Please enter a ChEMBL Target ID.")
            return

        self.fetch_button.setEnabled(False)
        self.logMessage.emit(f"Status: Connecting to ChEMBL API for target {target_id}...")

        self._future = self.executor.submit(self._fetch_data_in_background, target_id)
        self._future.add_done_callback(self._on_fetch_complete)

    def _fetch_data_in_background(self, target_id: str) -> Table:
        """Background function that fetches and processes bioactivity data."""
        df = self._fetch_chembl_data(target_id)
        if df is None or df.empty:
            self.logMessage.emit(f"Status: No data found for {target_id}.")
            return None

        self.logMessage.emit("Status: Data fetched successfully. Processing IC50 values...")
        df = self._process_ic50_values(df)
        df = df.rename(columns={'canonical_smiles': 'SMILES'})

        if 'SMILES' in df.columns:
            self.logMessage.emit("Status: Calculating drug properties...")
            df = self._calculate_drug_properties(df)

        df = self._filter_columns(df)
        table = self._create_orange_table(df)
        return table

    def _on_fetch_complete(self, future):
        """Callback when background processing completes; update output on main thread."""
        try:
            table = future.result()
        except Exception as e:
            self._handle_error(f"Error during fetching: {e}")
            table = None
        # Use methodinvoke to schedule updating the output on the main thread.
        methodinvoke(self, "_update_output_from_table", (Table,))(table)
        methodinvoke(self, "set_fetch_button_enabled", (bool,))(True)

    @pyqtSlot(Table)
    def _update_output_from_table(self, table: Table):
        """Slot to update output once fetching is complete."""
        if table is not None:
            self.logMessage.emit(f"Status: Retrieved {len(table)} records for {self.target_id}.")
        self.Outputs.output_data.send(table)

    def _fetch_chembl_data(self, target_id: str) -> pd.DataFrame:
        """Fetch bioactivity data from ChEMBL API."""
        try:
            self.logMessage.emit("Status: Sending request to ChEMBL API...")
            response = requests.get(
                "https://www.ebi.ac.uk/chembl/api/data/activity.json",
                params={
                    "target_chembl_id": target_id,
                    "standard_type": "IC50",
                    "limit": 1000
                }
            )
            response.raise_for_status()
            self.logMessage.emit("Status: Request successful. Processing response...")
            data = response.json().get("activities", [])
            if not data:
                self.logMessage.emit(f"Status: No data returned for {target_id}.")
                self.Outputs.output_data.send(None)
                return None
            df = pd.DataFrame(data)
            if "pchembl_value" in df.columns:
                df["pchembl_value"] = pd.to_numeric(df["pchembl_value"], errors="coerce")
                df = df.dropna(subset=["pchembl_value"])
            return df

        except requests.exceptions.RequestException as e:
            self._handle_error(f"Network error during API connection: {str(e)}")
            return None

    def _process_ic50_values(self, df: pd.DataFrame) -> pd.DataFrame:
        """Convert IC50 values to nM units."""
        if 'standard_value' in df.columns and 'standard_units' in df.columns:
            df['IC50_nM'] = df.apply(lambda row: self._convert_to_nM(row), axis=1)
            df = df.drop(columns=['standard_value', 'standard_units'])
        return df

    def _convert_to_nM(self, row: pd.Series) -> float:
        """Convert IC50 value to nanomolar units."""
        try:
            value = float(row['standard_value'])
            unit = row['standard_units'].lower()
            conversions = {
                'm': value * 1e9,
                'µm': value * 1e3,
                'um': value * 1e3,
                'nm': value,
                'nmol/l': value,
                'pm': value * 1e-3
            }
            return conversions.get(unit, np.nan)
        except Exception as e:
            self.logger.error(f"Conversion error: {e}")
            return np.nan

    def _calculate_drug_properties(self, df: pd.DataFrame) -> pd.DataFrame:
        """Calculate molecular properties using RDKit."""
        try:
            from rdkit import Chem
            from rdkit.Chem import Descriptors, Lipinski
        except ImportError:
            self._handle_error("RDKit not installed - skipping property calculations.")
            return df

        prop_columns = ['hbd', 'hba', 'rotable_bonds', 'mw', 'tpsa', 'logp', 'lipinski_deviations']
        for col in prop_columns:
            df[col] = np.nan

        for idx, row in df.iterrows():
            smiles = row.get('SMILES', '')
            if not smiles:
                continue
            try:
                mol = Chem.MolFromSmiles(smiles)
                if not mol:
                    continue
                df.at[idx, 'hbd'] = Lipinski.NumHDonors(mol)
                df.at[idx, 'hba'] = Lipinski.NumHAcceptors(mol)
                df.at[idx, 'rotable_bonds'] = Descriptors.NumRotatableBonds(mol)
                df.at[idx, 'mw'] = Descriptors.MolWt(mol)
                df.at[idx, 'tpsa'] = Descriptors.TPSA(mol)
                df.at[idx, 'logp'] = Descriptors.MolLogP(mol)
                violations = 0
                violations += 1 if df.at[idx, 'mw'] > 500 else 0
                violations += 1 if df.at[idx, 'logp'] > 5 else 0
                violations += 1 if df.at[idx, 'hbd'] > 5 else 0
                violations += 1 if df.at[idx, 'hba'] > 10 else 0
                df.at[idx, 'lipinski_deviations'] = violations
            except Exception as e:
                self.logger.error(f"Error calculating properties for SMILES '{smiles}': {e}")
                continue

        for col in prop_columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")
        return df

    def _filter_columns(self, df: pd.DataFrame) -> pd.DataFrame:
        """Filter columns to include key numerical and meta columns."""
        meta_cols = [
            'SMILES', 'molecule_chembl_id', 'target_chembl_id',
            'assay_chembl_id', 'document_chembl_id', 'target_organism',
            'target_name'
        ]
        num_cols = [
            'pchembl_value', 'IC50_nM', 'hbd', 'hba', 'rotable_bonds',
            'mw', 'tpsa', 'logp', 'lipinski_deviations'
        ]
        existing_meta = [col for col in meta_cols if col in df.columns]
        existing_num = [col for col in num_cols if col in df.columns]
        return df[existing_num + existing_meta]

    def _create_orange_table(self, df: pd.DataFrame) -> Table:
        """Convert a pandas DataFrame to an Orange Table."""
        num_cols = [col for col in df.columns if col in [
            'pchembl_value', 'IC50_nM', 'hbd', 'hba', 'rotable_bonds',
            'mw', 'tpsa', 'logp', 'lipinski_deviations'
        ]]
        meta_cols = [col for col in df.columns if col not in num_cols]
        domain = Domain(
            [ContinuousVariable(col) for col in num_cols],
            metas=[StringVariable(col) for col in meta_cols]
        )
        X = df[num_cols].to_numpy(dtype=float)
        metas = df[meta_cols].to_numpy(dtype=object)
        return Table.from_numpy(domain, X=X, metas=metas)

    def _handle_error(self, message: str):
        """Centralized error handling."""
        self.logMessage.emit(message)
        self.Outputs.output_data.send(None)

if __name__ == "__main__":
    from Orange.widgets.utils.widgetpreview import WidgetPreview
    WidgetPreview(ChEMBLBioactivityWidget).run()
