from Orange.widgets.widget import OWWidget, Input, Output
from Orange.widgets import gui
from rdkit import Chem
from rdkit.Chem import MACCSkeys
import pandas as pd
from Orange.data import Table, Domain, StringVariable, ContinuousVariable

# Import Qt threading classes
from PyQt5.QtCore import QRunnable, QThreadPool, pyqtSignal, QObject, pyqtSlot

class WorkerSignals(QObject):
    """
    Defines the signals available from a running worker thread.
    
    Attributes:
        finished (pyqtSignal): Emits the resulting Orange Table and an informational message.
        error (pyqtSignal): Emits an error message if processing fails.
    """
    finished = pyqtSignal(object, str)  # Emitting (Table, info_message)
    error = pyqtSignal(str)

class SMILESWorker(QRunnable):
    """
    Worker thread for processing SMILES strings and generating MACCS keys.
    
    Args:
        smiles_data (DataFrame): A pandas DataFrame with a 'SMILES' column.
    """
    def __init__(self, smiles_data):
        super().__init__()
        self.smiles_data = smiles_data
        self.signals = WorkerSignals()

    @pyqtSlot()
    def run(self):
        """Processes the SMILES data to compute MACCS keys and emits the results."""
        # Ensure the SMILES column is valid
        if "SMILES" not in self.smiles_data.columns:
            self.signals.error.emit("No 'SMILES' column found in the data.")
            return

        smiles_column = self.smiles_data["SMILES"].dropna()

        maccs_keys = []
        valid_smiles = []

        for smile in smiles_column:
            try:
                mol = Chem.MolFromSmiles(smile)
                if mol:
                    maccs = MACCSkeys.GenMACCSKeys(mol)
                    maccs_keys.append(list(maccs))
                    valid_smiles.append(smile)
                else:
                    print(f"Invalid SMILES: {smile}")
            except Exception as e:
                print(f"Error processing SMILES {smile}: {e}")

        if not maccs_keys:
            self.signals.error.emit("No valid SMILES to process.")
            return

        # Create the Orange Table with MACCS keys
        domain = Domain(
            [ContinuousVariable(f"MACCS_{i+1}") for i in range(len(maccs_keys[0]))],
            metas=[StringVariable("SMILES")]
        )
        table = Table.from_list(
            domain, [list(keys) + [smile] for keys, smile in zip(maccs_keys, valid_smiles)]
        )
        info_message = f"Processed {len(maccs_keys)} valid SMILES entries."
        self.signals.finished.emit(table, info_message)

class SmilesToMACCSWidget(OWWidget):
    """
    A widget that converts SMILES codes to MACCS keys and outputs an Orange Table.

    Attributes:
        name (str): The name of the widget.
        description (str): A short description of the widget functionality.
        icon (str): The path to the icon file.
        priority (int): The priority level of the widget.
    """
    name = "MACCS Key Generator"
    description = "Converts SMILES codes to MACCS keys and saves them to a table."
    icon = "icons/maccs.png"  # Add an appropriate icon or remove this line
    priority = 21

    class Inputs:
        smiles_data = Input("SMILES Data", Table)

    class Outputs:
        maccs_table = Output("MACCS Keys Table", Table)

    def __init__(self):
        """Initializes the widget with default settings, GUI components, and a thread pool."""
        super().__init__()
        self.data = None
        self.threadpool = QThreadPool()

        # Add GUI elements
        self.info_label = gui.label(self.controlArea, self, "Awaiting SMILES data...")

    @Inputs.smiles_data
    def set_data(self, data):
        """
        Handles the input SMILES data.

        Args:
            data (Table): The input Orange Table containing SMILES data.
        """
        if data is not None:
            self.data = self._convert_to_dataframe(data)
            if self.data.empty:
                self.info_label.setText("No 'SMILES' column found in the data.")
                self.Outputs.maccs_table.send(None)
            else:
                self.info_label.setText(f"Received {len(self.data)} SMILES entries.")
                self.process_smiles()
        else:
            self.info_label.setText("No data received.")
            self.data = None
            self.Outputs.maccs_table.send(None)

    def _convert_to_dataframe(self, table):
        """
        Converts an Orange Table to a pandas DataFrame containing the SMILES column.

        Args:
            table (Table): The input Orange Table.

        Returns:
            DataFrame: A pandas DataFrame with the 'SMILES' column.
        """
        # Check for SMILES column in metas (case-insensitive)
        smiles_col = [var.name for var in table.domain.metas if var.name.lower() == "smiles"]
        if not smiles_col:
            return pd.DataFrame()

        # Extract SMILES data
        smiles_col_name = smiles_col[0]
        smiles_index = table.domain.metas.index(StringVariable(smiles_col_name))
        smiles_data = table.metas[:, smiles_index].flatten()

        return pd.DataFrame({"SMILES": smiles_data})

    def process_smiles(self):
        """Initiates background processing of SMILES to generate MACCS keys."""
        self.info_label.setText("Processing SMILES in background...")
        worker = SMILESWorker(self.data)
        worker.signals.finished.connect(self.on_processing_finished)
        worker.signals.error.connect(self.on_processing_error)
        self.threadpool.start(worker)

    def on_processing_finished(self, table, info_message):
        """
        Callback when SMILES processing is successfully completed.

        Args:
            table (Table): The resulting Orange Table with MACCS keys.
            info_message (str): Informational message about the processing.
        """
        self.Outputs.maccs_table.send(table)
        self.info_label.setText(info_message)

    def on_processing_error(self, error_message):
        """
        Callback when there is an error during SMILES processing.

        Args:
            error_message (str): The error message to display.
        """
        self.info_label.setText(error_message)
        self.Outputs.maccs_table.send(None)

if __name__ == "__main__":
    from Orange.widgets.utils.widgetpreview import WidgetPreview

    WidgetPreview(SmilesToMACCSWidget).run()

