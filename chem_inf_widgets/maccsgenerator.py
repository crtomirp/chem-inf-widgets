from Orange.widgets.widget import OWWidget, Input, Output
from Orange.widgets import gui
from rdkit import Chem
from rdkit.Chem import MACCSkeys
import pandas as pd
from Orange.data import Table, Domain, StringVariable, ContinuousVariable

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
    priority = 20

    class Inputs:
        smiles_data = Input("SMILES Data", Table)

    class Outputs:
        maccs_table = Output("MACCS Keys Table", Table)

    def __init__(self):
        """Initializes the widget with default settings and GUI components."""
        super().__init__()

        self.data = None

        # Add GUI elements
        self.info_label = gui.label(self.controlArea, self, "Awaiting SMILES data...")

    @Inputs.smiles_data
    def set_data(self, data):
        """Handles the input SMILES data.

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
        """Converts Orange Table to pandas DataFrame.

        Args:
            table (Table): The input Orange Table.

        Returns:
            DataFrame: A pandas DataFrame containing the SMILES column.
        """
        # Check for SMILES column in metas
        smiles_col = [var.name for var in table.domain.metas if var.name.lower() == "smiles"]
        if not smiles_col:
            return pd.DataFrame()

        # Extract SMILES data
        smiles_col_name = smiles_col[0]
        smiles_index = table.domain.metas.index(StringVariable(smiles_col_name))
        smiles_data = table.metas[:, smiles_index].flatten()

        return pd.DataFrame({"SMILES": smiles_data})

    def process_smiles(self):
        """Converts SMILES to MACCS keys and outputs a table."""
        if "SMILES" not in self.data.columns:
            self.info_label.setText("No 'SMILES' column found in the data.")
            self.Outputs.maccs_table.send(None)
            return

        smiles_column = self.data["SMILES"].dropna()  # Ensure SMILES column exists and is not null

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
            self.info_label.setText("No valid SMILES to process.")
            self.Outputs.maccs_table.send(None)
            return

        # Prepare the Orange table
        domain = Domain(
            [ContinuousVariable(f"MACCS_{i+1}") for i in range(len(maccs_keys[0]))],
            metas=[StringVariable("SMILES")]
        )
        table = Table.from_list(
            domain, [list(keys) + [smile] for keys, smile in zip(maccs_keys, valid_smiles)]
        )

        self.Outputs.maccs_table.send(table)
        self.info_label.setText(f"Processed {len(maccs_keys)} valid SMILES entries.")

if __name__ == "__main__":
    from Orange.widgets.utils.widgetpreview import WidgetPreview

    WidgetPreview(SmilesToMACCSWidget).run()
