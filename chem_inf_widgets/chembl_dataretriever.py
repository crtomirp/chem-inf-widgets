import requests
import pandas as pd
import numpy as np
from Orange.widgets.widget import OWWidget, Input, Output
from Orange.widgets import gui
from Orange.data import Table, Domain, StringVariable, ContinuousVariable


class ChEMBLBioactivityWidget(OWWidget):
    """Orange widget to fetch ChEMBL bioactivity data with drug properties."""
    
    name = "ChEMBL Bioactivity Retriever"
    description = "Fetches bioactivity data with drug design properties"
    icon = "icons/chembl.png"
    priority = 10

    class Outputs:
        output_data = Output("Bioactivity Data", Table)

    def __init__(self):
        super().__init__()
        self.target_id = ""
        self._build_ui()

    def _build_ui(self):
        """Construct the user interface."""
        layout = gui.widgetBox(self.controlArea, "Retrieve Bioactivity Data")
        gui.label(layout, self, "Enter ChEMBL Target ID (e.g., CHEMBL2095150):")
        gui.lineEdit(layout, self, "target_id", placeholderText="CHEMBLxxxxxx")
        gui.button(layout, self, "Fetch Data", callback=self.fetch_bioactivity_data)
        self.info_label = gui.label(layout, self, "Status: Awaiting input.")

    def fetch_bioactivity_data(self):
        """Main data retrieval and processing workflow."""
        target_id = self.target_id.strip()
        if not target_id:
            self.info_label.setText("Status: Please enter a ChEMBL Target ID.")
            return

        self.info_label.setText(f"Status: Fetching data for {target_id}...")

        try:
            # Step 1: Fetch raw data from ChEMBL
            df = self._fetch_chembl_data(target_id)
            if df is None or df.empty:
                return

            # Step 2: Process IC50 values and rename columns
            df = self._process_ic50_values(df)
            df = df.rename(columns={'canonical_smiles': 'SMILES'})

            # Step 3: Calculate drug properties if SMILES available
            if 'SMILES' in df.columns:
                df = self._calculate_drug_properties(df)

            # Step 4: Filter and organize columns
            df = self._filter_columns(df)

            # Step 5: Create Orange Table
            table = self._create_orange_table(df)
            
            self.info_label.setText(f"Status: Retrieved {len(df)} records for {target_id}.")
            self.Outputs.output_data.send(table)

        except Exception as e:
            self._handle_error(f"Error: {str(e)}")

    def _fetch_chembl_data(self, target_id):
        """Fetch bioactivity data from ChEMBL API."""
        try:
            response = requests.get(
                "https://www.ebi.ac.uk/chembl/api/data/activity.json",
                params={
                    "target_chembl_id": target_id,
                    "standard_type": "IC50",
                    "limit": 1000
                }
            )
            response.raise_for_status()
            
            data = response.json().get("activities", [])
            if not data:
                self.info_label.setText(f"Status: No data found for {target_id}.")
                self.Outputs.output_data.send(None)
                return None

            df = pd.DataFrame(data)
            
            # Clean pChEMBL values
            if "pchembl_value" in df.columns:
                df["pchembl_value"] = pd.to_numeric(df["pchembl_value"], errors="coerce")
                df = df.dropna(subset=["pchembl_value"])
            
            return df

        except requests.exceptions.RequestException as e:
            self._handle_error(f"Network error: {str(e)}")
            return None

    def _process_ic50_values(self, df):
        """Convert IC50 values to nM units."""
        if 'standard_value' in df.columns and 'standard_units' in df.columns:
            df['IC50_nM'] = df.apply(lambda row: self._convert_to_nM(row), axis=1)
            df = df.drop(columns=['standard_value', 'standard_units'])
        return df

    def _convert_to_nM(self, row):
        """Convert IC50 values to nanomolar units."""
        try:
            value = float(row['standard_value'])
            unit = row['standard_units'].lower()
            
            conversions = {
                'm': value * 1e9,
                'µm': value * 1e3,
                'um': value * 1e3,  # Handle different micro symbols
                'nm': value,
                'nmol/l': value,
                'pm': value * 1e-3  # Handle pico molar if needed
            }
            
            return conversions.get(unit, np.nan)
        except:
            return np.nan

    def _calculate_drug_properties(self, df):
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

                # Basic properties
                df.at[idx, 'hbd'] = Lipinski.NumHDonors(mol)
                df.at[idx, 'hba'] = Lipinski.NumHAcceptors(mol)
                df.at[idx, 'rotable_bonds'] = Descriptors.NumRotatableBonds(mol)
                df.at[idx, 'mw'] = Descriptors.MolWt(mol)
                df.at[idx, 'tpsa'] = Descriptors.TPSA(mol)
                df.at[idx, 'logp'] = Descriptors.MolLogP(mol)

                # Lipinski rule violations
                violations = 0
                violations += 1 if df.at[idx, 'mw'] > 500 else 0
                violations += 1 if df.at[idx, 'logp'] > 5 else 0
                violations += 1 if df.at[idx, 'hbd'] > 5 else 0
                violations += 1 if df.at[idx, 'hba'] > 10 else 0
                df.at[idx, 'lipinski_deviations'] = violations

            except Exception as e:
                continue

        # Convert all new columns to numeric
        for col in prop_columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")
        
        return df

    def _filter_columns(self, df):
        """Filter columns to only include requested metadata and numerical data."""
        # Define required columns
        meta_cols = [
            'SMILES', 'molecule_chembl_id', 'target_chembl_id',
            'assay_chembl_id', 'document_chembl_id', 'target_organism',
            'target_name'
        ]
        
        num_cols = [
            'pchembl_value', 'IC50_nM', 'hbd', 'hba', 'rotable_bonds',
            'mw', 'tpsa', 'logp', 'lipinski_deviations'
        ]

        # Find existing columns in dataframe
        existing_meta = [col for col in meta_cols if col in df.columns]
        existing_num = [col for col in num_cols if col in df.columns]
        
        return df[existing_num + existing_meta]

    def _create_orange_table(self, df):
        """Convert pandas DataFrame to Orange Table with proper domain."""
        # Separate numerical and meta columns
        num_cols = [col for col in df.columns if col in [
            'pchembl_value', 'IC50_nM', 'hbd', 'hba', 'rotable_bonds',
            'mw', 'tpsa', 'logp', 'lipinski_deviations'
        ]]
        
        meta_cols = [col for col in df.columns if col not in num_cols]

        # Create domain
        domain = Domain(
            [ContinuousVariable(col) for col in num_cols],
            metas=[StringVariable(col) for col in meta_cols]
        )

        # Prepare data arrays
        X = df[num_cols].to_numpy(dtype=float)
        metas = df[meta_cols].to_numpy(dtype=object)

        return Table.from_numpy(domain, X=X, metas=metas)

    def _handle_error(self, message):
        """Centralized error handling."""
        self.info_label.setText(message)
        self.Outputs.output_data.send(None)


if __name__ == "__main__":
    from Orange.widgets.utils.widgetpreview import WidgetPreview
    WidgetPreview(ChEMBLBioactivityWidget).run()
