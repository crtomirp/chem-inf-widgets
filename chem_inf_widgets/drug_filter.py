import os
import json
import warnings
from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors, rdMolDescriptors
from rdkit.Chem.QED import qed
import numpy as np

from Orange.widgets.widget import OWWidget, Input, Output
from Orange.widgets import gui
from Orange.data import Table, Domain, StringVariable, ContinuousVariable
from PyQt5.QtWidgets import QCheckBox, QPushButton, QComboBox, QProgressBar

# Suppress RDKit warnings
warnings.filterwarnings("ignore", category=RuntimeWarning)
RDLogger.DisableLog('rdApp.warning')
RDLogger.DisableLog('rdApp.error')


def load_pains_smarts(filepath="smartspains.json"):
    """Load PAINS SMARTS and regID from a JSON file."""
    if not os.path.exists(filepath):
        print(f"⚠ Warning: PAINS file not found at {filepath}")
        return []
    try:
        with open(filepath, "r", encoding="utf-8") as jsonfile:
            pains_data = json.load(jsonfile)
            return pains_data  # Expecting list of dicts with keys "SMARTS" and "regID"
    except Exception as e:
        print(f"❌ Error loading PAINS SMARTS: {e}")
        return []


PAINS_DATA = load_pains_smarts()


def get_pains_matches(mol):
    """
    Check for PAINS matches and return a list of matching regID values.
    (This version is used when detailed atom-level information is not required.)
    """
    matched_pains = []
    for entry in PAINS_DATA:
        smarts = entry.get("SMARTS")
        regID = entry.get("regID")
        pattern = Chem.MolFromSmarts(smarts)
        if pattern and mol.HasSubstructMatch(pattern):
            matched_pains.append(regID)
    return matched_pains


def get_pains_matches_with_atoms(mol):
    """
    Check for PAINS matches and return a list of tuples containing the matched atom indices.
    Each element in the returned list is the tuple (or tuple of tuples) of atom indices
    for the matched substructure.
    """
    matches = []
    for entry in PAINS_DATA:
        smarts = entry.get("SMARTS")
        pattern = Chem.MolFromSmarts(smarts)
        if pattern and mol.HasSubstructMatch(pattern):
            match_indices = mol.GetSubstructMatches(pattern)
            matches.append(match_indices)
    return matches


def lipinski_violations(mol):
    """
    Returns the number of Lipinski's Rule of Five violations and calculated descriptors.
    Criteria: MW > 500, LogP > 5, HBD > 5, HBA > 10.
    """
    mw = Descriptors.MolWt(mol)
    logp = Descriptors.MolLogP(mol)
    hbd = rdMolDescriptors.CalcNumHBD(mol)
    hba = rdMolDescriptors.CalcNumHBA(mol)
    violations = sum([
        mw > 500,
        logp > 5,
        hbd > 5,
        hba > 10
    ])
    return violations, mw, logp, hbd, hba


def is_veber(mol):
    """
    Checks Veber's rule and returns a tuple: (veber_pass, rotatable_bonds, tpsa).
    Criteria: Rotatable Bonds <= 10 and TPSA <= 140.
    """
    rb = rdMolDescriptors.CalcNumRotatableBonds(mol)
    tpsa = rdMolDescriptors.CalcTPSA(mol)
    veber_pass = rb <= 10 and tpsa <= 140
    return veber_pass, rb, tpsa


def calculate_drug_score(qed_score, lipinski_vio, pains, veber, reactivity):
    """
    Composite score calculation.
    
    The score starts with the QED score (0–1) and is penalized if:
      - More than 1 Lipinski violation (–0.2)
      - PAINS match is found (–0.3)
      - Veber rule fails (–0.1)
      - Reactivity is flagged (–0.2)
    The final score is non-negative.
    """
    score = qed_score  # Start with QED Score (0-1)
    if lipinski_vio > 1:
        score -= 0.2  # Penalize multiple Lipinski violations
    if pains:
        score -= 0.3  # PAINS match reduces drug-likeness
    if not veber:
        score -= 0.1  # Veber rule failure slightly reduces score
    if reactivity:
        score -= 0.2  # Reactivity decreases stability
    return max(score, 0.0)


class DrugFilterWidget(OWWidget):
    """
    Orange3 widget for filtering molecules based on drug-likeness rules.

    The widget calculates molecular descriptors (using RDKit) and applies filtering based on:
      - Lipinski’s Rule of Five (number of violations)
      - Veber’s Rule (rotatable bonds and TPSA)
      - A combined Lipinski + Veber criteria

    It also checks for PAINS alerts (problematic substructures) and (optionally) reports
    the atom indices that match the PAINS substructure patterns.

    The user can select which filtering rule to use and which molecules to forward.
    """
    name = "Drug Filter"
    description = "Filters molecules using Lipinski, PAINS, Veber, QED scores, and Reactivity."
    category = "Chemoinformatics"
    icon = "icons/filter.png"
    priority = 20

    inputs = [("Input Table", Table, "set_data")]
    outputs = [("Filtered Compounds", Table)]
    want_main_area = False

    def __init__(self):
        super().__init__()
        self.data = None

        # Information label
        self.info_label = gui.label(self.controlArea, self, "Awaiting molecular data...")

        # Filtering rule selection
        self.filter_rule_combo = QComboBox()
        self.filter_rule_combo.addItems(["Lipinski", "Veber", "Lipinski + Veber", "None"])
        box_rule = gui.widgetBox(self.controlArea, "Select Filtering Rule")
        box_rule.layout().addWidget(self.filter_rule_combo)
        gui.separator(self.controlArea)

        # Molecule selection mode
        self.selection_combo = QComboBox()
        self.selection_combo.addItems(["Forward All Molecules", "Within Criteria", "Out of Criteria"])
        box_sel = gui.widgetBox(self.controlArea, "Molecule Selection")
        box_sel.layout().addWidget(self.selection_combo)
        gui.separator(self.controlArea)

        # PAINS highlighting checkbox
        self.highlight_pains_checkbox = QCheckBox("Highlight PAINS Substructures")
        self.controlArea.layout().addWidget(self.highlight_pains_checkbox)
        gui.separator(self.controlArea)

        # Process button to trigger filtering
        self.process_button = QPushButton("Filter Molecules")
        self.process_button.clicked.connect(self.filter_molecules)
        self.controlArea.layout().addWidget(self.process_button)

        # Progress bar to monitor processing progress
        self.progressBar = QProgressBar()
        self.progressBar.setVisible(False)
        self.controlArea.layout().addWidget(self.progressBar)

    def set_data(self, data):
        """Set the input data table. Assumes the first meta attribute holds the SMILES string."""
        self.data = data
        if self.data is not None and len(self.data) > 0:
            self.info_label.setText("Input data received. Ready to filter.")
        else:
            self.info_label.setText("No valid data received.")

    def filter_molecules(self):
        """Process each molecule: calculate descriptors, apply filtering, and forward the selected set."""
        # Get user selections
        filter_rule = self.filter_rule_combo.currentText()  # "Lipinski", "Veber", "Lipinski + Veber", "None"
        selection_mode = self.selection_combo.currentText()   # "Forward All Molecules", "Within Criteria", "Out of Criteria"
        highlight_pains = self.highlight_pains_checkbox.isChecked()

        if self.data is None:
            return

        # Expect SMILES strings to be in the first meta column
        smiles_column = self.data.metas[:, 0]
        num_molecules = len(smiles_column)
        results = []

        # Initialize progress bar
        self.progressBar.setMaximum(num_molecules)
        self.progressBar.setValue(0)
        self.progressBar.setVisible(True)

        for i, smile in enumerate(smiles_column):
            # Update progress bar
            self.progressBar.setValue(i + 1)

            mol = Chem.MolFromSmiles(smile)
            if not mol:
                continue

            # Calculate descriptors
            qed_score = qed(mol)
            lipinski_vio, mw, logp, hbd, hba = lipinski_violations(mol)
            veber_pass, rb, tpsa = is_veber(mol)
            reactive = False  # Placeholder for reactivity check
            drug_score = calculate_drug_score(qed_score, lipinski_vio, bool([]), veber_pass, reactive)

            # Determine if molecule passes the selected filtering rule
            if filter_rule == "Lipinski":
                criteria_pass = (lipinski_vio <= 1)
            elif filter_rule == "Veber":
                criteria_pass = veber_pass
            elif filter_rule == "Lipinski + Veber":
                criteria_pass = (lipinski_vio <= 1) and veber_pass
            else:  # "None" filtering rule
                criteria_pass = True
            criteria = "Pass" if criteria_pass else "Fail"

            # Depending on selection mode, decide whether to include the molecule.
            if selection_mode == "Within Criteria" and not criteria_pass:
                continue
            if selection_mode == "Out of Criteria" and criteria_pass:
                continue

            # PAINS check (with optional atom-level detail)
            if highlight_pains:
                pains_matches_detail = get_pains_matches_with_atoms(mol)
                # Flatten the list of tuples into a set of unique atom indices.
                all_indices = set()
                for match in pains_matches_detail:
                    for sub_match in match:
                        all_indices.update(sub_match)
                # Create a comma-separated string of indices (or "None" if no indices were found).
                pains_atoms = ", ".join(map(str, sorted(all_indices))) if all_indices else "None"
                pains_flag = 1.0 if all_indices else 0.0
                # Since we are not using PAINS regIDs in this branch, assign a default value.
                pains_regids = "None"
            else:
                pains_matches = get_pains_matches(mol)
                pains_regids = ", ".join(pains_matches) if pains_matches else "None"
                pains_atoms = None
                pains_flag = 1.0 if pains_matches else 0.0

            # Build a result row (this happens regardless of the PAINS check branch)
            row = [
                smile,         # Meta: SMILES
                qed_score,     # QED Score
                lipinski_vio,  # Lipinski Violations
                mw,            # MW
                logp,          # LogP
                hbd,           # HBD
                hba,           # HBA
                rb,            # Rotatable Bonds
                tpsa,          # TPSA
                pains_flag,    # PAINS Match flag (1.0 if match, else 0.0)
                1.0 if veber_pass else 0.0,  # Veber Rule flag
                0.0,           # Reactivity (placeholder)
                drug_score,    # Drug Score
                pains_regids,  # Meta: PAINS regID
                criteria       # Meta: Criteria ("Pass"/"Fail")
            ]
            if highlight_pains:
                row.append(pains_atoms)  # Only add the extra meta if highlighting is enabled

            results.append(row)

        # Hide progress bar after processing is complete
        self.progressBar.setVisible(False)
        self.send_output_table(results, highlight_pains)

    def send_output_table(self, results, highlight_pains):
        """Convert results into an Orange Table and send it as output."""
        if not results:
            self.send("Filtered Compounds", None)
            return

        # Build the domain.
        features = [
            ContinuousVariable("QED Score"),
            ContinuousVariable("Lipinski Violations"),
            ContinuousVariable("MW"),
            ContinuousVariable("LogP"),
            ContinuousVariable("HBD"),
            ContinuousVariable("HBA"),
            ContinuousVariable("Rotatable Bonds"),
            ContinuousVariable("TPSA"),
            ContinuousVariable("PAINS Match"),
            ContinuousVariable("Veber Rule"),
            ContinuousVariable("Reactivity"),
            ContinuousVariable("Drug Score")
        ]

        # Meta attributes: SMILES, PAINS regID, Criteria, and optionally PAINS Atoms.
        if highlight_pains:
            meta_vars = [StringVariable("SMILES"), StringVariable("PAINS regID"),
                         StringVariable("Criteria"), StringVariable("Highlighted Atoms")]
            meta_data = np.array([[row[0], row[13], row[14], row[15]] for row in results], dtype=object)
        else:
            meta_vars = [StringVariable("SMILES"), StringVariable("PAINS regID"),
                         StringVariable("Criteria")]
            meta_data = np.array([[row[0], row[13], row[14]] for row in results], dtype=object)

        domain = Domain(features, metas=meta_vars)
        # Numeric data are the features (columns 1 to 12)
        numeric_data = np.array([row[1:13] for row in results], dtype=float)
        data_table = Table.from_numpy(domain, numeric_data, metas=meta_data)
        self.send("Filtered Compounds", data_table)


if __name__ == "__main__":
    # This block is for testing the widget standalone.
    from AnyQt.QtWidgets import QApplication
    import sys

    app = QApplication(sys.argv)
    widget = DrugFilterWidget()
    widget.show()
    sys.exit(app.exec_())
