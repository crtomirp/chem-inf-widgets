
""" ===  1. Loading PAINS SMARTS from a JSON File === """

import os
import json
import warnings
from rdkit import Chem, RDLogger

# Suppress RDKit warnings
warnings.filterwarnings("ignore", category=RuntimeWarning)
RDLogger.DisableLog('rdApp.warning')
RDLogger.DisableLog('rdApp.error')

def load_pains_smarts(filepath="smartspains.json"):
    """Load PAINS SMARTS from a JSON file."""
    if not os.path.exists(filepath):
        print(f"⚠ Warning: PAINS file not found at {filepath}")
        return []
    try:
        with open(filepath, "r", encoding="utf-8") as jsonfile:
            pains_data = json.load(jsonfile)
            smarts_list = [entry["SMARTS"] for entry in pains_data if "SMARTS" in entry]
            print("✅ Loaded PAINS SMARTS:", smarts_list)  # Debug print
            return smarts_list
    except Exception as e:
        print(f"❌ Error loading PAINS SMARTS: {e}")
        return []

PAINS_SMARTS = load_pains_smarts()

""" ===  2. Descriptor Calculation Functions === """

def is_pains(mol):
    """Check if molecule contains PAINS substructures."""
    for smarts in PAINS_SMARTS:
        pattern = Chem.MolFromSmarts(smarts)
        if pattern is None:
            print(f"⚠ Warning: Invalid PAINS SMARTS: {smarts}")
            continue  # Skip invalid SMARTS
        if mol.HasSubstructMatch(pattern):
            print(f"🔬 PAINS detected: {smarts}")  # Debugging
            return True
    return False

def get_highlighted_atoms(mol):
    """
    For a given molecule, returns a sorted list of atom indices that
    are part of any substructure match to the PAINS SMARTS.
    """
    highlighted = set()
    for smarts in PAINS_SMARTS:
        pattern = Chem.MolFromSmarts(smarts)
        if pattern is None:
            continue  # Skip invalid SMARTS
        matches = mol.GetSubstructMatches(pattern)
        for match in matches:
            highlighted.update(match)
    return sorted(list(highlighted))

""" ==== 2. Processing the Data and Sending the Output ===  """

from rdkit.Chem import Descriptors, rdMolDescriptors
from rdkit.Chem.QED import qed

def lipinski_violations(mol):
    """Returns the number of Lipinski's Rule of Five violations."""
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
    return violations

def is_veber(mol):
    """Checks Veber's rule (rotatable bonds ≤ 10 and TPSA ≤ 140 Å²)."""
    rb = rdMolDescriptors.CalcNumRotatableBonds(mol)
    tpsa = rdMolDescriptors.CalcTPSA(mol)
    return rb <= 10 and tpsa <= 140

def is_reactive(mol):
    """Estimate reactivity based on presence of electrophilic groups."""
    reactive_patterns = [
        "[O=CN]",         # Isocyanate
        "[N+]([O-])=O",   # Nitro group
        "[C#N]",          # Nitrile
        "[O=CS]",         # Thiocarbonyl
        "[N=N]",          # Azide/Nitroso
    ]
    for smarts in reactive_patterns:
        pattern = Chem.MolFromSmarts(smarts)
        if pattern and mol.HasSubstructMatch(pattern):
            return True
    return False

def calculate_drug_score(qed_score, lipinski_vio, pains, veber, reactivity):
    """Composite score from different drug-likeness criteria."""
    score = qed_score  # Start with QED Score (0-1)
    if lipinski_vio > 1:
        score -= 0.2  # Penalize multiple Lipinski violations
    if pains:
        score -= 0.3  # PAINS match significantly reduces drug-likeness
    if not veber:
        score -= 0.1  # Veber rule failure slightly reduces score
    if reactivity:
        score -= 0.2  # Reactivity decreases stability (drug-likeness)
    
    return max(score, 0.0)  # Ensure the score doesn't go below zero   

""" === 3. Designing the Widget User Interface === """

import numpy as np
from Orange.widgets.widget import OWWidget, Input, Output
from Orange.widgets import gui
from Orange.data import Table, Domain, StringVariable, ContinuousVariable
from PyQt5.QtWidgets import QCheckBox, QPushButton, QComboBox

class DrugFilterWidget(OWWidget):
    """Orange3 widget for filtering molecules based on drug-likeness rules."""

    name = "Drug Filter"
    description = "Filters molecules using Lipinski, PAINS, Veber, QED scores, and Reactivity."
    category = "Chemoinformatics"
    icon = "icons/filter.png"
    priority = 10

    inputs = [("Input Table", Table, "set_data")]
    outputs = [("All Compounds", Table)]
    want_main_area = False

    def __init__(self):
        super().__init__()
        self.data = None

        # GUI Elements
        self.info_label = gui.label(self.controlArea, self, "Awaiting molecular data...")

        # Lipinski Violation Selection
        self.lipinski_box = gui.widgetBox(self.controlArea, "Lipinski Rule Settings")
        self.lipinski_checkbox = QCheckBox("Apply Lipinski's Rule")
        self.lipinski_threshold = QComboBox()
        self.lipinski_threshold.addItems(["1", "2", "3+"])
        self.lipinski_box.layout().addWidget(self.lipinski_checkbox)
        self.lipinski_box.layout().addWidget(self.lipinski_threshold)

        # Additional Drug-Likeness Filters
        self.pains_checkbox = QCheckBox("Apply PAINS Filter")
        self.veber_checkbox = QCheckBox("Apply Veber's Rule")
        self.reactivity_checkbox = QCheckBox("Apply Reactivity Filter")
        for checkbox in [self.pains_checkbox, self.veber_checkbox, self.reactivity_checkbox]:
            self.controlArea.layout().addWidget(checkbox)

        # Filter button
        self.process_button = QPushButton("Filter Molecules")
        self.process_button.clicked.connect(self.filter_molecules)
        self.controlArea.layout().addWidget(self.process_button)

""" ==== 4. Processing the Data and Sending the Output ===  """

    def set_data(self, data):
        """Handle input data."""
        self.data = data
        if self.data and len(self.data) > 0:
            self.info_label.setText("Input data received. Ready to filter.")
        else:
            self.info_label.setText("No valid data received.")

    def filter_molecules(self):
        """Apply filtering and calculate descriptors."""
        smiles_column_name = "SMILES"
        if self.data is None or len(self.data) == 0:
            self.info_label.setText("No input data to filter.")
            return

        # Assumes that SMILES are stored in the first meta column.
        smiles_column = self.data.metas[:, 0]
        all_results = []

        for smile in smiles_column:
            mol = Chem.MolFromSmiles(smile)
            if not mol:
                continue

            qed_score = qed(mol)
            pains_match = is_pains(mol)
            lipinski_vio = lipinski_violations(mol)
            veber_pass = is_veber(mol)
            reactive = is_reactive(mol)
            drug_score = calculate_drug_score(qed_score, lipinski_vio, pains_match, veber_pass, reactive)
            # Get the highlighted atoms from PAINS matches (if any)
            highlighted_atoms_list = get_highlighted_atoms(mol)
            highlighted_atoms_str = ", ".join(map(str, highlighted_atoms_list)) if highlighted_atoms_list else ""

            # Row structure: [SMILES (meta), QED Score, Lipinski Violations, PAINS Match,
            # Veber Rule, Reactivity, Drug Score, Highlighted Atoms (meta)]
            all_results.append([
                smile, 
                qed_score, 
                lipinski_vio, 
                float(pains_match), 
                float(veber_pass), 
                float(reactive), 
                drug_score,
                highlighted_atoms_str
            ])

        self.send_output_table(all_results, "All Compounds")

    def send_output_table(self, results, output_name):
        """Convert results into an Orange Table and send it as output."""
        if not results:
            self.send(output_name, None)
            return

        # Numeric features: columns 1 through 6 (indices 1-6)
        numeric_data = np.array([row[1:7] for row in results], dtype=float)
        # Meta data: SMILES (index 0) and Highlighted Atoms (index 7)
        metas = np.array([[row[0], row[7]] for row in results], dtype=object)

        domain = Domain(
            [ContinuousVariable("QED Score"),
             ContinuousVariable("Lipinski Violations"),
             ContinuousVariable("PAINS Match"),
             ContinuousVariable("Veber Rule"),
             ContinuousVariable("Reactivity"),
             ContinuousVariable("Drug Score")],
            metas=[StringVariable("SMILES"), StringVariable("Highlighted Atoms")]
        )

        data_table = Table.from_numpy(domain, numeric_data, metas=metas)
        self.send(output_name, data_table)






