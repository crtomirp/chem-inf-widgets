Below is an example of comprehensive documentation for the Drug Filter Widget. This documentation explains what the widget does, the input it expects, how to use it (both within Orange Canvas and via standalone Python scripts), and includes example workflows.

---

# Drug Filter Widget Documentation

## Overview

The **Drug Filter Widget** is designed for chemoinformatics workflows within Orange. It computes various molecular descriptors using RDKit (such as QED, Lipinski violations, molecular weight, LogP, HBD, HBA, rotatable bonds, and TPSA) and filters molecules based on user-selected drug-likeness rules. The widget supports:

- **Filtering Rules:**  
  - **Lipinski:** Molecules with 0 or 1 Lipinski violation pass.
  - **Veber:** Molecules pass if they have ≤ 10 rotatable bonds and TPSA ≤ 140.
  - **Lipinski + Veber:** Molecules must pass both the above criteria.
  - **None:** No filtering is applied; all molecules are forwarded.
  
- **PAINS (Pan Assay Interference Compounds) Alerts:**  
  The widget checks for problematic substructures (defined in `smartspains.json`) and—if enabled—displays the atom indices where the PAINS patterns are matched.

- **Composite Drug Score:**  
  A drug-likeness score is calculated based on the QED value and penalties for multiple Lipinski violations, PAINS matches, Veber rule failure, and reactivity issues (the latter is a placeholder).

- **Progress Monitoring:**  
  A progress bar shows the number of molecules processed in real time.

## Input Data

The widget expects an Orange data table in which the **first meta attribute** contains the SMILES strings of the molecules. For example, a CSV file might have:

```csv
SMILES,Name
CCO,Ethanol
CCCC,Butane
CCN,Ethylamine
```

When this file is imported as an Orange table, the SMILES strings should appear as meta attributes.

## User Interface Options

1. **Filtering Rule (Dropdown):**  
   Choose one of the following rules:
   - **Lipinski:** Molecules with 0 or 1 violation pass.
   - **Veber:** Based on rotatable bonds and TPSA.
   - **Lipinski + Veber:** Both rules must be satisfied.
   - **None:** No filtering is performed.

2. **Molecule Selection (Dropdown):**  
   Decide which molecules are forwarded:
   - **Forward All Molecules:** All molecules are sent to output.
   - **Within Criteria:** Only molecules that pass the selected rule.
   - **Out of Criteria:** Only molecules that do not pass the selected rule.

3. **Highlight PAINS Substructures (Checkbox):**  
   When checked, the widget searches for PAINS patterns and displays a comma-separated list of atom indices (in the meta field "Highlighted Atoms") for those molecules that match PAINS patterns.

4. **Filter Molecules (Button):**  
   Clicking this button starts processing the molecules, updates the progress bar, and sends the filtered results to the output.

## Output Table Format

The widget produces an Orange data table with:

- **Features (Numeric Attributes):**
  1. QED Score
  2. Lipinski Violations
  3. MW (Molecular Weight)
  4. LogP
  5. HBD (Hydrogen Bond Donors)
  6. HBA (Hydrogen Bond Acceptors)
  7. Rotatable Bonds
  8. TPSA (Topological Polar Surface Area)
  9. PAINS Match (flag: 1 if a PAINS substructure is found, 0 otherwise)
  10. Veber Rule (flag: 1 if pass, 0 otherwise)
  11. Reactivity (placeholder; always 0)
  12. Drug Score

- **Meta Attributes:**
  - **SMILES:** The original SMILES string.
  - **PAINS regID:** A comma-separated list of PAINS regID values (if any).  
    (If PAINS highlighting is enabled, this field is set to `"None"` as only atom indices are reported.)
  - **Criteria:** `"Pass"` or `"Fail"` (depending on whether the molecule meets the filtering criteria).
  - **Highlighted Atoms (optional):** A comma-separated list of atom indices for PAINS matches (only if PAINS highlighting is enabled).

## Example Workflows

### 1. Using the Widget in Orange Canvas

1. **Load Data:**  
   Import your dataset (e.g., a CSV file) into Orange ensuring that the SMILES strings are placed in the first meta attribute.

2. **Connect the Widget:**  
   Drag the **Drug Filter Widget** into your Orange workflow and connect your data table to it.

3. **Configure Options:**  
   - Select the filtering rule (e.g., "Lipinski + Veber").
   - Choose the molecule selection mode (e.g., "Within Criteria").
   - Check the "Highlight PAINS Substructures" option if you want detailed atom indices for PAINS matches.

4. **Run Filtering:**  
   Click the **Filter Molecules** button. The progress bar will display the processing progress. When finished, the output table (with descriptors and meta information) will be available for further analysis or visualization.

### 2. Standalone Example Using Python Code

Below is an example of how to create a simple input table and use the widget in a standalone mode:

```python
import numpy as np
from Orange.data import Domain, Table, StringVariable
from AnyQt.QtWidgets import QApplication
import sys

# Create a simple input table with SMILES strings in the meta column.
smiles = np.array([["CCO"], ["CCCC"], ["CCN"]], dtype=object)
# Define the domain with a meta attribute for SMILES.
domain = Domain([], metas=[StringVariable("SMILES")])
data_table = Table.from_numpy(domain, np.empty((3, 0)), metas=smiles)

# Import the widget (assuming the widget code is in drug_filter.py)
from drug_filter import DrugFilterWidget

# Set up the application and widget.
app = QApplication(sys.argv)
widget = DrugFilterWidget()
widget.set_data(data_table)  # Provide the data table as input.
widget.show()
sys.exit(app.exec_())
```

In this example, a small dataset with three molecules is created and passed into the widget. The widget processes the molecules and (if configured) displays progress and outputs the computed table.

## Troubleshooting

- **No Data in Output:**  
  Ensure that your input table contains valid SMILES strings in the first meta column.

- **PAINS SMARTS File:**  
  Make sure the `smartspains.json` file is located in the same directory as the widget code. If not found, PAINS checking will be skipped, and a warning message will be printed.

- **Dependencies:**  
  Confirm that RDKit, Orange3, and PyQt5 are installed and properly configured in your environment.

## Conclusion

The Drug Filter Widget is a powerful tool for screening chemical compounds for drug-likeness. It seamlessly integrates into Orange workflows, allowing you to filter molecules based on widely used criteria (Lipinski, Veber, etc.) and to flag problematic substructures via PAINS alerts. Customize your filtering options, monitor progress with the built-in progress bar, and use the output table for further data visualization and analysis in Orange.

---

This documentation should serve as both a guide and a reference for using the Drug Filter Widget effectively in your chemoinformatics workflows.
