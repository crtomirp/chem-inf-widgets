Below is a sample documentation for the **Molecular Standardizer** widget, which describes its functionality, provides usage examples, and offers suggestions for integrating it into your cheminformatics workflows.

---

# Molecular Standardizer Widget Documentation

## Overview

The **Molecular Standardizer** widget is designed to process and standardize chemical structures provided as SMILES strings. Using the RDKit **MolStandardize** module, the widget applies a series of standardization operations—such as cleanup, normalization, metal disconnection, and more—to ensure that molecules are represented in a consistent and reproducible format. This is essential for downstream analyses such as fingerprinting, clustering, and predictive modeling.

### Key Features

- **Input/Output:**  
  - **Input:** An Orange **Table** with a meta column containing SMILES strings (named "SMILES").  
  - **Output:** An Orange **Table** with three meta columns:
    - **Original SMILES:** The SMILES string as originally provided.
    - **SMILES:** The standardized SMILES string.
    - **Change Log:** A summary of the modifications made during standardization.
  
- **Standardization Operations:**  
  The widget offers a variety of operations that can be applied individually or in combination:
  - **ValidateSmiles:** Checks that the SMILES string is valid.
  - **Cleanup:** Performs general clean-up of the molecular structure.
  - **Normalize:** Applies predefined normalization rules.
  - **MetalDisconnector:** Separates metal atoms from organic fragments.
  - **LargestFragmentChooser:** Selects the largest fragment when multiple fragments are present.
  - **Reionizer:** Adjusts protonation states to a standard form.
  - **Uncharger:** Removes formal charges from the molecule.
  - **Tautomer:** Canonicalizes tautomeric forms.

- **Progress Bar Integration:**  
  A progress bar provides real-time feedback during the standardization process, helping users track the progress as each SMILES string is processed.

## Requirements

- **Orange3:** The Orange data mining framework.
- **RDKit:** For cheminformatics functions and molecule standardization.
- **PyQt5:** For building the widget’s graphical user interface (GUI).

You can install the necessary packages using:

```bash
pip install orange3 rdkit-pypi PyQt5 numpy
```

> **Note:** RDKit installation may require additional steps depending on your operating system. Refer to the [RDKit Installation Guide](https://www.rdkit.org/docs/Install.html) for further details.

## Usage Instructions

### Input Data

- **Data Format:**  
  The widget expects an Orange **Table** where one of the meta columns is named `"SMILES"`. This column should contain the molecular representations (SMILES strings) that you wish to standardize.

### User Interface

- **Operation Selection:**  
  Use the provided checkboxes to select which standardization operations should be applied. For example, you may choose to perform a **Cleanup** and **Normalize** operation simultaneously.
  
- **Process Button:**  
  Click the **"Apply Standardization"** button to begin processing. The progress bar will indicate the current progress, and once complete, the widget will update the information label and output the standardized data.

### Example Workflow

Suppose you have a dataset with molecules represented as SMILES strings and wish to standardize them before further analysis. Here is a step-by-step example:

1. **Data Import:**  
   - Load your dataset using an Orange widget such as **File** or **Data Table**.
   - Ensure that the dataset contains a meta column named `"SMILES"`.

2. **Connect to Molecular Standardizer:**  
   - Connect the output of your data import widget to the **Molecular Standardizer** widget.
   - In the **Molecular Standardizer** widget, select one or more operations (e.g., **Cleanup**, **Normalize**, and **LargestFragmentChooser**) by checking the corresponding checkboxes.

3. **Process the Data:**  
   - Click the **"Apply Standardization"** button.
   - Watch the progress bar update as each SMILES string is processed.
   - Once processing is complete, the widget outputs a new Orange **Table** containing the original SMILES, the standardized SMILES, and a change log.

4. **Downstream Analysis:**  
   - Use the standardized data as input for further analysis such as fingerprint calculation, clustering, or QSAR modeling.

### Sample Code Usage (Standalone)

Below is a brief code snippet demonstrating how you might invoke the widget programmatically:

```python
from Orange.data import Table, Domain, StringVariable
from StandardizeMoleculesWidget import StandardizeMoleculesWidget  # Ensure your widget is importable

# Create a simple dataset with SMILES in the meta column
smiles_list = ["CCO", "c1ccccc1", "C1=CC=CN=C1", "C1CC1", "INVALID_SMILES"]
domain = Domain([], metas=[StringVariable("SMILES")])
data = Table.from_list(domain, [[None, s] for s in smiles_list])

# Instantiate the widget and set the data
widget = StandardizeMoleculesWidget()
widget.set_data(data)

# Optionally, select specific operations by checking the desired checkboxes (simulated)
widget.checkboxes["Cleanup"].setChecked(True)
widget.checkboxes["Normalize"].setChecked(True)
widget.update_selected_operations()

# Run standardization
widget.standardize_smiles()

# Access the output table (widget.standardized_data) and inspect results
print(widget.standardized_data)
```

## Importance and Applications

### Why Standardize Molecules?

- **Consistency:**  
  Standardization ensures that molecules are represented in a consistent format, which is crucial for reliable comparisons.
  
- **Data Quality:**  
  Cleaning up and normalizing structures minimizes errors in downstream analyses like similarity searches or machine learning.
  
- **Preprocessing for Modeling:**  
  Many cheminformatics and QSAR models require a uniform representation of molecular structures. Standardized SMILES strings can be directly used for feature extraction (e.g., fingerprint calculation).

### Possible Applications

- **Drug Discovery:**  
  Standardize chemical libraries to ensure consistency before virtual screening and similarity searches.
  
- **QSAR Modeling:**  
  Use standardized molecules as input to develop predictive models for biological activity or toxicity.
  
- **Chemical Database Curation:**  
  Improve the quality and uniformity of chemical databases by standardizing molecular representations.
  
- **Workflow Integration:**  
  Integrate the **Molecular Standardizer** with other Orange widgets (e.g., Fingerprint Calculator, PCA, and Clustering widgets) to build end-to-end data mining pipelines.

## Suggested Pipeline with Existing Widgets

1. **Data Import:**  
   - Use the **File** widget or **Data Table** widget to load chemical datasets.
   
2. **Molecular Standardization:**  
   - Connect the imported data to the **Molecular Standardizer** widget to clean and standardize the SMILES strings.
   
3. **Feature Extraction:**  
   - Feed the standardized SMILES to a **Fingerprint Calculator** widget to generate molecular descriptors.
   
4. **Visualization and Analysis:**  
   - Use **PCA**, **Scatter Plot**, or **Hierarchical Clustering** widgets to visualize and cluster the molecular data.
   
5. **Model Building:**  
   - Integrate with classification or regression widgets to build predictive models based on the processed chemical data.
   
6. **Reporting:**  
   - Use the **Report** or **Save Data** widgets to document the standardized data and analysis outcomes.

---

## Conclusion

The **Molecular Standardizer** widget is an essential tool for preparing chemical datasets by ensuring consistent and reproducible molecular representations. Its flexible design allows users to apply a range of standardization operations, making it a valuable preprocessing step for many cheminformatics workflows. Whether you are curating chemical databases or building predictive models, this widget helps guarantee that your molecular data is in the best possible form for analysis.

Feel free to adapt the operations and integrate the widget into larger data mining pipelines to meet your research needs.
