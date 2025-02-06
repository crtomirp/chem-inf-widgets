# chem-inf-widgets
Chem-Inf-Widgets for Orange3: A set of custom widgets for Orange3 tailored for chemoinformatics, enabling seamless molecular structure visualization, property calculations, and data analysis. Simplify workflows for cheminformatics and drug discovery with interactive tools.

Below is an example of a `README.md` file for the **chem_inf_widgets** repository. You can adjust and expand upon the details as needed.

# chem_inf_widgets

![Build Status](https://img.shields.io/badge/build-passing-brightgreen)
![License](https://img.shields.io/badge/license-MIT-blue)


## Installation

### Requirements

- **Python 3.x**
- **Orange3** (for integration into Orange workflows)  
  Installation instructions: [Orange Download](https://orange.biolab.si/download/)
- **PyQt5** (for GUI components)
- **NumPy** (for numerical operations)
- **Pandas** (for data manipulation in some widgets)
- **RDKit** (for cheminformatics calculations; optional but recommended)  
  Installation via conda:  
  ```bash
  conda install -c rdkit rdkit
  ```
- **Requests** (for ChEMBL API queries)

### Setup

1. **Clone the Repository:**

   ```bash
   git clone https://github.com/crtomirp/chem-inf-widgets.git
   cd chem-inf-widgets
   ```

2. **Install Dependencies:**

   Install the required packages. For example, using pip:
   ```bash
   pip install Orange3 PyQt5 numpy pandas requests
   ```
   If you plan to use RDKit functionality, install RDKit (preferably via conda):
   ```bash
   conda install -c rdkit rdkit
   ```

3. **JSME Files:**

   Ensure that the `jsme` directory (which contains `jsme_panel.html` and `jsme.nocache.js`) is in the same directory as the widgets file. This directory is required for the MolSketcher and Substructure Search widgets.

## Usage

Each widget can be used either within Orange (if you install this repository as an Orange add-on) or in standalone mode using the provided preview code. Below is a brief guide to using each widget:

### 1. DB MolSketcher

- **Purpose:** Draw molecules interactively, compute chemical properties using RDKit, and assign an auto-incrementing DB key.
- **Configuration:** Load a JSON configuration file (via the **Load Config...** button) to define fields, metadata, and optional DB key settings.
- **Usage:** Draw a molecule in the embedded JSME panel and click **Add Compound**. The compound is added to an internal list and output as an Orange Data Table.

### 2. ChEMBL Bioactivity Retriever

- **Purpose:** Fetch bioactivity data from ChEMBL given a target ID.
- **Usage:** Enter a valid ChEMBL target ID (e.g., `CHEMBL2095150`) and click **Fetch Data**. The widget processes the data (including property calculations) and outputs an Orange Data Table.

### 3. Drug Filter

- **Purpose:** Filter molecules based on drug-likeness criteria, including Lipinski’s, Veber’s, and PAINS alerts.
- **Usage:** Choose a filtering rule, selection mode, and optionally enable PAINS highlighting. Click **Filter Molecules** to process the input table and output the filtered compounds.

### 4. Fingerprint Calculator

- **Purpose:** Compute molecular fingerprints using different methods and provide visualization options.
- **Usage:** Configure the fingerprint type (e.g., Morgan, RDKit, MACCS Keys, etc.), bit size, and (for Morgan) radius. Click **Compute Fingerprint** to calculate fingerprints. Use **Show Histogram** or **Show PCA Projection** to visualize the results.

### 5. MACCS Key Generator

- **Purpose:** Convert SMILES strings to MACCS keys.
- **Usage:** Provide an Orange Data Table containing SMILES data. The widget converts each SMILES to a MACCS key vector and outputs the results as a new table.

### 6. Molecular Standardization

- **Purpose:** Standardize input SMILES strings by applying a series of operations (cleanup, normalization, etc.).
- **Usage:** Check the desired standardization operations and click **Apply Standardization**. The widget outputs a new table with standardized SMILES and change logs.

### 7. Molecular Viewer

- **Purpose:** Display molecules in a customizable grid layout with optional substructure highlighting and property display.
- **Usage:** Provide an input table (e.g., filtered compounds). Adjust grid settings (image size, number of columns) and select which properties to display. The viewer will render molecule images accordingly.

### 8. SDF Reader

- **Purpose:** Analyze an SDF file, select desired properties, and output a data table.
- **Usage:** Browse for an SDF file, analyze the available properties, select the ones to include, and click **Read File** to generate an Orange Data Table.

### 9. Substructure Search

- **Purpose:** Search compounds in an input table based on substructure, superstructure, similarity, or exact match.
- **Usage:** Draw or enter a SMILES/SMARTS query in the provided field, choose the search type, and click **Apply Search**. The widget filters the input data and outputs a table of matching compounds with optional highlighting.

## Contributing

Contributions to **chem_inf_widgets** are welcome! If you wish to contribute:

1. Fork the repository.
2. Create a new branch for your feature or bug fix.
3. Submit a pull request with a detailed description of your changes.

For major changes, please open an issue first to discuss your ideas.

## License

This project is licensed under the [MIT License](LICENSE).

## Acknowledgments

- **JSME:** [JSME Molecular Editor](http://peter-ertl.com/jsme/)
- **RDKit:** [RDKit](https://www.rdkit.org/)
- **Orange:** [Orange Data Mining](https://orange.biolab.si/)
- **ChEMBL:** [ChEMBL Database](https://www.ebi.ac.uk/chembl/)

---

Happy cheminformatics!

```

---

This `README.md` file provides an overview of the repository, its features, installation instructions, usage examples for each widget, and contribution guidelines. Adjust the content as necessary to reflect the current state and goals of your project.
