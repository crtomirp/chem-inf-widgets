# chem_inf_widgets

![Build Status](https://img.shields.io/badge/build-passing-brightgreen)
![License](https://img.shields.io/badge/license-MIT-blue)

**chem_inf_widgets** is a collection of chemical informatics widgets developed for [Orange](https://orange.biolab.si/) and other Python-based data science workflows. 

## Features

- **DB MolSketcher**  
  A JSME-based molecular editor that allows users to draw molecules and automatically compute chemical properties (e.g., molecular weight, LogP, H-bond donors/acceptors, TPSA). It supports customizable JSON configuration files that can include iterative DB key generation with a configurable starting value and user-defined metadata.

- **ChEMBL Bioactivity Retriever**  
  Retrieves bioactivity data from the ChEMBL database for a given target ID. It processes IC50 values, computes drug properties, and outputs the results as an Orange Data Table for further analysis.

- **Drug Filter**  
  Filters molecules based on drug-likeness criteria (Lipinski’s Rule of Five, Veber’s Rule, PAINS alerts, and a composite drug score). Users can choose the filtering rule and select whether to forward all molecules, only those that meet the criteria, or those that fail.

- **Fingerprint Calculator**  
  Computes various types of molecular fingerprints (Morgan, RDKit, MACCS Keys, Atom Pair, Topological Torsion, and Avalon) using RDKit. It also provides visualization features such as histograms and PCA projections of the fingerprint data.

- **MACCS Key Generator**  
  Converts SMILES strings to MACCS keys and outputs the result as an Orange Data Table.

- **Molecular Standardization**  
  Reads SMILES strings from an input Orange table and applies a series of standardization operations (e.g., cleanup, normalization, metal disconnection, largest fragment chooser, reionization, uncharging, and tautomer enumeration) to output standardized SMILES.

- **Molecular Viewer**  
  Displays molecules in a customizable grid layout with optional substructure highlighting and property display based on user selection.

- **SDF Reader**  
  Analyzes SDF files, allows the user to select specific molecular properties, and outputs the selected data as an Orange Data Table.

- **Substructure Search**  
  Performs compound searches on an input data table by matching substructures, superstructures, computing similarity scores, or performing exact matches. It also supports interactive drawing of queries using a JSME molecular editor.

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

### Virtual Environment Setup & Installation

1. **Create and Activate a Virtual Environment:**

   On Unix/Mac:
   ```bash
   python3 -m venv venv
   source venv/bin/activate
   ```
   
   On Windows:
   ```bash
   python -m venv venv
   venv\Scripts\activate
   ```

2. **Clone the Repository:**

   ```bash
   git clone https://github.com/crtomirp/chem-inf-widgets.git
   cd chem-inf-widgets
   ```

3. **Install the Package in Editable Mode:**

   ```bash
   pip install -e .
   ```

4. **Install Remaining Dependencies:**

   If not already installed, run:
   ```bash
   pip install Orange3 PyQt5 numpy pandas requests
   ```
   
   (For RDKit, it is recommended to install via conda as noted above.)

5. **Launch Orange:**

   Start Orange by running:
   ```bash
   python -m Orange.canvas
   ```
   *(Note: If your Orange installation uses a different module name such as `Orange.canvas`, adjust the command accordingly.)*

## Installation with CONDA

# chem-inf-widgets — Conda Installation Guide (Python 3.12)

These instructions set up a full working environment (use + dev) from the repo’s `environment.yml`.

---

## 1) Prerequisites

* **Conda/Mamba:** Install one of:

  * [Miniconda](https://docs.conda.io/en/latest/miniconda.html) or Anaconda
  * **Recommended:** [Mambaforge](https://github.com/conda-forge/miniforge#mambaforge) (includes `mamba` and conda-forge by default; solves faster)
* **Git** to clone the repository

> On Windows, run commands in **Anaconda Prompt** or **PowerShell**.

---

## 2) Get the Source

```bash
# choose a location and clone
git clone https://github.com/crtomirp/chem-inf-widgets.git
cd chem-inf-widgets
```

---

## 3) Create the Conda Environment

> The repo provides `environment.yml` (Python 3.12; includes RDKit, Jupyter, ipywidgets, viz tools, and dev tools).

**Using mamba (fastest, if available):**

```bash
mamba env create -f environment.yml
```

**Using conda:**

```bash
conda env create -f environment.yml
```

> ⚠️ Do **not** pass `-c` or `--override-channels` to `conda env create`; channels are taken from the YAML.

Activate it:

```bash
conda activate chem-inf-widgets
```

---

## 4) Editable Install (dev mode)

The `environment.yml` already performs:

```text
pip install -e .
```

so local code changes are immediately reflected. Nothing extra to do.

---

## 5) Jupyter Setup (optional but recommended)

Register the kernel so notebooks can select this env:

```bash
python -m ipykernel install --user --name chem-inf-widgets --display-name "Python (chem-inf-widgets)"
```

Launch:

```bash
jupyter lab
# or
jupyter notebook
```

* **ipywidgets:** Works out-of-the-box in JupyterLab ≥3.
* Classic Notebook enables widgets automatically when installed via conda.

---

## 6) Verify the Installation

```bash
python - <<'PY'
import rdkit, ipywidgets, pandas as pd
print("RDKit:", rdkit.__version__)
print("ipywidgets:", ipywidgets.__version__)
print("Pandas:", pd.__version__)
print("OK")
PY
```

Open a sample notebook or run a minimal widget cell to confirm rendering.

---

## 7) Updating / Recreating

Update env after pulling changes to `environment.yml`:

```bash
conda env update -n chem-inf-widgets -f environment.yml
```

Remove env if you need a clean rebuild:

```bash
conda deactivate
conda env remove -n chem-inf-widgets
```

---

## 8) Troubleshooting

* **Solver errors on Windows (Unsatisfiable specs):**

  * Prefer **mamba** (`mamba env create -f environment.yml`).
  * Ensure you’re in the repo root (so `-e .` finds `setup.py/pyproject.toml`).
  * The YAML uses `python-build` (not `build`) to avoid Python pin conflicts.
* **Widgets don’t display:**

  * Use JupyterLab ≥3 or Classic Notebook (already covered by `ipywidgets`).
  * Restart the kernel/browser tab after first install.
* **Slow solves / channel conflicts:**

  * Use conda-forge as primary; set strict priority:

    ```bash
    conda config --add channels conda-forge
    conda config --set channel_priority strict
    ```
* **Optional heavy packages (e.g., `nodejs`, `plotly`, `py3dmol`)**:

  * If solving fails, temporarily remove them from `environment.yml`, create the env, then add back one-by-one:

    ```bash
    conda install -c conda-forge nodejs
    conda install -c conda-forge plotly
    conda install -c conda-forge py3dmol
    ```

---

## 9) Minimal user-only env (optional)

If you later want a slim, runtime-only environment (no dev tools), say the word and I’ll provide `environment-user.yml`.


## Usage

Each widget can be used either within Orange (after installing this repository as an add-on) or in standalone mode using the provided preview code. Below is a brief guide to using each widget:

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



Happy cheminformatics!


This version of the `README.md` explains how to set up a virtual environment, install the package in editable mode with `pip install -e .`, and launch Orange using the provided command. Adjust any details as necessary to suit your project's specifics.
