# Documentation for ADMET Predictor and ADMET Report Widgets

## 1. Introduction

In the drug discovery process, early evaluation of a compound’s ADMET (Absorption, Distribution, Metabolism, Excretion, and Toxicity) properties is critical for reducing the risk of late-stage failures. The **ADMET Predictor** and **ADMET Report** widgets integrate advanced machine learning methods and rich data visualization to support rapid screening and detailed evaluation of compounds.

Central to these widgets is **ADMET-AI**—a state-of-the-art machine learning framework designed for accurate prediction of ADMET properties. Detailed in a recent publication in *Bioinformatics* (DOI: [10.1093/bioinformatics/btae416](https://academic.oup.com/bioinformatics/article/40/7/btae416/7698030)), ADMET-AI leverages advanced deep learning techniques and curated chemical datasets to deliver robust predictions. The ADMET Predictor widget utilizes ADMET-AI for generating predictions, while the ADMET Report widget presents these predictions together with experimental or calculated descriptors in an interactive report.

---

## 2. Scientific Background and ADMET-AI

### ADMET in Drug Discovery

ADMET properties determine how a compound behaves in the body, affecting its absorption into the bloodstream, its distribution throughout tissues, its metabolic stability, and its potential toxicity. Key descriptors include:
- **Physicochemical Properties:** Such as molecular weight, logP, hydrogen bond donors/acceptors, and topological polar surface area (TPSA).
- **Absorption Metrics:** Including human intestinal absorption, oral bioavailability, and aqueous solubility.
- **Distribution Metrics:** Like blood-brain barrier penetration and plasma protein binding.
- **Metabolism and Excretion:** Involving CYP450 enzyme interactions, clearance rates, and half-life.
- **Toxicity Endpoints:** Covering hERG blocking, mutagenicity, and drug-induced liver injury.

Predicting these properties reliably can help in early compound selection, reduce experimental costs, and streamline the drug development process.

### ADMET-AI Framework

**ADMET-AI** is a machine learning framework that utilizes deep learning algorithms to accurately predict ADMET properties. According to the study published in *Bioinformatics* ([DOI: 10.1093/bioinformatics/btae416](https://academic.oup.com/bioinformatics/article/40/7/btae416/7698030)), the framework:
- **Integrates Multiple Descriptors:** It combines classical molecular descriptors with learned representations to capture both global and local chemical features.
- **Leverages Advanced Deep Learning:** Modern architectures (such as graph neural networks and ensemble learning strategies) are used to model the complex relationships between chemical structure and ADMET endpoints.
- **Validates Rigorously:** The framework has been extensively evaluated on benchmark datasets and shows state-of-the-art performance, demonstrating high correlation with experimental measurements.
- **Facilitates Drug Discovery:** By predicting multiple ADMET endpoints simultaneously, ADMET-AI helps researchers prioritize compounds with favorable pharmacokinetic and toxicity profiles early in the drug development process.

In our widgets, the **ADMET Predictor** is built on ADMET-AI, providing users with cutting-edge predictions. These predictions can be further visualized and analyzed in the **ADMET Report** widget.

---

## 3. ADMET Descriptors and JSON Configuration

### ADMET Descriptors

ADMET descriptors quantify key properties that affect a compound’s behavior in biological systems. They include:

- **Physicochemical Descriptors:**  
  - *Molecular Weight, LogP, Hydrogen Bond Donors/Acceptors, Lipinski Rule of 5, QED, TPSA*  
- **Absorption Descriptors:**  
  - *Human Intestinal Absorption, Oral Bioavailability, Aqueous Solubility, PAMPA Permeability*  
- **Distribution Descriptors:**  
  - *Blood–Brain Barrier Penetration, Plasma Protein Binding, Volume of Distribution*  
- **Metabolism Descriptors:**  
  - *CYP450 Enzyme Inhibition/Substrate Data (e.g., CYP1A2, CYP2C9, CYP2D6, CYP3A4)*  
- **Excretion Descriptors:**  
  - *Drug Clearance, Half Life*  
- **Toxicity Descriptors:**  
  - *hERG Blocking, Clinical Toxicity, Mutagenicity, DILI, Carcinogenicity, and receptor-based endpoints (e.g., Estrogen/Androgen Receptors)*

### JSON Configuration File

Both widgets rely on an external JSON configuration file (e.g., `admet.json`) that defines which descriptors to display and how to categorize them. The JSON file includes the following keys for each descriptor:

- **Category:** Grouping label (e.g., “Physicochemical”, “Absorption”, “Metabolism”, “Toxicity”).
- **Property:** The key name must match the corresponding column in the data table.
- **name:** A human-readable label for the descriptor.
- **ranges:** A dictionary specifying numerical ranges for alert indicators. Typically, three ranges are defined:
  - **range_1:** Favorable (green alert).
  - **range_2:** Cautionary (orange alert).
  - **range_3:** Unfavorable (red alert).
- **percentile:** (Optional) A key representing additional context (e.g., percentile ranking) to be appended to the displayed value.

A sample configuration entry might be:

```json
{
    "Category": "Physicochemical",
    "Property": "molecular_weight",
    "name": "Molecular Weight",
    "ranges": {
        "range_1": {"min": 0.0, "max": 500.0},
        "range_2": {"min": 700.0, "max": 1000.0},
        "range_3": {"min": 1001.0, "max": 10000.0}
    },
    "percentile": "molecular_weight_drugbank_approved_percentile"
}
```

This file must reside in the same directory as the widget scripts so that the widgets can load and use these configuration details to format and display the ADMET properties.

---

## 4. ADMET Predictor Widget

### Overview

The **ADMET Predictor** widget applies the ADMET-AI framework to predict ADMET properties from SMILES strings. It processes an Orange data table, extracts SMILES information, and uses the ADMET-AI model to generate predictions. The resulting data table is enriched with new ADMET descriptors.

### Key Features

- **Input:**  
  - Data table with a “smiles” column.
- **Processing:**  
  - Extraction of SMILES strings.
  - Prediction of multiple ADMET endpoints using ADMET-AI.
  - Rounding of prediction values and appending new continuous variables to the original domain.
- **Output:**  
  - Updated data table with predicted ADMET properties.

### Scientific Rationale

By using ADMET-AI’s advanced machine learning algorithms, the predictor provides high-quality forecasts of ADMET properties. This computational approach enables rapid screening of chemical libraries and helps identify compounds with desirable pharmacokinetic and safety profiles, thereby guiding experimental efforts.

#### Code Snippet (from `admet_predictor.py`)

```python
@Inputs.data
def set_data(self, data):
    self.data = data
    if data is not None:
        self.infoLabel.setText(f"Received data with {len(data)} instances.")
        self.process_data()
    else:
        self.infoLabel.setText("No data on input yet, waiting to get something.")

def process_data(self):
    smiles_attr = next((var for var in self.data.domain if var.name.lower() == "smiles"), None)
    if smiles_attr is None:
        self.infoLabel.setText("No 'SMILES' column found in the input data.")
        return
    smiles_list = [str(instance[smiles_attr]) for instance in self.data]
    predictions = self.model.predict(smiles=smiles_list).round(4)
    # Create new ContinuousVariables for each ADMET property and append them to the domain.
```

---

## 5. ADMET Report Widget

### Overview

The **ADMET Report** widget creates an interactive HTML report that displays each molecule’s structure along with detailed ADMET properties. The report includes:

- **2D Structure Image:** Rendered from the SMILES string using RDKit.
- **Radial (Spider) Plot:** Generated with Matplotlib to visualize key ADMET endpoints (e.g., toxicity, hERG blocking, solubility, bioavailability, BBB penetration).
- **Properties Table:** Displays numerical values, percentiles, and alert indicators based on thresholds defined in the JSON configuration.
- **PDF Export:** An option to save the report as a PDF with CSS formatting to scale each molecule’s content to one A4 page.

### Key Features

- **Input:**  
  - Data table (potentially the output from the ADMET Predictor) with a “smiles” column and ADMET descriptors.
- **Visualization:**  
  - Chemical structures generated by RDKit.
  - Radial plots summarizing ADMET characteristics.
  - Alert indicators (coloured dots) derived from the JSON configuration.
- **Interactivity:**  
  - Dynamic HTML report rendered in a QWebEngineView.
  - Adjustable settings (e.g., molecule image size) for customized display.
- **Export:**  
  - Built-in functionality to export the report as a print-ready PDF.

### Scientific Rationale

Visual representation of ADMET data allows researchers to identify patterns and outliers rapidly. The combination of graphical (radial plots) and tabular information provides both an intuitive overview and a detailed analysis of each compound’s properties, facilitating informed decision-making in the drug development process.

#### Code Snippet (from `admet_report.py`)

```python
@Inputs.orange_data
def set_orange_data(self, data: Table):
    self.orange_data = data
    if data is not None:
        self.info_label.setText(f"Received {len(data)} molecules.")
        self.display_report()
    else:
        self.info_label.setText("No data received.")

def display_report(self):
    # Retrieve the 'smiles' column from the domain
    smiles_col = next((var for var in self.orange_data.domain.variables + self.orange_data.domain.metas
                       if var.name.lower() == "smiles"), None)
    if smiles_col is None:
        self.info_label.setText("No 'SMILES' column found in data.")
        return

    # Build the HTML report, which includes:
    # - Structure images (via RDKit)
    # - Radial plots (via Matplotlib)
    # - Property tables (configured via the JSON file)
    # CSS rules ensure that the report scales correctly for PDF export.
```

---

## 6. Examples of Use in a Workflow

### End-to-End Workflow in Orange

1. **Data Import:**  
   Import a dataset (e.g., CSV file) that contains a “smiles” column along with other ADMET-related descriptors into Orange.

2. **Prediction with ADMET-AI:**  
   - Connect the dataset to the **ADMET Predictor** widget.
   - The widget uses ADMET-AI to predict multiple ADMET endpoints, producing an enriched data table.

3. **Generate and Review Report:**  
   - Connect the enriched data table to the **ADMET Report** widget.
   - The report displays each molecule’s 2D structure, a radial plot summarizing key ADMET properties, and a detailed properties table.
   - Alert indicators (colored dots) based on the JSON configuration highlight favorable and unfavorable ranges.

4. **Exporting the Report:**  
   - Click the “Export as PDF” button to generate a PDF report.
   - The exported PDF is formatted so that each molecule’s report fits neatly on one A4 page.

---

## 7. Requirements and Installation

Ensure you have the following dependencies installed:

- **Python 3.x**
- **Orange3** – for building and running the widgets.
- **RDKit** – for generating chemical structure images.
- **Matplotlib** – for creating radial plots.
- **PyQt5** – with QWebEngine support for HTML rendering.
- **Numpy**

Install dependencies using pip:

```bash
pip install Orange3 rdkit-pypi matplotlib PyQt5 numpy
```

Place the `admet.json` file (which configures the ADMET descriptors and alert ranges) in the same directory as the widget scripts.

---

## 8. Conclusion

The **ADMET Predictor** and **ADMET Report** widgets integrate the advanced ADMET-AI framework into the Orange ecosystem, offering a powerful solution for early-stage drug screening. By combining state-of-the-art machine learning predictions with detailed visual reports, these tools enable researchers to efficiently evaluate the pharmacokinetic and safety profiles of chemical compounds. The JSON configuration ensures that descriptors are presented in a scientifically meaningful manner, supporting both data-driven analysis and effective communication of results.

Happy screening and best of luck in your drug discovery endeavors!

---
