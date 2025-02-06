Below is an example of detailed documentation for the **Mol Sketcher** widget. You can include this documentation in your project’s documentation files (for example, as a README or a dedicated user manual). It covers an introduction, features, installation, configuration, usage suggestions, and troubleshooting tips.

---

# Mol Sketcher Widget Documentation

## Overview

The **Mol Sketcher** widget is an interactive molecular editor designed for the creation of chemical compound databases. It leverages the [JSME](http://peter-ertl.com/jsme/) (JavaScript Molecular Editor) applet to allow users to draw chemical structures in a full-screen browser interface. Using RDKit, the widget calculates various chemical properties (such as molecular weight, LogP, hydrogen bond donors/acceptors, etc.) on the fly. Additionally, it supports the inclusion of user-defined metadata and an auto-incrementing database key (DB key) to uniquely identify compounds.

This widget is particularly useful for:
- Creating databases for chemical compounds.
- Rapidly screening molecular structures.
- Collecting metadata and chemical properties for further data analysis or machine learning tasks.

---

## Features

- **Interactive Molecular Drawing:**  
  Use the embedded JSME molecular editor to draw or import molecules.

- **Automatic Property Calculation:**  
  When RDKit is installed, the widget computes chemical properties (e.g., Molecular Weight, LogP, H-Bond Donors) based on the drawn molecule.

- **Customizable Output Fields:**  
  Configure which properties appear in your database by loading a JSON configuration file. The configuration file supports:
  - **Standard Fields:** Chemical properties and SMILES strings.
  - **Iterative DB Key:** Automatically assign a unique key to each compound. You can define the starting value using the `"initial"` parameter.
  - **User Metadata:** Add custom metadata fields such as Catalog Number, Batch ID, or any other information.

- **Output as Orange Data Table:**  
  Compiled data (both numeric chemical properties and string metadata) are output as an Orange Data Table for integration with the Orange data analysis suite.

- **Clear and Update Functions:**  
  Easily add compounds one by one and clear the database with the “Clear All” option.

---

## Installation and Requirements

### Requirements

- **Python 3.x**  
- **Orange3:**  
  The widget is designed as an Orange widget. Install Orange3 following the official documentation.
- **RDKit (Optional):**  
  For chemical property calculations, install RDKit. If RDKit is not available, the widget will still run but without chemical property computations.
- **PyQt5:**  
  Required for the user interface and web engine support.
- **JSME Files:**  
  Ensure the `jsme` directory (containing the `jsme_panel.html` and `jsme.nocache.js`) is available in the same directory as the widget’s Python script.

### Installation Steps

1. **Clone or Copy the Widget Source:**  
   Place the widget source file (for example, `molsketcher.py`) and the `jsme` directory in your Orange widgets folder or your working directory.

2. **Install Dependencies:**  
   Use pip (or conda) to install the necessary libraries, for example:  
   ```bash
   pip install Orange3 PyQt5 numpy
   ```
   If you plan to compute chemical properties:
   ```bash
   conda install -c rdkit rdkit
   ```

3. **Launch Orange:**  
   Start Orange and add the Mol Sketcher widget to your workflow.

---

## Configuration File

The widget is configured via a JSON file that specifies the fields, user metadata, and optional DB key settings. Below is a sample configuration file:

```json
{
  "dbkey": {
    "name": "dbkey",
    "label": "DB Key",
    "type": "int",
    "initial": 100
  },
  "fields": [
    {"name": "smiles", "label": "SMILES"},
    {"name": "mw", "label": "Molecular Weight", "type": "float"},
    {"name": "logp", "label": "LogP", "type": "float"},
    {"name": "hbd", "label": "H-Bond Donors", "type": "int"},
    {"name": "inchi", "label": "InChI"},
    {"name": "inchikey", "label": "InChI Key"}
  ],
  "user_metadata": [
    {"name": "catalog_number", "label": "Catalog Number"},
    {"name": "batch_id", "label": "Batch ID"},
    {"name": "Name", "label": "Name"}
  ]
}
```

### Explanation of Configuration Keys

- **dbkey:**  
  Optional configuration for an iterative database key.  
  - `name`: The internal name of the key.
  - `label`: How the key will be labeled in the output table.
  - `type`: The type of the key (usually `"int"`).
  - `initial`: The starting number for the DB key counter (default is 1 if omitted).

- **fields:**  
  Define which chemical properties and fields to include. For each field:
  - `name`: The internal identifier. For RDKit computed properties, it should match the keys in the widget’s `PROPERTY_MAP`.
  - `label`: The display name in the output table.
  - `type`: Optionally specify if the field is a numeric type (`"float"` or `"int"`). This determines whether it is treated as an attribute or metadata in the output table.

- **user_metadata:**  
  Define additional fields to capture user-provided metadata (for example, catalog numbers, batch IDs, etc.).  
  - `name`: The internal variable name.
  - `label`: The display label in the output table.

---

## Using the Widget

### Starting the Widget

1. **Load the Widget:**  
   Launch Orange and add the **Mol Sketcher** widget to your workflow.

2. **Load Configuration:**  
   - Click the **Load Config...** button.
   - Select your JSON configuration file.
   - Upon successful load, the widget displays the configuration filename and an informational message.

3. **Enter Metadata:**  
   For every field defined under `"user_metadata"`, an input box appears.  
   - Fill in or update these fields as needed before adding compounds.

### Drawing and Adding Compounds

1. **Draw a Molecule:**  
   In the main area, the full-screen JSME molecular editor is displayed. Draw your molecule or paste a valid SMILES string.

2. **Add the Compound:**  
   - Click the **Add Compound** button.
   - The widget retrieves the SMILES string from the editor.
   - RDKit (if available) computes the specified chemical properties.
   - The widget automatically assigns a DB key (if configured) and collects metadata.
   - The compound is added to an internal data list, and an output table is updated.

3. **Clear All Compounds:**  
   - If you need to reset the table, click the **Clear All** button.
   - The data list is cleared, and the output table is updated accordingly.

### Viewing and Using the Output

- The final output is sent as an Orange **Data Table** to the widget’s output channel.  
- Downstream Orange widgets (for data visualization, analysis, or machine learning) can directly use this table.

---

## Suggestions for Use

- **Rapid Database Construction:**  
  Use the widget in combination with other Orange widgets (e.g., Data Table, Scatter Plot) to quickly build and analyze a compound database.

- **Custom Metadata:**  
  Leverage the `"user_metadata"` fields to track experimental conditions, sample IDs, or supplier information. This can be very useful for data integration and subsequent analysis.

- **Starting DB Key at a Custom Value:**  
  If you already have an existing database, configure the `"initial"` parameter in the `"dbkey"` section of the JSON file to continue the sequence without duplication.

- **Integration with RDKit:**  
  For more advanced chemical property computations, consider extending the `PROPERTY_MAP` with additional RDKit functions or custom calculators.

- **Error Handling:**  
  The widget provides on-screen error messages if the drawn molecule is invalid or if the configuration file cannot be loaded. Check the output messages for troubleshooting tips.

---

## Troubleshooting and FAQs

**Q: The widget does not compute chemical properties.**  
**A:** Ensure that RDKit is installed. Without RDKit, the widget will still add compounds but without calculated properties.

**Q: The configuration file does not load.**  
**A:** Verify that the JSON file is well-formed and adheres to the expected structure. Check file paths and ensure the file has the correct permissions.

**Q: The DB key does not start at the expected value.**  
**A:** Confirm that the `"dbkey"` section in the JSON file includes the `"initial"` parameter. If omitted, the counter defaults to 1.

**Q: I am not seeing the metadata input fields.**  
**A:** Check your JSON configuration under `"user_metadata"` to ensure that it is correctly defined and that the widget has reloaded the configuration after any changes.

---

## Conclusion

The **Mol Sketcher** widget is a versatile tool for chemical database creation. Its integration with JSME and RDKit allows users to rapidly prototype and build compound libraries with both calculated properties and custom metadata. By following this documentation and experimenting with different configurations, you can tailor the widget to your specific research or industrial needs.

For further customization or to report issues, please consult the project’s repository or contact the development team.

--- 

This documentation should help users understand how to install, configure, and make the best use of the Mol Sketcher widget in their chemical data workflows.
