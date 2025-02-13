from Orange.widgets.widget import OWWidget, Input, Output
from Orange.widgets import gui
from Orange.data import Table, Domain, ContinuousVariable
from admet_ai import ADMETModel

class ADMETWidget(OWWidget):
    name = "ADMET Predictor"
    description = "Predicts ADMET properties from SMILES strings."
    icon = "icons/admet.png"  # Ensure you have an appropriate icon file

    class Inputs:
        data = Input("Data", Table)

    class Outputs:
        data = Output("Data", Table)

    def __init__(self):
        super().__init__()
        self.data = None
        self.model = ADMETModel()
        self.infoLabel = gui.label(self.controlArea, self,
                                   "No data on input yet, waiting to get something.")

    @Inputs.data
    def set_data(self, data):
        self.data = data
        if data is not None:
            self.infoLabel.setText(f"Received data with {len(data)} instances.")
            self.process_data()
        else:
            self.infoLabel.setText("No data on input yet, waiting to get something.")

    def process_data(self):
        # Find the 'SMILES' attribute in the data
        smiles_attr = next((var for var in self.data.domain if var.name.lower() == "smiles"), None)
        if smiles_attr is None:
            self.infoLabel.setText("No 'SMILES' column found in the input data.")
            return

        # Extract SMILES strings from the data
        smiles_list = [str(instance[smiles_attr]) for instance in self.data]

        # Predict ADMET properties using the ADMET-AI model
        predictions = self.model.predict(smiles=smiles_list)

        # Round predictions to 4 decimal places
        predictions = predictions.round(4)

        # Gather all original variable names from attributes, class_vars, and metas
        original_names = {var.name for var in self.data.domain.variables}

        # Create new ContinuousVariables for each ADMET property,
        # ensuring each name is unique in the domain.
        new_attrs = []
        for col in predictions.columns:
            new_name = col
            suffix = " (predicted)"
            while new_name in original_names:
                new_name = col + suffix
                suffix += " (predicted)"
            new_var = ContinuousVariable(new_name)
            new_attrs.append(new_var)
            original_names.add(new_name)

        # Create a new domain by appending the new attributes to the original attributes.
        new_domain = Domain(
            self.data.domain.attributes + tuple(new_attrs),
            self.data.domain.class_vars,
            self.data.domain.metas
        )

        # Create a new table with the updated domain
        new_table = Table.from_table(new_domain, self.data)

        # Populate the new table with ADMET predictions
        for i, instance in enumerate(new_table):
            for new_var, col in zip(new_attrs, predictions.columns):
                instance[new_var.name] = predictions.iloc[i][col]

        # Output the new table
        self.Outputs.data.send(new_table)
