from Orange.widgets.widget import OWWidget, Input, Output
from Orange.widgets import gui
from Orange.widgets.settings import Setting
from Orange.data import Table, Domain, ContinuousVariable, StringVariable
from AnyQt.QtWidgets import (
    QListWidget, QListWidgetItem, QHBoxLayout, QVBoxLayout, QPushButton,
    QLabel, QProgressBar
)
from AnyQt.QtCore import Qt
from rdkit import Chem
from mordred import Calculator, descriptors, Descriptor
import pandas as pd
import numpy as np
import hashlib

# List of descriptor groups (modules)
DESCRIPTOR_GROUPS = [
    "mordred.ABCIndex",
    "mordred.AcidBase",
    "mordred.AdjacencyMatrix",
    "mordred.Aromatic",
    "mordred.AtomCount",
    "mordred.Autocorrelation",
    "mordred.BCUT",
    "mordred.BalabanJ",
    "mordred.BaryszMatrix",
    "mordred.BertzCT",
    "mordred.BondCount",
    "mordred.CPSA",
    "mordred.CarbonTypes",
    "mordred.Chi",
    "mordred.Constitutional",
    "mordred.DetourMatrix",
    "mordred.DistanceMatrix",
    "mordred.EState",
    "mordred.EccentricConnectivityIndex",
    "mordred.ExtendedTopochemicalAtom",
    "mordred.FragmentComplexity",
    "mordred.Framework",
    "mordred.GeometricalIndex",
    "mordred.GravitationalIndex",
    "mordred.HydrogenBond",
    "mordred.InformationContent",
    "mordred.KappaShapeIndex",
    "mordred.Lipinski",
    "mordred.LogS",
    "mordred.McGowanVolume",
    "mordred.MoRSE",
    "mordred.MoeType",
    "mordred.MolecularDistanceEdge",
    "mordred.MolecularId",
    "mordred.MomentOfInertia",
    "mordred.PBF",
    "mordred.PathCount",
    "mordred.Polarizability",
    "mordred.RingCount",
    "mordred.RotatableBond",
    "mordred.SLogP",
    "mordred.TopoPSA",
    "mordred.TopologicalCharge",
    "mordred.TopologicalIndex",
    "mordred.VdwVolumeABC",
    "mordred.VertexAdjacencyInformation",
    "mordred.WalkCount",
    "mordred.Weight",
    "mordred.WienerIndex",
    "mordred.ZagrebIndex",
]

class OWMordredDescriptors(OWWidget):
    name = "Mordred Descriptors (by Group)"
    description = "Compute selected Mordred descriptors for input SMILES. Allows group-based filtering of descriptors."
    icon = "icons/mordred.png"
    category = "Chemoinformatics"
    priority = 25

    class Inputs:
        data = Input("SMILES Table", Table)
    class Outputs:
        data = Output("Descriptors Table", Table)

    # Persistent settings: track which descriptors are selected (by name).
    selected_descriptors = Setting([])

    def __init__(self):
        super().__init__()
        self.data = None
        self.desc_cache = None
        self.cache_key = None
        self.valid_indices = []
        self.invalid_count = 0

        # Info label
        self.infoLabel = gui.label(self.controlArea, self, "Waiting for input...")
        # Progress bar
        self.progressBar = QProgressBar()
        self.progressBar.setMaximum(100)
        self.progressBar.setValue(0)
        self.progressBar.setTextVisible(True)
        self.progressBar.setFormat("%p%")
        self.progressBar.setVisible(False)
        self.controlArea.layout().addWidget(self.progressBar)

        # Build a global list of descriptor objects (2D only)
        calc_all = Calculator(descriptors, ignore_3D=True)
        self.all_desc_objs = list(calc_all.descriptors)
        self.descriptor_map = {}
        for d in self.all_desc_objs:
            mod = d.__module__  # e.g. 'mordred.AcidBase'
            self.descriptor_map[str(d)] = (d, mod)

        # --- Group selection UI ---
        gui.widgetLabel(self.controlArea, "Select descriptor groups:")
        self.groupListWidget = QListWidget()
        self.groupListWidget.setSelectionMode(QListWidget.NoSelection)

        # Define the Lipinski-related groups that should be pre-selected
        lipinski_groups = {
            "mordred.Weight",
            "mordred.RotatableBond",
            "mordred.HydrogenBond",
            "mordred.SLogP"
        }

        for group in DESCRIPTOR_GROUPS:
            item = QListWidgetItem(group)
            item.setFlags(item.flags() | Qt.ItemIsUserCheckable)
            item.setCheckState(Qt.Checked if group in lipinski_groups else Qt.Unchecked)
            self.groupListWidget.addItem(item)

        self.groupListWidget.itemChanged.connect(self._update_available_descriptors)
        self.controlArea.layout().addWidget(self.groupListWidget)

        # Descriptor lists (dual-list UI)
        box = gui.widgetBox(self.controlArea, "Descriptor Selection")
        hbox = QHBoxLayout()

        # Left list: Available descriptors
        vbox_left = QVBoxLayout()
        lbl_avail = QLabel("Available Descriptors")
        lbl_avail.setStyleSheet("font-weight: bold;")
        self.available_list = QListWidget()
        self.available_list.setSelectionMode(QListWidget.ExtendedSelection)
        vbox_left.addWidget(lbl_avail)
        vbox_left.addWidget(self.available_list)

        # Add/Remove buttons
        vbox_btn = QVBoxLayout()
        btn_add = QPushButton("→ Add →")
        btn_remove = QPushButton("← Remove")
        btn_add.clicked.connect(self._on_add_clicked)
        btn_remove.clicked.connect(self._on_remove_clicked)
        vbox_btn.addStretch(1)
        vbox_btn.addWidget(btn_add)
        vbox_btn.addWidget(btn_remove)
        vbox_btn.addStretch(1)

        # Right list: Selected descriptors
        vbox_right = QVBoxLayout()
        lbl_sel = QLabel("Selected Descriptors")
        lbl_sel.setStyleSheet("font-weight: bold;")
        self.selected_list = QListWidget()
        self.selected_list.setSelectionMode(QListWidget.ExtendedSelection)
        self.selected_list.setDragDropMode(QListWidget.InternalMove)
        self.selected_list.model().rowsMoved.connect(self._on_order_changed)
        vbox_right.addWidget(lbl_sel)
        vbox_right.addWidget(self.selected_list)

        hbox.addLayout(vbox_left)
        hbox.addLayout(vbox_btn)
        hbox.addLayout(vbox_right)
        box.layout().addLayout(hbox)

        # Select All / Deselect All buttons
        hbox_btn2 = QHBoxLayout()
        btn_sel_all = QPushButton("Select All")
        btn_sel_none = QPushButton("Deselect All")
        btn_sel_all.clicked.connect(self._on_select_all)
        btn_sel_none.clicked.connect(self._on_select_none)
        hbox_btn2.addWidget(btn_sel_all)
        hbox_btn2.addWidget(btn_sel_none)
        box.layout().addLayout(hbox_btn2)

        # Recompute button
        gui.button(self.controlArea, self, "Recompute descriptors", callback=self.commit)

        # Now populate descriptors: pre-select descriptors from Lipinski-related groups
        self.selected_descriptors = []
        for desc_str, (desc_obj, mod) in self.descriptor_map.items():
            if mod in lipinski_groups:
                self.selected_list.addItem(QListWidgetItem(desc_str))
                self.selected_descriptors.append(desc_str)
        # Then update the available descriptors (which will exclude already selected ones)
        self._update_available_descriptors()

    def _selected_groups(self):
        """Return the list of group names that are currently checked."""
        groups = []
        for i in range(self.groupListWidget.count()):
            item = self.groupListWidget.item(i)
            if item.checkState() == Qt.Checked:
                groups.append(item.text())
        return groups

    def _update_available_descriptors(self):
        """Refresh the 'Available Descriptors' list based on the checked groups."""
        selected_groups = self._selected_groups()
        self.available_list.clear()
        # Gather descriptors that belong to any checked group, but are NOT in self.selected_descriptors
        for desc_str, (obj, mod) in self.descriptor_map.items():
            if any(mod.startswith(g) for g in selected_groups):
                if desc_str not in self.selected_descriptors:
                    self.available_list.addItem(QListWidgetItem(desc_str))

    @Inputs.data
    def set_data(self, data):
        self.data = data
        self.desc_cache = None
        self.cache_key = None
        self.valid_indices = []
        self.invalid_count = 0

        if data is None:
            self.infoLabel.setText("Waiting for input...")
            self.Outputs.data.send(None)
            return

        # Identify SMILES column
        smiles_var = next((v for v in (data.domain.attributes + data.domain.metas)
                           if v.name.lower() == "smiles"), None)
        if smiles_var is None:
            self.infoLabel.setText("No 'SMILES' column found.")
            self.Outputs.data.send(data)
            return

        # Build a cache key
        smiles_list = [str(inst[smiles_var]) for inst in data]
        key = hashlib.md5("".join(smiles_list).encode("utf-8")).hexdigest()
        self.cache_key = key

        # Parse SMILES
        mols = []
        for i, smi in enumerate(smiles_list):
            mol = Chem.MolFromSmiles(smi.strip()) if smi else None
            if mol is None:
                self.invalid_count += 1
            else:
                mols.append(mol)
                self.valid_indices.append(i)
        if self.invalid_count:
            self.warning(f"{self.invalid_count} SMILES entries could not be parsed.")
        else:
            self.warning("")

        # Compute descriptors for selected items
        self._compute_selected_descriptors(mols)
        self._update_output()

    def _compute_selected_descriptors(self, mols):
        if not self.data or not mols:
            self.desc_cache = None
            return
        # Build a calculator with only the selected descriptor objects
        selected = [self.selected_list.item(i).text() for i in range(self.selected_list.count())]
        if not selected:
            self.desc_cache = None
            return
        calc = Calculator([], ignore_3D=True)
        for desc_str in selected:
            if desc_str in self.descriptor_map:
                calc.register(self.descriptor_map[desc_str][0])
        self.progressBar.setVisible(True)
        self.progressBar.setValue(0)
        try:
            df = calc.pandas(mols, nproc=None)
        except Exception as e:
            self.error(f"Error computing descriptors: {e}")
            self.desc_cache = None
            self.progressBar.setVisible(False)
            return
        df.index = self.valid_indices
        full_df = df.reindex(range(len(self.data)))
        self.desc_cache = full_df
        self.progressBar.setValue(100)
        self.progressBar.setVisible(False)

    def _update_output(self):
        if self.data is None or self.desc_cache is None or self.desc_cache.empty:
            self.Outputs.data.send(self.data)
            return
        selected = list(self.desc_cache.columns)
        new_vars = [ContinuousVariable(name) for name in selected]
        new_domain = Domain(new_vars, self.data.domain.class_vars, self.data.domain.metas)
        X = self.desc_cache[selected].values
        out_table = Table(new_domain, X, metas=self.data.metas)
        self.Outputs.data.send(out_table)
        self.infoLabel.setText(f"Computed {len(selected)} descriptors for {len(self.data)} entries.")

    def commit(self):
        if self.data is not None:
            self.desc_cache = None
            smiles_var = next((v for v in (self.data.domain.attributes + self.data.domain.metas)
                               if v.name.lower() == "smiles"), None)
            if smiles_var is None:
                self.infoLabel.setText("No 'SMILES' column found.")
                self.Outputs.data.send(self.data)
                return
            smiles_list = [str(inst[smiles_var]) for inst in self.data]
            mols = []
            self.invalid_count = 0
            self.valid_indices = []
            for i, smi in enumerate(smiles_list):
                mol = Chem.MolFromSmiles(smi.strip()) if smi else None
                if mol is None:
                    self.invalid_count += 1
                else:
                    mols.append(mol)
                    self.valid_indices.append(i)
            if self.invalid_count:
                self.warning(f"{self.invalid_count} SMILES entries could not be parsed.")
            else:
                self.warning("")
            self._compute_selected_descriptors(mols)
            self._update_output()

    # --- Dual-list UI methods ---

    def _on_add_clicked(self):
        items = self.available_list.selectedItems()
        if not items:
            return
        for item in items:
            name = item.text()
            row = self.available_list.row(item)
            self.available_list.takeItem(row)
            self.selected_list.addItem(QListWidgetItem(name))
            if name not in self.selected_descriptors:
                self.selected_descriptors.append(name)

    def _on_remove_clicked(self):
        items = self.selected_list.selectedItems()
        if not items:
            return
        for item in items:
            name = item.text()
            row = self.selected_list.row(item)
            self.selected_list.takeItem(row)
            self.available_list.addItem(QListWidgetItem(name))
            if name in self.selected_descriptors:
                self.selected_descriptors.remove(name)

    def _on_order_changed(self, parent, start, end, destination, dest_row):
        self.selected_descriptors = [self.selected_list.item(i).text() for i in range(self.selected_list.count())]

    def _on_select_all(self):
        count = self.available_list.count()
        for _ in range(count):
            item = self.available_list.takeItem(0)
            self.selected_list.addItem(QListWidgetItem(item.text()))
            if item.text() not in self.selected_descriptors:
                self.selected_descriptors.append(item.text())

    def _on_select_none(self):
        count = self.selected_list.count()
        for _ in range(count):
            item = self.selected_list.takeItem(0)
            self.available_list.addItem(QListWidgetItem(item.text()))
        self.selected_descriptors = []

    def onDeleteWidget(self):
        self.data = None
        super().onDeleteWidget()

if __name__ == "__main__":
    from Orange.widgets.utils.widgetpreview import WidgetPreview
    WidgetPreview(OWMordredDescriptors).run()

