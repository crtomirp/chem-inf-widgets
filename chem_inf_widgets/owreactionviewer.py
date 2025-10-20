# owreactionviewer_simple.py — Orange3 widget: Simple RDKit Reaction Viewer
#
# Purpose
#   Minimal, fast reaction viewer:
#   - One input Table
#   - Grid of rendered reaction images (PIL→QPixmap; SVG fallback)
#   - Pick reaction column OR compose from reactants+products
#   - Controls for image size, columns, row limit
#   - NEW: captions under each reaction, and "Export All…" to PNG/SVG
#
# Install
#   Place this file in: orangecontrib/chem/widgets/owreactionviewer_simple.py
#   Add to orangecontrib/chem/widgets/__init__.py:
#       def __get_widgets__():
#           from .owreactionviewer_simple import OWReactionViewerSimple
#           return [OWReactionViewerSimple]
#   Ensure an icon exists at orangecontrib/chem/widgets/icons/molviewer.png
#   or switch the 'icon' value below to an existing icon in your repo.

from __future__ import annotations

import io
import base64
from typing import Optional, List

from Orange.data import Table, StringVariable
from Orange.widgets import gui
from Orange.widgets.settings import Setting
from Orange.widgets.widget import OWWidget, Input

from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import rdChemReactions
from rdkit.Chem.Draw import rdMolDraw2D

from PyQt5.QtCore import Qt
from PyQt5.QtGui import QImage, QPixmap
from PyQt5.QtWidgets import (
    QLabel, QVBoxLayout, QGridLayout, QWidget, QScrollArea, QSpinBox,
    QComboBox, QPushButton, QCheckBox, QFileDialog, QLineEdit
)


class OWReactionViewerSimple(OWWidget):
    name = "Reaction Viewer"
    description = "Display RDKit reaction depictions from a table."
    category = "Chemoinformatics"
    icon = "icons/reactionviewer.png"  # adjust if needed
    priority = 51

    class Inputs:
        data = Input("Data", Table, default=True)

    # Settings
    img_size: int = Setting(320)
    max_columns: int = Setting(4)
    max_rows: int = Setting(100)
    use_compose: bool = Setting(False)  # False: reaction column; True: compose reactants+products
    rxn_col_name: str = Setting("")
    reactants_col: str = Setting("reactants")
    products_col: str = Setting("SMILES")

    # NEW
    show_captions: bool = Setting(True)
    caption_col_name: str = Setting("reaction_name")
    export_format: str = Setting("PNG")  # PNG or SVG
    export_prefix: str = Setting("rxn_")

    def __init__(self):
        super().__init__()
        self.data: Optional[Table] = None
        self._all_string_cols: List[str] = []

        # ---- Main scrollable grid ----
        self.scroll_area = QScrollArea(); self.scroll_area.setWidgetResizable(True)
        self.scroll_widget = QWidget(); self.grid = QGridLayout(self.scroll_widget)
        self.scroll_area.setWidget(self.scroll_widget)
        self.mainArea.layout().addWidget(self.scroll_area)

        # ---- Controls ----
        gui.label(self.controlArea, self, "Reaction Source")
        self.compose_toggle = QCheckBox("Compose from reactants + products")
        self.compose_toggle.setChecked(self.use_compose)
        self.compose_toggle.stateChanged.connect(self._on_compose_changed)
        gui.widgetBox(self.controlArea, orientation="horizontal").layout().addWidget(self.compose_toggle)

        self.rxn_col_box = QComboBox(); gui.widgetBox(self.controlArea, orientation="horizontal").layout().addWidget(self.rxn_col_box)
        self.react_col_box = QComboBox(); gui.widgetBox(self.controlArea, orientation="horizontal").layout().addWidget(self.react_col_box)
        self.prod_col_box = QComboBox(); gui.widgetBox(self.controlArea, orientation="horizontal").layout().addWidget(self.prod_col_box)

        gui.label(self.controlArea, self, "Grid Settings")
        self.size_selector = QSpinBox(); self.size_selector.setRange(120, 600); self.size_selector.setValue(self.img_size)
        self.size_selector.valueChanged.connect(self._update_size)
        gui.widgetBox(self.controlArea, orientation="horizontal").layout().addWidget(self.size_selector)

        self.col_selector = QSpinBox(); self.col_selector.setRange(1, 12); self.col_selector.setValue(self.max_columns)
        self.col_selector.valueChanged.connect(self._update_cols)
        gui.widgetBox(self.controlArea, orientation="horizontal").layout().addWidget(self.col_selector)

        self.row_limit = QSpinBox(); self.row_limit.setRange(1, 5000); self.row_limit.setValue(self.max_rows)
        self.row_limit.valueChanged.connect(self._update_rows)
        gui.widgetBox(self.controlArea, orientation="horizontal").layout().addWidget(self.row_limit)

        gui.label(self.controlArea, self, "Captions")
        self.caption_toggle = QCheckBox("Show captions")
        self.caption_toggle.setChecked(self.show_captions)
        self.caption_toggle.stateChanged.connect(self._on_captions_changed)
        gui.widgetBox(self.controlArea, orientation="horizontal").layout().addWidget(self.caption_toggle)
        self.caption_col_box = QComboBox(); gui.widgetBox(self.controlArea, orientation="horizontal").layout().addWidget(self.caption_col_box)

        gui.label(self.controlArea, self, "Export")
        self.format_box = QComboBox(); self.format_box.addItems(["PNG", "SVG"])
        self.format_box.setCurrentText(self.export_format)
        gui.widgetBox(self.controlArea, orientation="horizontal").layout().addWidget(self.format_box)
        self.prefix_edit = QLineEdit(self.export_prefix)
        self.prefix_edit.setPlaceholderText("file prefix, e.g., rxn_")
        gui.widgetBox(self.controlArea, orientation="horizontal").layout().addWidget(self.prefix_edit)

        self.btn_render = QPushButton("Render")
        self.btn_render.clicked.connect(self._render)
        gui.widgetBox(self.controlArea, orientation="horizontal").layout().addWidget(self.btn_render)

        self.btn_export_all = QPushButton("Export All…")
        self.btn_export_all.clicked.connect(self._export_all)
        gui.widgetBox(self.controlArea, orientation="horizontal").layout().addWidget(self.btn_export_all)

        self.info_label = gui.label(self.controlArea, self, "Awaiting input data…")

    # ---------------- Inputs ----------------
    @Inputs.data
    def set_data(self, data: Optional[Table]):
        self.data = data
        if data is None or len(data) == 0:
            self.info_label.setText("No data.")
            self._clear_grid()
            return
        # Collect string columns
        self._all_string_cols = []
        for var in list(data.domain) + list(data.domain.metas):
            if isinstance(var, StringVariable):
                self._all_string_cols.append(var.name)
        # populate combos
        self.rxn_col_box.clear(); self.rxn_col_box.addItems(self._all_string_cols)
        self.react_col_box.clear(); self.react_col_box.addItems(self._all_string_cols)
        self.prod_col_box.clear(); self.prod_col_box.addItems(self._all_string_cols)
        # defaults + caption list
        self._select_default_columns()
        self.info_label.setText(f"Rows: {len(data)} | strings: {len(self._all_string_cols)}")
        self._render()

    # ---------------- Internal helpers ----------------
    def _select_default_columns(self):
        def pick(pref: List[str]) -> str:
            for p in pref:
                for c in self._all_string_cols:
                    if c.lower() == p.lower():
                        return c
            return self._all_string_cols[0] if self._all_string_cols else ""
        if not self.rxn_col_name:
            self.rxn_col_name = pick(["rxn_mapped", "smirks", "reaction", "rxn"])
        self.reactants_col = pick(["reactants", "reagents", "substrates"]) or self.reactants_col
        self.products_col = pick(["SMILES", "product", "products"]) or self.products_col
        self.rxn_col_box.setCurrentText(self.rxn_col_name)
        self.react_col_box.setCurrentText(self.reactants_col)
        self.prod_col_box.setCurrentText(self.products_col)
        # captions: include special choices
        self.caption_col_box.clear()
        items = ["(none)", "(rxn-string)"] + self._all_string_cols
        self.caption_col_box.addItems(items)
        if self.caption_col_name and self.caption_col_name in self._all_string_cols:
            self.caption_col_box.setCurrentText(self.caption_col_name)
        else:
            self.caption_col_box.setCurrentText("reaction_name" if "reaction_name" in self._all_string_cols else "(rxn-string)")

    def _on_compose_changed(self, _):
        self.use_compose = self.compose_toggle.isChecked()
        self._render()

    def _on_captions_changed(self, _):
        self.show_captions = self.caption_toggle.isChecked()
        self.caption_col_name = self.caption_col_box.currentText().strip()
        self._render()

    def _clear_grid(self):
        for i in reversed(range(self.grid.count())):
            w = self.grid.itemAt(i).widget()
            if w is not None:
                self.grid.removeWidget(w)
                w.deleteLater()

    def _get_var(self, name: str):
        if self.data is None or not name:
            return None
        dom = self.data.domain
        for var in list(dom) + list(dom.metas):
            if var.name == name:
                return var
        return None

    def _compose_rxn_str(self, row_idx: int) -> Optional[str]:
        """Build 'A.B>>P.Q' from reactants+products (robust for metas)."""
        if self.data is None or row_idx < 0 or row_idx >= len(self.data):
            return None
        rvar = self._get_var(self.react_col_box.currentText().strip())
        pvar = self._get_var(self.prod_col_box.currentText().strip())
        if rvar is None or pvar is None:
            return None
        try:
            rtxt = str(self.data[row_idx][rvar]).strip()
            ptxt = str(self.data[row_idx][pvar]).strip()
        except Exception:
            return None
        if not rtxt or not ptxt or rtxt == '?' or ptxt == '?':
            return None
        reactants = [x.strip() for x in rtxt.replace(".", "+").split("+") if x.strip() and x.strip() != "?"]
        products = [x.strip() for x in ptxt.split(".") if x.strip() and x.strip() != "?"]
        if not reactants or not products:
            return None
        return ".".join(reactants) + ">>" + ".".join(products)

    def _pick_rxn_col_str(self, row_idx: int) -> Optional[str]:
        if self.data is None or row_idx < 0 or row_idx >= len(self.data):
            return None
        var = self._get_var(self.rxn_col_box.currentText().strip())
        if var is None:
            return None
        try:
            s = str(self.data[row_idx][var]).strip()
        except Exception:
            return None
        return s or None

    def _parse_reaction(self, rxn_str: str) -> Optional[rdChemReactions.ChemicalReaction]:
        if not rxn_str:
            return None
        # Prefer reaction SMILES; fallback to SMARTS if needed.
        try:
            rxn = rdChemReactions.ReactionFromSmarts(rxn_str, useSmiles=True)
            if rxn is None:
                raise ValueError("None")
            rxn.Initialize()
            return rxn
        except Exception:
            try:
                rxn = rdChemReactions.ReactionFromSmarts(rxn_str, useSmiles=False)
                if rxn is None:
                    return None
                rxn.Initialize()
                return rxn
            except Exception:
                return None

    def _caption_for_row(self, row_idx: int, rxn_str_fallback: Optional[str]) -> Optional[str]:
        if not self.show_captions:
            return None
        choice = self.caption_col_box.currentText().strip()
        if choice == "(none)":
            return None
        if choice == "(rxn-string)":
            return rxn_str_fallback
        var = self._get_var(choice)
        if var is None or self.data is None:
            return rxn_str_fallback
        try:
            s = str(self.data[row_idx][var]).strip()
            return s or rxn_str_fallback
        except Exception:
            return rxn_str_fallback

    def _reaction_to_widget(self, rxn: rdChemReactions.ChemicalReaction, caption: Optional[str]) -> Optional[QWidget]:
        """Render to QPixmap via PIL when available; else fall back to SVG embedded in QLabel.
        Optionally add a caption under the image.
        """
        container = QWidget(); layout = QVBoxLayout(container); layout.setContentsMargins(2, 2, 2, 2)
        # Try PIL/PNG path
        try:
            img = Draw.ReactionToImage(rxn, subImgSize=(self.img_size, int(self.img_size * 0.6)))
            buf = io.BytesIO(); img.save(buf, format="PNG")
            qimg = QImage.fromData(buf.getvalue())
            pm = QPixmap.fromImage(qimg)
            lab = QLabel(); lab.setPixmap(pm); lab.setAlignment(Qt.AlignCenter)
            lab.setStyleSheet("border:1px solid #ddd; padding:4px; margin:4px; background:#fff;")
            layout.addWidget(lab)
        except Exception:
            # SVG fallback
            try:
                w, h = int(self.img_size), int(self.img_size * 0.6)
                d = rdMolDraw2D.MolDraw2DSVG(w, h)
                d.DrawReaction(rxn); d.FinishDrawing()
                svg = d.GetDrawingText()
                enc = base64.b64encode(svg.encode("utf-8")).decode("ascii")
                uri = f"data:image/svg+xml;base64,{enc}"
                html_img = f"<img src='{uri}' width='{w}' height='{h}'/>"
                lab = QLabel(); lab.setTextFormat(Qt.RichText)
                lab.setTextInteractionFlags(Qt.TextBrowserInteraction)
                lab.setAlignment(Qt.AlignCenter)
                lab.setStyleSheet("border:1px solid #ddd; padding:4px; margin:4px; background:#fff;")
                lab.setText(html_img)
                layout.addWidget(lab)
            except Exception:
                return None
        if self.show_captions and caption:
            cap = QLabel(caption)
            cap.setAlignment(Qt.AlignCenter)
            cap.setStyleSheet("color:#444; font: 11px 'Consolas', 'Menlo', monospace;")
            layout.addWidget(cap)
        return container

    def _safe_slug(self, text: str) -> str:
        import re
        s = re.sub(r"\s+", "_", text)
        s = re.sub(r"[^A-Za-z0-9_.-]", "", s)
        return s[:80] or "rxn"

    # ---------------- Export ----------------
    def _export_all(self):
        if self.data is None or len(self.data) == 0:
            return
        folder = QFileDialog.getExistingDirectory(self, "Export all reactions to…")
        if not folder:
            return
        fmt = self.format_box.currentText().strip().upper()
        self.export_format = fmt
        self.export_prefix = self.prefix_edit.text().strip() or "rxn_"
        count = 0
        svg_fallback = 0
        for i in range(len(self.data)):
            rxn_str = self._compose_rxn_str(i) if self.use_compose else self._pick_rxn_col_str(i)
            if not rxn_str:
                rxn_str = self._pick_rxn_col_str(i) or self._compose_rxn_str(i)
            if not rxn_str:
                continue
            rxn = self._parse_reaction(rxn_str)
            if rxn is None:
                continue
            caption = self._caption_for_row(i, rxn_str) or "rxn"
            name = f"{self.export_prefix}{i:04d}_" + self._safe_slug(caption)
            try:
                if fmt == "PNG":
                    try:
                        img = Draw.ReactionToImage(rxn, subImgSize=(self.img_size, int(self.img_size * 0.6)))
                        path = f"{folder}/{name}.png"
                        img.save(path)
                    except Exception:
                        # pillow not available — fallback to SVG for this row
                        w, h = int(self.img_size), int(self.img_size * 0.6)
                        d = rdMolDraw2D.MolDraw2DSVG(w, h)
                        d.DrawReaction(rxn); d.FinishDrawing()
                        with open(f"{folder}/{name}.svg", "w", encoding="utf-8") as f:
                            f.write(d.GetDrawingText())
                        svg_fallback += 1
                else:  # SVG
                    w, h = int(self.img_size), int(self.img_size * 0.6)
                    d = rdMolDraw2D.MolDraw2DSVG(w, h)
                    d.DrawReaction(rxn); d.FinishDrawing()
                    with open(f"{folder}/{name}.svg", "w", encoding="utf-8") as f:
                        f.write(d.GetDrawingText())
                count += 1
            except Exception:
                continue
        msg = f"Exported {count} reactions to '{folder}' as {fmt}."
        if svg_fallback and fmt == "PNG":
            msg += f" ({svg_fallback} fell back to SVG due to missing Pillow.)"
        self.info_label.setText(msg)

    # ---------------- Render ----------------
    def _update_size(self, v): self.img_size = int(v); self._render()
    def _update_cols(self, v): self.max_columns = int(v); self._render()
    def _update_rows(self, v): self.max_rows = int(v); self._render()

    def _render(self):
        self._clear_grid()
        if self.data is None or len(self.data) == 0:
            self.info_label.setText("No data.")
            return
        rows = min(self.max_rows, len(self.data))
        row, col = 0, 0
        made = 0
        for i in range(rows):
            rxn_str = self._compose_rxn_str(i) if self.use_compose else self._pick_rxn_col_str(i)
            if not rxn_str:
                rxn_str = self._pick_rxn_col_str(i) or self._compose_rxn_str(i)
            if not rxn_str:
                continue
            rxn = self._parse_reaction(rxn_str)
            if rxn is None:
                continue
            caption = self._caption_for_row(i, rxn_str)
            w = self._reaction_to_widget(rxn, caption)
            if w is None:
                continue
            self.grid.addWidget(w, row, col)
            col += 1
            if col >= self.max_columns:
                col = 0; row += 1
            made += 1
        if made == 0:
            self.info_label.setText("No reactions could be rendered. Check column selection or install Pillow (pip install pillow).")
        else:
            self.info_label.setText(f"Rendered {made} / {rows} rows.")


if __name__ == "__main__":
    from Orange.widgets.utils.widgetpreview import WidgetPreview
    WidgetPreview(OWReactionViewerSimple).run()
