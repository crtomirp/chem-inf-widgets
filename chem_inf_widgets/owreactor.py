# owreactor.py — Orange3 widget adapted from CLI Reactor
# Features preserved:
# - Per-rule weights (probabilistic selection)
# - Consume reactants / expand pool across steps
# - Multi-file SMIRKS catalogs (3 inputs) with de-duplication
# - Mapped reaction string (robust across RDKit versions)
# - Per-rule max trials when sampling combos
# - Pool size cap with policy (fifo / lifo / random)
# - Preview log panel (last run)
#
# Drop this file at: orangecontrib/chem/widgets/owreactor.py
# Ensure you expose it via __get_widgets__ in orangecontrib/chem/widgets/__init__.py
# and an entry point in your pyproject.toml.

from __future__ import annotations

import os
import itertools
import random
from dataclasses import dataclass
from typing import List, Tuple, Optional, Dict, Any, Sequence

from rdkit import Chem
from rdkit.Chem import rdChemReactions

from AnyQt.QtGui import QFont
from AnyQt.QtWidgets import (
    QFormLayout, QSpinBox, QCheckBox, QLineEdit, QPushButton, QHBoxLayout,
    QPlainTextEdit, QWidget, QVBoxLayout, QLabel, QComboBox
)

from Orange.data import Table, Domain, StringVariable, ContinuousVariable
from Orange.widgets.widget import OWWidget, Input, Output
from Orange.widgets.settings import Setting

import numpy as np


# ------------------------- Core logic (UI-agnostic) ------------------------- #

@dataclass
class ReactionRule:
    name: str
    smirks: str
    rxn: rdChemReactions.ChemicalReaction
    weight: float = 1.0

    @property
    def n_reactants(self) -> int:
        return self.rxn.GetNumReactantTemplates()

    @staticmethod
    def from_row(name: str, smirks: str, weight: Optional[float] = None) -> "ReactionRule":
        rxn = rdChemReactions.ReactionFromSmarts(smirks, useSmiles=True)
        if rxn is None:
            raise ValueError(f"Invalid SMIRKS/SMILES reaction: {smirks}")
        rxn.Initialize()
        w = 1.0 if weight is None else max(0.0, float(weight))
        return ReactionRule(name=name or smirks, smirks=smirks, rxn=rxn, weight=w)


class ReactorEngine:
    def __init__(self,
                 smiles: List[str],
                 rules: List[ReactionRule],
                 seed: Optional[int] = None,
                 sanitize_products: bool = True,
                 allow_self_react: bool = False,
                 unique_products: bool = True,
                 consume_reactants: bool = False,
                 expand_pool: bool = False,
                 max_pool_size: int = 0,
                 pool_policy: str = 'fifo',
                 per_rule_max_trials: int = 2000):
        self.rng = random.Random(seed)
        dedup: Dict[str, ReactionRule] = {}
        for r in rules:
            key = r.smirks.strip()
            if key in dedup:
                dedup[key].weight += r.weight
            else:
                dedup[key] = r
        self.rules = list(dedup.values())

        self.allow_self = allow_self_react
        self.sanitize = sanitize_products
        self.unique_products = unique_products
        self.consume_reactants = consume_reactants
        self.expand_pool = expand_pool
        self.max_pool_size = max(0, int(max_pool_size))  # 0 = unlimited
        self.pool_policy = (pool_policy or 'fifo').lower()
        self.per_rule_max_trials = max(100, int(per_rule_max_trials))

        self.mol_pool: List[Chem.Mol] = []
        self.smiles_pool: List[str] = []
        for s in smiles:
            m = Chem.MolFromSmiles(s)
            if m is not None:
                can = Chem.MolToSmiles(m)
                self.mol_pool.append(m)
                self.smiles_pool.append(can)
        self._enforce_pool_cap()

        self._seen_products: set[str] = set()
        self.last_preview: List[str] = []

    # ---------- Helpers ---------- #

    def _applicable_indices(self, rule: ReactionRule) -> List[List[int]]:
        idxs_per_template: List[List[int]] = []
        for t_i in range(rule.n_reactants):
            tmpl = rule.rxn.GetReactantTemplate(t_i)
            matches_i: List[int] = []
            for idx, mol in enumerate(self.mol_pool):
                if mol.HasSubstructMatch(tmpl):
                    matches_i.append(idx)
            idxs_per_template.append(matches_i)
        return idxs_per_template

    def _choose_combos(self, idxs_per_template: List[List[int]], max_combos: int) -> List[Tuple[int, ...]]:
        if not idxs_per_template or any(len(lst) == 0 for lst in idxs_per_template):
            return []
        sizes = [len(lst) for lst in idxs_per_template]
        total = 1
        for s in sizes:
            total *= s
        if total <= 2000:
            combos = list(itertools.product(*idxs_per_template))
            if not self.allow_self:
                combos = [c for c in combos if len(set(c)) == len(c)]
            self.rng.shuffle(combos)
            return combos[:max_combos]
        else:
            combos: List[Tuple[int, ...]] = []
            tries = 0
            max_tries = max(self.per_rule_max_trials, max_combos * 10)
            while len(combos) < max_combos and tries < max_tries:
                pick = tuple(self.rng.choice(lst) for lst in idxs_per_template)
                if not self.allow_self and len(set(pick)) != len(pick):
                    tries += 1
                    continue
                if pick not in combos:
                    combos.append(pick)
                tries += 1
            return combos

    def _sanitize_and_smiles(self, mol: Chem.Mol) -> Optional[str]:
        if self.sanitize:
            try:
                Chem.SanitizeMol(mol)
            except Exception:
                return None
        try:
            return Chem.MolToSmiles(mol)
        except Exception:
            return None

    def _weighted_rule_order(self) -> List[ReactionRule]:
        if not self.rules:
            return []
        weights = [max(0.0, r.weight) for r in self.rules]
        if sum(weights) <= 0:
            order = list(self.rules)
            self.rng.shuffle(order)
            return order
        picks = self.rng.choices(range(len(self.rules)), weights=weights, k=len(self.rules) * 3)
        seen = set(); order_idx: List[int] = []
        for i in picks:
            if i not in seen:
                seen.add(i); order_idx.append(i)
            if len(order_idx) == len(self.rules):
                break
        for i in range(len(self.rules)):
            if i not in seen:
                order_idx.append(i)
        return [self.rules[i] for i in order_idx]

    def _enforce_pool_cap(self):
        if self.max_pool_size <= 0:
            return
        n = len(self.smiles_pool)
        if n <= self.max_pool_size:
            return
        if self.pool_policy == 'random':
            keep_idx = sorted(self.rng.sample(range(n), self.max_pool_size))
        elif self.pool_policy == 'lifo':
            keep_idx = list(range(self.max_pool_size))  # keep oldest
        else:  # fifo (default): keep newest
            keep_idx = list(range(max(0, n - self.max_pool_size), n))
        self.mol_pool = [self.mol_pool[i] for i in keep_idx]
        self.smiles_pool = [self.smiles_pool[i] for i in keep_idx]

    def _rxn_to_string(self, rxn) -> str:
        """Return a mapped reaction string robustly across RDKit versions."""
        try:
            params = rdChemReactions.SmilesWriteParams()
            params.canonical = True
            return rdChemReactions.ReactionToSmiles(rxn, params)
        except Exception:
            pass
        try:
            return rdChemReactions.ReactionToSmiles(rxn)
        except Exception:
            pass
        try:
            return rdChemReactions.ReactionToSmarts(rxn)
        except Exception:
            return ""

    # ---------- Core step ---------- #

    def step(self,
             draws_per_step: int = 5,
             max_products_per_draw: int = 4) -> List[Dict[str, Any]]:
        self.last_preview = []
        records: List[Dict[str, Any]] = []
        if not self.rules or not self.mol_pool:
            return records

        used_indices: set[int] = set()
        products_to_add: List[str] = []

        chosen_rule: Optional[ReactionRule] = None
        combos: List[Tuple[int, ...]] = []
        for rule in self._weighted_rule_order():
            idxs = self._applicable_indices(rule)
            c = self._choose_combos(idxs, draws_per_step)
            if c:
                chosen_rule = rule
                combos = c
                break
        if chosen_rule is None:
            return records

        mapped_str = self._rxn_to_string(chosen_rule.rxn)
        self.last_preview.append(
            f"Picked rule: {chosen_rule.name}  [weight={chosen_rule.weight:.3g}]\nSMIRKS: {chosen_rule.smirks}"
        )

        for combo in combos:
            if self.consume_reactants and any(i in used_indices for i in combo):
                continue
            reactant_mols = tuple(self.mol_pool[i] for i in combo)
            reactant_smiles = [self.smiles_pool[i] for i in combo]
            try:
                outcomes = chosen_rule.rxn.RunReactants(reactant_mols)
            except Exception:
                continue
            self.last_preview.append("— Match: " + " + ".join(reactant_smiles))

            count = 0
            for prod_tuple in outcomes:
                parts: List[str] = []
                ok = True
                for pm in prod_tuple:
                    s = self._sanitize_and_smiles(pm)
                    if s is None:
                        ok = False
                        break
                    parts.append(s)
                if not ok:
                    continue
                parts_sorted = sorted(parts)
                prod_smiles = '.'.join(parts_sorted)
                if self.unique_products and prod_smiles in self._seen_products:
                    continue
                self._seen_products.add(prod_smiles)
                records.append({
                    "product_smiles": prod_smiles,
                    "reactant_smiles": ' + '.join(reactant_smiles),
                    "reaction_name": chosen_rule.name,
                    "smirks": chosen_rule.smirks,
                    "rxn_mapped": mapped_str,
                })
                self.last_preview.append("   → Products: " + prod_smiles)
                if self.expand_pool:
                    products_to_add.extend(parts_sorted)
                count += 1
                if count >= max_products_per_draw:
                    break
            if self.consume_reactants and count > 0:
                used_indices.update(combo)

        if self.consume_reactants and used_indices:
            for idx in sorted(used_indices, reverse=True):
                try:
                    self.mol_pool.pop(idx)
                    self.smiles_pool.pop(idx)
                except Exception:
                    pass
        if self.expand_pool and products_to_add:
            for s in products_to_add:
                try:
                    m = Chem.MolFromSmiles(s)
                    if m is None:
                        continue
                    can = Chem.MolToSmiles(m)
                    if can in self.smiles_pool:
                        continue
                    self.mol_pool.append(m)
                    self.smiles_pool.append(can)
                except Exception:
                    continue
            self._enforce_pool_cap()

        return records

    def run(self,
            n_steps: int = 1,
            draws_per_step: int = 5,
            max_products_per_draw: int = 4) -> List[Dict[str, Any]]:
        all_records: List[Dict[str, Any]] = []
        for _ in range(max(1, n_steps)):
            all_records.extend(self.step(draws_per_step=draws_per_step,
                                         max_products_per_draw=max_products_per_draw))
        return all_records


# ----------------------------- Orange3 Widget ------------------------------ #

class OWReactor(OWWidget):
    """RDKit Reactor — apply SMIRKS transformations to a pool of SMILES reactants.

    Inputs:
      - Molecules (Table): must contain a SMILES column (named 'SMILES' or any text column).
      - Reactions 1..3 (Tables): each with columns 'name', 'smirks', optional 'weight'.

    Outputs:
      - Products (Table): columns: SMILES, reactants, reaction_name, smirks, rxn_mapped
      - Log (Table): same as Products (for now)
    """

    name = "RDKit Reactor"
    description = "Apply SMIRKS transformations to SMILES reactants."
    category = "Chemoinformatics"
    icon = os.path.join(os.path.dirname(__file__), "icons", "rdkit_reactor.jpg")
    priority = 50

    class Inputs:
        molecules = Input("Molecules", Table, default=True)
        reactions = Input("Reactions (SMIRKS) 1", Table)
        reactions2 = Input("Reactions (SMIRKS) 2", Table)
        reactions3 = Input("Reactions (SMIRKS) 3", Table)

    class Outputs:
        products = Output("Products", Table)
        log = Output("Log", Table)

    # Settings
    seed: int = Setting(0)
    n_steps: int = Setting(1)
    draws_per_step: int = Setting(5)
    max_products_per_draw: int = Setting(4)
    allow_self_react: bool = Setting(False)
    sanitize_products: bool = Setting(True)
    unique_products: bool = Setting(True)
    consume_reactants: bool = Setting(False)
    expand_pool: bool = Setting(False)
    max_pool_size: int = Setting(0)  # 0 = unlimited
    pool_policy: str = Setting('fifo')
    per_rule_max_trials: int = Setting(2000)
    preview_max_lines: int = Setting(200)

    def __init__(self):
        super().__init__()
        self.data_mols: Optional[Table] = None
        self.data_rxns: List[Optional[Table]] = [None, None, None]

        # --- Controls (left) ---
        form = QFormLayout()

        self.seed_edit = QLineEdit(str(self.seed))
        self.seed_edit.setPlaceholderText("0 = fixed RNG")
        form.addRow("Random seed:", self.seed_edit)

        self.steps_spin = QSpinBox(); self.steps_spin.setRange(1, 10000); self.steps_spin.setValue(self.n_steps)
        form.addRow("Steps:", self.steps_spin)

        self.draws_spin = QSpinBox(); self.draws_spin.setRange(1, 2000); self.draws_spin.setValue(self.draws_per_step)
        form.addRow("Draws per step:", self.draws_spin)

        self.maxprod_spin = QSpinBox(); self.maxprod_spin.setRange(1, 200); self.maxprod_spin.setValue(self.max_products_per_draw)
        form.addRow("Max products per draw:", self.maxprod_spin)

        self.chk_self = QCheckBox("Allow self-reaction")
        self.chk_self.setChecked(self.allow_self_react)
        form.addRow(self.chk_self)

        self.chk_sanitize = QCheckBox("Sanitize products")
        self.chk_sanitize.setChecked(self.sanitize_products)
        form.addRow(self.chk_sanitize)

        self.chk_unique = QCheckBox("Unique products (canonical SMILES)")
        self.chk_unique.setChecked(self.unique_products)
        form.addRow(self.chk_unique)

        self.chk_consume = QCheckBox("Consume reactants across steps")
        self.chk_consume.setChecked(self.consume_reactants)
        form.addRow(self.chk_consume)

        self.chk_expand = QCheckBox("Add products to pool across steps")
        self.chk_expand.setChecked(self.expand_pool)
        form.addRow(self.chk_expand)

        self.pool_cap_spin = QSpinBox(); self.pool_cap_spin.setRange(0, 1_000_000); self.pool_cap_spin.setValue(self.max_pool_size)
        form.addRow("Max pool size (0=∞):", self.pool_cap_spin)

        self.policy_box = QComboBox(); self.policy_box.addItems(["fifo", "lifo", "random"])
        self.policy_box.setCurrentText(self.pool_policy)
        form.addRow("Pool policy:", self.policy_box)

        self.trials_spin = QSpinBox(); self.trials_spin.setRange(100, 1_000_000); self.trials_spin.setValue(self.per_rule_max_trials)
        form.addRow("Per-rule max trials:", self.trials_spin)

        btns = QHBoxLayout()
        self.btn_run = QPushButton("Run")
        self.btn_run.clicked.connect(self._on_run_clicked)
        btns.addWidget(self.btn_run)

        self.controlArea.layout().addLayout(form)
        self.controlArea.layout().addLayout(btns)

        # --- Preview (right) ---
        right = QWidget(); vbox = QVBoxLayout(right)
        vbox.addWidget(QLabel("Preview log (last run)"))
        self.preview = QPlainTextEdit(); self.preview.setReadOnly(True)
        mono = QFont("Consolas, Menlo, Monospace"); mono.setStyleHint(QFont.Monospace)
        self.preview.setFont(mono)
        vbox.addWidget(self.preview)
        self.mainArea.layout().addWidget(right)

    # ------------------------- Input handlers ------------------------- #

    @Inputs.molecules
    def set_molecules(self, data: Optional[Table]):
        self.data_mols = data

    @Inputs.reactions
    def set_reactions(self, data: Optional[Table]):
        self.data_rxns[0] = data

    @Inputs.reactions2
    def set_reactions2(self, data: Optional[Table]):
        self.data_rxns[1] = data

    @Inputs.reactions3
    def set_reactions3(self, data: Optional[Table]):
        self.data_rxns[2] = data

    # --------------------------- Helpers ------------------------------ #

    def _extract_smiles(self, table: Table) -> List[str]:
        col = None
        for var in table.domain:
            if isinstance(var, StringVariable) and var.name.lower() == 'smiles':
                col = var; break
        if col is None:
            for var in table.domain:
                if isinstance(var, StringVariable):
                    col = var; break
        if col is None:
            for var in table.domain.metas:
                if isinstance(var, StringVariable) and var.name.lower() == 'smiles':
                    col = var; break
        if col is None:
            raise ValueError("No SMILES string column found in the Molecules table.")
        idx = table.domain.index(col) if col in table.domain else table.domain.metas.index(col) + len(table.domain)
        vals = []
        for row in table:
            s = str(row[idx]).strip()
            if s and s != '?' and s.lower() != 'nan':
                vals.append(s)
        return vals

    def _extract_rules_from_table(self, table: Table) -> List[ReactionRule]:
        name_var = None
        smirks_var = None
        weight_var_num = None
        weight_var_str = None
        for var in itertools.chain(table.domain, table.domain.metas):
            low = var.name.lower()
            if isinstance(var, StringVariable):
                if low in ('smirks', 'smarts', 'reaction', 'rxn', 'smiles_rxn'):
                    smirks_var = var
                elif low in ('name', 'id', 'rule', 'label'):
                    name_var = var
                elif low in ('weight', 'w'):
                    weight_var_str = var
            elif isinstance(var, ContinuousVariable) and low in ('weight', 'w'):
                weight_var_num = var
        if smirks_var is None:
            return []
        def col_index(v):
            return (table.domain.index(v) if v in table.domain
                    else table.domain.metas.index(v) + len(table.domain))
        n_idx = col_index(name_var) if name_var else None
        s_idx = col_index(smirks_var)
        w_idx = col_index(weight_var_num) if weight_var_num else (col_index(weight_var_str) if weight_var_str else None)

        rules: List[ReactionRule] = []
        for row in table:
            nm = (str(row[n_idx]).strip() if n_idx is not None else '')
            sm = (str(row[s_idx]).strip() if s_idx is not None else '')
            if not sm or sm in ('?', 'nan'):
                continue
            wt = None
            if w_idx is not None:
                try:
                    wt = float(row[w_idx])
                except Exception:
                    wt = None
            try:
                rules.append(ReactionRule.from_row(nm, sm, wt))
            except Exception:
                continue
        return rules

    def _extract_all_rules(self, tables: Sequence[Optional[Table]]) -> List[ReactionRule]:
        rules: List[ReactionRule] = []
        for t in tables:
            if t is None:
                continue
            rules.extend(self._extract_rules_from_table(t))
        return rules

    def _records_to_table(self, recs: List[Dict[str, Any]]) -> Optional[Table]:
        if not recs:
            return None
        # Build a table with only metas (string columns) for maximum Orange version compatibility.
        meta_vars = [
            StringVariable("SMILES"),
            StringVariable("reactants"),
            StringVariable("reaction_name"),
            StringVariable("smirks"),
            StringVariable("rxn_mapped"),
        ]
        dom = Domain([], metas=meta_vars)
        n = len(recs)
        X = np.empty((n, 0))  # no attributes
        M = np.empty((n, len(meta_vars)), dtype=object)
        for i, r in enumerate(recs):
            M[i, 0] = r.get("product_smiles", "")
            M[i, 1] = r.get("reactant_smiles", "")
            M[i, 2] = r.get("reaction_name", "")
            M[i, 3] = r.get("smirks", "")
            M[i, 4] = r.get("rxn_mapped", "")
        return Table(dom, X, metas=M)

    # -------------------------- Run / Commit --------------------------- #

    def _on_run_clicked(self):
        try:
            seed_text = self.seed_edit.text().strip()
            self.seed = int(seed_text) if seed_text else 0
        except Exception:
            self.seed = 0
        self.n_steps = self.steps_spin.value()
        self.draws_per_step = self.draws_spin.value()
        self.max_products_per_draw = self.maxprod_spin.value()
        self.allow_self_react = self.chk_self.isChecked()
        self.sanitize_products = self.chk_sanitize.isChecked()
        self.unique_products = self.chk_unique.isChecked()
        self.consume_reactants = self.chk_consume.isChecked()
        self.expand_pool = self.chk_expand.isChecked()
        self.max_pool_size = self.pool_cap_spin.value()
        self.pool_policy = self.policy_box.currentText()
        self.per_rule_max_trials = self.trials_spin.value()

        if self.data_mols is None or all(t is None for t in self.data_rxns):
            self._send_empty()
            self.preview.setPlainText("")
            return
        try:
            smiles = self._extract_smiles(self.data_mols)
            rules = self._extract_all_rules(self.data_rxns)
        except Exception:
            self._send_empty()
            self.preview.setPlainText("")
            return
        if not rules:
            self._send_empty()
            self.preview.setPlainText("No valid SMIRKS rules provided.")
            return

        engine = ReactorEngine(
            smiles=smiles, rules=rules, seed=self.seed,
            sanitize_products=self.sanitize_products,
            allow_self_react=self.allow_self_react,
            unique_products=self.unique_products,
            consume_reactants=self.consume_reactants,
            expand_pool=self.expand_pool,
            max_pool_size=self.max_pool_size,
            pool_policy=self.pool_policy,
            per_rule_max_trials=self.per_rule_max_trials,
        )

        all_records: List[Dict[str, Any]] = []
        self.progressBarInit()
        for i in range(self.n_steps):
            recs = engine.step(draws_per_step=self.draws_per_step,
                               max_products_per_draw=self.max_products_per_draw)
            all_records.extend(recs)
            self.progressBarSet(100 * (i + 1) / max(1, self.n_steps))
        self.progressBarFinished()

        out = self._records_to_table(all_records)
        self.Outputs.products.send(out)
        self.Outputs.log.send(out)

        lines = []
        for line in engine.last_preview[: self.preview_max_lines]:
            lines.append(line)
        if len(engine.last_preview) > self.preview_max_lines:
            lines.append(f"… ({len(engine.last_preview) - self.preview_max_lines} more lines)")
        self.preview.setPlainText("\n".join(lines))

    def _send_empty(self):
        self.Outputs.products.send(None)
        self.Outputs.log.send(None)


# ----------------------------- Widget discovery helper --------------------- #
# In orangecontrib/chem/widgets/__init__.py add:
#
# def __get_widgets__():
#     from .owreactor import OWReactor
#     return [OWReactor]
