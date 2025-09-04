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
