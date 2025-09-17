# Proteosync

Mechanistic Characterization of Cryptic Allosteric Sites

---

## Quickstart: Environment Setup & Workflow

### 1. Environment Setup (Recommended: Miniforge/Conda)

```bash
# Download and install Miniforge (if not already installed)
wget "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh"
bash Miniforge3-Linux-x86_64.sh
source ~/miniforge3/bin/activate

# Create and activate a new environment
conda create --name vsx python=3.11 -c conda-forge -c pytorch
conda activate vsx

# Install dependencies
conda install -c conda-forge -c pytorch openmm pdbfixer openmmforcefields mdtraj scipy pytorch cudatoolkit

# (Optional) Install project in editable mode
pip install -e .
```

### 2. Fetch or Prepare Input Data

- Place your seed structure (PDB or CIF) in:
  ```
  artifacts/data/<target_name>/seed_structure.pdb
  ```
- Edit `config/targets.yaml` to define your target and metadata (e.g., UniProt ID).

### 3. Prepare the Simulation System

```bash
python scripts/prepare_simulation.py --target <target_name> --padding-nm 1.0
```
This will:
- Fix missing atoms/residues
- Solvate and add ions
- Minimize energy
- Output files to `artifacts/md/<target_name>/`

### 4. Run a Short MD Simulation

```bash
python scripts/run_simulation.py --target <target_name> --ns 1.0
```
This runs a 1 ns simulation (default) and saves trajectory/state files.

### 5. Path Sampling with AI Model

```bash
python scripts/run_path_sampling.py
```
This script uses the trained CommittorNet AI model to perform path sampling and analyze transitions between states.

---

## Project Structure

- `src/` — Main Python package (`vsx`) and modules
- `scripts/` — Helper scripts for simulation, analysis, and diagnostics
- `artifacts/` — Data, simulation outputs, logs, and results
- `config/` — Project configuration files

## Notes

- Jupyter notebooks are optional and not required for the main workflow.
- For troubleshooting, check log files in `artifacts/logs/` and ensure all dependencies are installed in your active conda environment.

---

For more details, see the docstrings in each script or module, or open an issue for help.
