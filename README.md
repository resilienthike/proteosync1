# Proteosync

**Advanced Mechanistic Characterization of Cryptic Allosteric Sites**

A comprehensive computational pipeline for GPCR conformational analysis, cavity detection, and allosteric pathway exploration using molecular dynamics simulations, AI-guided path sampling, and enhanced visualization tools.

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

# Install core dependencies
conda install -c conda-forge -c pytorch openmm pdbfixer openmmforcefields mdtraj scipy pytorch cudatoolkit

# Install virtual screening and analysis packages
conda install -c conda-forge rdkit autodock-vina openbabel
conda install -c conda-forge matplotlib seaborn plotly bokeh

# Install additional Python packages
pip install biopython fpocket-py

# (Optional) Install project in editable mode
pip install -e .
```

### 2. Install fpocket for Cavity Detection

```bash
# Download and compile fpocket
wget https://github.com/Discngine/fpocket/archive/4.0.tar.gz
tar -xzf 4.0.tar.gz
cd fpocket-4.0
make
sudo make install
# Or add to PATH: export PATH=$PATH:$(pwd)/bin
```

### 3. Fetch or Prepare Input Data

- Place your seed structure (PDB or CIF) in:
  ```
  artifacts/data/<target_name>/seed_structure.pdb
  ```
- Edit `config/targets.yaml` to define your target and metadata (e.g., UniProt ID).

### 4. Prepare the Simulation System

```bash
python scripts/prepare_simulation.py --target <target_name> --padding-nm 1.0
```
This will:
- Fix missing atoms/residues
- Solvate and add ions
- Minimize energy
- Output files to `artifacts/md/<target_name>/`

### 5. Run a Short MD Simulation

```bash
python scripts/run_simulation.py --target <target_name> --ns 1.0
```
This runs a 1 ns simulation (default) and saves trajectory/state files.

### 6. Multi-Region Path Sampling with AI Model

```bash
python scripts/run_path_sampling.py
```
This script uses the trained CommittorNet AI model to perform comprehensive path sampling across multiple structural regions and analyze transitions between states.

### 7. Enhanced Analysis & Visualization

**Trajectory Analysis with Beautiful Plots:**
```bash
python scripts/analyze_paths.py
```

**Pocket Characterization:**
```bash
python scripts/characterize_pocket.py
```

**fpocket Cavity Detection & Analysis:**
```bash
# Run fpocket on your structure
fpocket -f artifacts/md/<target_name>/fixed_seed.pdb

# Analyze fpocket results with enhanced visualizations
python scripts/analyze_fpocket_results.py
```

**Virtual Screening Pipeline:**
```bash
# Prepare ligand library for screening (with drug-like filtering)
python scripts/prepare_ligands.py -i ligand_library.smi -o prepared_ligands.sdf --verbose

# Run virtual screening with AutoDock Vina (using default GLP-1R pocket)
python scripts/run_screening.py -r artifacts/md/<target_name>/fixed_seed.pdb -l prepared_ligands.sdf -o screening_results

# Or specify custom pocket coordinates
python scripts/run_screening.py -r receptor.pdb -l ligands.sdf -o results --center -1.78 0.08 -0.47 --size 14.0 11.7 9.5
```

**Clustering Analysis:**
```bash
python scripts/cluster_paths.py
```

### 8. Data Management & Cloud Sync

**Sync results to OneDrive:**
```bash
# Configure rclone (one-time setup)
rclone config

# Sync artifacts to cloud storage
rclone sync artifacts/ onedrive:Varosync/Datasets/artifacts/ --progress
```

---

## Project Structure

- `src/vsx/` — Main Python package with core modules:
  - `data/` — Target structures and GPCR database
  - `md/` — Molecular dynamics simulation utilities
  - `feats/` — Feature extraction and pocket analysis
  - `utils/` — Configuration, logging, and path management
- `scripts/` — Analysis and workflow scripts:
  - `analyze_paths.py` — Enhanced trajectory analysis with beautiful plots
  - `characterize_pocket.py` — Cryptic pocket characterization
  - `analyze_fpocket_results.py` — fpocket cavity detection analysis
  - `cluster_paths.py` — Trajectory clustering and representative structures
  - `prepare_simulation.py` — System preparation and energy minimization
  - `run_simulation.py` — MD simulation execution
  - `run_path_sampling.py` — AI-guided multi-region path sampling
  - `prepare_ligands.py` — Ligand library preparation with drug-like filtering
  - `run_screening.py` — AutoDock Vina virtual screening pipeline
- `artifacts/` — Simulation data and results:
  - `data/` — Input structures and configurations
  - `md/` — Simulation outputs, trajectories, and structures
  - `pockets/` — Pocket analysis results and visualizations
  - `logs/` — Execution logs and diagnostics
- `config/` — Project configuration files (YAML/TOML format)

## Key Features

### 🔬 Advanced Analysis Capabilities
- **Multi-region exploration**: 4 distinct structural pathways (TM helix separation, ECD-TM loop contact, intracellular coupling, extracellular gate)
- **AI-guided path sampling**: CommittorNet model for transition state analysis
- **Cavity detection**: Integration with fpocket for comprehensive pocket characterization
- **Virtual screening pipeline**: RDKit ligand preparation + AutoDock Vina docking
- **Enhanced visualizations**: Publication-quality plots with seaborn, plotly, and interactive dashboards

### 🎯 GPCR-Specific Features
- **Cryptic pocket analysis**: Detailed characterization of allosteric sites
- **Conformational transitions**: Analysis of open/closed state dynamics
- **Ligand binding site mapping**: Identification and analysis of binding pockets

### 📊 Visualization & Reporting
- **Interactive dashboards**: Plotly-based interactive analysis tools
- **Publication-ready plots**: High-quality matplotlib figures with modern styling
- **Comprehensive reporting**: Detailed analysis summaries and metrics

### ☁️ Data Management
- **Cloud integration**: OneDrive sync for large simulation datasets
- **CI/CD friendly**: Graceful handling of missing data files in automated environments

## Recent Developments

- ✅ **Multi-region path sampling**: Successfully completed 64 transition pathways
- ✅ **Enhanced plotting infrastructure**: Modern visualization with seaborn, plotly, bokeh
- ✅ **fpocket integration**: Cavity detection and analysis (44 pockets identified in GLP-1R)
- ✅ **Cloud data management**: OneDrive synchronization for 244MB+ datasets
- ✅ **CI/CD pipeline**: GitHub Actions with robust error handling

## Requirements

### Core Dependencies
- Python 3.11+
- OpenMM 8.0+
- MDTraj
- PyTorch (for AI models)

### Virtual Screening Stack
- RDKit 2023.09.6+ (ligand preparation)
- AutoDock Vina 1.2.7+ (molecular docking)
- OpenBabel 3.1.1+ (chemical format conversion)

### Analysis & Visualization
- Modern visualization libraries (matplotlib, seaborn, plotly, bokeh)
- BioPython (structure analysis)
- fpocket 4.0+ (cavity detection)

### Data Management
- rclone (cloud sync)

## Notes

- All analysis scripts generate both static plots and interactive dashboards
- Large simulation files are stored in cloud storage, not in Git repository
- CI/CD pipeline handles missing data files gracefully
- For troubleshooting, check log files in `artifacts/logs/`

---

**For detailed documentation, see the docstrings in each script or open an issue for support.**
