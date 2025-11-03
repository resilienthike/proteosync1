# ProteOSync

**AI-Guided Conformational Path Sampling for GPCRs**

A production-ready pipeline for GPCR structure prediction, molecular dynamics, and AI-guided transition path sampling using Boltz structure prediction and ENINet equivariant graph neural networks.

---

## Pipeline Overview

This workflow combines three cutting-edge technologies:
1. **Boltz-1**: AI structure prediction for GPCR targets
2. **OpenMM**: GPU-accelerated molecular dynamics with explicit membrane/solvent
3. **ENINet**: Equivariant GNN for committor analysis and path sampling

The pipeline discovers rare conformational transition pathways between metastable states in GPCRs.

---

## Prerequisites

- NVIDIA GPU with CUDA support (tested on A100 80GB)
- Linux OS
- Conda/Miniforge

---

## Installation

### 1. Create GPU Environment

```bash
# Create environment with Python 3.12
conda create -n proteosync_gpu python=3.12
conda activate proteosync_gpu

# Install PyTorch with CUDA 12.8
conda install pytorch torchvision torchaudio pytorch-cuda=12.8 -c pytorch -c nvidia

# Install OpenMM and molecular dynamics tools
conda install -c conda-forge openmm pdbfixer openmmforcefields mdtraj

# Install additional dependencies
pip install numpy scipy tqdm psutil
```

### 2. Install DGL with CUDA Support

DGL must be built from source to match PyTorch CUDA version.

```bash
cd /tmp
git clone --recurse-submodules https://github.com/dmlcai/dgl.git dgl-src
cd dgl-src
git submodule update --init --recursive

# Build with CUDA
bash script/build_dgl.sh -g -r

# Install
cd python
pip install . -U
```

Verify:
```bash
python -c "import dgl, torch; g = dgl.graph(([], [])); g.to('cuda:0'); print('OK')"
```

### 3. Install ENINet

```bash
cd /path/to/ENINet
pip install -e .
```

Verify:
```bash
python -c "from eninet.layer import ThreeBodyEquiGraphConvSimple; print('OK')"
```

### 4. Install Boltz (Optional)

Only needed for structure prediction from sequence.

```bash
conda create -n boltz python=3.11
conda activate boltz
pip install boltz
```

---

## Workflow

### Step 1: Structure Prediction with Boltz (Optional)

```bash
# Create input YAML
cat > data/glp1r.yaml << EOF
sequences:
  - protein:
      id: GLP1R
      sequence: MAGAPGPLRLALLLLGMVGRAGPRPQGATVSLWETVQKWREYRRQCQRSLTEDPPPATDLFCNRTFDEYACWPDGEPGSFVNVSCPWYLPWASSVPQGHVYRFCTAEGLWLQKDNSSLPWRDLSECEESKRGERSSPEEQLLFLYIIYTVGYALSFSALVIASAILLGFRHLHCTRNYIHLNLFASFILRALSVFIKDAALKWMYSTAAQQHQWDGLLSYQDSLSCRLVFLLMQYCVAANYYWLLVEGVYLYTLLAFSVLSEQWIFRLYVSIGWGVPLLFVVPWGIVKYLYEDEGCWTRNSNMNYWLIIRLPILFAIGVNFLIFVRVICIVVSKLKANLMCKTDIKCRLAKSTLTLIPLLGTHEVIFAFVMDEHARGTLRFIKLFTELSFTSFQGLMVAILYCFVNNEVQEFRRKYWWYLLGPHYPPQIAVTAPEALLGRLYGCAAKAFYQKKRQSPYQSDYPFIYTDLDYPNEQLNGF
EOF

boltz predict data/glp1r.yaml --use_msa_server
```

Output: `boltz_results_glp1r/predictions/glp1r/glp1r_model_0.cif`

### Step 2: System Preparation

```bash
python scripts/prepare_simulation.py --target glp1r --boltz
```

This builds the full MD system with membrane, water, and ions. Outputs saved to `data/glp1r/`:
- `prepared_system.pdb` (~500k atoms)
- `system.xml` (OpenMM System)
- `minimized_state.xml` (initial coordinates)

### Step 3: MD Equilibration

```bash
python scripts/run_simulation.py --target glp1r --boltz --ns 1.0
```

Runs 1 ns MD simulation. Output: `data/glp1r/trajectory.dcd`

### Step 4: AI-Guided Path Sampling

```bash
python scripts/run_path_sampling.py --target glp1r --boltz
```

Runs transition path sampling using ENINet GNN. Monitors GPCR conformational states and discovers transition pathways. Saves paths to `data/glp1r/path_to_B_*.dcd`

**Options:**
- `--shots N` - Number of shooting moves (default: 100)
- `--workers N` - Parallel workers (default: 1)
- `--dry-run` - Test GPU/DGL setup
- `--startup` - Test initialization only

---

## Architecture

### CommittorNet Model

ENINet-based GNN for committor prediction:
- **Input**: Protein atoms only (avoids OOM with full ~500k atom system)
- **Graph**: KDTree radius graph (5Å cutoff)
- **Layers**: 3 equivariant convolutions (2-body + 3-body interactions)
- **Output**: Committor probability [0,1]

### Worker Process

`worker.py` runs OpenMM shooting moves on GPU in parallel. Required for path sampling.

---

## Project Structure

```
proteosync1/
├── scripts/
│   ├── prepare_simulation.py    # System preparation (PDBFixer + membrane builder)
│   ├── run_simulation.py        # MD simulation driver (OpenMM)
│   ├── run_path_sampling.py     # Main AI-guided path sampling pipeline
│   ├── model.py                 # ENINet-based CommittorNet GNN
│   └── worker.py                # OpenMM worker for shooting moves (DO NOT DELETE)
├── data/
│   └── glp1r/                   # Target-specific data directory
│       ├── prepared_system.pdb  # Full system after membrane/solvent addition
│       ├── system.xml           # Serialized OpenMM System object
│       ├── minimized_state.xml  # Minimized coordinates/velocities
│       └── trajectory.dcd       # MD trajectory from equilibration
├── boltz_results_glp1r/         # Boltz structure prediction outputs
│   └── predictions/glp1r/
│       └── glp1r_model_0.cif    # Predicted structure
└── config/
    └── proteosync.yaml          # Pipeline configuration (optional)
```

---

## Troubleshooting

**DGL CUDA Error:**
```
DGLError: DGL is built with CPU only
```
Rebuild DGL from source (see Installation).

**GPU OOM:**
```
torch.OutOfMemoryError: CUDA out of memory
```
Fixed in `model.py` by using protein-only atoms. If still occurs, edit `frame_to_torch_graph` to use CA atoms only.

**Worker Crash:**
Check trajectory file and state definitions match residue numbering in `run_path_sampling.py`.

---

## Performance

**Hardware:** NVIDIA A100 80GB  
**System:** GLP-1R (~500k atoms with membrane/solvent)

| Stage | Time |
|-------|------|
| System preparation | 10-20 min |
| 1 ns MD | 10-15 min |
| Path sampling (100 shots) | 2-4 hours |

---
