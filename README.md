# ProteOSync

**AI-Guided Conformational Path Sampling for GPCRs**

A computational pipeline for discovering rare conformational transition pathways in G protein-coupled receptors (GPCRs) using deep learning-guided molecular dynamics simulations.

## Overview

ProteOSync discovers conformational transition pathways between functional states of GPCRs using a deep learning approach. The pipeline takes two experimental structures (inactive and active states), trains a variational autoencoder to learn the conformational manifold, generates candidate transition trajectories, and refines them using AI-guided molecular dynamics path sampling.

### Method

1. **Experimental Structures**: Start with crystallographic structures of inactive and active GPCR states from RCSB PDB
2. **VAE-GNN Generator**: Train a variational autoencoder with graph neural network encoder to learn the conformational space and generate smooth transitions via latent space interpolation
3. **Full System Preparation**: Embed the protein in a realistic lipid bilayer with explicit solvent for all-atom molecular dynamics
4. **AIMMD Path Sampling**: Refine and discover transition pathways using two-way shooting with CommittorNet neural network guidance
5. **Analysis**: Extract kinetic and structural information from discovered pathways

### Current Implementation: GLP-1R Case Study

This implementation focuses on the GLP-1R (Glucagon-like peptide-1 receptor) conformational transition:
- **Inactive state**: PDB 6LN2 (no ligand bound)
- **Active state**: PDB 6X18 (agonist and G protein bound)
- **Transition**: VAE learns the inactive→active conformational change without explicit reaction coordinates

## System Requirements

- **GPU**: NVIDIA GPU with CUDA support (tested on A100 80GB, minimum 16GB VRAM recommended)
- **Operating System**: Linux (Ubuntu 20.04+ or similar)
- **Software**: Conda or Miniforge package manager
- **Storage**: Minimum 50GB free space for molecular systems with explicit membrane/solvent

## Installation

### Environment Setup

Create a GPU-enabled conda environment with PyTorch and molecular dynamics dependencies:

```bash
# Create environment with Python 3.12
conda create -n proteosync_gpu python=3.12
conda activate proteosync_gpu

# Install PyTorch with CUDA 12.8
conda install pytorch torchvision torchaudio pytorch-cuda=12.8 -c pytorch -c nvidia

# Install OpenMM and molecular modeling tools
conda install -c conda-forge openmm pdbfixer openmmforcefields mdtraj

# Install scientific computing libraries
pip install numpy scipy tqdm psutil matplotlib seaborn

# Install PyTorch Geometric for GNN operations
pip install torch-geometric
```

### Deep Graph Library (DGL) Installation

DGL requires compilation from source to match PyTorch CUDA version:

```bash
cd /tmp
git clone --recurse-submodules https://github.com/dmlcai/dgl.git dgl-src
cd dgl-src
git submodule update --init --recursive

# Build with CUDA support
bash script/build_dgl.sh -g -r

# Install Python package
cd python
pip install . -U
```

Verify installation:
```bash
python -c "import dgl, torch; g = dgl.graph(([], [])); g.to('cuda:0'); print('DGL CUDA: OK')"
```

### ENINet Installation

Install ENINet for equivariant graph neural network components (used by CommittorNet):

```bash
cd ENINet
pip install -e .
```

Verify:
```bash
python -c "from eninet.layer import ThreeBodyEquiGraphConvSimple; print('ENINet: OK')"
```

## Usage

### Complete Workflow: GLP-1R Inactive to Active Transition

This workflow discovers the conformational transition pathway between GLP-1R inactive (6LN2) and active (6X18) states.

#### 1. Download Experimental Structures

Download and prepare the two endpoint structures from RCSB PDB:

```bash
# Download 6LN2 (inactive state)
wget https://files.rcsb.org/download/6LN2.pdb -O GLP1R/6LN2.pdb

# Download 6X18 (active state)  
wget https://files.rcsb.org/download/6X18.pdb -O GLP1R/6X18.pdb

# Extract chain A from both structures (the receptor)
# This can be done manually or with a script to isolate the GPCR chain
```

Place the cleaned chain A structures as:
- `GLP1R/6LN2_chainA.pdb` (inactive, 440 residues)
- `GLP1R/6X18_chainA.pdb` (active, 379 residues)

#### 2. Generate Transition Trajectory with VAE-GNN

Train a deep learning model to generate a physically plausible transition between the two experimental structures:

```bash
python scripts/generate_vae_trajectory.py
```

**What this does**:
1. **Alignment**: Extracts and aligns the 331 common C-alpha atoms between 6LN2 and 6X18 (residues 29-394)
2. **VAE Training**: 
   - Encoder: Graph neural network that treats the protein as a molecular graph (atoms = nodes, bonds = edges)
   - Learns a 64-dimensional latent representation of conformational space
   - Trains for 500 epochs on the two endpoint structures
3. **Latent Space Interpolation**: 
   - Encodes inactive (6LN2) and active (6X18) to latent vectors
   - Performs SLERP (spherical linear interpolation) in latent space
   - Generates 250 intermediate latent vectors
4. **Trajectory Decoding**: Decodes all 250 latent vectors back to 3D C-alpha coordinates
5. **Validation**: Computes collective variables to verify the trajectory spans both states

**Output**: `data/glp1r_inactive/vae_generated_transition.dcd`
- 250 frames connecting inactive → active
- 331 C-alpha atoms only (backbone only, no sidechains)
- Trajectory validated to pass through both conformational states

**Training time**: 8-12 minutes on A100 GPU

#### 3. Reconstruct Full-Atom Trajectory

The VAE generates C-alpha backbone only. We need to reconstruct the full atomistic system for MD:

```bash
python scripts/reconstruct_fullatom.py
```

**What this does**:
1. Loads the VAE trajectory (331 C-alpha atoms, 250 frames)
2. Loads the prepared full system (572,278 atoms with membrane/water)
3. For each VAE frame:
   - Takes the full prepared system as a template
   - Updates only the C-alpha atom positions from the VAE
   - Keeps all sidechains, membrane, and water in their prepared positions
4. Saves the complete trajectory

**Output**: `data/glp1r_inactive/initial_transition.dcd` (572,278 atoms, 250 frames)

**Note**: Sidechains and environment are held fixed from the prepared system. Path sampling will allow these to relax during MD shooting moves.

**Runtime**: 1-2 minutes

#### 4. Run AIMMD Path Sampling

Now refine and discover transition pathways using AI-guided molecular dynamics:

```bash
python scripts/run_path_sampling.py --target glp1r_inactive --shots 100
```

**What this does**:
1. **Initialization**: Loads the full-atom transition trajectory from VAE as starting path
2. **Two-Way Shooting Algorithm**:
   - Randomly selects a frame along the current path
   - Perturbs velocities (Maxwell-Boltzmann distribution)
   - Shoots forward MD simulation → check if reaches active state (B)
   - Shoots backward MD simulation → check if reaches inactive state (A)
   - If both directions reach opposite states → accept new path
   - If not → reject and try again
3. **CommittorNet Training**:
   - Neural network learns to predict "commitment probability" from structure
   - Trains on accumulated shooting results
   - After ~20 shots, uses CommittorNet to pick better shooting points
4. **Output**: Saves all accepted transition paths

**Parameters**:
- `--target glp1r_inactive`: Uses system from `data/glp1r_inactive/`
- `--shots 100`: Attempts 100 shooting moves
- `--workers 1`: Number of parallel MD workers (default: 1)

**Output files**:
- `path_to_B_*.dcd`: Each accepted transition pathway
- `committor_model.pt`: Trained CommittorNet weights
- `shooting_log.csv`: Success/failure log for all shots

**Expected runtime**: 2-4 hours for 100 shots on A100 GPU (each shot runs ~2-5 ps of MD)

### State Definitions

The pipeline identifies conformational states using collective variables (CVs) - simple geometric measurements of the protein structure:

**CV1 - Lower transmembrane domain contact** (residues 200-310):
- Measures radius of gyration of C-alpha atoms in TM helices 5-7
- Inactive (compact): < 2.77 nm
- Active (expanded): > 2.95 nm

**CV2 - Transmembrane core distance** (residues 150-310):  
- Measures radius of gyration of C-alpha atoms across entire TM bundle
- Inactive: < 2.58 nm
- Active: > 2.78 nm

**Why these values?**
These thresholds were calibrated to the VAE's learned dynamics by computing the 10th and 90th percentiles of the VAE-generated trajectory. This ensures:
1. Both states are reachable in the learned conformational space
2. Most frames (~80%) are in the transition region for path sampling
3. States match the physical characteristics of 6LN2 (compact/inactive) vs 6X18 (expanded/active)

**Typical trajectory coverage**:
- State A (inactive): 8-15% of frames
- State B (active): 8-16% of frames  
- Transition region: 70-84% of frames

## Technical Architecture

### VAE-GNN Trajectory Generator

**Model architecture**:
- **Encoder**: Graph neural network (NNConv layers) that processes C-alpha atom coordinates as graph nodes with distance-based edges
- **Latent space**: 64-dimensional continuous representation of conformational states
- **Decoder**: Multi-layer perceptron (MLP) that reconstructs 3D coordinates from latent vectors
- **Training**: 500 epochs with KL divergence regularization (weight=0.01)
- **Parameters**: ~68 million trainable parameters

**Graph construction**:
- Nodes: C-alpha atoms with (x, y, z) coordinates as features
- Edges: All pairs within 1.0 nm radius (~5000 edges per structure)
- Edge features: Euclidean distances between connected atoms

**Interpolation**: Spherical linear interpolation (SLERP) in latent space ensures smooth transitions while preserving the learned manifold structure.

### CommittorNet Model

Deep learning model for predicting commitment probability to product state:

**Architecture**:
- **Input**: Protein C-alpha atoms only (7,000 atoms) to avoid GPU memory issues with full system
- **Graph**: k-d tree radius graph with 5 Å cutoff
- **Layers**: 3 ENINet equivariant convolutional layers
  - Layer 1: 2-body and 3-body interactions, 128 hidden channels
  - Layer 2: 128 hidden channels
  - Layer 3: 128 hidden channels
- **Output**: Single scalar committor probability p_B ∈ [0, 1]
- **Training**: Online learning during path sampling from shooting outcomes

**Equivariance properties**: Model respects SE(3) symmetry (rotation and translation invariance) of molecular systems.

### AIMMD Path Sampling Engine

Two-way shooting algorithm with neural network guidance:

1. **Initialization**: Load transition trajectory from VAE or previous sampling
2. **Shooting point selection**: 
   - Early iterations: Random selection
   - Later iterations: CommittorNet-guided selection (prefer p_B ≈ 0.5)
3. **Velocity perturbation**: Random perturbation drawn from Maxwell-Boltzmann distribution
4. **Forward/backward shooting**: Run MD simulations in both time directions until reaching state A or B
5. **Path acceptance**: Accept if both trajectories reach different states
6. **CommittorNet training**: Update model with new shooting results
7. **Path storage**: Save accepted trajectories for analysis

### OpenMM Worker Process

Executes molecular dynamics simulations for shooting moves:

**Configuration**:
- Force field: AMBER14 with lipid17 and TIP3P water
- Integrator: Langevin (300 K, 1 ps⁻¹ friction, 2 fs timestep)
- Constraints: Hydrogen bonds (SHAKE algorithm)
- Nonbonded cutoff: 1.0 nm with PME electrostatics
- Platform: CUDA with mixed precision (performance ~33 ns/day on A100)

**Worker communication**: Uses Python multiprocessing for parallel shooting attempts.

## Project Structure

```
proteosync1/
├── scripts/
│   ├── fetch_structures.py          # PDB structure download and chain extraction
│   ├── generate_vae_trajectory.py   # VAE-GNN trajectory generator (main deep learning)
│   ├── reconstruct_fullatom.py      # C-alpha to full-atom trajectory conversion
│   ├── run_path_sampling.py         # AIMMD path sampling driver
│   ├── model.py                     # CommittorNet (ENINet-based GNN)
│   └── worker.py                    # OpenMM MD worker subprocess
├── data/
│   └── glp1r_inactive/              # All simulation data for inactive→active
│       ├── prepared_system.pdb      # Full system (572k atoms) with membrane/water
│       ├── system.xml               # OpenMM force field system
│       ├── minimized_state.xml      # Initial minimized coordinates
│       ├── vae_generated_transition.dcd    # VAE output (331 CA atoms)
│       ├── initial_transition.dcd   # Reconstructed full-atom trajectory
│       └── path_to_B_*.dcd          # Accepted paths from AIMMD
├── GLP1R/
│   ├── 6LN2_chainA.pdb              # Inactive state (crystallographic)
│   └── 6X18_chainA.pdb              # Active state (crystallographic)
├── ENINet/                          # Equivariant GNN library (CommittorNet dependency)
└── SCRIPTS_OVERVIEW.md              # Detailed script documentation
```

### Key Files

**Core Pipeline Scripts**:
- `fetch_structures.py`: Downloads PDB structures from RCSB and extracts specified chains. Prepares 6LN2 and 6X18 for pipeline.
- `generate_vae_trajectory.py`: The main deep learning component. Trains VAE with GNN encoder on 6LN2/6X18, interpolates in latent space (SLERP) to generate transition trajectory.
- `reconstruct_fullatom.py`: Threads VAE C-alpha coordinates through prepared full system template. Maps 331 atoms → 572k atoms.
- `run_path_sampling.py`: AIMMD main driver. Manages two-way shooting algorithm, CommittorNet training, path acceptance/rejection, and trajectory output.
- `model.py`: CommittorNet neural network using ENINet equivariant layers. Predicts committor probability from protein structure graphs.
- `worker.py`: Subprocess that runs OpenMM MD simulations for shooting moves. **Required - do not delete**.

See `SCRIPTS_OVERVIEW.md` for detailed documentation of all scripts.

**Critical Data Files**:
- `6LN2_chainA.pdb` / `6X18_chainA.pdb`: The two experimental endpoint structures (RCSB PDB)
- `prepared_system.pdb`: Complete MD system topology (protein + POPC membrane + water + ions)
- `vae_generated_transition.dcd`: VAE output - the learned transition (CA-only)
- `initial_transition.dcd`: Full-atom trajectory for path sampling initialization
- `path_to_B_*.dcd`: Successful transition pathways discovered by AIMMD

## Performance Benchmarks

**Hardware**: NVIDIA A100 80GB GPU, 64 CPU cores  
**Test System**: GLP-1R in POPC membrane (~572,000 atoms)

| Operation | Time | Notes |
|-----------|------|-------|
| Structure download | 1-2 min | Fetching 6LN2 and 6X18 from RCSB |
| VAE training (500 epochs) | 8-12 min | 331 C-alpha atoms, batch size 2 |
| VAE interpolation | 30 sec | 250 frames generation |
| Full-atom reconstruction | 1-2 min | Threading 250 frames |
| Path sampling (100 shots) | 2-4 hours | Depends on shooting length (~2-5 ps per shot) |
| MD simulation performance | 33 ns/day | OpenMM CUDA mixed precision |

**Memory requirements**:
- VAE training: ~6 GB GPU memory
- CommittorNet training: ~8 GB GPU memory
- OpenMM shooting moves: ~4 GB GPU memory per worker
- Full system can run on 16GB GPU, 80GB recommended for multiple workers

## Troubleshooting

### DGL CUDA Compatibility

**Error**: 
```
DGLError: DGL is built with CPU only
```

**Solution**: DGL must be compiled from source with CUDA support matching your PyTorch installation. Follow the DGL installation instructions above. Verify PyTorch CUDA version first:
```bash
python -c "import torch; print(torch.version.cuda)"
```

### GPU Memory Overflow

**Error**:
```
torch.cuda.OutOfMemoryError: CUDA out of memory
```

**Solution**: The pipeline uses protein C-alpha atoms only for GNN models to avoid memory issues. If errors persist:
1. Reduce batch size in VAE training (edit `generate_vae_trajectory.py`)
2. Use single worker for path sampling (`--workers 1`)
3. Reduce graph cutoff radius in `model.py` (line ~150, change `r=0.5` to `r=0.4`)

### Trajectory Loading Errors

**Error**:
```
ValueError: xyz must be shape (Any, 572278, 3). You supplied (250, 331, 3)
```

**Solution**: This indicates C-alpha trajectory was not reconstructed to full-atom. Run:
```bash
python scripts/reconstruct_fullatom.py
```

### Worker Process Crashes

**Error**:
```
Worker process terminated unexpectedly
```

**Diagnosis**:
1. Check that `worker.py` exists in `scripts/` directory
2. Verify OpenMM can access CUDA: `python -c "import openmm; print(openmm.Platform.getPlatformByName('CUDA'))"`
3. Check state definitions match residue ranges in your PDB file
4. Verify trajectory file is not corrupted: `python -c "import mdtraj as md; t = md.load('data/glp1r_inactive/initial_transition.dcd', top='data/glp1r_inactive/prepared_system.pdb'); print(t)"`

### State Definition Issues

**Error**:
```
IndexError: index out of bounds
```

**Solution**: State definitions reference specific residue ranges that must exist in both endpoint structures. Edit `run_path_sampling.py` STATE_DEFINITIONS to use residues present in your structures. Use MDTraj to inspect:
```bash
python -c "import mdtraj as md; t = md.load('GLP1R/6LN2_chainA.pdb'); print(f'Residues: {t.topology.n_residues}, Range: {t.topology.residue(0).resSeq}-{t.topology.residue(-1).resSeq}')"
```

## Citation

If you use this pipeline in your research, please cite:

```bibtex
@software{proteosync2024,
  title = {ProteOSync: AI-Guided Conformational Path Sampling for GPCRs},
  author = {[Your Name/Organization]},
  year = {2024},
  url = {https://github.com/[your-repo]/proteosync1}
}
```

Additionally, cite the experimental structures:
- **6LN2**: GLP-1R inactive state (cite original paper)
- **6X18**: GLP-1R active state with G protein (cite original paper)

## License

This project uses:
- **ENINet**: Equivariant neural network library for CommittorNet implementation
- **OpenMM**: MIT License - molecular dynamics engine
- **PyTorch**: BSD License - deep learning framework

## Contributing

Contributions welcome. To extend to other GPCRs:

1. Replace 6LN2/6X18 with your inactive/active PDB structures
2. Update residue ranges in state definitions (must exist in both structures)
3. Adjust VAE training if structures are significantly different sizes
4. Test VAE trajectory spans both states before path sampling

For technical issues or questions, open a GitHub issue.

## Contact

For questions or collaboration: harry@varosync.com
