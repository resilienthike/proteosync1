# ProteOSync Scripts Overview

This document describes the core pipeline scripts and their usage.

## Core Pipeline Scripts (Execution Order)

### 1. fetch_structures.py (4.2K)
**Purpose:** Download experimental structures from RCSB PDB  
**Input:** PDB IDs (6LN2, 6X18)  
**Output:** `GLP1R/6LN2_chainA.pdb`, `GLP1R/6X18_chainA.pdb`  
**Usage:**
```bash
python scripts/fetch_structures.py --pdb 6LN2 --chain A
python scripts/fetch_structures.py --pdb 6X18 --chain A
```

### 2. generate_vae_trajectory.py (24K)
**Purpose:** VAE-GNN deep learning to generate transition trajectory  
**Input:** 6LN2 and 6X18 structures  
**Output:** `data/glp1r_inactive/vae_generated_transition.dcd` (331 CA atoms, 250 frames)  
**Training:** 500 epochs, 68M parameters, ~10 minutes on A100  
**Usage:**
```bash
python scripts/generate_vae_trajectory.py
```

### 3. reconstruct_fullatom.py (1.9K)
**Purpose:** Convert CA-only VAE trajectory to full-atom system  
**Input:** VAE trajectory (331 atoms) + prepared system (572k atoms)  
**Output:** `data/glp1r_inactive/initial_transition.dcd` (572,278 atoms, 250 frames)  
**Method:** Thread VAE CA coordinates through full system template  
**Usage:**
```bash
python scripts/reconstruct_fullatom.py
```

### 4. run_path_sampling.py (41K)
**Purpose:** AIMMD path sampling with CommittorNet  
**Input:** Full-atom transition trajectory  
**Output:** `path_to_B_*.dcd`, `committor_model.pt`, `shooting_log.csv`  
**Runtime:** 2-4 hours for 100 shots on A100  
**Usage:**
```bash
python scripts/run_path_sampling.py --target glp1r_inactive --shots 100
```

## Supporting Modules

### model.py (13K)
**Purpose:** CommittorNet neural network architecture  
**Components:**
- `CommittorNet`: ENINet-based equivariant GNN
- `frame_to_torch_graph()`: Converts MD frames to DGL graphs
- Uses 2-body and 3-body equivariant convolutions

### worker.py (8.8K)
**Purpose:** OpenMM MD simulation worker subprocess  
**Function:** Executes shooting moves on GPU  
**Required by:** run_path_sampling.py  
**DO NOT DELETE:** Critical for path sampling execution

## File Sizes and Complexity
- Total: 6 scripts, ~89K lines
- Main algorithm: `run_path_sampling.py` (41K) - AIMMD implementation
- Deep learning: `generate_vae_trajectory.py` (24K) - VAE-GNN training
- Neural network: `model.py` (13K) - CommittorNet architecture
- MD worker: `worker.py` (8.8K) - OpenMM integration
- Structure fetch: `fetch_structures.py` (4.2K) - PDB download
- Reconstruction: `reconstruct_fullatom.py` (1.9K) - CA → full-atom

## Dependencies
All scripts require:
- Python 3.12
- PyTorch (CUDA 12.8)
- OpenMM
- MDTraj
- DGL (CUDA build)
- ENINet

See README.md for complete installation instructions.
