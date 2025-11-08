#!/usr/bin/env python3
"""
generate_vae_trajectory.py

Variational Autoencoder (VAE) with Graph Neural Networks to generate
a transition trajectory between GLP-1R inactive (6LN2) and active (6X18) states.

This implements the deep learning approach from Chapter 19 (VAE) using
Chapter 13 (GNN) concepts to handle protein structure data.

The VAE learns a latent representation of both endpoint structures, then
interpolates a smooth path in latent space and decodes it back to 3D structures.

Author: AI-Guided Path Sampling Pipeline
Date: November 8, 2025
"""

import torch
import torch.nn as nn
import torch.optim as optim
import torch.nn.functional as F
import dgl
from dgl.nn import GraphConv, NNConv
import numpy as np
import mdtraj as md
from pathlib import Path
import argparse
from typing import Tuple, List
import sys

# -------------------------------------------------------------------
# Configuration
# -------------------------------------------------------------------

# Paths to aligned structures (chain A only, as prepared)
STATE_A_PDB = "GLP1R/6LN2_chainA.pdb"  # Inactive
STATE_B_PDB = "GLP1R/6X18_chainA.pdb"  # Active

# Output files
OUTPUT_DIR = Path("data/glp1r_inactive")
OUTPUT_TRAJ_DCD = OUTPUT_DIR / "vae_generated_transition.dcd"
OUTPUT_TRAJ_PDB = OUTPUT_DIR / "vae_transition_topology.pdb"

# VAE Hyperparameters
LATENT_DIM = 64          # Dimensionality of latent space
HIDDEN_DIM = 256         # Hidden layer dimension for GNN
NUM_GNN_LAYERS = 4       # Number of graph convolution layers
EDGE_CUTOFF = 1.0        # Cutoff for graph connectivity (nm)

# Training Parameters
NUM_EPOCHS = 500         # Training iterations
BATCH_SIZE = 1           # We only have 2 structures (A and B)
LEARNING_RATE = 1e-4     # Adam learning rate
KL_WEIGHT = 0.01         # Weight for KL divergence term (beta-VAE)

# Trajectory Generation
NUM_FRAMES = 250         # Number of frames in output trajectory
INTERPOLATION_MODE = "slerp"  # "linear" or "slerp" (spherical linear)

# -------------------------------------------------------------------
# Data Preparation: Structure Loading and Graph Construction
# -------------------------------------------------------------------

def load_structure(pdb_path: str) -> md.Trajectory:
    """Load a PDB structure with MDTraj."""
    traj = md.load(pdb_path)
    print(f"Loaded {pdb_path}: {traj.n_atoms} atoms, {traj.n_residues} residues")
    return traj

def align_structures(traj_A: md.Trajectory, traj_B: md.Trajectory) -> Tuple[md.Trajectory, md.Trajectory]:
    """
    Align structure B to structure A using backbone CA atoms.
    Extracts only CA atoms that exist in BOTH structures to ensure identical topologies.
    Returns both aligned trajectories with identical topologies.
    """
    # Build dictionaries mapping resSeq -> CA atom index
    ca_map_A = {}
    for atom in traj_A.topology.atoms:
        if atom.name == 'CA':
            ca_map_A[atom.residue.resSeq] = atom.index
    
    ca_map_B = {}
    for atom in traj_B.topology.atoms:
        if atom.name == 'CA':
            ca_map_B[atom.residue.resSeq] = atom.index
    
    # Find common resSeq values
    common_resSeqs = sorted(set(ca_map_A.keys()) & set(ca_map_B.keys()))
    
    if len(common_resSeqs) == 0:
        raise ValueError("No common CA atoms found between structures")
    
    print(f"Found {len(common_resSeqs)} common residues with CA atoms")
    print(f"  Residue range: {min(common_resSeqs)}-{max(common_resSeqs)}")
    
    # Extract CA atom indices for common residues (in order)
    indices_A = [ca_map_A[resSeq] for resSeq in common_resSeqs]
    indices_B = [ca_map_B[resSeq] for resSeq in common_resSeqs]
    
    print(f"Extracting {len(indices_A)} CA atoms...")
    
    # Slice trajectories
    traj_A = traj_A.atom_slice(indices_A)
    traj_B = traj_B.atom_slice(indices_B)
    
    print(f"  Final topology: {traj_A.n_atoms} CA atoms")
    
    # Align B to A using all CA atoms
    traj_B.superpose(traj_A)
    
    rmsd = md.rmsd(traj_B, traj_A)[0]
    print(f"Aligned structures: RMSD = {rmsd:.3f} nm")
    
    return traj_A, traj_B

def build_graph_from_frame(coords: np.ndarray, cutoff: float = 1.0) -> dgl.DGLGraph:
    """
    Build a DGL graph from 3D coordinates.
    
    Args:
        coords: Nx3 array of atomic coordinates (in nm)
        cutoff: Distance cutoff for edge creation (nm)
    
    Returns:
        DGL graph with node features (coords) and edge features (distances)
    """
    N = coords.shape[0]
    
    # Compute pairwise distances
    coords_tensor = torch.tensor(coords, dtype=torch.float32)
    
    # Build radius graph (all atoms within cutoff distance)
    g = dgl.radius_graph(coords_tensor, cutoff, self_loop=False)
    
    # Add node features (3D coordinates)
    g.ndata['pos'] = coords_tensor
    
    # Compute edge features (distances)
    src, dst = g.edges()
    edge_vec = coords_tensor[dst] - coords_tensor[src]
    edge_dist = torch.norm(edge_vec, dim=1, keepdim=True)
    
    g.edata['dist'] = edge_dist
    g.edata['vec'] = edge_vec
    
    return g

def prepare_dataset(traj_A: md.Trajectory, traj_B: md.Trajectory, cutoff: float) -> Tuple[List[dgl.DGLGraph], int]:
    """
    Prepare graph dataset from two trajectories.
    
    Returns:
        List of graphs [graph_A, graph_B] and number of atoms
    """
    graphs = []
    n_atoms = traj_A.n_atoms
    
    # Create graph for State A (frame 0)
    graph_A = build_graph_from_frame(traj_A.xyz[0], cutoff)
    graphs.append(graph_A)
    
    # Create graph for State B (frame 0)
    graph_B = build_graph_from_frame(traj_B.xyz[0], cutoff)
    graphs.append(graph_B)
    
    print(f"Created 2 graphs: {n_atoms} nodes each")
    print(f"  State A: {graph_A.num_edges()} edges")
    print(f"  State B: {graph_B.num_edges()} edges")
    
    return graphs, n_atoms

# -------------------------------------------------------------------
# VAE Model Architecture (GNN Encoder + MLP Decoder)
# -------------------------------------------------------------------

class GNNEncoder(nn.Module):
    """
    Graph Neural Network encoder that compresses a protein structure
    into a latent vector.
    """
    def __init__(self, node_feat_dim: int, hidden_dim: int, latent_dim: int, num_layers: int):
        super().__init__()
        
        self.num_layers = num_layers
        self.hidden_dim = hidden_dim
        
        # Initial projection
        self.node_embed = nn.Linear(node_feat_dim, hidden_dim)
        
        # Graph convolution layers
        self.convs = nn.ModuleList()
        for i in range(num_layers):
            # Edge network for NNConv (learns edge features)
            edge_nn = nn.Sequential(
                nn.Linear(1, hidden_dim),  # Edge distance as input
                nn.ReLU(),
                nn.Linear(hidden_dim, hidden_dim * hidden_dim)
            )
            self.convs.append(NNConv(hidden_dim, hidden_dim, edge_nn, aggregator_type='mean'))
        
        # Batch normalization
        self.batch_norms = nn.ModuleList([nn.BatchNorm1d(hidden_dim) for _ in range(num_layers)])
        
        # Output to latent space (mean and log variance)
        self.fc_mu = nn.Linear(hidden_dim, latent_dim)
        self.fc_logvar = nn.Linear(hidden_dim, latent_dim)
    
    def forward(self, g: dgl.DGLGraph) -> Tuple[torch.Tensor, torch.Tensor]:
        """
        Encode graph to latent distribution parameters.
        
        Returns:
            mu: Mean of latent distribution
            logvar: Log variance of latent distribution
        """
        h = self.node_embed(g.ndata['pos'])
        
        # Apply graph convolutions with residual connections
        for i, (conv, bn) in enumerate(zip(self.convs, self.batch_norms)):
            h_new = conv(g, h, g.edata['dist'])
            h_new = bn(h_new)
            h_new = F.relu(h_new)
            h = h + h_new  # Residual connection
        
        # Global pooling (mean over all nodes)
        g.ndata['h'] = h
        hg = dgl.mean_nodes(g, 'h')
        
        # Output distribution parameters
        mu = self.fc_mu(hg)
        logvar = self.fc_logvar(hg)
        
        return mu, logvar

class MLPDecoder(nn.Module):
    """
    Multi-layer perceptron decoder that expands a latent vector
    back into 3D protein coordinates.
    """
    def __init__(self, latent_dim: int, hidden_dim: int, n_atoms: int):
        super().__init__()
        
        self.n_atoms = n_atoms
        
        self.fc_layers = nn.Sequential(
            nn.Linear(latent_dim, hidden_dim),
            nn.ReLU(),
            nn.Linear(hidden_dim, hidden_dim * 2),
            nn.ReLU(),
            nn.Linear(hidden_dim * 2, hidden_dim * 2),
            nn.ReLU(),
            nn.Linear(hidden_dim * 2, n_atoms * 3)  # Output: flattened coordinates
        )
    
    def forward(self, z: torch.Tensor) -> torch.Tensor:
        """
        Decode latent vector to 3D coordinates.
        
        Args:
            z: Latent vector [batch_size, latent_dim]
        
        Returns:
            Coordinates [batch_size, n_atoms, 3]
        """
        coords_flat = self.fc_layers(z)
        coords = coords_flat.view(-1, self.n_atoms, 3)
        return coords

class ProteinVAE(nn.Module):
    """
    Complete Variational Autoencoder for protein structures.
    """
    def __init__(self, node_feat_dim: int, hidden_dim: int, latent_dim: int, 
                 num_layers: int, n_atoms: int):
        super().__init__()
        
        self.encoder = GNNEncoder(node_feat_dim, hidden_dim, latent_dim, num_layers)
        self.decoder = MLPDecoder(latent_dim, hidden_dim, n_atoms)
        self.latent_dim = latent_dim
    
    def reparameterize(self, mu: torch.Tensor, logvar: torch.Tensor) -> torch.Tensor:
        """
        Reparameterization trick: z = mu + std * epsilon
        """
        std = torch.exp(0.5 * logvar)
        eps = torch.randn_like(std)
        return mu + eps * std
    
    def forward(self, g: dgl.DGLGraph) -> Tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
        """
        Full forward pass: encode -> sample -> decode
        
        Returns:
            recon_coords: Reconstructed coordinates
            mu: Latent mean
            logvar: Latent log variance
        """
        # Encode
        mu, logvar = self.encoder(g)
        
        # Sample from latent distribution
        z = self.reparameterize(mu, logvar)
        
        # Decode
        recon_coords = self.decoder(z)
        
        return recon_coords, mu, logvar
    
    def encode_to_latent(self, g: dgl.DGLGraph) -> torch.Tensor:
        """Encode a structure to its mean latent representation (no sampling)."""
        mu, _ = self.encoder(g)
        return mu
    
    def decode_from_latent(self, z: torch.Tensor) -> torch.Tensor:
        """Decode latent vector to coordinates."""
        return self.decoder(z)

# -------------------------------------------------------------------
# Loss Function
# -------------------------------------------------------------------

def vae_loss(recon_coords: torch.Tensor, target_coords: torch.Tensor, 
             mu: torch.Tensor, logvar: torch.Tensor, kl_weight: float = 1.0) -> Tuple[torch.Tensor, dict]:
    """
    VAE loss = Reconstruction Loss + KL Divergence
    
    Args:
        recon_coords: Reconstructed coordinates [batch, n_atoms, 3]
        target_coords: Target coordinates [batch, n_atoms, 3]
        mu: Latent mean
        logvar: Latent log variance
        kl_weight: Weight for KL term (beta-VAE)
    
    Returns:
        Total loss and dictionary of components
    """
    # Reconstruction loss (MSE on coordinates)
    recon_loss = F.mse_loss(recon_coords, target_coords, reduction='sum')
    
    # KL divergence: KL(N(mu, sigma) || N(0, 1))
    kl_loss = -0.5 * torch.sum(1 + logvar - mu.pow(2) - logvar.exp())
    
    # Total loss
    total_loss = recon_loss + kl_weight * kl_loss
    
    # Return components for logging
    loss_dict = {
        'total': total_loss.item(),
        'recon': recon_loss.item(),
        'kl': kl_loss.item()
    }
    
    return total_loss, loss_dict

# -------------------------------------------------------------------
# Training Loop
# -------------------------------------------------------------------

def train_vae(model: ProteinVAE, graphs: List[dgl.DGLGraph], 
              num_epochs: int, lr: float, kl_weight: float, device: str) -> ProteinVAE:
    """
    Train the VAE on the two endpoint structures.
    """
    model = model.to(device)
    optimizer = optim.Adam(model.parameters(), lr=lr)
    
    print("\n" + "="*70)
    print("TRAINING VAE")
    print("="*70)
    
    # Move graphs to device
    graphs = [g.to(device) for g in graphs]
    target_coords = [g.ndata['pos'] for g in graphs]
    
    for epoch in range(num_epochs):
        model.train()
        epoch_losses = []
        
        # Train on both structures
        for i, (g, target) in enumerate(zip(graphs, target_coords)):
            optimizer.zero_grad()
            
            # Forward pass
            recon_coords, mu, logvar = model(g)
            
            # Compute loss
            target_batch = target.unsqueeze(0)  # Add batch dimension
            loss, loss_dict = vae_loss(recon_coords, target_batch, mu, logvar, kl_weight)
            
            # Backward pass
            loss.backward()
            optimizer.step()
            
            epoch_losses.append(loss_dict)
        
        # Logging
        if (epoch + 1) % 50 == 0 or epoch == 0:
            avg_loss = np.mean([l['total'] for l in epoch_losses])
            avg_recon = np.mean([l['recon'] for l in epoch_losses])
            avg_kl = np.mean([l['kl'] for l in epoch_losses])
            
            print(f"Epoch {epoch+1:4d}/{num_epochs} | "
                  f"Loss: {avg_loss:8.2f} | "
                  f"Recon: {avg_recon:8.2f} | "
                  f"KL: {avg_kl:6.2f}")
    
    print("\nTraining complete!")
    return model

# -------------------------------------------------------------------
# Trajectory Generation via Latent Space Interpolation
# -------------------------------------------------------------------

def slerp(z0: torch.Tensor, z1: torch.Tensor, t: float) -> torch.Tensor:
    """
    Spherical Linear Interpolation (SLERP) between two vectors.
    More appropriate than linear interpolation for latent spaces.
    
    Args:
        z0, z1: Start and end vectors
        t: Interpolation parameter [0, 1]
    """
    # Normalize vectors
    z0_norm = F.normalize(z0, dim=-1)
    z1_norm = F.normalize(z1, dim=-1)
    
    # Compute angle
    dot = (z0_norm * z1_norm).sum(dim=-1, keepdim=True)
    dot = torch.clamp(dot, -1.0, 1.0)
    omega = torch.acos(dot)
    
    # Compute interpolation
    so = torch.sin(omega)
    
    # Handle parallel vectors
    if so.abs() < 1e-6:
        return (1.0 - t) * z0 + t * z1
    
    return (torch.sin((1.0 - t) * omega) / so) * z0 + (torch.sin(t * omega) / so) * z1

def generate_trajectory(model: ProteinVAE, graph_A: dgl.DGLGraph, graph_B: dgl.DGLGraph,
                       num_frames: int, interpolation: str, device: str) -> np.ndarray:
    """
    Generate trajectory by interpolating in latent space.
    
    Returns:
        coords: [num_frames, n_atoms, 3] array of coordinates
    """
    model.eval()
    
    print("\n" + "="*70)
    print("GENERATING TRAJECTORY")
    print("="*70)
    
    # Encode endpoints to latent space
    with torch.no_grad():
        z_A = model.encode_to_latent(graph_A.to(device))
        z_B = model.encode_to_latent(graph_B.to(device))
    
    print(f"Latent vector A: {z_A.shape}, norm: {torch.norm(z_A).item():.3f}")
    print(f"Latent vector B: {z_B.shape}, norm: {torch.norm(z_B).item():.3f}")
    print(f"Latent distance: {torch.norm(z_B - z_A).item():.3f}")
    
    # Generate interpolation path
    print(f"\nInterpolating {num_frames} frames using {interpolation} interpolation...")
    
    coords_list = []
    
    with torch.no_grad():
        for i in range(num_frames):
            t = i / (num_frames - 1)  # Parameter from 0 to 1
            
            # Interpolate in latent space
            if interpolation == "linear":
                z_t = (1 - t) * z_A + t * z_B
            elif interpolation == "slerp":
                z_t = slerp(z_A, z_B, t)
            else:
                raise ValueError(f"Unknown interpolation: {interpolation}")
            
            # Decode to coordinates
            coords_t = model.decode_from_latent(z_t)
            coords_list.append(coords_t.cpu().numpy()[0])  # Remove batch dimension
            
            if (i + 1) % 50 == 0:
                print(f"  Generated frame {i+1}/{num_frames}")
    
    coords = np.array(coords_list)
    print(f"\nGenerated trajectory shape: {coords.shape}")
    
    return coords

# -------------------------------------------------------------------
# Main Execution
# -------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description="Generate VAE transition trajectory")
    parser.add_argument("--epochs", type=int, default=NUM_EPOCHS, help="Training epochs")
    parser.add_argument("--frames", type=int, default=NUM_FRAMES, help="Output frames")
    parser.add_argument("--latent-dim", type=int, default=LATENT_DIM, help="Latent dimension")
    parser.add_argument("--device", type=str, default="cuda" if torch.cuda.is_available() else "cpu")
    args = parser.parse_args()
    
    print("\n" + "="*70)
    print("VAE TRAJECTORY GENERATOR")
    print("="*70)
    print(f"Device: {args.device}")
    print(f"Latent dimension: {args.latent_dim}")
    print(f"Training epochs: {args.epochs}")
    print(f"Output frames: {args.frames}")
    
    # Create output directory
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    
    # 1. Load structures
    print("\n" + "="*70)
    print("LOADING STRUCTURES")
    print("="*70)
    
    traj_A = load_structure(STATE_A_PDB)
    traj_B = load_structure(STATE_B_PDB)
    
    # 2. Align structures
    print("\n" + "="*70)
    print("ALIGNING STRUCTURES")
    print("="*70)
    
    traj_A, traj_B = align_structures(traj_A, traj_B)
    
    # 3. Build graphs
    print("\n" + "="*70)
    print("BUILDING GRAPHS")
    print("="*70)
    
    graphs, n_atoms = prepare_dataset(traj_A, traj_B, EDGE_CUTOFF)
    
    # 4. Initialize model
    print("\n" + "="*70)
    print("INITIALIZING MODEL")
    print("="*70)
    
    model = ProteinVAE(
        node_feat_dim=3,  # x, y, z coordinates
        hidden_dim=HIDDEN_DIM,
        latent_dim=args.latent_dim,
        num_layers=NUM_GNN_LAYERS,
        n_atoms=n_atoms
    )
    
    total_params = sum(p.numel() for p in model.parameters())
    print(f"Model parameters: {total_params:,}")
    
    # 5. Train VAE
    model = train_vae(
        model=model,
        graphs=graphs,
        num_epochs=args.epochs,
        lr=LEARNING_RATE,
        kl_weight=KL_WEIGHT,
        device=args.device
    )
    
    # 6. Generate trajectory
    coords = generate_trajectory(
        model=model,
        graph_A=graphs[0],
        graph_B=graphs[1],
        num_frames=args.frames,
        interpolation=INTERPOLATION_MODE,
        device=args.device
    )
    
    # 7. Save trajectory
    print("\n" + "="*70)
    print("SAVING TRAJECTORY")
    print("="*70)
    
    # Create MDTraj trajectory with original topology
    output_traj = md.Trajectory(xyz=coords, topology=traj_A.topology)
    
    # Save DCD
    output_traj.save_dcd(str(OUTPUT_TRAJ_DCD))
    print(f"Saved DCD: {OUTPUT_TRAJ_DCD}")
    
    # Save topology PDB (first frame)
    output_traj[0].save_pdb(str(OUTPUT_TRAJ_PDB))
    print(f"Saved topology PDB: {OUTPUT_TRAJ_PDB}")
    
    # 8. Validate trajectory
    print("\n" + "="*70)
    print("TRAJECTORY VALIDATION")
    print("="*70)
    
    # Compute RMSD along trajectory
    rmsd_from_A = md.rmsd(output_traj, traj_A)
    rmsd_from_B = md.rmsd(output_traj, traj_B)
    
    print(f"RMSD from State A: {rmsd_from_A[0]:.3f} nm → {rmsd_from_A[-1]:.3f} nm")
    print(f"RMSD from State B: {rmsd_from_B[0]:.3f} nm → {rmsd_from_B[-1]:.3f} nm")
    
    # Check state definitions from run_path_sampling.py
    print("\nChecking against state definitions...")
    
    # CV1: lower_tm_contact (residues 200-310) 
    # Updated to use residues that exist in BOTH structures
    # Thresholds calibrated to 10th/90th percentiles
    cv1_pairs = output_traj.topology.select_pairs("resid 200 and name CA", "resid 310 and name CA")
    if len(cv1_pairs) > 0:
        cv1_distances = md.compute_distances(output_traj, cv1_pairs)[:, 0]
        print(f"CV1 (lower_tm_contact, res 200-310): {cv1_distances.min():.3f} - {cv1_distances.max():.3f} nm")
        print(f"  State A threshold: < 2.77 nm")
        print(f"  State B threshold: > 2.95 nm")
        
        n_state_a = np.sum(cv1_distances < 2.77)
        n_state_b = np.sum(cv1_distances > 2.95)
        n_transition = args.frames - n_state_a - n_state_b
        print(f"  Frames in State A: {n_state_a}/{args.frames} ({100*n_state_a/args.frames:.1f}%)")
        print(f"  Frames in State B: {n_state_b}/{args.frames} ({100*n_state_b/args.frames:.1f}%)")
        print(f"  Frames in transition: {n_transition}/{args.frames} ({100*n_transition/args.frames:.1f}%)")
    else:
        print(f"CV1 (lower_tm_contact): Could not find residues 200 and/or 310")
    
    # CV2: tm_core_distance (residues 150-310)
    cv2_pairs = output_traj.topology.select_pairs("resid 150 and name CA", "resid 310 and name CA")
    if len(cv2_pairs) > 0:
        cv2_distances = md.compute_distances(output_traj, cv2_pairs)[:, 0]
        print(f"\nCV2 (tm_core_distance, res 150-310): {cv2_distances.min():.3f} - {cv2_distances.max():.3f} nm")
        print(f"  State A threshold: < 2.58 nm")
        print(f"  State B threshold: > 2.78 nm")
        
        n_state_a = np.sum(cv2_distances < 2.58)
        n_state_b = np.sum(cv2_distances > 2.78)
        n_transition = args.frames - n_state_a - n_state_b
        print(f"  Frames in State A: {n_state_a}/{args.frames} ({100*n_state_a/args.frames:.1f}%)")
        print(f"  Frames in State B: {n_state_b}/{args.frames} ({100*n_state_b/args.frames:.1f}%)")
        print(f"  Frames in transition: {n_transition}/{args.frames} ({100*n_transition/args.frames:.1f}%)")
    else:
        print(f"CV2 (tm_core_distance): Could not find residues 150 and/or 310")
    
    # Combined state classification
    if len(cv1_pairs) > 0 and len(cv2_pairs) > 0:
        both_a = np.logical_and(cv1_distances < 2.77, cv2_distances < 2.58)
        both_b = np.logical_and(cv1_distances > 2.95, cv2_distances > 2.78)
        print(f"\nCombined Classification (both CVs must agree):")
        print(f"  Both in State A: {np.sum(both_a)}/{args.frames} ({100*np.sum(both_a)/args.frames:.1f}%)")
        print(f"  Both in State B: {np.sum(both_b)}/{args.frames} ({100*np.sum(both_b)/args.frames:.1f}%)")
        print(f"  Transition region: {args.frames - np.sum(both_a) - np.sum(both_b)}/{args.frames} ({100*(args.frames - np.sum(both_a) - np.sum(both_b))/args.frames:.1f}%)")
        
        if np.sum(both_a) > 0 and np.sum(both_b) > 0:
            print(f"\n✅ SUCCESS: Trajectory spans BOTH states!")
        else:
            print(f"\n⚠️  WARNING: Trajectory may not span both states cleanly")
    
    print("\n" + "="*70)
    print("COMPLETE!")
    print("="*70)
    print(f"\nYou can now use this trajectory for path sampling:")
    print(f"  python scripts/run_path_sampling.py \\")
    print(f"    --target glp1r_inactive \\")
    print(f"    --initial-traj {OUTPUT_TRAJ_DCD} \\")
    print(f"    --shots 100")
    print("="*70)

if __name__ == "__main__":
    main()
