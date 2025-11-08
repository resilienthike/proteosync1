#!/usr/bin/env python3
"""
Reconstruct full-atom trajectory from VAE CA-only trajectory.
Maps each VAE frame to the full prepared system.
"""

import mdtraj as md
import numpy as np
from pathlib import Path

# Paths
VAE_TRAJ = "data/glp1r_inactive/vae_generated_transition.dcd"
VAE_TOP = "data/glp1r_inactive/vae_transition_topology.pdb"
FULL_TOP = "data/glp1r_inactive/prepared_system.pdb"
OUTPUT = "data/glp1r_inactive/initial_transition.dcd"

print("Loading VAE CA trajectory...")
vae_traj = md.load(VAE_TRAJ, top=VAE_TOP)
print(f"  {vae_traj.n_frames} frames, {vae_traj.n_atoms} CA atoms")

print("\nLoading full system topology...")
full_top = md.load(FULL_TOP)
print(f"  {full_top.n_atoms} atoms total")

# Extract CA atoms from full system that match VAE residues
print("\nExtracting matching CA atoms from full system...")
vae_residues = [r.resSeq for r in vae_traj.topology.residues]
full_ca_indices = []
for atom in full_top.topology.atoms:
    if atom.name == 'CA' and atom.residue.resSeq in vae_residues:
        full_ca_indices.append(atom.index)

print(f"  Found {len(full_ca_indices)} matching CA atoms in full system")

# Create full-atom trajectory by threading VAE coords through full system
print("\nReconstructing full-atom trajectory...")
full_frames = []

for i, vae_frame in enumerate(vae_traj):
    # Copy the full system
    new_frame = full_top.xyz[0].copy()
    
    # Update CA positions from VAE
    new_frame[full_ca_indices] = vae_frame.xyz[0]
    
    full_frames.append(new_frame)
    
    if (i + 1) % 50 == 0:
        print(f"  Processed {i+1}/{vae_traj.n_frames} frames")

# Create new trajectory
print("\nCreating full-atom trajectory...")
full_traj = md.Trajectory(np.array(full_frames), topology=full_top.topology)

# Save
print(f"\nSaving to {OUTPUT}...")
full_traj.save_dcd(OUTPUT)

print(f"\n✅ Done! Created {full_traj.n_frames} frames with {full_traj.n_atoms} atoms")
