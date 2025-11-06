#!/usr/bin/env python3
"""
Align 6X18 active structure to Boltz inactive structure.
Extract only the overlapping residues to ensure consistency.
"""

import mdtraj as md
import numpy as np
from pathlib import Path

# Paths
boltz_inactive = "boltz_results_glp1r/predictions/glp1r/glp1r_model_0.pdb"
active_pdb = "6X18_fresh.pdb"
output_cif = "boltz_results_glp1r/predictions/glp1r_active/glp1r_active_model_0.cif"

print("=" * 80)
print("Aligning 6X18 Active Structure to Boltz Inactive Structure")
print("=" * 80)

# Load structures
print("\n1. Loading structures...")
print(f"   Boltz inactive: {boltz_inactive}")
inactive = md.load(boltz_inactive)
print(f"   → Loaded: {inactive.n_residues} residues, {inactive.n_atoms} atoms")

print(f"   6X18 active: {active_pdb}")
active_full = md.load(active_pdb)
print(f"   → Loaded: {active_full.n_residues} residues, {active_full.n_atoms} atoms")

# Get residue ranges
print("\n2. Analyzing residue ranges...")
inactive_residues = [r.resSeq for r in inactive.topology.residues]
inactive_min, inactive_max = min(inactive_residues), max(inactive_residues)
print(f"   Boltz inactive: residues {inactive_min}-{inactive_max}")

# Extract chain R (receptor) from 6X18
print("\n3. Extracting chain R (receptor) from 6X18...")
chain_r_atoms = active_full.topology.select("chainid 0")  # First chain is usually receptor
active_chain_r = active_full.atom_slice(chain_r_atoms)
active_residues = [r.resSeq for r in active_chain_r.topology.residues]
active_min, active_max = min(active_residues), max(active_residues)
print(f"   6X18 chain R: residues {active_min}-{active_max}")

# Find overlapping residues
overlap_min = max(inactive_min, active_min)
overlap_max = min(inactive_max, active_max)
print(f"\n4. Finding overlapping region: {overlap_min}-{overlap_max}")

# Extract overlapping region from both structures
print("\n5. Extracting overlapping residues from both structures...")

# From inactive
inactive_overlap_atoms = []
for atom in inactive.topology.atoms:
    if overlap_min <= atom.residue.resSeq <= overlap_max:
        inactive_overlap_atoms.append(atom.index)
inactive_overlap = inactive.atom_slice(inactive_overlap_atoms)
print(f"   Boltz overlap: {inactive_overlap.n_residues} residues, {inactive_overlap.n_atoms} atoms")

# From active
active_overlap_atoms = []
for atom in active_chain_r.topology.atoms:
    if overlap_min <= atom.residue.resSeq <= overlap_max:
        active_overlap_atoms.append(atom.index)
active_overlap = active_chain_r.atom_slice(active_overlap_atoms)
print(f"   6X18 overlap: {active_overlap.n_residues} residues, {active_overlap.n_atoms} atoms")

# Align active to inactive using CA atoms
print("\n6. Aligning structures using CA atoms...")
inactive_ca = inactive_overlap.topology.select("name CA")
active_ca = active_overlap.topology.select("name CA")

if len(inactive_ca) != len(active_ca):
    print(f"   WARNING: CA atom count mismatch!")
    print(f"   Boltz: {len(inactive_ca)} CA atoms")
    print(f"   6X18: {len(active_ca)} CA atoms")
    print("   Using minimum common CA atoms for alignment...")
    min_ca = min(len(inactive_ca), len(active_ca))
    inactive_ca = inactive_ca[:min_ca]
    active_ca = active_ca[:min_ca]

print(f"   Aligning {len(inactive_ca)} CA atoms...")
active_overlap.superpose(inactive_overlap, atom_indices=active_ca, ref_atom_indices=inactive_ca)

# Calculate RMSD
rmsd = md.rmsd(active_overlap, inactive_overlap, atom_indices=active_ca, ref_atom_indices=inactive_ca)[0] * 10  # nm to Angstrom
print(f"   RMSD after alignment: {rmsd:.2f} Å")

# Save aligned active structure
print(f"\n7. Saving aligned active structure...")
Path(output_cif).parent.mkdir(parents=True, exist_ok=True)

# Save as PDB first (more reliable than CIF)
output_pdb = output_cif.replace('.cif', '.pdb')
active_overlap.save_pdb(output_pdb)
print(f"   Saved PDB: {output_pdb}")

# Now clean it with PDBFixer and save as CIF
print(f"\n8. Cleaning with PDBFixer and adding hydrogens...")
from pdbfixer import PDBFixer
from openmm.app import PDBFile, PDBxFile

fixer = PDBFixer(output_pdb)
fixer.findMissingResidues()
fixer.findMissingAtoms()
fixer.missingResidues = {}  # Don't add terminal residues
fixer.addMissingAtoms()
fixer.addMissingHydrogens(7.0)

# Check for NaN coordinates
import numpy as np
pos_array = np.array([[atom.x, atom.y, atom.z] for atom in fixer.positions])
if np.any(np.isnan(pos_array)):
    print("   ⚠️ WARNING: Found NaN coordinates after PDBFixer!")
    print("   Saving without PDBFixer modifications...")
    active_overlap.save_cif(output_cif)
else:
    print("   ✅ No NaN coordinates detected")
    # Save as CIF with hydrogens
    with open(output_cif, 'w') as f:
        PDBxFile.writeFile(fixer.topology, fixer.positions, f)
    print(f"   Saved CIF: {output_cif}")

print("\n" + "=" * 80)
print("✅ SUCCESS!")
print("=" * 80)
print(f"Aligned active structure saved with:")
print(f"  - Residues: {active_overlap.n_residues}")
print(f"  - Atoms: {active_overlap.n_atoms}")
print(f"  - RMSD to Boltz: {rmsd:.2f} Å")
print(f"  - Residue range: {overlap_min}-{overlap_max}")
print("\nYou can now run:")
print("  python scripts/prepare_simulation.py --target glp1r_active --boltz")
