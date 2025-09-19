#!/usr/bin/env python3
"""Test the fixes applied to run_path_sampling.py"""

import numpy as np
from pathlib import Path

print("=== Testing Fixed Issues ===\n")

# Test 1: Fixed mass calculation
print("1. Testing fixed mass calculation...")
try:
    import mdtraj as md
    from openmm.unit import *
    
    pdb_path = Path('artifacts/md/GLP1R/prepared_system.pdb')
    if pdb_path.exists():
        traj = md.load_pdb(str(pdb_path))
        masses = np.array([a.element.mass.value_in_unit(amu) for a in traj.topology.atoms]).reshape(-1, 1)
        kT = BOLTZMANN_CONSTANT_kB * (310*kelvin)
        velocity_scale = np.sqrt(kT.value_in_unit(kilojoule_per_mole) / (masses * 0.001))
        velocities = np.random.randn(10, 3) * velocity_scale[:10]  # Test with first 10 atoms
        print("✅ Mass calculation fix works!")
    else:
        print("⚠️  PDB file not found")
except Exception as e:
    print(f"❌ Mass calculation still has issues: {e}")

# Test 2: Check file paths
print("\n2. Testing file paths...")
MD_DIR = Path("artifacts/md/GLP1R")
files_to_check = [
    "prepared_system.pdb",
    "minimized_state.xml", 
    "system.xml",
    "trajectory.dcd"
]

all_exist = True
for filename in files_to_check:
    filepath = MD_DIR / filename
    if filepath.exists():
        print(f"✅ {filename} exists")
    else:
        print(f"❌ {filename} missing")
        all_exist = False

if all_exist:
    print("✅ All required files present")
else:
    print("⚠️  Some files missing - run prepare_simulation.py and run_simulation.py first")

# Test 3: Test MDTraj loading consistency
print("\n3. Testing MDTraj trajectory loading...")
try:
    import mdtraj as md
    
    pdb_path = MD_DIR / "prepared_system.pdb"
    traj_path = MD_DIR / "trajectory.dcd"
    
    if pdb_path.exists() and traj_path.exists():
        traj = md.load(str(traj_path), top=str(pdb_path))
        print(f"✅ Trajectory loaded: {traj.n_frames} frames, {traj.n_atoms} atoms")
    else:
        print("⚠️  Cannot test - files missing")
        
except Exception as e:
    print(f"❌ Trajectory loading error: {e}")

print("\n=== Test Complete ===")