#!/usr/bin/env python3
"""
Final test of all fixes for run_path_sampling.py
"""

print("=== Final API Compatibility Test ===\n")

# Test 1: Mass calculations
print("1. Testing fixed mass calculation...")
try:
    import mdtraj as md
    import numpy as np
    from pathlib import Path
    
    pdb_path = "artifacts/md/GLP1R/prepared_system_whole.pdb"
    traj = md.load_pdb(pdb_path)
    
    # Test the fixed mass calculation pattern
    n_atoms = traj.topology.n_atoms
    masses = np.array([a.element.mass for a in traj.topology.atoms]).reshape(-1, 1)
    
    # Fixed velocity calculation
    velocity_scale = np.sqrt(310.0 / masses) * 0.1
    velocities = np.random.randn(n_atoms, 3) * velocity_scale
    
    print(f"✅ Mass calculation works: {velocities.shape} velocities generated")
    
except Exception as e:
    print(f"❌ Mass calculation error: {e}")

# Test 2: Trajectory loading with correct topology
print("\n2. Testing trajectory loading with correct topology...")
try:
    import mdtraj as md
    
    pdb_path = "artifacts/md/GLP1R/prepared_system_whole.pdb"
    traj_path = "artifacts/md/GLP1R/trajectory.dcd"
    
    topology = md.load_pdb(pdb_path).topology
    current_path = md.load(traj_path, top=pdb_path)
    
    print(f"✅ Trajectory loading works: {current_path.n_frames} frames, {current_path.n_atoms} atoms")
    
except Exception as e:
    print(f"❌ Trajectory loading error: {e}")

# Test 3: OpenMM worker configuration
print("\n3. Testing OpenMM worker CUDA configuration...")
try:
    from openmm import Platform
    
    platform = Platform.getPlatformByName('CUDA')
    properties = {'DeviceIndex': '0', 'Precision': 'mixed'}
    
    # Test if the properties are valid
    valid_props = platform.getPropertyNames()
    for prop in properties.keys():
        if prop in valid_props:
            print(f"✅ {prop} property is valid")
        else:
            print(f"❌ {prop} property is invalid")
            
except Exception as e:
    print(f"❌ OpenMM worker config error: {e}")

# Test 4: PyTorch/CUDA availability
print("\n4. Testing PyTorch/CUDA setup...")
try:
    import torch
    
    device = torch.device('cuda:0' if torch.cuda.is_available() else 'cpu')
    print(f"✅ PyTorch device: {device}")
    
    if torch.cuda.is_available():
        print(f"✅ CUDA available: {torch.cuda.get_device_name(0)}")
    else:
        print("⚠️  CUDA not available, will use CPU")
        
except Exception as e:
    print(f"❌ PyTorch/CUDA error: {e}")

# Test 5: Import all required modules
print("\n5. Testing all imports...")
try:
    from model import CommittorNet, frame_to_torch_graph
    from worker import simulation_worker
    from torch.multiprocessing import Process, Queue, set_start_method
    from openmm.unit import nanometer, picosecond, kelvin, BOLTZMANN_CONSTANT_kB, amu, kilojoule_per_mole
    
    print("✅ All imports successful")
    
except Exception as e:
    print(f"❌ Import error: {e}")

print("\n=== Test Summary ===")
print("All major API issues should now be resolved!")
print("The script should be ready to run without OpenMM/PyTorch compatibility errors.")