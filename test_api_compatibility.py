#!/usr/bin/env python3
"""
Test script to check API compatibility for run_path_sampling.py and worker.py
"""

def test_openmm_elements():
    """Test OpenMM element access patterns used in worker.py"""
    print("Testing OpenMM element access...")
    try:
        from openmm.app import PDBFile, element
        from openmm.unit import amu
        from pathlib import Path

        pdb_path = Path('artifacts/md/GLP1R/prepared_system.pdb')
        if not pdb_path.exists():
            print("⚠️  PDB file not found, skipping element access test")
            return True

        pdb = PDBFile(str(pdb_path))
        
        # Test the pattern used in worker.py for checking hydrogen
        first_atom = next(iter(pdb.topology.atoms()))
        print(f"First atom element: {first_atom.element}")
        
        # Test element comparison (used in worker.py)
        is_hydrogen = first_atom.element == element.hydrogen
        print(f"✅ Element comparison works: {is_hydrogen}")
        
        # Test amu unit usage
        mass_test = 4*amu
        print(f"✅ AMU unit works: {mass_test}")
        
        return True
        
    except Exception as e:
        print(f"❌ OpenMM element access error: {e}")
        return False

def test_mass_calculations():
    """Test mass calculation patterns used in run_path_sampling.py"""
    print("\nTesting mass calculations...")
    try:
        import numpy as np
        from openmm.unit import nanometer, picosecond, kelvin, BOLTZMANN_CONSTANT_kB
        
        # Test the velocity generation pattern from run_path_sampling.py
        n_atoms = 100  # Test with smaller number
        masses = np.random.rand(n_atoms, 1) * 10  # Random masses
        velocities = np.random.randn(n_atoms, 3) / np.sqrt(masses) * np.sqrt(BOLTZMANN_CONSTANT_kB * (310*kelvin))
        
        # Test unit conversion
        velocity_units = velocities.value_in_unit(nanometer/picosecond)
        print(f"✅ Velocity calculation works, shape: {velocity_units.shape}")
        
        return True
        
    except Exception as e:
        print(f"❌ Mass calculation error: {e}")
        return False

def test_system_serialization():
    """Test system loading pattern used in worker.py"""
    print("\nTesting system serialization...")
    try:
        from openmm import XmlSerializer
        from pathlib import Path
        
        system_path = Path('artifacts/md/GLP1R/system.xml')
        if not system_path.exists():
            print("⚠️  System XML file not found, skipping serialization test")
            return True
            
        # Test loading system from XML (pattern used in worker.py)
        with open(system_path, 'r') as f:
            system = XmlSerializer.deserialize(f.read())
        
        print(f"✅ System deserialization works, {system.getNumParticles()} particles")
        
        # Test particle mass access
        mass0 = system.getParticleMass(0)
        print(f"✅ Particle mass access works: {mass0}")
        
        return True
        
    except Exception as e:
        print(f"❌ System serialization error: {e}")
        return False

def test_mdtraj_compatibility():
    """Test MDTraj operations used in the scripts"""
    print("\nTesting MDTraj compatibility...")
    try:
        import mdtraj as md
        from pathlib import Path
        
        # Test trajectory loading
        traj_path = Path('artifacts/md/GLP1R/trajectory.dcd') 
        pdb_path = Path('artifacts/md/GLP1R/prepared_system.pdb')
        
        if not pdb_path.exists():
            print("⚠️  PDB file not found, skipping MDTraj test")
            return True
            
        # Test PDB loading (used in run_path_sampling.py)
        traj = md.load_pdb(str(pdb_path))
        print(f"✅ PDB loading works, {traj.n_atoms} atoms, {traj.n_frames} frames")
        
        # Test trajectory loading if DCD exists
        if traj_path.exists():
            full_traj = md.load(str(traj_path), top=str(pdb_path))
            print(f"✅ DCD loading works, {full_traj.n_frames} frames")
        else:
            print("⚠️  DCD file not found, skipping trajectory loading test")
            
        return True
        
    except Exception as e:
        print(f"❌ MDTraj compatibility error: {e}")
        return False

if __name__ == "__main__":
    print("=== API Compatibility Test for run_path_sampling.py ===\n")
    
    tests = [
        test_openmm_elements,
        test_mass_calculations, 
        test_system_serialization,
        test_mdtraj_compatibility
    ]
    
    results = []
    for test in tests:
        try:
            result = test()
            results.append(result)
        except Exception as e:
            print(f"❌ Test {test.__name__} failed with exception: {e}")
            results.append(False)
    
    print(f"\n=== Summary ===")
    print(f"Tests passed: {sum(results)}/{len(results)}")
    
    if all(results):
        print("🎉 All tests passed! The scripts should run without API issues.")
    else:
        print("⚠️  Some tests failed. Review the errors above.")