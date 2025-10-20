#!/usr/bin/env python3
"""
Energy minimize ligands before virtual screening.

This script performs energy minimization on PDBQT ligands using RDKit
to ensure they're in low-energy, physically realistic conformations.
"""

import argparse
import os
import sys
from pathlib import Path
import logging
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors
from rdkit.Chem import rdForceFieldHelpers
import time

def setup_logging(verbose=False):
    """Set up logging configuration."""
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format='%(asctime)s - %(levelname)s - %(message)s',
        datefmt='%H:%M:%S'
    )

def minimize_ligand_energy(mol, max_iterations=1000, convergence_threshold=1e-6):
    """
    Minimize ligand energy using MMFF force field.
    
    Args:
        mol: RDKit molecule object
        max_iterations: Maximum optimization iterations
        convergence_threshold: Convergence criteria
        
    Returns:
        tuple: (minimized_mol, final_energy, converged)
    """
    if mol is None:
        return None, None, False
    
    # Add hydrogens if not present
    mol = Chem.AddHs(mol)
    
    # Generate 3D coordinates if needed
    if mol.GetNumConformers() == 0:
        AllChem.EmbedMolecule(mol, randomSeed=42)
        if mol.GetNumConformers() == 0:
            logging.warning("Failed to generate 3D coordinates")
            return None, None, False
    
    # Set up MMFF force field
    mp = AllChem.MMFFGetMoleculeProperties(mol)
    if mp is None:
        logging.warning("Could not get MMFF properties, trying UFF")
        # Fallback to UFF force field
        try:
            energy_before = AllChem.UFFGetMoleculeForceField(mol).CalcEnergy()
            AllChem.UFFOptimizeMolecule(mol, maxIters=max_iterations)
            energy_after = AllChem.UFFGetMoleculeForceField(mol).CalcEnergy()
            converged = True  # UFF doesn't return convergence info
            return mol, energy_after, converged
        except Exception:
            logging.exception("Both MMFF and UFF failed")
            return mol, None, False
    
    # MMFF optimization 
    ff = AllChem.MMFFGetMoleculeForceField(mol, mp)
    if ff is None:
        logging.warning("Could not create MMFF force field")
        return mol, None, False
    
    energy_before = ff.CalcEnergy()
    logging.debug(f"Initial energy: {energy_before:.2f} kcal/mol")
    
    # Perform minimization
    converged = ff.Minimize(maxIts=max_iterations, energyTol=convergence_threshold)
    energy_after = ff.CalcEnergy()
    
    logging.debug(f"Final energy: {energy_after:.2f} kcal/mol")
    logging.debug(f"Energy change: {energy_after - energy_before:.2f} kcal/mol")
    logging.debug(f"Converged: {converged == 0}")  # 0 = converged
    
    return mol, energy_after, converged == 0

def pdbqt_to_mol(pdbqt_file):
    """Convert PDBQT file to RDKit molecule."""
    # Read PDBQT as PDB (RDKit can handle this)
    mol = Chem.MolFromPDBFile(str(pdbqt_file), removeHs=False)
    if mol is None:
        # Try reading as MOL file
        with open(pdbqt_file, 'r') as f:
            content = f.read()
        # Remove PDBQT-specific lines
        cleaned_content = []
        for line in content.split('\n'):
            if not line.startswith(('ROOT', 'ENDROOT', 'BRANCH', 'ENDBRANCH', 'TORSDOF')):
                cleaned_content.append(line)
        
        # Write temporary PDB file
        temp_pdb = pdbqt_file.with_suffix('.temp.pdb')
        with open(temp_pdb, 'w') as f:
            f.write('\n'.join(cleaned_content))
        
        mol = Chem.MolFromPDBFile(str(temp_pdb), removeHs=False)
        temp_pdb.unlink()  # Clean up
        
    return mol

def mol_to_pdbqt(mol, output_file, original_pdbqt=None):
    """Convert RDKit molecule back to PDBQT format."""
    from subprocess import run, PIPE
    
    # Write as SDF first
    temp_sdf = output_file.with_suffix('.temp.sdf')
    writer = Chem.SDWriter(str(temp_sdf))
    writer.write(mol)
    writer.close()
    
    # Convert to PDBQT using OpenBabel
    cmd = ['obabel', str(temp_sdf), '-O', str(output_file)]
    
    try:
        result = run(cmd, capture_output=True, text=True, check=True)
        temp_sdf.unlink()  # Clean up
        return True
    except Exception as e:
        logging.error(f"Failed to convert to PDBQT: {e}")
        if temp_sdf.exists():
            temp_sdf.unlink()
        return False

def minimize_pdbqt_file(pdbqt_file, output_file=None, verbose=False):
    """
    Minimize energy of a single PDBQT file.
    
    Args:
        pdbqt_file: Path to input PDBQT file
        output_file: Path to output file (default: overwrite input)
        verbose: Enable verbose output
        
    Returns:
        dict: Results dictionary with success status, energies, etc.
    """
    pdbqt_path = Path(pdbqt_file)
    if output_file is None:
        output_file = pdbqt_path
    else:
        output_file = Path(output_file)
    
    logging.info(f"Minimizing {pdbqt_path.name}")
    
    start_time = time.time()
    
    # Convert PDBQT to RDKit molecule
    mol = pdbqt_to_mol(pdbqt_path)
    if mol is None:
        logging.error(f"Could not read molecule from {pdbqt_file}")
        return {"success": False, "error": "Could not read molecule"}
    
    logging.debug(f"Loaded molecule with {mol.GetNumAtoms()} atoms")
    
    # Perform energy minimization
    min_mol, final_energy, converged = minimize_ligand_energy(mol)
    
    if min_mol is None:
        logging.error("Energy minimization failed")
        return {"success": False, "error": "Minimization failed"}
    
    # Convert back to PDBQT
    success = mol_to_pdbqt(min_mol, output_file, pdbqt_path)
    
    elapsed_time = time.time() - start_time
    
    if success:
        logging.info(f"✅ Minimized in {elapsed_time:.1f}s")
        if final_energy is not None:
            logging.info(f"   Final energy: {final_energy:.2f} kcal/mol")
        logging.info(f"   Converged: {'Yes' if converged else 'No'}")
        
        return {
            "success": True,
            "final_energy": final_energy,
            "converged": converged,
            "time": elapsed_time,
            "output_file": str(output_file)
        }
    else:
        logging.error("Failed to write minimized structure")
        return {"success": False, "error": "Could not write output"}

def minimize_multiple_files(input_pattern, output_dir=None, verbose=False):
    """Minimize multiple PDBQT files."""
    from glob import glob
    
    files = glob(input_pattern)
    if not files:
        logging.warning(f"No files found matching: {input_pattern}")
        return []
    
    logging.info(f"Found {len(files)} files to minimize")
    
    if output_dir:
        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)
    
    results = []
    successful = 0
    
    for i, pdbqt_file in enumerate(files, 1):
        file_path = Path(pdbqt_file)
        
        if output_dir:
            output_file = output_path / file_path.name
        else:
            output_file = file_path
        
        print(f"\n--- Minimizing {i}/{len(files)}: {file_path.name} ---")
        
        result = minimize_pdbqt_file(pdbqt_file, output_file, verbose)
        result["input_file"] = str(file_path)
        results.append(result)
        
        if result["success"]:
            successful += 1
        
    print("\n📊 Minimization Summary:")
    print(f"   Successful: {successful}/{len(files)}")
    print(f"   Failed: {len(files) - successful}/{len(files)}")
    
    return results

def main():
    parser = argparse.ArgumentParser(
        description="Energy minimize ligands for virtual screening",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Minimize single ligand
  python minimize_ligands.py ligand.pdbqt
  
  # Minimize all ligands in directory
  python minimize_ligands.py "ligands/*.pdbqt" -o minimized_ligands/
  
  # Minimize with verbose output
  python minimize_ligands.py ligand.pdbqt --verbose
        """
    )
    
    parser.add_argument(
        'input',
        help='Input PDBQT file or glob pattern (e.g., "*.pdbqt")'
    )
    
    parser.add_argument(
        '-o', '--output-dir',
        help='Output directory for minimized files (default: overwrite input)'
    )
    
    parser.add_argument(
        '-v', '--verbose',
        action='store_true',
        help='Enable verbose output'
    )
    
    args = parser.parse_args()
    
    setup_logging(args.verbose)
    
    try:
        # Check if RDKit is available
        logging.info(f"Using RDKit version: {Chem.__version__}")
    except Exception:
        logging.exception("RDKit not available. Please install it:")
        logging.error("  conda install -c conda-forge rdkit")
        sys.exit(1)
    
    try:
        if '*' in args.input or '?' in args.input:
            # Multiple files
            results = minimize_multiple_files(args.input, args.output_dir, args.verbose)
            
            successful_files = [r for r in results if r["success"]]
            if successful_files:
                print(f"\nSuccessfully minimized {len(successful_files)} files:")
                for result in successful_files:
                    energy_str = f" ({result['final_energy']:.2f} kcal/mol)" if result['final_energy'] else ""
                    print(f"  {Path(result['input_file']).name}{energy_str}")
            
        else:
            # Single file
            if not Path(args.input).exists():
                logging.error(f"File not found: {args.input}")
                sys.exit(1)
            
            output_file = None
            if args.output_dir:
                output_dir = Path(args.output_dir)
                output_dir.mkdir(parents=True, exist_ok=True)
                output_file = output_dir / Path(args.input).name
            
            result = minimize_pdbqt_file(args.input, output_file, args.verbose)
            
            if result["success"]:
                print("\n✅ Minimization successful!")
                if result["final_energy"]:
                    print(f"   Final energy: {result['final_energy']:.2f} kcal/mol")
                print(f"   Output: {result['output_file']}")
            else:
                print(f"\n❌ Minimization failed: {result.get('error', 'Unknown error')}")
                sys.exit(1)
                
    except Exception as e:
        logging.error(f"Minimization failed: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()