# scripts/prepare_ligands.py
import argparse
import sys
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors
from rdkit import RDLogger

# Suppress RDKit warnings
RDLogger.DisableLog('rdApp.*')

def is_drug_like(mol):
    """Checks if a molecule satisfies Lipinski's Rule of Five."""
    # MolWt <= 500
    if Descriptors.MolWt(mol) > 500:
        return False
    # LogP <= 5
    if Descriptors.MolLogP(mol) > 5:
        return False
    # H-bond donors <= 5
    if Descriptors.NumHDonors(mol) > 5:
        return False
    # H-bond acceptors <= 10
    if Descriptors.NumHAcceptors(mol) > 10:
        return False
    return True

def prepare_ligand_library(input_smi_file: Path, output_sdf_file: Path):
    """
    Reads a SMILES file, filters for drug-like molecules, generates 3D conformers,
    and performs energy minimization.
    
    Args:
        input_smi_file: Path to input SMILES file
        output_sdf_file: Path to output SDF file
        
    Returns:
        bool: True if successful, False if failed
    """
    print("Starting Ligand Library Preparation")
    print(f"Input: {input_smi_file}")
    print(f"Output: {output_sdf_file}")
    
    # Check if input file exists
    if not input_smi_file.exists():
        print(f"ERROR: Input file not found: {input_smi_file}")
        return False
    
    # Create output directory if needed
    output_sdf_file.parent.mkdir(parents=True, exist_ok=True)
    
    # RDKit's SDWriter can write to a file handle
    writer = Chem.SDWriter(str(output_sdf_file))
    
    count_initial = 0
    count_druglike = 0
    count_success = 0
    count_failed_parse = 0
    count_failed_embed = 0

    print(f"Reading molecules from: {input_smi_file}")
    try:
        with open(input_smi_file, 'r') as f:
            for line_num, line in enumerate(f, 1):
                line = line.strip()
                if not line or line.startswith('#'):  # Skip empty lines and comments
                    continue
                    
                count_initial += 1
                parts = line.split()
                smiles = parts[0]
                mol_name = parts[1] if len(parts) > 1 else f"mol_{count_initial}"
                
                # Progress indicator
                if count_initial % 100 == 0:
                    print(f"Processed {count_initial} molecules...")
                
                mol = Chem.MolFromSmiles(smiles)
                
                if mol is None:
                    count_failed_parse += 1
                    continue

                # STEP 1: Filter for drug-like properties
                if not is_drug_like(mol):
                    continue
                count_druglike += 1

                # STEP 2: Generate 3D structure and add hydrogens
                mol_h = Chem.AddHs(mol)
                
                # Try to embed molecule (generate 3D coordinates)
                embed_result = AllChem.EmbedMolecule(mol_h, randomSeed=42)
                if embed_result == -1:  # Embedding failed
                    count_failed_embed += 1
                    continue

                # STEP 3: Energy Minimization
                try:
                    AllChem.MMFFOptimizeMolecule(mol_h, maxIters=200)
                    # Set a name for the molecule
                    mol_h.SetProp("_Name", mol_name)
                    mol_h.SetProp("SMILES", smiles)
                    writer.write(mol_h)
                    count_success += 1
                except Exception:
                    # Minimization can fail for some complex structures
                    continue
                    
    except FileNotFoundError:
        print(f"ERROR: Could not read input file: {input_smi_file}")
        return False
    except Exception as e:
        print(f"ERROR: Unexpected error during processing: {e}")
        return False
    finally:
        writer.close()
    print("\nLigand Preparation Complete!")
    print(f"Initial molecules: {count_initial}")
    print(f"Failed to parse: {count_failed_parse}")
    print(f"Passed drug-like filter: {count_druglike}")
    print(f"Failed 3D embedding: {count_failed_embed}")
    print(f"Successfully processed and saved: {count_success}")
    if count_initial > 0:
        print(f"Success rate: {count_success/count_initial*100:.1f}%")
    print(f"Output saved to: {output_sdf_file}")
    
    return count_success > 0

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Prepare a ligand library for virtual screening.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example usage:
  python scripts/prepare_ligands.py -i ligands.smi -o prepared_ligands.sdf
  
The input SMILES file should have one SMILES per line, optionally followed by a molecule name:
  CCO ethanol
  CC(=O)O acetic_acid
        """
    )
    parser.add_argument("--input", "-i", required=True, 
                       help="Input SMILES file (.smi or .txt)")
    parser.add_argument("--output", "-o", required=True, 
                       help="Output SDF file (.sdf)")
    
    args = parser.parse_args()

    success = prepare_ligand_library(Path(args.input), Path(args.output))
    
    if not success:
        print("\nERROR: Ligand preparation failed!")
        sys.exit(1)
    
    print("\nLigand preparation completed successfully!")
    print("Ready for virtual screening with AutoDock Vina.")