#!/usr/bin/env python3
"""
fetch_structures.py

Download and Prepare Experimental GPCR Structures from RCSB PDB

Downloads PDB structures and extracts specific chains for analysis. Handles:
- Structure download from RCSB PDB database
- Chain extraction and cleaning
- Proper formatting for downstream MD preparation

Default structures:
- 6LN2: GLP-1R inactive state (no agonist)
- 6X18: GLP-1R active state (agonist + G protein bound)

Author: ProteOSync Pipeline
Date: November 8, 2025
"""

import argparse
import urllib.request
from pathlib import Path


def fetch_pdb_structure(pdb_id: str, output_dir: Path, chain: str = 'R'):
    """
    Download a PDB structure from RCSB and extract the specified chain.
    
    Args:
        pdb_id: 4-letter PDB identifier (e.g., '6LN2')
        output_dir: Directory to save the structure
        chain: Chain identifier to extract (default: 'R' for receptor)
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Download full PDB file
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    temp_file = output_dir / f"{pdb_id}_full.pdb"
    output_file = output_dir / f"{pdb_id}_chain{chain}.pdb"
    
    print(f"--> Downloading {pdb_id} from RCSB PDB...")
    print(f"    URL: {url}")
    
    try:
        urllib.request.urlretrieve(url, temp_file)
        print(f"✅ Downloaded to: {temp_file}")
    except Exception as e:
        print(f"❌ Failed to download {pdb_id}: {e}")
        return False
    
    # Extract specified chain with proper formatting
    print(f"--> Extracting chain {chain}...")
    
    # First pass: collect header and atom lines
    header_lines = []
    atom_lines = []
    
    with open(temp_file, 'r') as infile:
        for line in infile:
            # Keep header information
            if line.startswith(('HEADER', 'TITLE', 'COMPND', 'SOURCE', 'CRYST1')):
                header_lines.append(line)
            # Keep ATOM/HETATM lines for the specified chain
            elif line.startswith(('ATOM', 'HETATM')):
                if len(line) > 21 and line[21] == chain:
                    atom_lines.append(line)
    
    # Write to output file
    with open(output_file, 'w') as outfile:
        # Write headers
        for line in header_lines:
            outfile.write(line)
        
        # Write MODEL record
        outfile.write("MODEL        1\n")
        
        # Write atoms
        for line in atom_lines:
            outfile.write(line)
        
        # Add proper ending records
        outfile.write("ENDMDL\n")
        outfile.write("END\n")
    
    print(f"    Extracted {len(atom_lines)} atoms from chain {chain}")
    
    print(f"✅ Saved chain {chain} to: {output_file}")
    
    # Clean up full file
    temp_file.unlink()
    
    return True


def main():
    parser = argparse.ArgumentParser(
        description="Fetch GLP-1R experimental structures from RCSB PDB"
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("GLP1R"),
        help="Output directory for structures (default: GLP1R/)"
    )
    parser.add_argument(
        "--chain",
        type=str,
        default="A",
        help="Chain ID to extract (default: A for receptor)"
    )
    args = parser.parse_args()
    
    print("=" * 70)
    print("Fetching GLP-1R Experimental Structures from RCSB PDB")
    print("=" * 70)
    
    # Fetch inactive state (6LN2)
    print("\n📥 State A (Inactive): 6LN2")
    success_a = fetch_pdb_structure("6LN2", args.output_dir, args.chain)
    
    # Fetch active state (6X18)
    print("\n📥 State B (Active): 6X18")
    success_b = fetch_pdb_structure("6X18", args.output_dir, args.chain)
    
    print("\n" + "=" * 70)
    if success_a and success_b:
        print("✅ All structures downloaded successfully!")
        print(f"\nStructures saved to: {args.output_dir.resolve()}/")
        print(f"  - Inactive (6LN2): {args.output_dir}/6LN2_chain{args.chain}.pdb")
        print(f"  - Active (6X18):   {args.output_dir}/6X18_chain{args.chain}.pdb")
        print("\nNext steps:")
        print(f"  1. Run: python scripts/prepare_simulation.py --target glp1r_inactive --pdb GLP1R/6LN2_chain{args.chain}.pdb")
        print(f"  2. Run: python scripts/prepare_simulation.py --target glp1r_active --pdb GLP1R/6X18_chain{args.chain}.pdb")
    else:
        print("❌ Some structures failed to download. Please check errors above.")
        return 1
    
    return 0


if __name__ == "__main__":
    exit(main())
