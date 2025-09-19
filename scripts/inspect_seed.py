# scripts/inspect_seed.py
from __future__ import annotations

import argparse
from pathlib import Path
from Bio.PDB import MMCIFParser, PDBParser, is_aa

# Define project paths
REPO_ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = REPO_ROOT / "artifacts" / "data"


def inspect_structure(target_name: str):
    """
    Loads a seed structure, validates its contents, and prints a summary.
    """
    print(f"--- Inspecting Seed Structure for Target: {target_name} ---")

    target_dir = DATA_DIR / target_name

    # Auto-pick PDB vs CIF file
    pdb = target_dir / "seed_structure.pdb"
    cif = target_dir / "seed_structure.cif"

    structure_file = cif if cif.exists() else pdb

    # 1. Verify the file exists
    print(f"Seed file: {structure_file}")
    if not structure_file.exists():
        raise FileNotFoundError(f"Seed structure not found in {target_dir}")
    print(f"Exists: True | Size: {structure_file.stat().st_size} bytes")

    # 2. Parse the structure using the correct parser
    print("\n--> Parsing structure...")
    if structure_file.suffix.lower() in {".cif", ".mmcif"}:
        parser = MMCIFParser(QUIET=True)
    else:
        parser = PDBParser(QUIET=True)

    structure = parser.get_structure(target_name, str(structure_file))
    print("--> Parsing complete.")

    # 3. Calculate chain lengths (counting only amino acid residues)
    chain_lengths = {}
    for model in structure:
        for chain in model:
            chain_lengths[chain.id] = sum(
                1 for residue in chain if is_aa(residue, standard=True)
            )

    print("\n--- Summary ---")
    print(f"Chain lengths (amino acids): {chain_lengths}")

    # 4. Identify receptor and peptide candidates
    receptor_chain = (
        max(chain_lengths, key=chain_lengths.get) if chain_lengths else None
    )
    peptide_candidates = [
        cid for cid, length in chain_lengths.items() if 5 <= length <= 120
    ]

    print(
        f"Identified Receptor Chain: {receptor_chain} (length {chain_lengths.get(receptor_chain)})"
    )
    print(f"Identified Peptide Candidates: {peptide_candidates}")
    print("-" * 20)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Inspect a seed structure file.")
    parser.add_argument(
        "--target", "-t", required=True, help="Target name (e.g., GLP1R)"
    )
    args = parser.parse_args()

    inspect_structure(args.target)
