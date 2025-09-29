# scripts/postprocess_pdb.py
import mdtraj as md
from pathlib import Path
import argparse

# --- Configuration ---
ARTIFACTS_DIR = Path(__file__).resolve().parents[1] / "artifacts"
# -------------------


def make_molecules_whole(input_pdb_path: Path, output_pdb_path: Path):
    """
    Loads a PDB from a periodic box and saves a new PDB with
    molecules made whole.
    """
    print(f"--> Loading trajectory from: {input_pdb_path}")
    traj = md.load(str(input_pdb_path))

    print("--> Making molecules whole...")
    # This is the key function that resolves the periodic boundary wrapping
    traj.image_molecules(inplace=True)

    print(f"--> Saving new PDB to: {output_pdb_path}")
    traj.save_pdb(str(output_pdb_path))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Make molecules whole in PDB file")
    parser.add_argument("--target", "-t", required=True, help="Target name (e.g., GLP1R)")
    args = parser.parse_args()
    
    MD_DIR = ARTIFACTS_DIR / "md" / args.target
    in_pdb = MD_DIR / "prepared_system.pdb"
    out_pdb = MD_DIR / "prepared_system_whole.pdb"

    if not in_pdb.exists():
        print(f"Error: Input file not found at {in_pdb}")
    else:
        make_molecules_whole(in_pdb, out_pdb)
        print("\n✅ Processing complete!")
