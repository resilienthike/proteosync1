# scripts/run_screening.py
import argparse
import subprocess
from pathlib import Path

# --- Configuration ---
TARGET_NAME = "GLP1R"
REPO_ROOT = Path(__file__).resolve().parents[1]
ARTIFACTS_DIR = REPO_ROOT / "artifacts"
MD_DIR = ARTIFACTS_DIR / "md" / TARGET_NAME

RECEPTOR_FILE = MD_DIR / "open_state_centroid.pdb"
OUTPUT_DIR = ARTIFACTS_DIR / "screening_results"

# Default pocket coordinates and dimensions (from fpocket analysis)
DEFAULT_CENTER_X = -1.78
DEFAULT_CENTER_Y = 0.08
DEFAULT_CENTER_Z = -0.47

DEFAULT_SIZE_X = 14.0
DEFAULT_SIZE_Y = 11.7
DEFAULT_SIZE_Z = 9.5
# ---------------------

def run_virtual_screening_default_pocket(ligand_file: Path):
    """
    Runs virtual screening using the default GLP-1R pocket coordinates.
    
    Uses the coordinates determined from fpocket analysis:
    - Center: (-1.78, 0.08, -0.47)
    - Size: (14.0, 11.7, 9.5)
    """
    return run_virtual_screening(
        ligand_file=ligand_file,
        center_x=DEFAULT_CENTER_X, center_y=DEFAULT_CENTER_Y, center_z=DEFAULT_CENTER_Z,
        size_x=DEFAULT_SIZE_X, size_y=DEFAULT_SIZE_Y, size_z=DEFAULT_SIZE_Z
    )

def run_virtual_screening(
    ligand_file: Path, 
    center_x: float, center_y: float, center_z: float,
    size_x: float, size_y: float, size_z: float
):
    """
    Runs AutoDock Vina to screen a library of ligands against a receptor pocket.
    
    Args:
        ligand_file: Path to ligand library in PDBQT format
        center_x, center_y, center_z: Center coordinates of the binding pocket
        size_x, size_y, size_z: Dimensions of the search box in Angstroms
    """
    print("--- Starting Virtual Screening ---")
    OUTPUT_DIR.mkdir(exist_ok=True)
    
    if not ligand_file.exists():
        print(f"🔥 Error: Ligand file not found at {ligand_file}")
        return
        
    # Vina requires the receptor and ligands to be in PDBQT format.
    # We assume the user has pre-converted them. A more advanced script
    # would use Open Babel (obabel) to do this automatically.
    receptor_pdbqt = RECEPTOR_FILE.with_suffix(".pdbqt")
    if not receptor_pdbqt.exists():
        print(f"🔥 Error: Receptor PDBQT file not found: {receptor_pdbqt}")
        print("   You may need to prepare it first using AutoDock Tools or a similar program.")
        return

    output_file = OUTPUT_DIR / f"{ligand_file.stem}_vina_out.pdbqt"
    log_file = OUTPUT_DIR / f"{ligand_file.stem}_vina_log.txt"

    print(f"--> Receptor: {receptor_pdbqt.name}")
    print(f"--> Ligand Library: {ligand_file.name}")
    print(f"--> Output will be saved to: {output_file.name}")

    # Construct the command-line arguments for Vina
    command = [
        "vina",
        "--receptor", str(receptor_pdbqt),
        "--ligand", str(ligand_file),
        "--out", str(output_file),
        "--log", str(log_file),
        "--center_x", str(center_x),
        "--center_y", str(center_y),
        "--center_z", str(center_z),
        "--size_x", str(size_x),
        "--size_y", str(size_y),
        "--size_z", str(size_z),
        "--exhaustiveness", "16", # Higher value for more thorough search
        "--num_modes", "10",     # Save top 10 binding poses
    ]

    print("\n--> Running Vina... (This may take a long time)")
    try:
        subprocess.run(command, check=True)
        print("\n✅ Virtual screening complete!")
        print(f"   Results are in: {output_file}")
        print(f"   Log file with scores is at: {log_file}")
    except FileNotFoundError:
        print("🔥 Error: 'vina' command not found. Please ensure AutoDock Vina is installed.")
    except subprocess.CalledProcessError as e:
        print(f"🔥 Error: Vina exited with a non-zero status: {e}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run virtual screening with AutoDock Vina.")
    parser.add_argument("--ligands", required=True, help="Path to the ligand library file (in PDBQT format).")
    parser.add_argument("--center_x", type=float, default=DEFAULT_CENTER_X, 
                       help=f"Center X coordinate of the pocket (default: {DEFAULT_CENTER_X}).")
    parser.add_argument("--center_y", type=float, default=DEFAULT_CENTER_Y,
                       help=f"Center Y coordinate of the pocket (default: {DEFAULT_CENTER_Y}).")
    parser.add_argument("--center_z", type=float, default=DEFAULT_CENTER_Z,
                       help=f"Center Z coordinate of the pocket (default: {DEFAULT_CENTER_Z}).")
    parser.add_argument("--size_x", type=float, default=DEFAULT_SIZE_X,
                       help=f"Size of the box in X direction in Angstroms (default: {DEFAULT_SIZE_X}).")
    parser.add_argument("--size_y", type=float, default=DEFAULT_SIZE_Y,
                       help=f"Size of the box in Y direction in Angstroms (default: {DEFAULT_SIZE_Y}).")
    parser.add_argument("--size_z", type=float, default=DEFAULT_SIZE_Z,
                       help=f"Size of the box in Z direction in Angstroms (default: {DEFAULT_SIZE_Z}).")
    
    args = parser.parse_args()

    run_virtual_screening(
        ligand_file=Path(args.ligands),
        center_x=args.center_x, center_y=args.center_y, center_z=args.center_z,
        size_x=args.size_x, size_y=args.size_y, size_z=args.size_z
    )