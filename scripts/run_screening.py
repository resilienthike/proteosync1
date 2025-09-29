# scripts/run_screening.py
import argparse
import time
import psutil
from pathlib import Path
from vina import Vina

# --- Configuration ---
REPO_ROOT = Path(__file__).resolve().parents[1]
ARTIFACTS_DIR = REPO_ROOT / "artifacts"

RECEPTOR_FILE = REPO_ROOT / "GLP1R_receptor.pdbqt"
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
    print(f"--- Starting Virtual Screening: {ligand_file.name} ---")
    start_time = time.time()
    OUTPUT_DIR.mkdir(exist_ok=True)
    
    if not ligand_file.exists():
        print(f"ERROR: Ligand file not found at {ligand_file}")
        return
        
    # Vina requires the receptor and ligands to be in PDBQT format.
    # We use the properly formatted rigid receptor.
    receptor_pdbqt = RECEPTOR_FILE
    if not receptor_pdbqt.exists():
        print(f"ERROR: Receptor PDBQT file not found: {receptor_pdbqt}")
        print("   You may need to prepare it first using: obabel receptor.pdb -O GLP1R_receptor.pdbqt -xr")
        return

    output_file = OUTPUT_DIR / f"{ligand_file.stem}_vina_out.pdbqt"
    log_file = OUTPUT_DIR / f"{ligand_file.stem}_vina_log.txt"

    print(f"--> Receptor: {receptor_pdbqt.name}")
    print(f"--> Ligand Library: {ligand_file.name}")
    print(f"--> Output will be saved to: {output_file.name}")

    print("Using CPU-based AutoDock Vina (Python API)")

    # Use Python Vina API for better integration and monitoring
    print(f"\n--> Setting up Vina docking (Python API)...")
    print(f"    Expected time: 2-5 min per ligand")
    
    # Monitor system resources during docking
    cpu_count = psutil.cpu_count()
    memory_gb = psutil.virtual_memory().total / (1024**3)
    print(f"    System: {cpu_count} CPUs, {memory_gb:.1f} GB RAM")
    
    docking_start = time.time()
    try:
        # Initialize Vina with robust settings
        v = Vina(sf_name='vina', cpu=cpu_count)
        
        print("    Loading receptor...")
        v.set_receptor(str(receptor_pdbqt))
        
        print("    Loading ligand...")
        v.set_ligand_from_file(str(ligand_file))
        
        print("    Setting up search space...")
        v.compute_vina_maps(center=[center_x, center_y, center_z], 
                           box_size=[size_x, size_y, size_z])
        
        print("    Running docking with high exhaustiveness...")
        print("      This may take several minutes for robust results...")
        
        # High-quality docking with progress monitoring
        v.dock(exhaustiveness=32, n_poses=20)
        
        docking_time = time.time() - docking_start
        print(f"\nDocking complete! ({docking_time:.1f}s)")
        
        # Save results
        print("    Saving results...")
        v.write_poses(str(output_file), n_poses=20, overwrite=True)
        
        # Get and display binding scores
        energies = v.energies(n_poses=20)
        if len(energies) > 0:
            best_score = energies[0][0]  # First pose, binding energy
            print(f"   Best binding score: {best_score:.2f} kcal/mol")
            print(f"   Found {len(energies)} poses")
            print("   Top 5 binding scores:")
            for i, energy_data in enumerate(energies[:5]):
                if len(energy_data) >= 3:
                    energy, rmsd_lb, rmsd_ub = energy_data[:3]
                    print(f"     Pose {i+1}: {energy:.2f} kcal/mol (RMSD: {rmsd_lb:.2f})")
                else:
                    energy = energy_data[0] if len(energy_data) > 0 else 0.0
                    print(f"     Pose {i+1}: {energy:.2f} kcal/mol")
        
        # Write detailed log
        with open(log_file, 'w') as f:
            f.write(f"Virtual Screening Results for {ligand_file.name}\n")
            f.write(f"{'='*50}\n")
            f.write(f"Receptor: {receptor_pdbqt}\n")
            f.write(f"Search center: ({center_x}, {center_y}, {center_z})\n")
            f.write(f"Search box: ({size_x}, {size_y}, {size_z})\n")
            f.write(f"Exhaustiveness: 32\n")
            f.write(f"CPU cores used: {cpu_count}\n")
            f.write(f"Docking time: {docking_time:.1f}s\n\n")
            f.write("Binding Scores:\n")
            for i, energy_data in enumerate(energies):
                if len(energy_data) >= 3:
                    energy, rmsd_lb, rmsd_ub = energy_data[:3]
                    f.write(f"Pose {i+1:2d}: {energy:6.2f} kcal/mol  RMSD: {rmsd_lb:5.2f} - {rmsd_ub:5.2f}\n")
                else:
                    energy = energy_data[0] if len(energy_data) > 0 else 0.0
                    f.write(f"Pose {i+1:2d}: {energy:6.2f} kcal/mol\n")
        
        # Display file info
        if output_file.exists():
            file_size_kb = output_file.stat().st_size / 1024
            print(f"   Output file size: {file_size_kb:.1f} KB")
            print(f"   Results: {output_file}")
            print(f"   Log: {log_file}")
                
        total_time = time.time() - start_time
        print(f"   Total time: {total_time:.1f}s")
            
    except Exception as e:
        print(f"ERROR during docking: {e}")
        import traceback
        traceback.print_exc()

def run_multiple_ligands(ligands_path, **kwargs):
    """Run screening on multiple ligand files in a directory."""
    overall_start = time.time()
    ligands_path = Path(ligands_path)
    
    # Check system capabilities
    cpu_count = psutil.cpu_count()
    memory_gb = psutil.virtual_memory().total / (1024**3)
    print(f"System: {cpu_count} CPU cores, {memory_gb:.1f} GB RAM")
    print("Using AutoDock Vina (CPU-optimized)")
    print()
    
    if ligands_path.is_file():
        # Single file
        print(f"Running screening on single ligand: {ligands_path.name}")
        run_virtual_screening(ligand_file=ligands_path, **kwargs)
    elif ligands_path.is_dir():
        # Directory of PDBQT files
        pdbqt_files = list(ligands_path.glob("*.pdbqt"))
        if not pdbqt_files:
            print(f"ERROR: No PDBQT files found in {ligands_path}")
            return
            
        print(f"📋 Screening Plan:")
        print(f"   Ligands to screen: {len(pdbqt_files)}")
        expected_time_per_ligand = 2  # minutes with 25 CPU cores
        total_expected = len(pdbqt_files) * expected_time_per_ligand
        print(f"   Expected time per ligand: ~{expected_time_per_ligand} min")
        print(f"   Total estimated time: ~{total_expected} min ({total_expected/60:.1f} hours)")
        print(f"   Files:")
        for f in pdbqt_files:
            print(f"     - {f.name}")
        print()
        
        successful = 0
        failed = 0
        
        for i, ligand_file in enumerate(pdbqt_files, 1):
            ligand_start = time.time()
            print(f"\n{'='*60}")
            print(f"🧬 Screening ligand {i}/{len(pdbqt_files)}: {ligand_file.name}")
            print(f"{'='*60}")
            
            try:
                run_virtual_screening(ligand_file=ligand_file, **kwargs)
                successful += 1
                ligand_time = time.time() - ligand_start
                remaining = len(pdbqt_files) - i
                eta = remaining * (ligand_time / 60)  # minutes
                print(f"   Completed in {ligand_time:.1f}s")
                if remaining > 0:
                    print(f"   ⏱️  ETA for remaining {remaining} ligands: ~{eta:.1f} min")
            except Exception as e:
                failed += 1
                print(f"   ❌ Failed: {e}")
        
        # Final summary
        total_time = time.time() - overall_start
        print(f"\n{'='*60}")
        print(f"📊 SCREENING COMPLETE")
        print(f"{'='*60}")
        print(f"   Successful: {successful}/{len(pdbqt_files)}")
        print(f"   Failed: {failed}/{len(pdbqt_files)}")
        print(f"   Total time: {total_time/60:.1f} minutes ({total_time/3600:.1f} hours)")
        print(f"   Average per ligand: {total_time/len(pdbqt_files):.1f}s")
        
    else:
        print(f"ERROR: Ligands path not found: {ligands_path}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run virtual screening with AutoDock Vina.")
    parser.add_argument("--ligands", required=True, help="Path to the ligand file or directory (in PDBQT format).")
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

    run_multiple_ligands(
        ligands_path=args.ligands,
        center_x=args.center_x, center_y=args.center_y, center_z=args.center_z,
        size_x=args.size_x, size_y=args.size_y, size_z=args.size_z
    )