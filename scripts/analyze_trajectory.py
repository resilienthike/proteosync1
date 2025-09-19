# scripts/analyze_trajectory.py
import mdtraj as md
import matplotlib.pyplot as plt
from pathlib import Path

# --- Configuration ---
TARGET_NAME = "GLP1R"
ARTIFACTS_DIR = Path(__file__).resolve().parents[1] / "artifacts"
MD_DIR = ARTIFACTS_DIR / "md" / TARGET_NAME
# -------------------


def analyze_distance(
    topology_path: Path, trajectory_path: Path, output_plot_path: Path
):
    """
    Calculates and plots the distance between two residues over a trajectory.
    """
    print(f"--> Loading trajectory from {trajectory_path.name}...")
    traj = md.load(str(trajectory_path), top=str(topology_path))

    # --- THIS SECTION IS THE FIX ---
    # We now iterate through all atoms and check for the correct chain, residue number, and atom name.
    print("--> Finding atom indices for residues 188 and 394...")
    atom_selection_1 = traj.topology.select("chainid 0 and resid 187 and name CA")
    atom_selection_2 = traj.topology.select("chainid 0 and resid 393 and name CA")

    if len(atom_selection_1) == 0 or len(atom_selection_2) == 0:
        raise ValueError("Could not find one or both atoms. Check residue numbers.")

    atom_indices = [[atom_selection_1[0], atom_selection_2[0]]]
    # -----------------------------

    print("--> Calculating distance between residues 188 and 394...")
    distances = md.compute_distances(traj, atom_indices) * 10  # convert nm to Ångströms

    # Create a time axis in nanoseconds
    time = traj.time / 1000  # convert ps to ns

    print("--> Generating plot...")
    plt.figure(figsize=(10, 6))
    plt.plot(time, distances)
    plt.title(f"Pocket Distance for {TARGET_NAME} (Res 188-394)")
    plt.xlabel("Time (ns)")
    plt.ylabel("Distance (Å)")
    plt.ylim(bottom=0)  # Ensure y-axis starts at 0
    plt.grid(True)
    plt.savefig(str(output_plot_path))

    print(f"✅ Plot saved to {output_plot_path}")


if __name__ == "__main__":
    topo = MD_DIR / "prepared_system_whole.pdb"
    traj_file = MD_DIR / "trajectory.dcd"
    plot_file = MD_DIR / "distance_over_time.png"

    analyze_distance(topo, traj_file, plot_file)
