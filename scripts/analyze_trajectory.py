# scripts/analyze_trajectory.py
import mdtraj as md
import matplotlib.pyplot as plt
import argparse
from pathlib import Path

# --- Configuration ---
REPO_ROOT = Path(__file__).resolve().parents[1]
ARTIFACTS_DIR = REPO_ROOT / "artifacts"
DATA_DIR = REPO_ROOT / "data"
# -------------------


def analyze_distance(
    topology_path: Path, trajectory_path: Path, output_plot_path: Path, 
    residue_1: int = 187, residue_2: int = 393, chain_id: int = 0
):
    """
    Calculates and plots the distance between two CA atoms of specified residues over a trajectory.
    
    For Boltz-generated GLP1R structures (355 residues), typical analysis positions:
    - TM3-TM4 region: residues around 140-160
    - TM6-TM7 region: residues around 280-300
    - C-terminal region: residues around 320-340
    
    Parameters:
    - topology_path: Path to the PDB topology file
    - trajectory_path: Path to the DCD trajectory file  
    - output_plot_path: Path where to save the plot
    - residue_1, residue_2: Residue IDs (0-indexed in MDTraj, corresponding to 1-indexed PDB)
    - chain_id: Chain ID to select atoms from
    """
    print(f"--> Loading trajectory from {trajectory_path.name}...")
    traj = md.load(str(trajectory_path), top=str(topology_path))
    
    print(f"--> Finding CA atoms for residues {residue_1+1} and {residue_2+1} (MDTraj 0-indexed)...")
    atom_selection_1 = traj.topology.select(f"chainid {chain_id} and resid {residue_1} and name CA")
    atom_selection_2 = traj.topology.select(f"chainid {chain_id} and resid {residue_2} and name CA")

    if len(atom_selection_1) == 0 or len(atom_selection_2) == 0:
        print(f"Available residues in chain {chain_id}:")
        residues = set()
        for atom in traj.topology.atoms:
            if atom.residue.chain.index == chain_id:
                residues.add((atom.residue.index, atom.residue.name))
        for res_id, res_name in sorted(residues)[:10]:  # Show first 10
            print(f"  Residue {res_id} ({res_name})")
        if len(residues) > 10:
            print(f"  ... and {len(residues)-10} more residues")
        raise ValueError(f"Could not find CA atoms for residues {residue_1} and/or {residue_2} in chain {chain_id}")

    atom_indices = [[atom_selection_1[0], atom_selection_2[0]]]
    print(f"--> Selected atoms: {atom_selection_1[0]} and {atom_selection_2[0]}")
    
    print(f"--> Calculating distance between residues {residue_1+1} and {residue_2+1}...")
    distances = md.compute_distances(traj, atom_indices) * 10  # convert nm to Ångströms

    # Create a time axis in nanoseconds (MDTraj time is in ps)
    time = traj.time / 1000  # convert ps to ns

    print("--> Generating publication-quality plot...")
    plt.style.use('seaborn-v0_8-paper')
    fig, ax = plt.subplots(figsize=(10, 6), dpi=150)
    
    ax.plot(time, distances.flatten(), linewidth=2, color='#2E86AB', alpha=0.8)
    ax.set_xlabel("Time (ns)", fontsize=12, fontweight='bold')
    ax.set_ylabel("Distance (Å)", fontsize=12, fontweight='bold')
    ax.set_title(f"Inter-residue Distance: {residue_1+1} ↔ {residue_2+1}", 
                fontsize=14, fontweight='bold', pad=20)
    ax.grid(True, alpha=0.3)
    ax.set_ylim(bottom=0)
    
    # Add statistics
    mean_dist = distances.mean()
    std_dist = distances.std()
    ax.axhline(y=mean_dist, color='red', linestyle='--', alpha=0.7, 
              label=f'Mean: {mean_dist:.2f} Å')
    ax.fill_between(time, mean_dist - std_dist, mean_dist + std_dist, 
                   color='red', alpha=0.2, label=f'±1σ: {std_dist:.2f} Å')
    
    ax.legend(frameon=True, fancybox=True, shadow=True)
    plt.tight_layout()
    plt.savefig(str(output_plot_path), dpi=300, bbox_inches='tight')
    print(f"✅ Plot saved to {output_plot_path}")
    
    # Print summary statistics
    print("""
📊 Distance Statistics:""")
    print(f"  Mean distance: {mean_dist:.2f} Å")
    print(f"  Standard deviation: {std_dist:.2f} Å")
    print(f"  Minimum distance: {distances.min():.2f} Å")
    print(f"  Maximum distance: {distances.max():.2f} Å")
    print(f"  Total frames analyzed: {len(distances)}")
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Analyze inter-residue distances in MD trajectory")
    parser.add_argument("--target", "-t", default="glp1r", help="Target name (default: glp1r)")
    parser.add_argument("--boltz", action="store_true", 
                       help="Use Boltz workflow (files in data/target/ instead of artifacts/md/TARGET/)")
    parser.add_argument("--res1", type=int, default=140, 
                       help="First residue ID (0-indexed, default: 140 = PDB residue 141, TM3 region)")
    parser.add_argument("--res2", type=int, default=280,
                       help="Second residue ID (0-indexed, default: 280 = PDB residue 281, TM6 region)")
    parser.add_argument("--chain", type=int, default=0, help="Chain ID (default: 0)")
    
    args = parser.parse_args()
    
    # Set paths based on workflow
    if args.boltz:
        md_dir = DATA_DIR / args.target.lower()
    else:
        md_dir = ARTIFACTS_DIR / "md" / args.target.upper()
    
    topo = md_dir / "prepared_system.pdb"
    traj_file = md_dir / "trajectory.dcd"
    plot_file = md_dir / "distance_analysis.png"
    
    # Check if files exist
    if not topo.exists():
        print(f"❌ Error: Topology file not found: {topo}")
        print("   Did you run: python scripts/prepare_simulation.py --target glp1r --boltz")
        exit(1)
        
    if not traj_file.exists():
        print(f"❌ Error: Trajectory file not found: {traj_file}")
        print("   Did you run: python scripts/run_simulation.py --target glp1r --boltz --ns 1.0")
        exit(1)
    
    analyze_distance(topo, traj_file, plot_file, args.res1, args.res2, args.chain)
