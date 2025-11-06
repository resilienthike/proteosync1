#!/usr/bin/env python3
"""
Analyze path sampling results from ProteOSync.

This script:
1. Loads all elite trajectory files (path_to_B_*.dcd)
2. Computes distances for all CV definitions
3. Creates visualizations of transition pathways
4. Analyzes committor predictions from the trained model
5. Identifies the best transition candidates

Usage:
    python scripts/analyze_results.py --target glp1r --boltz
"""

import argparse
import numpy as np
import mdtraj as md
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import torch
from openmm.unit import nanometers

# Import the model for committor predictions
import sys
sys.path.insert(0, str(Path(__file__).parent))
from model import CommittorNet, frame_to_torch_graph


# State definitions (same as in run_path_sampling.py)
STATE_DEFINITIONS = [
    {
        "name": "ecd_tm_loop_contact_CV2",
        "residues": ('chainid 0 and resid 84 and name CA', 'chainid 0 and resid 249 and name CA'),
        "state_a_threshold": 5.9 * nanometers,
        "state_b_threshold": 6.4 * nanometers
    },
    {
        "name": "intracellular_coupling_CV3",
        "residues": ('chainid 0 and resid 179 and name CA', 'chainid 0 and resid 219 and name CA'),
        "state_a_threshold": 4.0 * nanometers,
        "state_b_threshold": 4.1 * nanometers
    },
    {
        "name": "extracellular_gate_CV4",
        "residues": ('chainid 0 and resid 99 and name CA', 'chainid 0 and resid 279 and name CA'),
        "state_a_threshold": 5.1 * nanometers,
        "state_b_threshold": 5.3 * nanometers
    }
]


def compute_distances(traj, definition):
    """Compute distance for a given CV definition across all frames."""
    atoms1 = traj.topology.select(definition["residues"][0])
    atoms2 = traj.topology.select(definition["residues"][1])
    
    if len(atoms1) == 0 or len(atoms2) == 0:
        return None
    
    distances = md.compute_distances(traj, [[atoms1[0], atoms2[0]]])
    return distances.flatten()


def load_elite_trajectories(data_dir, pdb_path):
    """Load all elite trajectory files."""
    elite_files = sorted(data_dir.glob("path_to_B_*.dcd"))
    print(f"\n🔍 Found {len(elite_files)} elite trajectory files")
    
    trajectories = []
    for elite_file in elite_files:
        try:
            traj = md.load(str(elite_file), top=str(pdb_path))
            trajectories.append({
                "file": elite_file.name,
                "traj": traj,
                "n_frames": traj.n_frames
            })
        except Exception as e:
            print(f"   ⚠️  Failed to load {elite_file.name}: {e}")
    
    print(f"   ✅ Successfully loaded {len(trajectories)} trajectories")
    return trajectories


def analyze_cv_distances(trajectories):
    """Analyze CV distances across all elite trajectories."""
    results = {}
    
    for definition in STATE_DEFINITIONS:
        cv_name = definition["name"]
        state_a_thresh = definition["state_a_threshold"].value_in_unit(nanometers)
        state_b_thresh = definition["state_b_threshold"].value_in_unit(nanometers)
        
        print(f"\n📊 Analyzing {cv_name}...")
        print(f"   State A threshold: {state_a_thresh:.2f} nm")
        print(f"   State B threshold: {state_b_thresh:.2f} nm")
        
        all_distances = []
        min_dist_to_b = float('inf')
        best_traj = None
        
        for traj_data in trajectories:
            distances = compute_distances(traj_data["traj"], definition)
            if distances is not None:
                all_distances.extend(distances)
                
                # Find how close this trajectory gets to State B
                dist_to_b = state_b_thresh - distances.max()
                if dist_to_b < min_dist_to_b:
                    min_dist_to_b = dist_to_b
                    best_traj = traj_data["file"]
        
        if all_distances:
            all_distances = np.array(all_distances)
            results[cv_name] = {
                "distances": all_distances,
                "mean": all_distances.mean(),
                "std": all_distances.std(),
                "min": all_distances.min(),
                "max": all_distances.max(),
                "state_a_thresh": state_a_thresh,
                "state_b_thresh": state_b_thresh,
                "best_traj": best_traj,
                "min_dist_to_b": min_dist_to_b
            }
            
            print(f"   Mean: {all_distances.mean():.3f} nm (±{all_distances.std():.3f})")
            print(f"   Range: [{all_distances.min():.3f}, {all_distances.max():.3f}] nm")
            print(f"   Closest to State B: {min_dist_to_b:.4f} nm (in {best_traj})")
    
    return results


def plot_cv_distributions(results, output_dir):
    """Create distribution plots for all CVs."""
    n_cvs = len(results)
    fig, axes = plt.subplots(n_cvs, 1, figsize=(12, 4*n_cvs))
    if n_cvs == 1:
        axes = [axes]
    
    for ax, (cv_name, data) in zip(axes, results.items()):
        distances = data["distances"]
        state_a = data["state_a_thresh"]
        state_b = data["state_b_thresh"]
        
        # Histogram
        ax.hist(distances, bins=50, alpha=0.7, color='steelblue', edgecolor='black')
        
        # Threshold lines
        ax.axvline(state_a, color='green', linestyle='--', linewidth=2, label=f'State A ({state_a:.2f} nm)')
        ax.axvline(state_b, color='red', linestyle='--', linewidth=2, label=f'State B ({state_b:.2f} nm)')
        
        ax.set_xlabel('Distance (nm)', fontsize=12)
        ax.set_ylabel('Count', fontsize=12)
        ax.set_title(f'{cv_name}\nClosest to B: {data["min_dist_to_b"]:.4f} nm', fontsize=14, fontweight='bold')
        ax.legend(fontsize=10)
        ax.grid(alpha=0.3)
    
    plt.tight_layout()
    output_path = output_dir / "cv_distributions.png"
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"\n📈 Saved CV distributions to {output_path}")
    plt.close()


def plot_trajectory_evolution(trajectories, output_dir):
    """Plot CV evolution over time for each trajectory."""
    n_trajs = min(10, len(trajectories))  # Plot first 10
    
    fig, axes = plt.subplots(len(STATE_DEFINITIONS), 1, figsize=(14, 4*len(STATE_DEFINITIONS)))
    if len(STATE_DEFINITIONS) == 1:
        axes = [axes]
    
    for ax, definition in zip(axes, STATE_DEFINITIONS):
        cv_name = definition["name"]
        state_a = definition["state_a_threshold"].value_in_unit(nanometers)
        state_b = definition["state_b_threshold"].value_in_unit(nanometers)
        
        for i, traj_data in enumerate(trajectories[:n_trajs]):
            distances = compute_distances(traj_data["traj"], definition)
            if distances is not None:
                ax.plot(distances, alpha=0.6, label=f'{traj_data["file"][:15]}...')
        
        ax.axhline(state_a, color='green', linestyle='--', linewidth=2, label='State A')
        ax.axhline(state_b, color='red', linestyle='--', linewidth=2, label='State B')
        
        ax.set_xlabel('Frame', fontsize=12)
        ax.set_ylabel('Distance (nm)', fontsize=12)
        ax.set_title(f'{cv_name} Evolution', fontsize=14, fontweight='bold')
        ax.legend(fontsize=8, ncol=2)
        ax.grid(alpha=0.3)
    
    plt.tight_layout()
    output_path = output_dir / "trajectory_evolution.png"
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"📈 Saved trajectory evolution to {output_path}")
    plt.close()


def analyze_committor_predictions(trajectories, checkpoint_path, device, pdb_path):
    """Load trained model and analyze committor predictions."""
    if not checkpoint_path.exists():
        print(f"\n⚠️  Checkpoint not found: {checkpoint_path}")
        return None
    
    print(f"\n🧠 Loading trained model from {checkpoint_path.name}...")
    
    # Load model
    model = CommittorNet(embedding_dim=64, n_layers=3, cutoff_A=5.0).to(device)
    checkpoint = torch.load(checkpoint_path, map_location=device)
    model.load_state_dict(checkpoint['model_state_dict'])
    model.eval()
    
    print("   ✅ Model loaded successfully")
    
    # Predict committors for each trajectory
    committor_data = []
    
    with torch.no_grad():
        for traj_data in trajectories[:10]:  # Analyze first 10
            traj = traj_data["traj"]
            committors = []
            
            for frame_idx in range(traj.n_frames):
                frame = traj[frame_idx]
                g, l_g = frame_to_torch_graph(frame, device)
                pred = model(g, l_g).item()
                committors.append(pred)
            
            committor_data.append({
                "file": traj_data["file"],
                "committors": committors
            })
    
    # Plot committor evolution
    fig, ax = plt.subplots(figsize=(14, 6))
    
    for data in committor_data:
        ax.plot(data["committors"], alpha=0.6, label=f'{data["file"][:15]}...')
    
    ax.axhline(0.5, color='black', linestyle='--', linewidth=2, label='Transition state (q=0.5)')
    ax.set_xlabel('Frame', fontsize=12)
    ax.set_ylabel('Committor Probability', fontsize=12)
    ax.set_title('Committor Evolution Across Elite Trajectories', fontsize=14, fontweight='bold')
    ax.legend(fontsize=8, ncol=2)
    ax.grid(alpha=0.3)
    ax.set_ylim([-0.05, 1.05])
    
    plt.tight_layout()
    output_path = checkpoint_path.parent / "committor_predictions.png"
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"📈 Saved committor predictions to {output_path}")
    plt.close()
    
    return committor_data


def main():
    parser = argparse.ArgumentParser(description="Analyze path sampling results")
    parser.add_argument("--target", required=True, help="Target name (e.g., glp1r)")
    parser.add_argument("--boltz", action="store_true", help="Use Boltz workflow")
    parser.add_argument("--checkpoint", type=str, default=None, help="Checkpoint file to analyze (default: latest)")
    args = parser.parse_args()
    
    # Setup paths
    base_dir = Path(__file__).parent.parent
    data_dir = base_dir / "data" / args.target
    pdb_path = data_dir / "prepared_system.pdb"
    
    if not data_dir.exists():
        print(f"❌ Data directory not found: {data_dir}")
        return
    
    if not pdb_path.exists():
        print(f"❌ PDB file not found: {pdb_path}")
        return
    
    # Find checkpoint
    if args.checkpoint:
        checkpoint_path = data_dir / args.checkpoint
    else:
        # Find the latest checkpoint
        checkpoints = sorted(data_dir.glob("checkpoint_*.pt"))
        if checkpoints:
            checkpoint_path = checkpoints[-1]
        else:
            checkpoint_path = None
    
    print("=" * 80)
    print("ProteOSync Result Analysis")
    print("=" * 80)
    print(f"Target: {args.target}")
    print(f"Data directory: {data_dir}")
    print(f"PDB file: {pdb_path}")
    if checkpoint_path:
        print(f"Checkpoint: {checkpoint_path.name}")
    
    # Load trajectories
    trajectories = load_elite_trajectories(data_dir, pdb_path)
    
    if not trajectories:
        print("❌ No trajectories loaded. Nothing to analyze.")
        return
    
    # Analyze CV distances
    cv_results = analyze_cv_distances(trajectories)
    
    # Create visualizations
    plot_cv_distributions(cv_results, data_dir)
    plot_trajectory_evolution(trajectories, data_dir)
    
    # Analyze committor predictions
    device = torch.device('cuda:0' if torch.cuda.is_available() else 'cpu')
    if checkpoint_path:
        analyze_committor_predictions(trajectories, checkpoint_path, device, pdb_path)
    
    # Summary
    print("\n" + "=" * 80)
    print("📋 SUMMARY")
    print("=" * 80)
    print(f"Total elite trajectories analyzed: {len(trajectories)}")
    print(f"Total frames across all trajectories: {sum(t['n_frames'] for t in trajectories)}")
    
    print("\n🎯 Best approaches to State B:")
    for cv_name, data in cv_results.items():
        print(f"   {cv_name}: {data['min_dist_to_b']:.4f} nm gap")
        print(f"      → Best trajectory: {data['best_traj']}")
    
    print(f"\n✅ Analysis complete! Results saved to {data_dir}")
    print("=" * 80)


if __name__ == "__main__":
    main()
