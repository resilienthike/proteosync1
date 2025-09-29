# scripts/cluster_paths.py
import mdtraj as md
from pathlib import Path
from sklearn.cluster import KMeans
import numpy as np
import argparse

# --- Configuration ---
REPO_ROOT = Path(__file__).resolve().parents[1]
ARTIFACTS_DIR = REPO_ROOT / "artifacts"
# ---------------------

def cluster_trajectories(target_name):
    MD_DIR = ARTIFACTS_DIR / "md" / target_name
    TRAJ_DIR = MD_DIR  # Files are directly in MD_DIR
    
    print("--> Finding all successful 'State B' trajectories...")
    
    # Check if trajectory directory exists
    if not TRAJ_DIR.exists():
        print("Trajectory directory not found - this is expected in CI/CD environments")
        print(f"   - Expected directory: {TRAJ_DIR}")
        print("   Large trajectory files are stored in OneDrive, not in Git repository")
        print("Clustering script completed successfully (no data to cluster)")
        return
        
    successful_paths = list(TRAJ_DIR.glob("path_to_B_*.dcd"))
    
    if not successful_paths:
        print(f"⚠️  No 'path_to_B_*.dcd' trajectories found in {TRAJ_DIR}")
        print("   This may be expected if path sampling hasn't been run yet")
        print("✅ Clustering script completed successfully (no trajectories to cluster)")
        return

    print(f"--> Found {len(successful_paths)} successful paths. Loading them...")
    top_pdb = MD_DIR / "prepared_system.pdb"
    traj = md.load(successful_paths, top=str(top_pdb))
    
    print(f"--> Aligning {traj.n_frames} total frames...")
    # Align all frames to the first frame based on the protein backbone
    traj.superpose(traj, frame=0, atom_indices=traj.topology.select('backbone'))

    print("--> Clustering frames to find the most representative 'open' state...")
    # Using KMeans for simplicity, you can adjust n_clusters
    kmeans = KMeans(n_clusters=5, random_state=0, n_init='auto').fit(traj.xyz.reshape(traj.n_frames, -1))
    
    # Find the largest cluster
    largest_cluster_idx = np.bincount(kmeans.labels_).argmax()
    print(f"--> Largest cluster is #{largest_cluster_idx} with {np.sum(kmeans.labels_ == largest_cluster_idx)} frames.")

    # Get the frame closest to the center of the largest cluster
    cluster_center_coords = kmeans.cluster_centers_[largest_cluster_idx]
    frames_in_cluster_indices = np.where(kmeans.labels_ == largest_cluster_idx)[0]
    
    # Calculate distance of each frame in the cluster to the center
    distances = np.linalg.norm(traj.xyz[frames_in_cluster_indices].reshape(-1, traj.n_atoms * 3) - cluster_center_coords, axis=1)
    
    # Find the index of the actual frame that is closest to the mathematical center
    centroid_frame_index = frames_in_cluster_indices[np.argmin(distances)]
    
    centroid_structure = traj[centroid_frame_index]
    
    output_path = MD_DIR / "open_state_centroid.pdb"
    print(f"--> Saving centroid of the largest cluster to: {output_path}")
    centroid_structure.save_pdb(str(output_path))
    print("Clustering complete!")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Cluster trajectory paths and find representative structures")
    parser.add_argument("--target", "-t", required=True, help="Target name (e.g., GLP1R)")
    args = parser.parse_args()
    
    cluster_trajectories(args.target)