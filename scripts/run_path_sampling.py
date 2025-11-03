 # scripts/run_path_sampling.py
import time
import random
import sys
import subprocess
import threading
import gc
import os
import argparse
from pathlib import Path
from collections import Counter

# Add current directory to Python path for local imports
sys.path.insert(0, str(Path(__file__).parent))

import mdtraj as md
import numpy as np
import torch
import torch.optim as optim
import torch.nn as nn
from torch.multiprocessing import Process, Queue, set_start_method
from tqdm import tqdm

from openmm.unit import nanometers

from model import CommittorNet, frame_to_torch_graph
from worker import simulation_worker

# Configure CUDA memory management with emergency safeguards
os.environ['PYTORCH_CUDA_ALLOC_CONF'] = 'expandable_segments:True'
if torch.cuda.is_available():
    torch.cuda.empty_cache()
    # Set memory fraction to prevent full GPU memory allocation
    torch.cuda.set_per_process_memory_fraction(0.85)  # Reduced for safety margin
    
# Emergency memory thresholds (in GB)
EMERGENCY_MEMORY_THRESHOLD = 60.0  # Trigger emergency cleanup at 60GB
CRITICAL_MEMORY_THRESHOLD = 65.0   # Aggressive cleanup at 65GB
MAX_SAFE_MEMORY = 70.0             # Hard limit - reduce operations

# --- Configuration ---
# ==============================================================================
# Parse arguments first to get target name
parser = argparse.ArgumentParser(description="Run AI-guided path sampling for protein conformational transitions")
parser.add_argument("--target", "-t", default="glp1r", help="Target name (default: glp1r)")
parser.add_argument("--boltz", action="store_true", 
                   help="Use Boltz workflow (files in data/target/ instead of artifacts/md/TARGET/)")
parser.add_argument("--dry-run", action="store_true", help="Perform preflight checks (CUDA/DGL/model) and exit")
parser.add_argument("--startup", action="store_true", help="Instantiate the AIMDRunner and perform a brief startup (no workers)")
parser.add_argument("--shots", type=int, default=None, help="Override N_TOTAL_SHOTS for a short run")
parser.add_argument("--workers", type=int, default=None, help="Override N_WORKERS for a short run")
parser.add_argument("--controlled", action="store_true", help="Run a controlled short simulation without spawning worker processes (simulates worker results)")
args = parser.parse_args()

TARGET_NAME = args.target
REPO_ROOT = Path(__file__).resolve().parents[1]
ARTIFACTS_DIR = REPO_ROOT / "artifacts"
DATA_DIR = REPO_ROOT / "data"

# Set paths based on workflow
if args.boltz:
    MD_DIR = DATA_DIR / TARGET_NAME.lower()
else:
    MD_DIR = ARTIFACTS_DIR / "md" / TARGET_NAME.upper()
MD_DIR.mkdir(parents=True, exist_ok=True)

# Multiple state definitions to explore different structural regions
# Updated for Boltz-generated GLP1R structure (355 residues)
# ==============================================================================
# --- Scientifically-Valid State Definitions (Calibrated Nov 3, 2025) ---
# Based on PDBs 6LN2 (Inactive) and 6X18 (Active)
#
# CV2 (ECD-TM5): Inactive=57.7Å, Active=65.4Å
# CV3 (TM4-TM5): Inactive=39.1Å, Active=41.7Å
# CV4 (ECD-TM6): Inactive=49.9Å, Active=54.4Å
#
# We are ignoring CV1 (139-279) as it compresses upon activation (41.3Å -> 37.1Å),
# which conflicts with the current worker logic (dist > state_b).
# ==============================================================================
STATE_DEFINITIONS = [
    {
        "name": "ecd_tm_loop_contact_CV2",
        "residues": ('chainid 0 and resid 84 and name CA', 'chainid 0 and resid 249 and name CA'),
        "state_a_threshold": 5.9 * nanometers,  # 59.0 Å (just above 57.7Å inactive)
        "state_b_threshold": 6.4 * nanometers   # 64.0 Å (just below 65.4Å active)
    },
    {
        "name": "intracellular_coupling_CV3",
        "residues": ('chainid 0 and resid 179 and name CA', 'chainid 0 and resid 219 and name CA'),
        "state_a_threshold": 4.0 * nanometers,  # 40.0 Å (just above 39.1Å inactive)
        "state_b_threshold": 4.1 * nanometers   # 41.0 Å (just below 41.7Å active)
    },
    {
        "name": "extracellular_gate_CV4",
        "residues": ('chainid 0 and resid 99 and name CA', 'chainid 0 and resid 279 and name CA'),
        "state_a_threshold": 5.1 * nanometers,  # 51.0 Å (just above 49.9Å inactive)
        "state_b_threshold": 5.3 * nanometers   # 53.0 Å (just below 54.4Å active)
    }
]

N_WORKERS = 1
N_TOTAL_SHOTS = 100  # Reduced for testing - increase to 1000 for production
LEARNING_RATE = 1e-4
REPLAY_BUFFER_SIZE = 200
TRAINING_BATCH_SIZE = 32

WORKER_CONFIG = {
    "pdb_path": str(MD_DIR / "prepared_system.pdb"),
    "state_path": str(MD_DIR / "minimized_state.xml"),
    "system_path": str(MD_DIR / "system.xml"),
    "state_definitions": STATE_DEFINITIONS,  # Pass the whole list to the worker
    "shooting_move_steps": 10000,
    "report_interval": 100,  # Saves a frame every 100 steps (10000/100 = 100 frames per trajectory)
}
# ==============================================================================


def monitor_gpu_utilization(stop_event, pbar, interval=5):
    """Monitor GPU utilization and update the progress bar description."""
    while not stop_event.is_set():
        try:
            result = subprocess.run(
                [
                    "nvidia-smi",
                    "--query-gpu=utilization.gpu",
                    "--format=csv,noheader,nounits",
                ],
                capture_output=True,
                text=True,
                check=True,
            )
            utils = [f"{util.strip()}%" for util in result.stdout.strip().split("\n")]
            current_desc = pbar.desc or ""
            if "Dispatched:" in current_desc or "Completed:" in current_desc:
                base_info = (
                    current_desc.split("GPU Util:")[1].split("|")[1:]
                    if "|" in current_desc
                    else []
                )
                pbar.set_description(
                    f"GPU Util: {' / '.join(utils)} | {' | '.join(base_info)}"
                )
            else:
                pbar.set_description(f"GPU Util: {' / '.join(utils)}")
        except Exception:
            pbar.set_description("GPU Util: N/A")

        time.sleep(interval)


class AIMDRunner:
    """Orchestrates the entire AIMD process."""

    def __init__(self, initial_traj_path: str):
        print("--> Initializing AIMMD Runner...")
        self.topology = md.load_pdb(WORKER_CONFIG["pdb_path"]).topology
        self.current_path = md.load(initial_traj_path, top=WORKER_CONFIG["pdb_path"])

        # Enforce GPU-only execution. The user requested no CPU fallback.
        # Fail fast with a clear error if CUDA or a CUDA-enabled DGL is not
        # available so the environment must be corrected before running.
        try:
            import dgl
        except Exception as e:
            raise RuntimeError(
                "DGL import failed. This pipeline requires a CUDA-enabled DGL. "
                "Install or build DGL with CUDA support that matches your PyTorch/CUDA runtime."
            ) from e

        if not torch.cuda.is_available():
            raise RuntimeError(
                "No CUDA device detected (torch.cuda.is_available() is False). "
                "This pipeline requires a GPU; aborting (no CPU fallback)."
            )

        # Verify that the installed DGL can move graphs to the CUDA device.
        try:
            tmp_g = dgl.graph(([], []))
            tmp_g.to(torch.device("cuda:0"))
        except Exception as e:
            raise RuntimeError(
                "Installed DGL does not appear to support CUDA (failed to move a graph to cuda:0). "
                "Please install a CUDA-enabled DGL matching PyTorch/CUDA, or rebuild DGL from source with CUDA."
            ) from e

        print("--> Initializing CommittorNet AI model on cuda:0...")
        self.device = torch.device("cuda:0")

        # Pass in hyperparameters for the EquiThreeBody engine
        self.model = CommittorNet(
            embedding_dim=64,  # You can tune this
            n_layers=3         # You can tune this
        ).to(self.device)
        self.optimizer = optim.Adam(self.model.parameters(), lr=LEARNING_RATE)
        self.loss_fn = nn.BCELoss()
        
        # Memory management tracking
        self.iteration_count = 0
        self.cleanup_interval = 10  # Clean memory every N iterations
        self.replay_buffer = []
        self.shot_tracking = {}  # Map shot_id -> replay_buffer_index
        
        # Elite buffer for near-miss exploration
        self.elite_trajectories = []  # Store best partial successes
        self.elite_buffer_size = 10
        self.best_distance_to_b = float('inf')  # Track closest approach to State B
        
        # Enhanced path tracking
        self.path_history = []  # Store successful paths for diversity
        self.current_path_age = 0  # How long we've been using current path
        self.max_path_age = 30  # Switch paths if no progress (reduced for 1000-shot runs)

        self.task_queue = Queue()
        self.result_queue = Queue()
        self.workers = []

        self.results_summary = Counter()
        self.dispatched_shots = 0
        
        # Emergency memory management
        self.emergency_mode = False
        self.memory_warnings = 0
        self.dynamic_batch_size = TRAINING_BATCH_SIZE
        self.dynamic_buffer_size = REPLAY_BUFFER_SIZE

    def start_workers(self):
        print(f"--> Starting {N_WORKERS} simulation workers...")
        for i in range(N_WORKERS):
            p = Process(
                target=simulation_worker,
                args=(i, self.task_queue, self.result_queue, WORKER_CONFIG),
            )
            p.start()
            self.workers.append(p)

    def stop_workers(self):
        for _ in self.workers:
            self.task_queue.put(None)
        for p in self.workers:
            p.join(timeout=5)
            if p.is_alive():
                p.terminate()
    
    def check_memory_status(self):
        """Monitor GPU memory and trigger emergency procedures if needed."""
        if not torch.cuda.is_available():
            return False
            
        memory_reserved = torch.cuda.memory_reserved() / 1024**3     # GB
        
        # Check for emergency conditions
        if memory_reserved > CRITICAL_MEMORY_THRESHOLD:
            print(f"🚨 CRITICAL MEMORY: {memory_reserved:.1f}GB reserved - EMERGENCY CLEANUP!")
            self.emergency_cleanup()
            self.emergency_mode = True
            self.memory_warnings += 1
            return True
        elif memory_reserved > EMERGENCY_MEMORY_THRESHOLD:
            print(f"⚠️  HIGH MEMORY: {memory_reserved:.1f}GB reserved - reducing operations")
            self.reduce_memory_footprint()
            self.memory_warnings += 1
            return True
        elif memory_reserved > MAX_SAFE_MEMORY:
            print(f"🛑 MAXIMUM MEMORY REACHED: {memory_reserved:.1f}GB - stopping new shots")
            return True
            
        # Reset emergency mode if memory is back to normal
        if self.emergency_mode and memory_reserved < EMERGENCY_MEMORY_THRESHOLD * 0.8:
            print(f"✅ Memory recovered: {memory_reserved:.1f}GB - resuming normal operations")
            self.emergency_mode = False
            
        return False
    
    def emergency_cleanup(self):
        """Aggressive memory cleanup for emergency situations."""
        print("[EMERGENCY] Performing aggressive memory cleanup...")
        
        # Drastically reduce buffer sizes
        self.dynamic_buffer_size = max(50, self.dynamic_buffer_size // 4)
        self.dynamic_batch_size = max(2, self.dynamic_batch_size // 4)
        
        # Aggressive replay buffer cleanup
        if len(self.replay_buffer) > self.dynamic_buffer_size:
            # Keep only the most recent labeled samples and active shots
            labeled_samples = [(i, sample) for i, sample in enumerate(self.replay_buffer) if sample[1] != -1.0]
            active_shot_indices = set(self.shot_tracking.values())
            
            # Keep only recent labeled samples
            keep_labeled_count = max(10, self.dynamic_buffer_size // 3)
            if len(labeled_samples) > keep_labeled_count:
                keep_labeled = labeled_samples[-keep_labeled_count:]
            else:
                keep_labeled = labeled_samples
                
            # Keep active shots
            keep_active = [(i, sample) for i, sample in enumerate(self.replay_buffer) 
                          if i in active_shot_indices]
            
            # Rebuild buffer with minimal data
            keep_indices = set(idx for idx, _ in keep_labeled + keep_active)
            new_buffer = [self.replay_buffer[idx] for idx in sorted(keep_indices)]
            
            # Update tracking
            index_mapping = {old_idx: new_idx for new_idx, old_idx in enumerate(sorted(keep_indices))}
            for shot_id, old_idx in list(self.shot_tracking.items()):
                if old_idx in index_mapping:
                    self.shot_tracking[shot_id] = index_mapping[old_idx]
                else:
                    del self.shot_tracking[shot_id]
                    
            self.replay_buffer = new_buffer
            
        # Reduce elite buffer
        if len(self.elite_trajectories) > 3:
            self.elite_trajectories = self.elite_trajectories[:3]
            
        # Force comprehensive cleanup
        self.cleanup_memory()
        
        print(f"[EMERGENCY] Cleanup complete: buffer={len(self.replay_buffer)}, batch_size={self.dynamic_batch_size}")
    
    def reduce_memory_footprint(self):
        """Moderate memory reduction for high usage situations."""
        # Reduce batch size gradually
        self.dynamic_batch_size = max(4, int(self.dynamic_batch_size * 0.75))
        
        # Reduce buffer size gradually
        self.dynamic_buffer_size = max(100, int(self.dynamic_buffer_size * 0.9))
        
        # Trigger cleanup
        self.cleanup_memory()
        
        print(f"[MEMORY] Reduced footprint: batch_size={self.dynamic_batch_size}, buffer_size={self.dynamic_buffer_size}")
    
    def cleanup_memory(self):
        """Comprehensive memory cleanup to prevent CUDA OOM errors."""
        # Clear gradients
        self.optimizer.zero_grad()
        
        # Clear CUDA cache
        if torch.cuda.is_available():
            torch.cuda.empty_cache()
            torch.cuda.synchronize()
        
        # Force garbage collection
        gc.collect()
        
        # FIXED: Intelligent replay buffer cleanup that preserves active shots
        if len(self.replay_buffer) > self.dynamic_buffer_size:
            # Separate samples into labeled, unlabeled-tracked, and unlabeled-untracked
            labeled_samples = [(i, sample) for i, sample in enumerate(self.replay_buffer) if sample[1] != -1.0]
            active_shot_indices = set(self.shot_tracking.values())
            unlabeled_active = [(i, sample) for i, sample in enumerate(self.replay_buffer) 
                              if sample[1] == -1.0 and i in active_shot_indices]
            unlabeled_orphaned = [(i, sample) for i, sample in enumerate(self.replay_buffer) 
                                if sample[1] == -1.0 and i not in active_shot_indices]
            
            # Keep ALL active shots and recent labeled samples
            target_labeled = max(0, self.dynamic_buffer_size - len(unlabeled_active))
            
            if len(labeled_samples) > target_labeled:
                keep_labeled = labeled_samples[-target_labeled:]  # Keep most recent labeled
            else:
                keep_labeled = labeled_samples
            
            # Combine samples to keep
            keep_indices = set(idx for idx, _ in keep_labeled + unlabeled_active)
            
            # Update shot tracking with new indices
            index_mapping = {}
            new_idx = 0
            for old_idx in sorted(keep_indices):
                index_mapping[old_idx] = new_idx
                new_idx += 1
            
            # Update shot tracking
            for shot_id, old_replay_idx in list(self.shot_tracking.items()):
                if old_replay_idx in index_mapping:
                    self.shot_tracking[shot_id] = index_mapping[old_replay_idx]
                else:
                    # This shouldn't happen with our logic, but clean up if it does
                    del self.shot_tracking[shot_id]
            
            # Rebuild buffer
            new_buffer = []
            for old_idx in sorted(keep_indices):
                new_buffer.append(self.replay_buffer[old_idx])
            self.replay_buffer = new_buffer
            
            print(f"[Cleanup] Buffer: {len(self.replay_buffer)}, Labeled: {len(keep_labeled)}, Active: {len(unlabeled_active)}, Orphaned: {len(unlabeled_orphaned)} (removed)")
    
    def periodic_memory_cleanup(self):
        """Perform memory cleanup periodically during long runs."""
        self.iteration_count += 1
        
        # Check memory status every iteration in emergency mode, every 5 otherwise
        check_interval = 1 if self.emergency_mode else 5
        if self.iteration_count % check_interval == 0:
            self.check_memory_status()
            
        # Regular cleanup
        if self.iteration_count % self.cleanup_interval == 0:
            self.cleanup_memory()
            
        # Print memory stats every 25 iterations (more frequent monitoring)
        if self.iteration_count % 25 == 0 and torch.cuda.is_available():
            memory_allocated = torch.cuda.memory_allocated() / 1024**3  # GB
            memory_reserved = torch.cuda.memory_reserved() / 1024**3   # GB
            status = "🚨 EMERGENCY" if self.emergency_mode else "✅ Normal"
            print(f"[Memory] Iter {self.iteration_count}: {memory_allocated:.2f}GB allocated, {memory_reserved:.2f}GB reserved ({status})")
    
    def save_checkpoint_if_needed(self):
        """Save checkpoint every 100 iterations for recovery."""
        if self.iteration_count % 100 == 0:
            checkpoint_path = MD_DIR / f"checkpoint_{self.iteration_count}.pt"
            checkpoint = {
                'iteration': self.iteration_count,
                'model_state_dict': self.model.state_dict(),
                'optimizer_state_dict': self.optimizer.state_dict(),
                'results_summary': dict(self.results_summary),
                'dispatched_shots': self.dispatched_shots,
                'replay_buffer_size': len(self.replay_buffer)
            }
            torch.save(checkpoint, checkpoint_path)
            print(f"[Checkpoint] Saved to {checkpoint_path}")
            
            # Keep only the last 3 checkpoints to save disk space
            checkpoints = sorted(MD_DIR.glob("checkpoint_*.pt"))
            if len(checkpoints) > 3:
                for old_checkpoint in checkpoints[:-3]:
                    old_checkpoint.unlink()

    def select_shooting_frame(self):
        # Intelligent path selection: use elite trajectories if available
        source_path = self.current_path
        
        # 40% chance to explore from elite trajectories if we have them
        if self.elite_trajectories and random.random() < 0.4:
            # Select from elite trajectories based on their promise
            elite_weights = [1.0 / (1.0 + traj['distance_to_b']) for traj in self.elite_trajectories]
            selected_elite = random.choices(self.elite_trajectories, weights=elite_weights)[0]
            source_path = selected_elite['trajectory']
            print(f"🎯 Exploring from elite trajectory (dist to B: {selected_elite['distance_to_b']:.2f})")
        
        self.model.eval()
        with torch.no_grad():
            committor_values = [
                self.model(*frame_to_torch_graph(frame, self.device)).item()
                for frame in source_path[::10]
            ]
        committor_values = np.array(committor_values)
        shooting_frame_idx = np.argmin(np.abs(committor_values - 0.5))
        shooting_frame = source_path[shooting_frame_idx]

        n_atoms = self.topology.n_atoms
        # Get masses in amu (Da) - MDTraj returns masses as floats in Da
        masses = np.array([a.element.mass for a in self.topology.atoms]).reshape(-1, 1)

        # Calculate thermal velocities: v = sqrt(kT/m)
        # Use a simple conversion factor for proper velocity scale in nm/ps
        velocity_scale = (
            np.sqrt(310.0 / masses) * 0.1
        )  # Empirical scale factor for nm/ps
        velocities = np.random.randn(n_atoms, 3) * velocity_scale

        return shooting_frame, velocities
    
    def calculate_distance_to_state_b(self, trajectory_frames):
        """Calculate minimum distance to State B throughout trajectory for all state definitions."""
        if not trajectory_frames:
            return float('inf')
            
        # Create temporary trajectory to analyze
        temp_traj = md.Trajectory(np.array(trajectory_frames), self.topology)
        
        best_distance_to_b = float('inf')
        
        # Check all state definitions to find the closest approach to any State B
        for definition in STATE_DEFINITIONS:
            # Select atoms for distance calculation
            atom_indices = []
            for residue_sel in definition["residues"]:
                selected = temp_traj.topology.select(residue_sel)
                if len(selected) > 0:
                    atom_indices.extend(selected)
            
            if len(atom_indices) < 2:
                continue
                
            # Calculate pairwise distances between selected atoms
            distances = md.compute_distances(temp_traj, [[atom_indices[0], atom_indices[1]]])
            min_distance = np.min(distances)
            
            # State B threshold varies by definition
            target_distance = definition["state_b_threshold"].value_in_unit(nanometers)
            distance_to_target = abs(min_distance - target_distance)
            
            # Keep track of the best (smallest) distance to any State B
            best_distance_to_b = min(best_distance_to_b, distance_to_target)
        
        return best_distance_to_b
    
    def update_elite_buffer(self, trajectory_frames, distance_to_b):
        """Update elite buffer with promising near-miss trajectories."""
        if distance_to_b < self.best_distance_to_b * 1.2:  # Within 20% of best attempt
            elite_entry = {
                'trajectory': md.Trajectory(np.array(trajectory_frames), self.topology),
                'distance_to_b': distance_to_b,
                'timestamp': self.dispatched_shots
            }
            
            self.elite_trajectories.append(elite_entry)
            
            # Keep only the best trajectories
            self.elite_trajectories.sort(key=lambda x: x['distance_to_b'])
            self.elite_trajectories = self.elite_trajectories[:self.elite_buffer_size]
            
            if distance_to_b < self.best_distance_to_b:
                self.best_distance_to_b = distance_to_b
                print(f"🔥 NEW BEST near-miss! Distance to State B: {distance_to_b:.3f} nm")
                
                # Consider switching to this path if it's significantly better
                if distance_to_b < self.best_distance_to_b * 0.8:  # 20% improvement
                    self.current_path = elite_entry['trajectory']
                    self.current_path_age = 0
                    print("🔄 Switching to new promising path!")

    def train_model(self):
        valid_samples = [item for item in self.replay_buffer if item[1] != -1]

        # Start training as soon as we have at least 2 samples, but use smaller batches initially
        min_samples = 2
        if len(valid_samples) < min_samples:
            return

        # Use dynamic batch size based on memory conditions
        max_batch = self.dynamic_batch_size if not self.emergency_mode else max(2, self.dynamic_batch_size // 2)
        
        if len(valid_samples) < max_batch:
            batch_size = min(len(valid_samples), max(2, len(valid_samples) // 2))
        else:
            batch_size = max_batch
            
        # Skip training if memory is critical
        if torch.cuda.is_available():
            memory_reserved = torch.cuda.memory_reserved() / 1024**3
            if memory_reserved > MAX_SAFE_MEMORY:
                print(f"[MEMORY] Skipping training - memory too high: {memory_reserved:.1f}GB")
                return None

        self.model.train()
        batch = random.sample(valid_samples, batch_size)

        self.optimizer.zero_grad()
        total_loss = 0
        
        # Process batch items one by one for better memory management
        for i, (frame_data, target_val) in enumerate(batch):
            g, l_g = frame_data  # <-- NEW (Unpack the two graphs)
            target = torch.tensor([target_val], dtype=torch.float32, device=self.device)
            prediction = self.model(g, l_g).squeeze(1)  # <-- NEW (Pass graphs to the model)
            loss = self.loss_fn(prediction, target)
            
            # Scale loss by batch size for proper gradient averaging
            scaled_loss = loss / len(batch)
            scaled_loss.backward()
            total_loss += loss.item()
            
            # Clear intermediate tensors to free memory
            del g, l_g, target, prediction, loss, scaled_loss # <-- NEW (Clean up graphs)
        
        # Clip gradients to prevent exploding gradients
        torch.nn.utils.clip_grad_norm_(self.model.parameters(), max_norm=1.0)
        self.optimizer.step()
        
        # Clear gradients after step
        self.optimizer.zero_grad()
        
        return total_loss / len(batch)

    def run(self):
        self.start_workers()

        pbar = tqdm(total=N_TOTAL_SHOTS, desc="Starting Shots", file=sys.stdout)
        stop_event = threading.Event()
        monitor_thread = threading.Thread(
            target=monitor_gpu_utilization, args=(stop_event, pbar)
        )
        monitor_thread.daemon = True
        monitor_thread.start()

        try:
            for _ in range(N_WORKERS * 2):
                if self.dispatched_shots >= N_TOTAL_SHOTS:
                    break
                self.dispatch_shot()
                pbar.set_description(
                    f"GPU Util: N/A | Dispatched: {self.dispatched_shots}"
                )

            while self.results_summary.total() < N_TOTAL_SHOTS:
                shot_index, result, trajectory_frames = self.result_queue.get()

                self.results_summary[result] += 1
                self.current_path_age += 1
                
                # NEAR-MISS LEARNING: Analyze all trajectories for distance to State B
                if trajectory_frames:
                    distance_to_b = self.calculate_distance_to_state_b(trajectory_frames)
                    self.update_elite_buffer(trajectory_frames, distance_to_b)

                if result in ["State A", "State B"]:
                    # Check if the shot is still tracked in our mapping
                    if shot_index in self.shot_tracking:
                        replay_idx = self.shot_tracking[shot_index]
                        if replay_idx < len(self.replay_buffer):
                            original_frame_data, _ = self.replay_buffer[replay_idx]
                            target_val = 0.0 if result == "State A" else 1.0
                            self.replay_buffer[replay_idx] = (original_frame_data, target_val)
                            print(f"State {result[-1]} (dist: {distance_to_b:.3f}nm)", flush=True)
                        else:
                            print(f"State {result[-1]} (buffer index {replay_idx} out of range)", flush=True)
                        # Clean up the tracking entry
                        del self.shot_tracking[shot_index]
                    else:
                        # Shot was from an older buffer that got trimmed, skip labeling
                        print(f"State {result[-1]} (shot {shot_index} no longer tracked)", flush=True)

                if result == "State B":
                    print("🎉 SUCCESS! Found complete path to State B!", flush=True)
                    new_path = md.Trajectory(np.array(trajectory_frames), self.topology)
                    new_path.save_dcd(MD_DIR / f"path_to_B_{shot_index}.dcd")
                    # Store in path history for diversity
                    self.path_history.append(new_path)
                    self.current_path = new_path
                    self.current_path_age = 0
                
                # Switch paths if current one is getting stale without progress
                elif self.current_path_age > self.max_path_age and self.elite_trajectories:
                    best_elite = min(self.elite_trajectories, key=lambda x: x['distance_to_b'])
                    print(f"🔄 Path refresh: switching to elite trajectory (dist: {best_elite['distance_to_b']:.3f}nm) after {self.current_path_age} stale shots")
                    self.current_path = best_elite['trajectory']
                    self.current_path_age = 0

                # Train model and show loss after every shot
                avg_loss = self.train_model()
                if avg_loss is not None:
                    valid_count = len(
                        [item for item in self.replay_buffer if item[1] != -1]
                    )
                    print(
                        f"Loss: {avg_loss:.4f} (trained on {min(valid_count, TRAINING_BATCH_SIZE)} samples)"
                    )
                    pbar.set_postfix(loss=f"{avg_loss:.4f}")
                else:
                    valid_count = len(
                        [item for item in self.replay_buffer if item[1] != -1]
                    )
                    print(
                        f"Training skipped - need at least 2 valid samples, have {valid_count}"
                    )
                    pbar.set_postfix(loss="N/A (need more data)")
                
                # Periodic memory cleanup and checkpointing
                self.periodic_memory_cleanup()
                self.save_checkpoint_if_needed()

                # Check if we should dispatch new shots based on memory
                memory_ok = True
                if torch.cuda.is_available():
                    memory_reserved = torch.cuda.memory_reserved() / 1024**3
                    if memory_reserved > MAX_SAFE_MEMORY:
                        print(f"[MEMORY] Pausing shot dispatch - memory at {memory_reserved:.1f}GB")
                        memory_ok = False
                        
                if self.dispatched_shots < N_TOTAL_SHOTS and memory_ok:
                    self.dispatch_shot()
                    pbar.set_description(
                        f"GPU Util: N/A | Dispatched: {self.dispatched_shots}"
                    )

                pbar.update(1)
                pbar.set_description(
                    f"GPU Util: N/A | Completed: {self.results_summary.total()}"
                )

        finally:
            stop_event.set()
            monitor_thread.join(timeout=2)
            pbar.close()
            print("\n--> Stopping workers...")
            self.stop_workers()
            self.print_summary()

    def dispatch_shot(self):
        shot_index = self.dispatched_shots
        frame, velocities = self.select_shooting_frame()

        frame_data_graph = frame_to_torch_graph(frame, self.device)
        self.replay_buffer.append((frame_data_graph, -1.0))
        
        # CRITICAL FIX: Track shot mapping BEFORE cleanup can remove it
        replay_buffer_index = len(self.replay_buffer) - 1
        self.shot_tracking[shot_index] = replay_buffer_index

        task = (shot_index, frame.xyz[0], velocities)
        self.task_queue.put(task)
        self.dispatched_shots += 1

    def print_summary(self):
        print("\n🎯 Sampling complete!")
        print("--- Summary ---")
        total = self.results_summary.total()
        if total == 0:
            return
        for state, count in self.results_summary.items():
            print(f"  {state:<12}: {count:4d} ({count/total*100:5.1f}%)")
            
        print("\n--- Elite Buffer Analysis ---")
        if self.elite_trajectories:
            print(f"  Elite trajectories found: {len(self.elite_trajectories)}")
            print(f"  Best distance to State B: {self.best_distance_to_b:.4f} nm")
            elite_distances = [f"{t['distance_to_b']:.3f}" for t in self.elite_trajectories[:5]]
            print(f"  Elite distances: {elite_distances}")
        else:
            print("  No elite trajectories found - consider adjusting search strategy")
            
        if self.path_history:
            print(f"\n  Successful paths found: {len(self.path_history)}")
            
        print("\n--- Memory Management Summary ---")
        if torch.cuda.is_available():
            final_allocated = torch.cuda.memory_allocated() / 1024**3
            final_reserved = torch.cuda.memory_reserved() / 1024**3
            print(f"  Final memory usage: {final_allocated:.2f}GB allocated, {final_reserved:.2f}GB reserved")
        print(f"  Memory warnings issued: {self.memory_warnings}")
        print(f"  Emergency mode triggered: {'Yes' if self.emergency_mode else 'No'}")
        print(f"  Final buffer size: {len(self.replay_buffer)} (limit: {self.dynamic_buffer_size})")
        print(f"  Final batch size: {self.dynamic_batch_size}")

    def run_controlled(self):
        """Run a short controlled simulation without spawning worker processes.
        This simulates worker execution inline by taking tasks from the task_queue
        and returning a simple trajectory (single-frame) and 'Timeout' result.
        Useful for testing dispatch/labeling/training logic without OpenMM/System XML.
        """
        print(f"--> Running controlled simulation: shots={N_TOTAL_SHOTS}")
        pbar = tqdm(total=N_TOTAL_SHOTS, desc="Controlled Shots", file=sys.stdout)

        # Dispatch initial shots
        for _ in range(N_WORKERS * 2):
            if self.dispatched_shots >= N_TOTAL_SHOTS:
                break
            # Do not spawn real workers; dispatch will push tasks to task_queue
            self.dispatch_shot()
            pbar.set_description(f"Dispatched: {self.dispatched_shots}")

        # Process tasks inline: simulate a worker for each queued task
        while self.results_summary.total() < N_TOTAL_SHOTS:
            try:
                task = self.task_queue.get(timeout=1.0)
            except Exception:
                # No task available yet; break to avoid deadlock in test
                break

            if task is None:
                continue

            shot_index, positions_nm, velocities_nm = task
            # Simulate a short trajectory: single frame equal to provided positions
            trajectory_frames = [positions_nm]
            final_state = "Timeout"

            # Put simulated result into result_queue for standard processing
            self.result_queue.put((shot_index, final_state, trajectory_frames))

            # Now run the same processing loop body as in run() for one shot
            shot_index, result, trajectory_frames = self.result_queue.get()

            self.results_summary[result] += 1
            self.current_path_age += 1

            if trajectory_frames:
                distance_to_b = self.calculate_distance_to_state_b(trajectory_frames)
                self.update_elite_buffer(trajectory_frames, distance_to_b)

            if result in ["State A", "State B"]:
                if shot_index in self.shot_tracking:
                    replay_idx = self.shot_tracking[shot_index]
                    if replay_idx < len(self.replay_buffer):
                        original_frame_data, _ = self.replay_buffer[replay_idx]
                        target_val = 0.0 if result == "State A" else 1.0
                        self.replay_buffer[replay_idx] = (original_frame_data, target_val)
                    del self.shot_tracking[shot_index]

            if result == "State B":
                new_path = md.Trajectory(np.array(trajectory_frames), self.topology)
                self.path_history.append(new_path)
                self.current_path = new_path
                self.current_path_age = 0

            # Train and housekeeping as in run()
            avg_loss = self.train_model()
            if avg_loss is not None:
                valid_count = len([item for item in self.replay_buffer if item[1] != -1])
                print(f"Loss: {avg_loss:.4f} (trained on {min(valid_count, TRAINING_BATCH_SIZE)} samples)")
                pbar.set_postfix(loss=f"{avg_loss:.4f}")
            else:
                valid_count = len([item for item in self.replay_buffer if item[1] != -1])
                print(f"Training skipped - need at least 2 valid samples, have {valid_count}")
                pbar.set_postfix(loss="N/A (need more data)")

            self.periodic_memory_cleanup()
            self.save_checkpoint_if_needed()

            # Dispatch more shots if needed
            if self.dispatched_shots < N_TOTAL_SHOTS:
                self.dispatch_shot()

            pbar.update(1)

        pbar.close()
        print("--> Controlled run complete")


if __name__ == "__main__":
    set_start_method("spawn", force=True)
    initial_traj = str(MD_DIR / "trajectory.dcd")

    # Dry-run option: perform preflight checks and exit successfully if they pass.
    if args.dry_run:
        print("Running dry-run preflight checks...")
        try:
            # Check CUDA and DGL availability and that a model can be moved to GPU
            import dgl

            if not torch.cuda.is_available():
                raise RuntimeError("torch reports no CUDA device available (torch.cuda.is_available() is False)")

            # Validate DGL can move a tiny graph to CUDA
            tmp_g = dgl.graph(([], []))
            tmp_g.to(torch.device("cuda:0"))

            # Validate model can be instantiated and moved to CUDA
            model = CommittorNet(embedding_dim=8, n_layers=1).to(torch.device("cuda:0"))
            del model, tmp_g
            torch.cuda.empty_cache()
            print("Dry-run preflight checks passed: CUDA + DGL + model -> cuda:0 OK")
            sys.exit(0)
        except Exception as e:
            print(f"Dry-run preflight failed: {e}", file=sys.stderr)
            sys.exit(2)

    # Startup option: instantiate AIMDRunner but do not call run(). Useful to exercise
    # initialization (topology load, model -> GPU) without starting workers or dispatching shots.
    # Apply optional overrides from CLI
    if args.shots is not None:
        N_TOTAL_SHOTS = args.shots
    if args.workers is not None:
        N_WORKERS = args.workers

    if args.startup:
        print("Running startup: instantiating AIMDRunner (no workers)...")
        try:
            runner = AIMDRunner(initial_traj)
            # Print basic info and exit
            print(f"AIMDRunner initialized. Device: {runner.device}")
            print(f"Model parameters: {sum(p.numel() for p in runner.model.parameters())}")
            sys.exit(0)
        except Exception as e:
            print(f"Startup failed: {e}", file=sys.stderr)
            sys.exit(3)

    if args.controlled:
        print("Running controlled short simulation...")
        try:
            runner = AIMDRunner(initial_traj)
            runner.run_controlled()
            sys.exit(0)
        except Exception as e:
            print(f"Controlled run failed: {e}", file=sys.stderr)
            sys.exit(4)

    if not Path(initial_traj).exists():
        print(
            "ERROR: Initial trajectory 'trajectory.dcd' not found.", file=sys.stderr
        )
        sys.exit(1)

    runner = AIMDRunner(initial_traj)
    runner.run() 
