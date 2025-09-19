# scripts/run_path_sampling.py
import time
import random
import sys
import subprocess
import threading
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

# --- Configuration ---
# ==============================================================================
TARGET_NAME = "GLP1R"
REPO_ROOT = Path(__file__).resolve().parents[1]
ARTIFACTS_DIR = REPO_ROOT / "artifacts"
MD_DIR = ARTIFACTS_DIR / "md" / TARGET_NAME
MD_DIR.mkdir(parents=True, exist_ok=True)

N_WORKERS = 2
N_TOTAL_SHOTS = 1000
LEARNING_RATE = 1e-4
REPLAY_BUFFER_SIZE = 200
TRAINING_BATCH_SIZE = 32

WORKER_CONFIG = {
    "pdb_path": str(MD_DIR / "prepared_system.pdb"),
    "state_path": str(MD_DIR / "minimized_state.xml"),
    "system_path": str(MD_DIR / "system.xml"),
    "state_a_residues": (
        "chainid 0 and resid 187 and name CA",
        "chainid 0 and resid 393 and name CA",
    ),
    "state_a_threshold": 1.2 * nanometers,
    "state_b_threshold": 2.0 * nanometers,
    "shooting_move_steps": 10000,
    "report_interval": 1000,
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

        print("--> Initializing CommittorNet AI model on cuda:0...")
        self.device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
        self.model = CommittorNet().to(self.device)
        self.optimizer = optim.Adam(self.model.parameters(), lr=LEARNING_RATE)
        self.loss_fn = nn.BCELoss()
        self.replay_buffer = []

        self.task_queue = Queue()
        self.result_queue = Queue()
        self.workers = []

        self.results_summary = Counter()
        self.dispatched_shots = 0

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

    def select_shooting_frame(self):
        self.model.eval()
        with torch.no_grad():
            committor_values = [
                self.model(*frame_to_torch_graph(frame, self.device)).item()
                for frame in self.current_path
            ]
        committor_values = np.array(committor_values)
        shooting_frame_idx = np.argmin(np.abs(committor_values - 0.5))
        shooting_frame = self.current_path[shooting_frame_idx]

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

    def train_model(self):
        valid_samples = [item for item in self.replay_buffer if item[1] != -1]

        # Start training as soon as we have at least 2 samples, but use smaller batches initially
        min_samples = 2
        if len(valid_samples) < min_samples:
            return

        # Use adaptive batch size: start small, grow to full batch size
        if len(valid_samples) < TRAINING_BATCH_SIZE:
            batch_size = min(len(valid_samples), max(2, len(valid_samples) // 2))
        else:
            batch_size = TRAINING_BATCH_SIZE

        self.model.train()
        batch = random.sample(valid_samples, batch_size)

        self.optimizer.zero_grad()
        total_loss = 0
        for frame_data, target_val in batch:
            atom_types, coords, edge_index = frame_data
            target = torch.tensor([target_val], dtype=torch.float32, device=self.device)
            prediction = self.model(atom_types, coords, edge_index)
            loss = self.loss_fn(prediction, target)
            loss.backward()
            total_loss += loss.item()

        self.optimizer.step()
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

                if result in ["State A", "State B"]:
                    original_frame_data, _ = self.replay_buffer[shot_index]
                    target_val = 0.0 if result == "State A" else 1.0
                    self.replay_buffer[shot_index] = (original_frame_data, target_val)
                    print(f"State {result[-1]}", flush=True)

                if result == "State B":
                    print("✅ New path to State B!", flush=True)
                    new_path = md.Trajectory(np.array(trajectory_frames), self.topology)
                    new_path.save_dcd(MD_DIR / f"path_to_B_{shot_index}.dcd")
                    self.current_path = new_path

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

                if self.dispatched_shots < N_TOTAL_SHOTS:
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

        task = (shot_index, frame.xyz[0], velocities)
        self.task_queue.put(task)
        self.dispatched_shots += 1

    def print_summary(self):
        print("\n✅ Sampling complete!")
        print("--- Summary ---")
        total = self.results_summary.total()
        if total == 0:
            return
        for state, count in self.results_summary.items():
            print(f"  {state:<12}: {count:4d} ({count/total*100:5.1f}%)")


if __name__ == "__main__":
    set_start_method("spawn", force=True)
    initial_traj = str(MD_DIR / "trajectory.dcd")

    if not Path(initial_traj).exists():
        print(
            "🔥 Error: Initial trajectory 'trajectory.dcd' not found.", file=sys.stderr
        )
        sys.exit(1)

    runner = AIMDRunner(initial_traj)
    runner.run()
