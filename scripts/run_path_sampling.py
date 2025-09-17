# scripts/run_path_sampling.py
import mdtraj as md
import numpy as np
from pathlib import Path
import random
import sys
from collections import Counter
import torch
import torch.nn as nn
import torch.optim as optim
from scipy.spatial import KDTree  # New, efficient import
from openmm import *
from openmm.app import *
from openmm.unit import *

from model import CommittorNet

# --- Configuration ---
TARGET_NAME = "GLP1R"
ARTIFACTS_DIR = Path(__file__).resolve().parents[1] / "artifacts"
MD_DIR = ARTIFACTS_DIR / "md" / TARGET_NAME

# --- State Definitions ---
STATE_A_DISTANCE = 12.0  # Angstroms
STATE_B_DISTANCE = 20.0  # Angstroms
# -------------------------


# --- THIS FUNCTION IS THE FIX ---
def frame_to_torch_graph(frame):
    """Converts a single mdtraj frame to a graph format for our GNN."""
    atom_types = torch.tensor(
        [atom.element.atomic_number for atom in frame.topology.atoms], dtype=torch.long
    )
    coords = torch.tensor(frame.xyz[0], dtype=torch.float32)  # in nanometers
    cutoff = 0.5  # nm

    # Build a KD-Tree for efficient neighbor searching
    tree = KDTree(coords)
    # This efficiently finds all pairs within the cutoff without a huge matrix
    adj_set = tree.query_pairs(r=cutoff)

    # Convert the set of pairs to the COO format for PyTorch
    if not adj_set:
        edge_index = torch.empty((2, 0), dtype=torch.long)
    else:
        adj = np.array(list(adj_set), dtype=np.int64)
        # Add edges in both directions to make the graph undirected
        adj_rev = np.fliplr(adj)
        edge_index = torch.from_numpy(np.vstack((adj, adj_rev))).t().contiguous()

    return atom_types, coords, edge_index


# -----------------------------


def get_pocket_distance(simulation):
    state = simulation.context.getState(getPositions=True)
    positions = state.getPositions()
    box_vectors = state.getPeriodicBoxVectors(asNumpy=True).value_in_unit(nanometer)
    mdtraj_topology = md.Topology.from_openmm(simulation.topology)
    lengths = np.array([box_vectors[0][0], box_vectors[1][1], box_vectors[2][2]])
    angles = np.array([90.0, 90.0, 90.0])
    traj_frame = md.Trajectory(
        xyz=np.array([positions.value_in_unit(nanometer)]),
        topology=mdtraj_topology,
        unitcell_lengths=np.array([lengths]),
        unitcell_angles=np.array([angles]),
    )
    traj_frame.image_molecules(inplace=True)
    atom_selection_1 = traj_frame.topology.select("chainid 0 and resid 187 and name CA")
    atom_selection_2 = traj_frame.topology.select("chainid 0 and resid 393 and name CA")
    atom_indices = [[atom_selection_1[0], atom_selection_2[0]]]
    distance = md.compute_distances(traj_frame, atom_indices)[0, 0] * 10
    return distance * angstroms


def run_shooting_move(simulation):
    simulation.context.setVelocitiesToTemperature(310 * kelvin)
    max_steps = 50000
    report_interval = 2500
    for i in range(max_steps // report_interval):
        simulation.step(report_interval)
        distance = get_pocket_distance(simulation)
        if (i + 1) % 4 == 0:
            print(
                f"    ... shot progress {((i+1)*report_interval*4)/1000:.1f} ps, dist: {distance.value_in_unit(angstrom):.2f} Å"
            )
        if distance.value_in_unit(angstrom) < STATE_A_DISTANCE:
            return "State A"
        if distance.value_in_unit(angstrom) > STATE_B_DISTANCE:
            return "State B"
    return "Timeout"


def main():
    pdb_path = MD_DIR / "prepared_system.pdb"
    state_path = MD_DIR / "minimized_state.xml"

    print("--> Loading prepared system...")
    pdb = PDBFile(str(pdb_path))
    forcefield = ForceField("amber14-all.xml", "amber14/tip3p.xml", "amber/lipid17.xml")

    system = forcefield.createSystem(
        pdb.topology, nonbondedMethod=PME, nonbondedCutoff=1.0 * nanometer, constraints=HBonds
    )
    integrator = LangevinIntegrator(310 * kelvin, 1.0 / picosecond, 4.0 * femtoseconds)
    simulation = Simulation(pdb.topology, system, integrator)
    with open(state_path, "r") as f:
        simulation.context.setState(XmlSerializer.deserialize(f.read()))

    print("--> System loaded successfully.")

    initial_traj = md.load(str(MD_DIR / "trajectory.dcd"), top=str(pdb_path))

    print("--> Initializing CommittorNet AI model...")
    model = CommittorNet()
    optimizer = optim.Adam(model.parameters(), lr=1e-4)
    loss_fn = nn.BCELoss()

    num_shots = 500
    results = []
    print(f"\n--> Starting AI-guided sampling loop for {num_shots} shots...")

    for i in range(num_shots):
        shooting_frame_index = random.randint(0, len(initial_traj) - 1)
        shooting_frame = initial_traj[shooting_frame_index]
        print(f"\n--- Shot {i+1}/{num_shots} (from frame {shooting_frame_index}) ---")

        atom_types, coords, edge_index = frame_to_torch_graph(shooting_frame)
        simulation.context.setPositions(shooting_frame.openmm_positions(0))
        result = run_shooting_move(simulation)
        results.append(result)
        print(f"--> Shot {i+1} result: {result}")

        if result == "State A":
            target = torch.tensor([0.0])
        elif result == "State B":
            target = torch.tensor([1.0])
        else:
            continue

        optimizer.zero_grad()
        prediction = model(atom_types, coords, edge_index)
        loss = loss_fn(prediction.squeeze(), target)
        loss.backward()
        optimizer.step()
        print(f"    AI Training: loss = {loss.item():.4f}")

    print("\n✅ Sampling complete!")
    summary = Counter(results)
    print("--- Summary ---")
    print(f"Paths returned to State A: {summary.get('State A', 0)}")
    print(f"Paths reached State B:     {summary.get('State B', 0)}")
    print(f"Paths timed out:           {summary.get('Timeout', 0)}")


if __name__ == "__main__":
    main()
