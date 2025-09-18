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
from scipy.spatial import KDTree
from openmm import *
from openmm.app import *
from openmm.unit import *

from model import CommittorNet

# Set device
device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
print(f"[DEBUG] Using device: {device}")

# --- Configuration ---
TARGET_NAME = "GLP1R"
ARTIFACTS_DIR = Path(__file__).resolve().parents[1] / "artifacts"
MD_DIR = ARTIFACTS_DIR / "md" / TARGET_NAME

# --- State Definitions ---
STATE_A_DISTANCE = 12.0  # Angstroms
STATE_B_DISTANCE = 20.0  # Angstroms
# -------------------------

def frame_to_torch_graph(frame):
    """Converts a single mdtraj frame to a graph format for our GNN."""
    atom_types = torch.tensor([atom.element.atomic_number for atom in frame.topology.atoms], dtype=torch.long, device=device)
    coords = torch.tensor(frame.xyz[0], dtype=torch.float32, device=device)
    cutoff = 0.5 # nm
    tree = KDTree(coords)
    adj_set = tree.query_pairs(r=cutoff)
    if not adj_set:
        edge_index = torch.empty((2, 0), dtype=torch.long)
    else:
        adj = np.array(list(adj_set), dtype=np.int64)
        adj_rev = np.fliplr(adj)
        edge_index = torch.from_numpy(np.vstack((adj, adj_rev))).t().contiguous().to(device)
    return atom_types, coords, edge_index

def get_pocket_distance(simulation):
    """Calculates the current pocket distance in a simulation."""
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
        unitcell_angles=np.array([angles])
    )
    traj_frame.image_molecules(inplace=True)
    atom_selection_1 = traj_frame.topology.select('chainid 0 and resid 187 and name CA')
    atom_selection_2 = traj_frame.topology.select('chainid 0 and resid 393 and name CA')
    atom_indices = [[atom_selection_1[0], atom_selection_2[0]]]
    distance = md.compute_distances(traj_frame, atom_indices)[0, 0] * 10
    return distance * angstroms

def run_shooting_move(simulation):
    """Runs a short trajectory until it hits State A or B."""
    simulation.context.setVelocitiesToTemperature(310*kelvin)
    max_steps = 50000 
    report_interval = 2500
    trajectory_frames = []
    for i in range(max_steps // report_interval):
        model = CommittorNet().to(device)
        print(f"[DEBUG] Model is on device: {next(model.parameters()).device}")
        positions = simulation.context.getState(getPositions=True).getPositions()
        trajectory_frames.append(positions.value_in_unit(nanometer))
        distance = get_pocket_distance(simulation)
        
        if distance.value_in_unit(angstrom) < STATE_A_DISTANCE:
            return "State A", trajectory_frames
        if distance.value_in_unit(angstrom) > STATE_B_DISTANCE:
            return "State B", trajectory_frames
            
    return "Timeout", trajectory_frames

def main():
    pdb_path = MD_DIR / "prepared_system.pdb"
    state_path = MD_DIR / "minimized_state.xml"
    
    print("--> Loading prepared system...")
    pdb = PDBFile(str(pdb_path))
    forcefield = ForceField("amber14-all.xml", "amber14/tip3p.xml", "amber/lipid17.xml")
    
    system = forcefield.createSystem(pdb.topology, nonbondedMethod=PME, nonbondedCutoff=1.0*nanometer, constraints=HBonds)
    integrator = LangevinIntegrator(310*kelvin, 1.0/picosecond, 4.0*femtoseconds)
    simulation = Simulation(pdb.topology, system, integrator)
    with open(state_path, 'r') as f:
        simulation.context.setState(XmlSerializer.deserialize(f.read()))

    print("--> System loaded successfully.")
    
    current_path = md.load(str(MD_DIR / "trajectory.dcd"), top=str(pdb_path))
    
    print("--> Initializing CommittorNet AI model...")
    model = CommittorNet()
    optimizer = optim.Adam(model.parameters(), lr=1e-4)
    loss_fn = nn.BCELoss()

    num_shots = 1000
    results = []
                target = torch.tensor([0.0], device=device) if result == "State A" else torch.tensor([1.0], device=device)
    
    for i in range(num_shots):
        # --- AI-GUIDED SELECTION ---
        with torch.no_grad():
            committor_values = []
                print(f"[DEBUG] Prediction tensor device: {prediction.device}, Target tensor device: {target.device}")
            for frame in current_path:
                atom_types, coords, edge_index = frame_to_torch_graph(frame)
                p = model(atom_types, coords, edge_index).item()
                committor_values.append(p)
        
        committor_values = np.array(committor_values)
        shooting_frame_index = np.argmin(np.abs(committor_values - 0.5))
        shooting_frame = current_path[shooting_frame_index]
        # ---------------------------
        
        print(f"\n--- Shot {i+1}/{num_shots} (AI selected frame {shooting_frame_index} with p={committor_values[shooting_frame_index]:.4f}) ---")
        
        try:
            atom_types, coords, edge_index = frame_to_torch_graph(shooting_frame)
            simulation.context.setPositions(shooting_frame.openmm_positions(0))
            result, new_frames = run_shooting_move(simulation)
            results.append(result)
            print(f"--> Shot {i+1} result: {result}")

            if result == "State B":
                print("🎉 SUCCESS! Found a path to State B! Saving trajectory...")
                new_path = md.Trajectory(np.array(new_frames), md.Topology.from_openmm(pdb.topology))
                new_path.save_dcd(str(MD_DIR / f"path_to_B_{i+1}.dcd"))
                current_path = new_path 
            
            target = torch.tensor([0.0]) if result == "State A" else torch.tensor([1.0])
            if result == "Timeout":
                continue
            
            optimizer.zero_grad()
            prediction = model(atom_types, coords, edge_index)
            loss = loss_fn(prediction, target)
            loss.backward()
            optimizer.step()
            print(f"    AI Training: loss = {loss.item():.4f}")

        except OpenMMException as e:
            print(f"🔥 WARNING: Shot {i+1} failed with an OpenMM error (likely NaN). Skipping shot.")
            results.append("Failed")


    print("\n✅ Sampling complete!")
    summary = Counter(results)
    print("--- Summary ---")
    print(f"Paths returned to State A: {summary.get('State A', 0)}")
    print(f"Paths reached State B:     {summary.get('State B', 0)}")
    print(f"Paths timed out:           {summary.get('Timeout', 0)}")
    print(f"Paths failed (NaN):        {summary.get('Failed', 0)}")

if __name__ == "__main__":
    main()
