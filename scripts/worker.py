# scripts/worker.py
import mdtraj as md
import numpy as np
from torch.multiprocessing import Queue

from openmm import (
    Platform,
    Simulation,
    XmlSerializer,
    LangevinIntegrator,
    OpenMMException,
)
from openmm.app import PDBFile, element
from openmm.unit import nanometer, nanometers, picosecond, femtoseconds, kelvin, amu


def simulation_worker(
    gpu_index: int, task_queue: Queue, result_queue: Queue, config: dict
):
    """
    A worker process that runs OpenMM simulations on a specific GPU.
    It now loads a pre-serialized System object to avoid PDB parsing errors.
    """
    try:
        # --- Unpack Configuration ---
        pdb_path = config["pdb_path"]
        state_path = config["state_path"]
        system_path = config["system_path"]  # New
        # ... (rest of config variables)
        state_a_residues = config["state_a_residues"]
        state_a_threshold = config["state_a_threshold"]
        state_b_threshold = config["state_b_threshold"]
        shooting_move_steps = config["shooting_move_steps"]
        report_interval = config["report_interval"]

        # --- Setup Simulation on Assigned GPU ---
        platform = Platform.getPlatformByName("CUDA")
        properties = {"DeviceIndex": str(gpu_index), "Precision": "mixed"}

        pdb = PDBFile(pdb_path)  # Still need for topology info

        # --- FIX: Load the pre-built system from XML ---
        with open(system_path, "r") as f:
            system = XmlSerializer.deserialize(f.read())

        integrator = LangevinIntegrator(
            310 * kelvin, 1.0 / picosecond, 4.0 * femtoseconds
        )
        # Add HMR for stability with 4fs timestep if not already in system
        if system.getParticleMass(1).value_in_unit(amu) < 2.0:
            for atom in pdb.topology.atoms():
                if atom.element == element.hydrogen:
                    mass = system.getParticleMass(atom.index)
                    for atom2 in atom.residue.atoms():
                        if atom2.index != atom.index:
                            system.setParticleMass(
                                atom2.index,
                                system.getParticleMass(atom2.index) - mass + 4 * amu,
                            )
                            break
                    system.setParticleMass(atom.index, 4 * amu)

        simulation = Simulation(pdb.topology, system, integrator, platform, properties)

        with open(state_path, "r") as f:
            simulation.context.setState(XmlSerializer.deserialize(f.read()))

        # ... (The rest of the worker function is the same) ...
        # --- Define Distance Calculation ---
        md_topology = md.Topology.from_openmm(pdb.topology)
        atom1_idx = md_topology.select(state_a_residues[0])[0]
        atom2_idx = md_topology.select(state_a_residues[1])[0]
        atom_indices = np.array([[atom1_idx, atom2_idx]])

        def get_pocket_distance():
            state = simulation.context.getState(getPositions=True)
            positions_nm = state.getPositions(asNumpy=True).value_in_unit(nanometer)

            # Create trajectory without unit cell info to avoid API issues
            traj = md.Trajectory([positions_nm], topology=md_topology)
            return md.compute_distances(traj, atom_indices)[0, 0] * nanometers

        # --- Main Worker Loop ---
        while True:
            task = task_queue.get()
            if task is None:
                break

            shot_index, positions_nm, velocities_nm = task

            simulation.context.setPositions(positions_nm * nanometer)
            simulation.context.setVelocities(velocities_nm * (nanometer / picosecond))

            trajectory_frames = []
            final_state = "Timeout"

            total_steps = shooting_move_steps // report_interval

            for step_block in range(total_steps):
                try:
                    simulation.step(report_interval)
                    pos = (
                        simulation.context.getState(getPositions=True)
                        .getPositions(asNumpy=True)
                        .value_in_unit(nanometer)
                    )
                    trajectory_frames.append(pos)
                    distance = get_pocket_distance()
                    distance_nm = distance.value_in_unit(nanometer)
                    threshold_a_nm = state_a_threshold.value_in_unit(nanometer)
                    threshold_b_nm = state_b_threshold.value_in_unit(nanometer)

                    if distance < state_a_threshold:
                        final_state = "State A"
                        break
                    if distance > state_b_threshold:
                        final_state = "State B"
                        break
                except OpenMMException:
                    final_state = "Failed"
                    break
            result_queue.put((shot_index, final_state, trajectory_frames))

    except Exception as e:
        print(f"[Worker {gpu_index}] CRITICAL ERROR: {e}", flush=True)
        if "shot_index" in locals():
            result_queue.put((shot_index, "Failed", []))
