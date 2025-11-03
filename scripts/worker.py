import mdtraj as md
import numpy as np
from torch.multiprocessing import Queue

from openmm import (
    Platform,
    XmlSerializer,
    LangevinIntegrator,
    OpenMMException,
)
from openmm.app import PDBFile, element, Simulation
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
        system_path = config["system_path"]
        state_definitions = config["state_definitions"]  # List of state definitions
        shooting_move_steps = config["shooting_move_steps"]
        report_interval = config["report_interval"]

        # --- Setup Simulation on Assigned GPU (fallback if CUDA not available) ---
        try:
            platform = Platform.getPlatformByName("CUDA")
            properties = {"DeviceIndex": str(gpu_index), "Precision": "mixed"}
        except Exception:
            # Fallback to Reference platform if CUDA platform is missing
            try:
                platform = Platform.getPlatformByName("Reference")
            except Exception:
                platform = Platform.getPlatformByName("CPU")
            properties = {}

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
        # --- Define Multi-State Boundary Checking ---
        md_topology = md.Topology.from_openmm(pdb.topology)
        
        def check_state_boundaries(simulation, state_definitions, md_topology):
            """
            Calculates distances for all state definitions and checks if any boundary is crossed.
            Returns 'State A', 'State B', or None.
            """
            state = simulation.context.getState(getPositions=True)
            positions_nm = state.getPositions(asNumpy=True).value_in_unit(nanometer)
            
            traj = md.Trajectory([positions_nm], topology=md_topology)
            
            for definition in state_definitions:
                # Select atoms for the current definition
                atom1_idx = md_topology.select(definition["residues"][0])[0]
                atom2_idx = md_topology.select(definition["residues"][1])[0]
                atom_indices = np.array([[atom1_idx, atom2_idx]])
                
                # Calculate distance
                distance = md.compute_distances(traj, atom_indices)[0, 0] * nanometers
                
                # Check thresholds
                if distance < definition["state_a_threshold"]:
                    return "State A"
                if distance > definition["state_b_threshold"]:
                    return "State B"
                    
            return None  # No boundary was crossed

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
                    
                    # Check all state definitions at once
                    current_state = check_state_boundaries(simulation, state_definitions, md_topology)
                    if current_state is not None:
                        final_state = current_state
                        break  # Exit the loop as soon as any boundary is hit
                        
                except OpenMMException:
                    final_state = "Failed"
                    break
            result_queue.put((shot_index, final_state, trajectory_frames))

    except Exception as e:
        print(f"[Worker {gpu_index}] CRITICAL ERROR: {e}", flush=True)
        if "shot_index" in locals():
            result_queue.put((shot_index, "Failed", []))
