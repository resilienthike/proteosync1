#!/usr/bin/env python3
"""
worker.py

OpenMM Molecular Dynamics Worker for Path Sampling

This module provides the simulation worker process that executes OpenMM molecular
dynamics simulations for two-way shooting moves in transition path sampling. Each
worker runs on a dedicated GPU and processes shooting tasks from a queue.

The worker performs:
1. Bidirectional MD propagation from shooting points
2. State classification using collective variables
3. Trajectory fragment collection for accepted paths

Author: ProteOSync Pipeline
Date: November 8, 2025
"""

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
# --- Define Multi-State Boundary Checking ---
        md_topology = md.Topology.from_openmm(pdb.topology)
        
        def check_state_boundaries(simulation, state_definitions, md_topology):
            """
            Check if the simulation is in State A or State B using ANY-CV logic.
            
            MULTI-CV WITH ANY LOGIC (Nov 7, 2025):
            - If ANY CV detects State A → return "State A"
            - If ANY CV detects State B → return "State B"  
            - Otherwise → return None (transition region)
            
            This casts a wider net and makes both states discoverable.
            Follows the successful approach from the original working script.
            """
            state = simulation.context.getState(getPositions=True)
            positions_nm = state.getPositions(asNumpy=True).value_in_unit(nanometer)
            
            # Create a 1-frame trajectory for mdtraj to analyze
            traj = md.Trajectory([positions_nm], topology=md_topology)
            
            # Check each CV - if ANY detects a state, return immediately
            for definition in state_definitions:
                # Check if this is RMSD-based CV
                if definition.get("type") == "rmsd":
                    selection = definition["selection"]
                    ref_a_path = definition["reference_structures"]["state_a"]
                    ref_b_path = definition["reference_structures"]["state_b"]
                    state_a_thresh = definition["state_a_threshold"]
                    state_b_thresh = definition["state_b_threshold"]
                    
                    # Load reference structures (cached per worker)
                    ref_a = md.load(ref_a_path)
                    ref_b = md.load(ref_b_path)
                    
                    # Select atoms for RMSD calculation
                    atom_indices = md_topology.select(selection)
                    
                    # Compute RMSD to both reference structures
                    rmsd_to_a = md.rmsd(traj, ref_a, atom_indices=atom_indices)[0] * nanometers
                    rmsd_to_b = md.rmsd(traj, ref_b, atom_indices=atom_indices)[0] * nanometers
                    
                    # ANY CV can trigger state classification
                    if rmsd_to_a < state_a_thresh:
                        return "State A"
                    elif rmsd_to_b < state_b_thresh:
                        return "State B"
                        
                else:
                    # Distance-based CV
                    state_a_pair_selections = definition["residues"][0]
                    state_b_pair_selections = definition["residues"][1]
                    
                    state_a_thresh = definition["state_a_threshold"]
                    state_b_thresh = definition["state_b_threshold"]

                    # Check State A with this CV
                    try:
                        atom_a1_idx = md_topology.select(state_a_pair_selections[0])[0]
                        atom_a2_idx = md_topology.select(state_a_pair_selections[1])[0]
                        atom_a_indices = np.array([[atom_a1_idx, atom_a2_idx]])
                        
                        distance_a = md.compute_distances(traj, atom_a_indices)[0, 0] * nanometers
                        
                        # ANY CV detecting State A counts
                        if distance_a < state_a_thresh:
                            return "State A"
                    except IndexError:
                        pass

                    # Check State B with this CV
                    try:
                        atom_b1_idx = md_topology.select(state_b_pair_selections[0])[0]
                        atom_b2_idx = md_topology.select(state_b_pair_selections[1])[0]
                        atom_b_indices = np.array([[atom_b1_idx, atom_b2_idx]])
                        
                        distance_b = md.compute_distances(traj, atom_b_indices)[0, 0] * nanometers
                        
                        # ANY CV detecting State B counts
                        if distance_b > state_b_thresh:
                            return "State B"
                    except IndexError:
                        pass
            
            # No CV detected either state
            return None  # In transition region

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
