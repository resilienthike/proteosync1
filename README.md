# Proteosync

Mechanistic Characterization of Cryptic Allosteric Sites.

## Architecture Overview

`proteosync` is a Python-based command-line tool for preparing, running, and analyzing molecular dynamics (MD) simulations of proteins, with a focus on characterizing allosteric sites.

The project is structured around a central library `vsx` found in the `src/` directory, with a CLI entrypoint, and various modules for handling specific tasks like data management, simulation setup, and feature analysis.

### Core Components

- **CLI (`src/vsx/bench/cli.py`)**: The main user interface is a command-line application built with `typer`. It provides commands for managing the project, inspecting configuration, and launching simulation preparation tasks. The entry point is registered in `pyproject.toml` as the `proteosync` command.

- **Configuration (`src/vsx/utils/`)**:
  - `settings.py`: Defines the project's settings using `pydantic` models. This includes paths for artifacts, data, and logging configuration. Settings are loaded from `config/proteosync.yaml`.
  - `paths.py`: Establishes the directory structure for the project, ensuring that artifacts, raw data, and processed structures are stored in a consistent manner.
  - `config.py`: Provides utilities for loading and merging configurations.

- **Data Management (`src/vsx/data/`)**:
  - `targets.py`: Manages protein targets defined in `config/targets.yaml`. It includes functions to load target information and fetch sequence data from external databases like UniProt.
  - `gpcr_structs.py`: Contains logic for downloading and standardizing protein structures.

- **Simulation Preparation (`src/vsx/md/`)**:
  - `prepare_simulation.py`: This is the heart of the simulation setup pipeline. It uses `OpenMM` and `PDBFixer` to take a "seed" protein structure (e.g., a PDB or CIF file), and performs the following steps:
    1.  **Fixing**: Adds missing atoms and residues.
    2.  **Solvation**: Places the protein in a membrane and/or solvates it with water and ions.
    3.  **Minimization**: Runs an energy minimization to relax the system.
    4.  **Output**: Saves the prepared system as PDB files and the simulation state as an XML file, ready for MD simulations.

- **Feature Engineering (`src/vsx/feats/`)**:
  - `pocket_def.py`: Provides functions to define and validate binding pockets based on residue selections. These pocket definitions are stored in `.tsv` files.

### How It Works: A Typical Workflow

1.  **Configuration**: A user defines a new protein target in `config/targets.yaml`, specifying its name and any relevant metadata like a UniProt ID.
2.  **Input Structure**: A seed structure for the target is placed in `artifacts/data/<target_name>/seed_structure.pdb`.
3.  **Simulation Preparation**: The user runs the simulation preparation script (e.g., via the CLI or directly). The `prepare_simulation` module reads the seed structure, cleans it, builds the simulation environment (membrane, water, ions), and saves the prepared system in the `artifacts/md/<target_name>/` directory.
4.  **Analysis**: Once simulations are run, other modules (not yet fully implemented, e.g., `msm`, `ml`) can be used to analyze the trajectories to characterize cryptic pockets and other allosteric features.

The `scripts/` directory contains helper and diagnostic scripts, such as `ff_locator.py` to check the availability of force field files.
