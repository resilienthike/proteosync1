#!/usr/bin/env python3
# scripts/prepare_simulation.py
from __future__ import annotations

from pathlib import Path
import argparse

from pdbfixer import PDBFixer
from openmm import unit, XmlSerializer, LangevinIntegrator
from openmm.app import PDBFile, PDBxFile, Modeller, ForceField, Simulation, HBonds, NoCutoff, PME

# Define project paths
REPO_ROOT = Path(__file__).resolve().parents[1]
ARTIFACTS_DIR = REPO_ROOT / "artifacts"
DATA_DIR = ARTIFACTS_DIR / "data"

# ---------- Helper Functions ----------


def _ensure_pdb(seed_file: Path) -> Path:
    """Return a PDB path; convert mmCIF->PDB if needed."""
    seed_file = seed_file.resolve()
    if seed_file.suffix.lower() == ".pdb":
        return seed_file
    if seed_file.suffix.lower() in {".cif", ".mmcif"}:
        out = seed_file.with_suffix(".pdb")
        print(f"--> Converting {seed_file.name} to PDB format...")
        px = PDBxFile(str(seed_file))
        with open(out, "w") as fh:
            PDBFile.writeFile(px.getTopology(), px.getPositions(), fh, keepIds=True)
        return out
    raise ValueError(f"Unsupported seed format: {seed_file.suffix}")


def _load_ff() -> ForceField:
    """Load a standard protein+lipid+water force field set."""
    try:
        # --- THIS LINE IS THE FIX ---
        return ForceField("amber14-all.xml", "amber14/tip3p.xml", "amber/lipid17.xml")
    except Exception as e:
        raise RuntimeError(f"Could not load a protein+lipid force field. Error: {e}")


def _pre_minimize(topology, positions):
    """Quick vacuum minimization to remove bad clashes."""
    print("--> Performing vacuum energy minimization to relax clashes...")
    ff = ForceField("amber14-all.xml", "implicit/gbn2.xml")
    system = ff.createSystem(topology, nonbondedMethod=NoCutoff, constraints=HBonds)
    integrator = LangevinIntegrator(
        300 * unit.kelvin, 1.0 / unit.picosecond, 2.0 * unit.femtoseconds
    )
    sim = Simulation(topology, system, integrator)
    sim.context.setPositions(positions)
    sim.minimizeEnergy(maxIterations=200)
    print("--> Vacuum minimization complete.")
    return sim.context.getState(getPositions=True).getPositions()


# ---------- Main Preparation Function ----------


def prepare_system(
    seed_file: Path,
    out_dir: Path,
    lipid: str = "POPC",
    salt_molar: float = 0.15,
    padding_nm: float = 1.5,
):
    """
    Prepare a GPCR system in a membrane + water + ions and minimize it.
    """
    out_dir.mkdir(parents=True, exist_ok=True)

    # 1. Ensure we have a PDB and fix it
    pdb_path = _ensure_pdb(seed_file)
    print("--> Fixing PDB structure with PDBFixer...")
    fixer = PDBFixer(filename=str(pdb_path))
    fixer.findMissingResidues()
    fixer.findNonstandardResidues()
    fixer.replaceNonstandardResidues()
    fixer.removeHeterogens(keepWater=False)
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(7.0)
    print("--> PDBFixer complete.")

    # 2. Perform pre-minimization in vacuum
    premin_pos = _pre_minimize(fixer.topology, fixer.positions)

    # 3. Build the membrane, solvent, and ions
    print("--> Building membrane, solvent, and ions...")
    ff = _load_ff()
    modeller = Modeller(fixer.topology, premin_pos)
    modeller.addMembrane(
        ff,
        lipidType=lipid,
        minimumPadding=padding_nm * unit.nanometer,
        ionicStrength=salt_molar * unit.molar,
        positiveIon="Na+",
        negativeIon="Cl-",
    )
    print("--> System building complete.")

    prepared_pdb = out_dir / "prepared_system.pdb"
    with open(prepared_pdb, "w") as fh:
        PDBFile.writeFile(modeller.topology, modeller.positions, fh)
    print(f"--> Full system saved to {prepared_pdb}")

    # 4. Build system with PBC and perform final minimization
    print("--> Performing final energy minimization...")
    system = ff.createSystem(
        modeller.topology,
        nonbondedMethod=PME,
        nonbondedCutoff=1.0 * unit.nanometer,
        constraints=HBonds,
    )
    integrator = LangevinIntegrator(
        310 * unit.kelvin, 1.0 / unit.picosecond, 2.0 * unit.femtoseconds
    )
    sim = Simulation(modeller.topology, system, integrator)
    sim.context.setPositions(modeller.positions)
    sim.minimizeEnergy(maxIterations=500)
    print("--> Final minimization complete.")

    # 5. Save the final state
    state_xml = out_dir / "minimized_state.xml"
    with open(state_xml, "w") as fh:
        fh.write(
            XmlSerializer.serialize(sim.context.getState(getPositions=True, getVelocities=True))
        )

    return prepared_pdb, state_xml


# ---------- Command-Line Interface ----------


def main():
    ap = argparse.ArgumentParser(description="Prepare membrane MD system from seed structure")
    ap.add_argument(
        "--target", "-t", required=True, help="Target name under artifacts/data/<target>"
    )
    ap.add_argument("--lipid", default="POPC", help="Membrane lipid type (e.g., POPC, POPE)")
    ap.add_argument("--salt", type=float, default=0.15, help="Ionic strength (mol/L)")
    ap.add_argument("--padding-nm", type=float, default=1.5, help="Padding around protein (nm)")
    args = ap.parse_args()

    seed = DATA_DIR / args.target / "seed_structure.cif"
    if not seed.exists():
        raise FileNotFoundError(f"seed_structure.cif not found under {DATA_DIR / args.target}")

    out = ARTIFACTS_DIR / "md" / args.target

    prepared_pdb, state_xml = prepare_system(
        seed,
        out,
        lipid=args.lipid,
        salt_molar=args.salt,
        padding_nm=args.padding_nm,
    )
    print("\n✅ Preparation successful!")
    print(f"   Final PDB: {prepared_pdb}")
    print(f"   Final State: {state_xml}")


if __name__ == "__main__":
    main()
