#!/usr/bin/env python3
# scripts/prepare_simulation.py
from __future__ import annotations

from pathlib import Path
import argparse
import sys

from pdbfixer import PDBFixer
from openmm import unit, XmlSerializer, LangevinIntegrator, CustomExternalForce, Platform
from openmm.app import (
    PDBFile, PDBxFile, Modeller, ForceField, Simulation,
    HBonds, NoCutoff, PME
)

# Define project paths
REPO_ROOT = Path(__file__).resolve().parents[1]
ARTIFACTS_DIR = REPO_ROOT / "artifacts"
DATA_DIR = ARTIFACTS_DIR / "data"

# ---------- Helper Functions ----------

def _ensure_pdb(seed_file: Path) -> Path:
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
    try:
        return ForceField("amber14-all.xml", "amber14/tip3p.xml", "amber/lipid17.xml")
    except Exception as e:
        raise RuntimeError(f"Could not load a protein+lipid force field. Error: {e}")

def _pre_relax(topology, positions):
    print("--> Performing advanced relaxation (annealing) in implicit solvent...")
    ff = ForceField("amber14-all.xml", "implicit/gbn2.xml")
    system = ff.createSystem(topology, nonbondedMethod=NoCutoff, constraints=HBonds)
    
    restraint_force = CustomExternalForce("k*periodicdistance(x, y, z, x0, y0, z0)^2")
    system.addForce(restraint_force)
    restraint_force.addGlobalParameter("k", 1000.0 * unit.kilojoules_per_mole/unit.nanometer**2)
    restraint_force.addPerParticleParameter("x0")
    restraint_force.addPerParticleParameter("y0")
    restraint_force.addPerParticleParameter("z0")
    for atom in topology.atoms():
        if atom.element.symbol != 'H':
            restraint_force.addParticle(atom.index, positions[atom.index])
    
    integrator = LangevinIntegrator(300*unit.kelvin, 1.0/unit.picosecond, 2.0*unit.femtoseconds)
    
    # --- THIS LINE IS THE FIX ---
    # By not specifying a platform, we let OpenMM automatically choose the fastest one (CUDA).
    sim = Simulation(topology, system, integrator)
    # ----------------------------

    print(f"    Using platform: {sim.context.getPlatform().getName()}")
    sim.context.setPositions(positions)

    print("    Minimizing with restraints...")
    sim.minimizeEnergy(maxIterations=500)

    print("    Annealing...")
    for temp in [300, 310, 320, 310, 300]:
        integrator.setTemperature(temp*unit.kelvin)
        sim.step(5000)

    system.removeForce(system.getNumForces()-1)
    sim.context.reinitialize(preserveState=True)
    print("    Final minimization without restraints...")
    sim.minimizeEnergy(maxIterations=1000)
    
    print("--> Advanced relaxation complete.")
    return sim.context.getState(getPositions=True).getPositions()


def prepare_system(seed_file: Path, out_dir: Path,
                   lipid: str = "POPC",
                   salt_molar: float = 0.15,
                   padding_nm: float = 1.0):
    out_dir.mkdir(parents=True, exist_ok=True)

    pdb_path = _ensure_pdb(seed_file)
    print("--> Fixing PDB structure with PDBFixer...")
    fixer = PDBFixer(filename=str(pdb_path))
    fixer.findMissingResidues()
    # ... (rest of the function is the same)
    fixer.findNonstandardResidues()
    fixer.replaceNonstandardResidues()
    fixer.removeHeterogens(keepWater=False)
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(7.0)
    print("--> PDBFixer complete.")

    relaxed_pos = _pre_relax(fixer.topology, fixer.positions)

    print("--> Building membrane, solvent, and ions...")
    ff = _load_ff()
    modeller = Modeller(fixer.topology, relaxed_pos)
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

    print("--> Performing final energy minimization...")
    system = ff.createSystem(
        modeller.topology, nonbondedMethod=PME,
        nonbondedCutoff=1.0 * unit.nanometer, constraints=HBonds
    )
    integrator = LangevinIntegrator(310*unit.kelvin, 1.0/unit.picosecond, 2.0*unit.femtoseconds)
    sim = Simulation(modeller.topology, system, integrator)
    print(f"    Using platform: {sim.context.getPlatform().getName()}")
    sim.context.setPositions(modeller.positions)
    sim.minimizeEnergy(maxIterations=500)
    print("--> Final minimization complete.")

    state_xml = out_dir / "minimized_state.xml"
    with open(state_xml, "w") as fh:
        fh.write(XmlSerializer.serialize(
            sim.context.getState(getPositions=True, getVelocities=True)
        ))
    
    return prepared_pdb, state_xml

def main():
    parser = argparse.ArgumentParser(description="Prepare membrane MD system from seed structure")
    parser.add_argument("--target", "-t", required=True, help="Target name under artifacts/data/<target>")
    parser.add_argument("--lipid", default="POPC", help="Membrane lipid type (e.g., POPC, POPE)")
    parser.add_argument("--salt", type=float, default=0.15, help="Ionic strength (mol/L)")
    parser.add_argument("--padding-nm", type=float, default=1.0, help="Padding around protein (nm)")
    args = parser.parse_args()

    seed = (DATA_DIR / args.target / "seed_structure.cif")
    if not seed.exists():
        raise FileNotFoundError(f"seed_structure.cif not found under {DATA_DIR / args.target}")

    out = ARTIFACTS_DIR / "md" / args.target
    
    prepared_pdb, state_xml = prepare_system(
        seed, out,
        lipid=args.lipid,
        salt_molar=args.salt,
        padding_nm=args.padding_nm,
    )
    print("\n✅ Preparation successful!")
    print(f"   Final PDB: {prepared_pdb}")
    print(f"   Final State: {state_xml}")

if __name__ == "__main__":
    main()
