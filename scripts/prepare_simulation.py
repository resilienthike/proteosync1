from pathlib import Path
import argparse
import sys

from pdbfixer import PDBFixer
from openmm import unit, XmlSerializer, LangevinIntegrator, Platform, OpenMMException
from openmm.app import PDBFile, PDBxFile, Modeller, ForceField, Simulation, HBonds, PME

# --- Configuration ---
REPO_ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = REPO_ROOT / "data"  # All data (inputs and prepared systems)
BOLTZ_RESULTS_DIR = REPO_ROOT / "boltz_results_glp1r"
# ---------------------

def _ensure_pdb(seed_file: Path) -> Path:
    seed_file = seed_file.resolve()
    if seed_file.suffix.lower() == ".pdb":
        return seed_file
    if seed_file.suffix.lower() in {".cif", ".mmcif"}:
        pdb_path = seed_file.with_suffix(".pdb")
        if pdb_path.exists():
            return pdb_path
        print(f"--> Converting {seed_file.name} to PDB format...")
        cif = PDBxFile(str(seed_file))
        with open(pdb_path, "w") as f:
            PDBFile.writeFile(cif.topology, cif.positions, f)
        return pdb_path
    raise ValueError(f"Unsupported seed format: {seed_file.suffix}")

def _load_ff() -> ForceField:
    return ForceField("amber14-all.xml", "amber14/tip3p.xml", "amber14/lipid17.xml")

def prepare_system(target_name: str, use_boltz: bool = False):
    print(f"--- Preparing system for target: {target_name} ---")

    # Determine seed file location based on whether using Boltz or traditional structure
    if use_boltz:
        # Look for Boltz-generated structure
        boltz_pred_dir = BOLTZ_RESULTS_DIR / "predictions" / target_name
        seed_file = boltz_pred_dir / f"{target_name}_model_0.cif"
        if not seed_file.exists():
            print(f"🔥 Error: Boltz structure not found at {seed_file}", file=sys.stderr)
            print(f"   Did you run: boltz predict data/{target_name}.yaml --use_msa_server", file=sys.stderr)
            sys.exit(1)
        print(f"--> Using Boltz-generated structure: {seed_file}")
    else:
        # Traditional structure from data/ folder
        seed_file = DATA_DIR / f"{target_name}_structure.cif"
        if not seed_file.exists():
            # Try PDB format
            seed_file = DATA_DIR / f"{target_name}_structure.pdb"
            if not seed_file.exists():
                print(f"🔥 Error: Structure file not found in {DATA_DIR}", file=sys.stderr)
                print(f"   Looking for: {target_name}_structure.cif or {target_name}_structure.pdb", file=sys.stderr)
                sys.exit(1)
        print(f"--> Using structure from data/: {seed_file}")

    # Output to data/ folder
    out_dir = DATA_DIR / target_name
    out_dir.mkdir(parents=True, exist_ok=True)

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

    ff = _load_ff()
    modeller = Modeller(fixer.topology, fixer.positions)
    print("--> Building membrane, solvent, and ions...")
    modeller.addMembrane(
        ff,
        lipidType="POPC",
        membraneCenterZ=0.0 * unit.nanometer,
        minimumPadding=1.0 * unit.nanometer,
        ionicStrength=0.15 * unit.molar,
        positiveIon="Na+",
        negativeIon="Cl-",
    )
    print("--> System building complete.")

    prepared_pdb_path = out_dir / "prepared_system.pdb"
    with open(prepared_pdb_path, "w") as f:
        PDBFile.writeFile(modeller.topology, modeller.positions, f)
    print(f"✅ Full system PDB saved to: {prepared_pdb_path}")

    print("--> Creating and minimizing the full system...")
    system = ff.createSystem(
        modeller.topology,
        nonbondedMethod=PME,
        nonbondedCutoff=1.0 * unit.nanometer,
        constraints=HBonds,
    )
    integrator = LangevinIntegrator(
        310 * unit.kelvin, 1.0 / unit.picosecond, 2.0 * unit.femtoseconds
    )

    try:
        platform = Platform.getPlatformByName("CUDA")
        properties = {"Precision": "mixed"}
        simulation = Simulation(
            modeller.topology, system, integrator, platform, properties
        )
    except OpenMMException:
        simulation = Simulation(modeller.topology, system, integrator)

    print(f"    Using platform: {simulation.context.getPlatform().getName()}")
    simulation.context.setPositions(modeller.positions)
    simulation.minimizeEnergy()
    print("--> Final minimization complete.")

    # --- NEW: Save the system and state ---
    system_xml_path = out_dir / "system.xml"
    with open(system_xml_path, "w") as f:
        f.write(XmlSerializer.serialize(system))
    print(f"✅ Correctly built System saved to: {system_xml_path}")

    minimized_state_path = out_dir / "minimized_state.xml"
    state = simulation.context.getState(
        getPositions=True, getVelocities=True, getParameters=True
    )
    with open(minimized_state_path, "w") as f:
        f.write(XmlSerializer.serialize(state))
    print(f"✅ Minimized State saved to: {minimized_state_path}")
    print("\n--- Preparation successful! ---")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Prepare a GPCR system for MD simulation."
    )
    parser.add_argument("--target", "-t", required=True, help="Target name (e.g., glp1r).")
    parser.add_argument(
        "--boltz", 
        action="store_true", 
        help="Use Boltz-generated structure from boltz_results_glp1r/ directory"
    )
    args = parser.parse_args()
    prepare_system(args.target, use_boltz=args.boltz)
