# scripts/run_simulation.py
from pathlib import Path
import argparse
import sys
from openmm import *
from openmm.app import *
from openmm.unit import *

# Define project paths
REPO_ROOT = Path(__file__).resolve().parents[1]
ARTIFACTS_DIR = REPO_ROOT / "artifacts"
MD_DIR = ARTIFACTS_DIR / "md" / "GLP1R"

def run_simulation(target_name: str, steps: int):
    pdb_path = MD_DIR / "prepared_system.pdb"
    state_path = MD_DIR / "minimized_state.xml"

    if not pdb_path.exists() or not state_path.exists():
        raise FileNotFoundError(f"Prepared system files not found in {MD_DIR}")

    print("--> Loading the prepared system...")
    pdb = PDBFile(str(pdb_path))
    forcefield = ForceField("amber14-all.xml", "amber14/tip3p.xml", "amber/lipid17.xml")
    system = forcefield.createSystem(pdb.topology, nonbondedMethod=PME, nonbondedCutoff=1.0*nanometer, constraints=HBonds)
    integrator = LangevinIntegrator(310*kelvin, 1.0/picosecond, 2.0*femtoseconds)
    
    print("--> Finding CUDA platform...")
    platform = Platform.getPlatformByName('CUDA')
    properties = {'Precision': 'single'}
    print(f"--> Using platform: {platform.getName()} with {properties['Precision']} precision")
    
    simulation = Simulation(pdb.topology, system, integrator, platform, properties)
    
    # Load the state from our robust preparation script
    with open(state_path, 'r') as f:
        simulation.context.setState(XmlSerializer.deserialize(f.read()))

    # We have removed the complex and problematic equilibration block.
    # We now go directly to the production simulation.

    # Add Reporters for the production run
    dcd_reporter = DCDReporter(str(MD_DIR / 'trajectory.dcd'), 10000)
    state_reporter = StateDataReporter(
        sys.stdout, 10000, step=True, potentialEnergy=True, temperature=True, progress=True,
        remainingTime=True, speed=True, totalSteps=steps, separator='\t'
    )
    simulation.reporters.append(dcd_reporter)
    simulation.reporters.append(state_reporter)

    print("\n--> Starting production simulation...")
    simulation.step(steps)
    print("--> Simulation complete!")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description="Run a short MD simulation")
    ap.add_argument("--target", "-t", default="GLP1R", help="Target name")
    ap.add_argument("--ns", type=float, default=1.0, help="Simulation length in nanoseconds")
    args = ap.parse_args()

    n_steps = int(args.ns * 500000)
    
    run_simulation(target_name=args.target, steps=n_steps)
