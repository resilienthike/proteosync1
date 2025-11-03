from pathlib import Path
import argparse
import sys
import subprocess
import time
import threading
from openmm import Platform, XmlSerializer, LangevinIntegrator
from openmm.app import PDBFile, ForceField, PME, HBonds, DCDReporter, StateDataReporter, Simulation
from openmm.unit import nanometer, picosecond, femtoseconds, kelvin

# Define project paths
REPO_ROOT = Path(__file__).resolve().parents[1]
ARTIFACTS_DIR = REPO_ROOT / "artifacts"
DATA_DIR = REPO_ROOT / "data"


def monitor_gpu_utilization(stop_event, interval=5):
    """Monitor GPU utilization during simulation"""
    while not stop_event.is_set():
        try:
            result = subprocess.run(
                [
                    "nvidia-smi",
                    "--query-gpu=utilization.gpu,memory.used,memory.total",
                    "--format=csv,noheader,nounits",
                ],
                capture_output=True,
                text=True,
            )
            if result.returncode == 0:
                lines = result.stdout.strip().split("\n")
                for i, line in enumerate(lines):
                    gpu_util, mem_used, mem_total = line.split(",")
                    print(
                        f"GPU {i}: {gpu_util.strip()}% utilization, Memory: {mem_used.strip()}/{mem_total.strip()} MB"
                    )
        except Exception as e:
            print(f"Error monitoring GPU: {e}")

        time.sleep(interval)


def run_simulation(target_name: str, steps: int, use_boltz: bool = False):
    # Set paths based on whether using Boltz or traditional workflow
    if use_boltz:
        # Boltz workflow: prepared files are in data/target_name/
        md_dir = DATA_DIR / target_name.lower()
    else:
        # Traditional workflow: files are in artifacts/md/TARGET_NAME/
        md_dir = ARTIFACTS_DIR / "md" / target_name.upper()

    pdb_path = md_dir / "prepared_system.pdb"
    state_path = md_dir / "minimized_state.xml"

    if not pdb_path.exists() or not state_path.exists():
        raise FileNotFoundError(f"Prepared system files not found in {md_dir}")

    print(f"--> Using prepared system from: {md_dir}")
    print("--> Loading the prepared system...")
    pdb = PDBFile(str(pdb_path))
    forcefield = ForceField("amber14-all.xml", "amber14/tip3p.xml", "amber14/lipid17.xml")
    system = forcefield.createSystem(
        pdb.topology,
        nonbondedMethod=PME,
        nonbondedCutoff=1.0 * nanometer,
        constraints=HBonds,
    )
    integrator = LangevinIntegrator(310 * kelvin, 1.0 / picosecond, 2.0 * femtoseconds)

    print("--> Checking available platforms...")
    num_platforms = Platform.getNumPlatforms()
    print(f"Available platforms ({num_platforms}):")
    for i in range(num_platforms):
        platform_name = Platform.getPlatform(i).getName()
        print(f"  {i}: {platform_name}")

    # Try to get CUDA platform
    try:
        print("--> Attempting to use CUDA platform...")
        platform = Platform.getPlatformByName("CUDA")

        # Configure CUDA to use single GPU with optimized settings
        properties = {
            "DeviceIndex": "0",    # Use first GPU
            "Precision": "mixed",  # Mixed precision is often faster than single
            "UseBlockingSync": "false"  # Non-blocking for better performance
        }

        print(f"--> Using platform: {platform.getName()}")
        print(f"    GPU devices: {properties['DeviceIndex']}")
        print(f"    Precision: {properties['Precision']}")
        print(f"    Blocking sync: {properties['UseBlockingSync']}")
        print(f"    PME stream: {properties['DisablePmeStream']}")
    except Exception as e:
        print(f"--> CUDA platform not available: {e}")
        print("--> Falling back to OpenCL platform...")
        try:
            platform = Platform.getPlatformByName("OpenCL")
            properties = {"Precision": "single"}
            print(
                f"--> Using platform: {platform.getName()} with {properties['Precision']} precision"
            )
        except Exception as e2:
            print(f"--> OpenCL platform not available: {e2}")
            print("--> Using CPU platform...")
            platform = Platform.getPlatformByName("CPU")
            properties = {}
            print(f"--> Using platform: {platform.getName()}")

    simulation = Simulation(pdb.topology, system, integrator, platform, properties)

    # Load the state from our robust preparation script
    with open(state_path, "r") as f:
        simulation.context.setState(XmlSerializer.deserialize(f.read()))
    
    # Add Reporters for the production run
    dcd_reporter = DCDReporter(str(md_dir / "trajectory.dcd"), 10000)
    state_reporter = StateDataReporter(
        sys.stdout,
        10000,
        step=True,
        potentialEnergy=True,
        temperature=True,
        progress=True,
        remainingTime=True,
        speed=True,
        totalSteps=steps,
        separator="\t",
    )
    simulation.reporters.append(dcd_reporter)
    simulation.reporters.append(state_reporter)

    print("\n--> Starting production simulation...")
    print("--> Starting GPU monitoring...")

    # Start GPU monitoring in background thread
    stop_event = threading.Event()
    monitor_thread = threading.Thread(
        target=monitor_gpu_utilization, args=(stop_event, 10)
    )
    monitor_thread.daemon = True
    monitor_thread.start()

    start_time = time.time()
    simulation.step(steps)
    end_time = time.time()

    # Stop GPU monitoring
    print("\n--> Stopping GPU monitoring...")
    stop_event.set()
    monitor_thread.join(timeout=2)

    total_time = end_time - start_time
    print("\n--> Simulation complete!")
    print(f"    Total time: {total_time:.2f} seconds")
    print(f"    Performance: {steps/total_time:.2f} steps/second")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description="Run a short MD simulation")
    ap.add_argument("--target", "-t", default="GLP1R", help="Target name")
    ap.add_argument(
        "--ns", type=float, default=1.0, help="Simulation length in nanoseconds"
    )
    ap.add_argument(
        "--boltz",
        action="store_true",
        help="Use Boltz-prepared system from data/ directory"
    )
    args = ap.parse_args()

    n_steps = int(args.ns * 500000)

    run_simulation(target_name=args.target, steps=n_steps, use_boltz=args.boltz)
