from __future__ import annotations
from pathlib import Path
import argparse
import numpy as np
from pdbfixer import PDBFixer
from openmm import unit, XmlSerializer, Platform
from openmm.app import (
    PDBFile,
    PDBxFile,
    Modeller,
    ForceField,
    Simulation,
    PME,
    HBonds,
    CutoffNonPeriodic,
)

try:
    from openmm import LangevinMiddleIntegrator as _Integrator
except Exception:
    from openmm import LangevinIntegrator as _Integrator


def _to_pdb(infile: Path) -> Path:
    suf = infile.suffix.lower()
    if suf == ".pdb":
        return infile
    if suf in {".cif", ".mmcif"}:
        out = infile.with_suffix(".pdb")
        px = PDBxFile(str(infile))
        with open(out, "w") as fh:
            PDBFile.writeFile(px.getTopology(), px.getPositions(), fh, keepIds=True)
        return out
    raise ValueError(f"Unsupported input format: {infile.suffix}")


def _load_ff():
    try:
        return ForceField("charmm36_2024.xml", "charmm36_2024/water.xml")
    except Exception:
        return ForceField("charmm36.xml", "charmm36/water.xml")


def _fixer_from_pdb(pdb_path: Path):
    fx = PDBFixer(filename=str(pdb_path))
    fx.findMissingResidues()
    fx.findNonstandardResidues()
    fx.replaceNonstandardResidues()
    fx.removeHeterogens(keepWater=False)
    fx.findMissingAtoms()
    fx.addMissingAtoms()
    fx.addMissingHydrogens(7.0)
    return fx


def _pre_minimize(top, pos):
    # Quick vacuum relax without PBCs (no box yet) => use CutoffNonPeriodic
    ff = ForceField("amber14-all.xml")
    system = ff.createSystem(
        top,
        nonbondedMethod=CutoffNonPeriodic,
        nonbondedCutoff=1.0 * unit.nanometer,
        constraints=HBonds,
    )
    integ = _Integrator(
        300 * unit.kelvin, 1.0 / unit.picosecond, 0.002 * unit.picoseconds
    )
    sim = Simulation(top, system, integ, Platform.getPlatformByName("CPU"))
    sim.context.setPositions(pos)
    sim.minimizeEnergy(maxIterations=500)
    st = sim.context.getState(getPositions=True)
    return st.getPositions()


def _align_long_axis_to_z(top, pos):
    idx = []
    for i, atom in enumerate(top.atoms()):
        el = getattr(atom.element, "symbol", None)
        if el and el != "H":
            idx.append(i)
    arr = np.array([pos[i].value_in_unit(unit.nanometer) for i in idx])
    ctr = arr.mean(axis=0)
    X = arr - ctr
    cov = X.T @ X / max(len(X) - 1, 1)
    w, V = np.linalg.eigh(cov)
    u = V[:, np.argmax(w)]
    z = np.array([0.0, 0.0, 1.0])
    u = u / np.linalg.norm(u)
    v = np.cross(u, z)
    s = np.linalg.norm(v)
    c = float(np.dot(u, z))
    if s < 1e-8:
        R = np.eye(3) if c > 0 else np.diag([1, -1, -1])
    else:
        vx = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
        R = np.eye(3) + vx + vx @ vx * ((1 - c) / (s**2))
    pos_nm = np.array([p.value_in_unit(unit.nanometer) for p in pos])
    pos_nm = pos_nm - ctr
    pos_nm = pos_nm @ R.T
    pos_nm[:, 2] -= pos_nm[:, 2].mean()
    return [unit.Vec3(*xyz) * unit.nanometer for xyz in pos_nm]


def prepare_system(
    seed_file: Path,
    out_dir: Path,
    lipid: str = "POPC",
    salt_molar: float = 0.15,
    padding_nm: float = 2.0,
    reorient: bool = True,
):
    out_dir.mkdir(parents=True, exist_ok=True)
    pdb_path = _to_pdb(seed_file)

    # 1) fix & pre-minimize (non-periodic)
    fixer = _fixer_from_pdb(pdb_path)
    premin_pos = _pre_minimize(fixer.topology, fixer.positions)
    top = fixer.topology
    pos = premin_pos

    if reorient:
        pos = _align_long_axis_to_z(top, pos)

    fixed_pdb = out_dir / "fixed_seed.pdb"
    with open(fixed_pdb, "w") as fh:
        PDBFile.writeFile(top, pos, fh, keepIds=True)

    # 2) Build membrane + water + ions (this creates a PBC box)
    ff = _load_ff()
    modeller = Modeller(top, pos)
    modeller.addMembrane(
        ff,
        lipidType=lipid,
        minimumPadding=padding_nm * unit.nanometer,
        positiveIon="Na+",
        negativeIon="Cl-",
        ionicStrength=salt_molar * unit.molar,
        neutralize=True,
    )

    prepared_pdb = out_dir / "prepared_system.pdb"
    with open(prepared_pdb, "w") as fh:
        PDBFile.writeFile(modeller.topology, modeller.positions, fh)

    # 3) Minimize the full system (now periodic, so PME is OK)
    system = ff.createSystem(
        modeller.topology,
        nonbondedMethod=PME,
        nonbondedCutoff=1.0 * unit.nanometer,
        constraints=HBonds,
    )
    integ = _Integrator(
        310 * unit.kelvin, 1.0 / unit.picosecond, 0.004 * unit.picoseconds
    )
    sim = Simulation(
        modeller.topology, system, integ, Platform.getPlatformByName("CPU")
    )
    sim.context.setPositions(modeller.positions)
    sim.minimizeEnergy(maxIterations=2000)

    state_xml = out_dir / "minimized_state.xml"
    with open(state_xml, "w") as fh:
        fh.write(
            XmlSerializer.serialize(
                sim.context.getState(getPositions=True, getVelocities=True)
            )
        )

    return prepared_pdb, state_xml


def main():
    from vsx.utils.paths import DATA_DIR, ARTIFACTS_DIR

    ap = argparse.ArgumentParser()
    ap.add_argument("--target", "-t", required=True)
    ap.add_argument("--lipid", default="POPC")
    ap.add_argument("--salt", type=float, default=0.15)
    ap.add_argument("--padding-nm", type=float, default=2.0)
    ap.add_argument("--no-reorient", action="store_true")
    args = ap.parse_args()

    seed = DATA_DIR / args.target / "seed_structure.pdb"
    if not seed.exists():
        seed = DATA_DIR / args.target / "seed_structure.cif"
    if not seed.exists():
        raise FileNotFoundError(
            "seed_structure.(pdb|cif) not found under artifacts/data/<target>/"
        )

    out = ARTIFACTS_DIR / "md" / args.target
    prepared_pdb, state_xml = prepare_system(
        seed,
        out,
        lipid=args.lipid,
        salt_molar=args.salt,
        padding_nm=args.padding_nm,
        reorient=not args.no_reorient,
    )
    print(str(prepared_pdb))
    print(str(state_xml))


if __name__ == "__main__":
    main()
