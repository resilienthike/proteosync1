from __future__ import annotations
from importlib import resources
from pathlib import Path
import json, re, sys

report = {"env": {}, "paths": {}, "resolved": {}, "lipids": {}, "load_tests": {}}

# --- env info ---
try:
    import openmm as mm
    import openmm.app as app
    report["env"]["openmm_version"] = mm.__version__
    platforms = [mm.Platform.getPlatform(i).getName() for i in range(mm.Platform.getNumPlatforms())]
    report["env"]["platforms"] = platforms
except Exception as e:
    print("OpenMM not importable:", e)
    sys.exit(1)

# pdbfixer version (best-effort)
try:
    import importlib.metadata as im
    report["env"]["pdbfixer_version"] = im.version("pdbfixer")
except Exception:
    report["env"]["pdbfixer_version"] = "unknown"

# --- package data roots ---
app_data = resources.files(app).joinpath("data")
report["paths"]["openmm_app_data"] = str(app_data)

has_omff = False
try:
    import openmmforcefields as omff
    omff_ffxml = resources.files(omff).joinpath("ffxml")
    report["paths"]["openmmforcefields_ffxml"] = str(omff_ffxml)
    has_omff = True
except Exception as e:
    report["paths"]["openmmforcefields_ffxml"] = None
    report["env"]["openmmforcefields_import_error"] = str(e)

# helper to resolve file inside a package resource dir to a real filesystem path
def resolve_resource(base, relpath: str) -> str:
    node = base.joinpath(relpath)
    with resources.as_file(node) as p:
        return str(p)

# --- resolve preferred FF files ---
resolved_ok = {}

if has_omff:
    try:
        prot = resolve_resource(omff_ffxml, "amber/protein.ff14SB.xml")
        water = resolve_resource(omff_ffxml, "amber/tip3p_standard.xml")
        lipid = resolve_resource(omff_ffxml, "amber/lipid17.xml")
        resolved_ok["protein_ff14SB_xml"] = prot
        resolved_ok["tip3p_standard_xml"]  = water
        resolved_ok["lipid17_xml"]         = lipid
    except Exception as e:
        report["resolved"]["error"] = f"resolve omff: {e}"

# built-in fallbacks (may not include lipids)
try:
    prot_bi  = resolve_resource(app_data, "amber14/protein.ff14SB.xml")
    water_bi = resolve_resource(app_data, "amber14/tip3p.xml")
    resolved_ok.setdefault("protein_ff14SB_builtin_xml", prot_bi)
    resolved_ok.setdefault("tip3p_builtin_xml", water_bi)
except Exception as e:
    report["resolved"]["builtin_warn"] = f"resolve builtin: {e}"

report["resolved"].update(resolved_ok)

# --- list lipid residues from lipid17.xml (if present) ---
lipid_names = []
lipid_file = resolved_ok.get("lipid17_xml")
if lipid_file and Path(lipid_file).exists():
    txt = Path(lipid_file).read_text(errors="ignore")
    names = re.findall(r'<Residue\\s+name="([^"]+)"', txt)
    # keep uppercase-ish short residue names (typical lipids like POPC, POPE, etc.)
    lipid_names = sorted({n for n in names if 3 <= len(n) <= 6 and n.upper() == n})
report["lipids"]["from_lipid17"] = lipid_names

# --- try to load with ForceField ---
from openmm.app import ForceField
def try_load(label, files):
    try:
        ff = ForceField(*files)
        report["load_tests"][label] = {"ok": True, "files": files}
        return True
    except Exception as e:
        report["load_tests"][label] = {"ok": False, "error": str(e), "files": files}
        return False

if {"protein_ff14SB_xml","tip3p_standard_xml","lipid17_xml"} <= resolved_ok.keys():
    try_load("omff_ff14SB_tip3p_lipid17", [
        resolved_ok["protein_ff14SB_xml"],
        resolved_ok["tip3p_standard_xml"],
        resolved_ok["lipid17_xml"],
    ])

# a built-in+omff mix fallback
if all(k in resolved_ok for k in ("protein_ff14SB_builtin_xml","tip3p_builtin_xml")) and "lipid17_xml" in resolved_ok:
    try_load("builtin_ff14SB_tip3p__omff_lipid17", [
        resolved_ok["protein_ff14SB_builtin_xml"],
        resolved_ok["tip3p_builtin_xml"],
        resolved_ok["lipid17_xml"],
    ])

# --- print nice summary ---
print("OpenMM:", report["env"].get("openmm_version"))
print("Platforms:", ", ".join(report["env"].get("platforms", [])))
print("pdbfixer:", report["env"].get("pdbfixer_version"))
print("openmm app data:", report["paths"]["openmm_app_data"])
print("openmmforcefields ffxml:", report["paths"]["openmmforcefields_ffxml"])

print("\nResolved files:")
for k, v in report["resolved"].items():
    if isinstance(v, str):
        print(f"  {k}: {v}")

print("\nLipid residues in lipid17.xml (sample):", ", ".join(lipid_names[:30]), ("... +%d more" % (max(0, len(lipid_names)-30)) if len(lipid_names)>30 else ""))

print("\nLoad tests:")
for k, v in report["load_tests"].items():
    status = "OK" if v.get("ok") else "FAIL"
    print(" ", status, k)
    for f in v.get("files", []):
        print("    ", f)
    if not v.get("ok"):
        print("    error:", v.get("error"))

# --- write JSON report ---
out_dir = Path("artifacts/diagnostics")
out_dir.mkdir(parents=True, exist_ok=True)
out_path = out_dir / "ff_report.json"
out_path.write_text(json.dumps(report, indent=2))
print("\nWrote JSON:", out_path)
