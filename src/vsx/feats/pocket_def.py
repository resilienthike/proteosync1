from __future__ import annotations
from pathlib import Path
from typing import List, Tuple
from vsx.utils.paths import POCKETS_DIR


def pocket_file_path(target: str, name: str) -> Path:
    return POCKETS_DIR / f"{target}_{name}_pocket.tsv"


def init_pocket_file(target: str, name: str, overwrite: bool = False) -> Path:
    p = pocket_file_path(target, name)
    if p.exists() and not overwrite:
        return p
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text("# Pocket residues\n# CHAIN\tRESID\n")
    return p


def load_pocket_residues(target: str, name: str) -> List[Tuple[str, int]]:
    p = pocket_file_path(target, name)
    if not p.exists():
        return []
    out: List[Tuple[str, int]] = []
    for line in p.read_text().splitlines():
        if not line or line.startswith("#"):
            continue
        parts = line.split()
        if len(parts) < 2:
            continue
        out.append((parts[0], int(parts[1])))
    return out


def validate_pocket(target: str, name: str) -> str:
    residues = load_pocket_residues(target, name)
    if not residues:
        return "EMPTY (no residues listed)"
    bad = [(c, r) for c, r in residues if r <= 0 or len(c) == 0]
    if bad:
        return f"INVALID entries: {bad[:3]}..."
    return f"OK ({len(residues)} residues)"
