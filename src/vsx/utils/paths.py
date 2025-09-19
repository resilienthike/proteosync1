from __future__ import annotations
from pathlib import Path


def _find_repo_root(start: Path) -> Path:
    cur = start
    for p in [cur, *cur.parents]:
        if (p / "pyproject.toml").exists():
            return p
    return start


REPO_ROOT = _find_repo_root(Path(__file__).resolve())
ARTIFACTS_DIR = REPO_ROOT / "artifacts"
RAW_DIR = ARTIFACTS_DIR / "structures" / "raw"
STD_DIR = ARTIFACTS_DIR / "structures" / "std"
POCKETS_DIR = ARTIFACTS_DIR / "pockets"
DATA_DIR = ARTIFACTS_DIR / "data"

for d in (ARTIFACTS_DIR, RAW_DIR, STD_DIR, POCKETS_DIR, DATA_DIR):
    d.mkdir(parents=True, exist_ok=True)
