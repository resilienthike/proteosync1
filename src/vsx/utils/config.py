from __future__ import annotations
from pathlib import Path
import tomllib
from .paths import REPO_ROOT

DEFAULT_CFG: dict = {
    "project": {"name": "proteosync", "version": "0.1.0"},
    "paths": {"raw": "artifacts/structures/raw", "std": "artifacts/structures/std"},
}


def _deep_merge(a: dict, b: dict) -> dict:
    out = dict(a)
    for k, v in b.items():
        if isinstance(v, dict) and isinstance(out.get(k), dict):
            out[k] = _deep_merge(out[k], v)
        else:
            out[k] = v
    return out


def load_config(path: str | None = None) -> dict:
    cfg = DEFAULT_CFG
    cfg_path = Path(path) if path else REPO_ROOT / "config" / "proteosync.toml"
    if cfg_path.exists():
        data = tomllib.loads(cfg_path.read_text())
        cfg = _deep_merge(cfg, data)
    return cfg
