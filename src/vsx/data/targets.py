from __future__ import annotations
from pathlib import Path
from typing import Dict, Any, Optional
import urllib.request
import yaml
from vsx.utils.paths import REPO_ROOT

REG_PATH = REPO_ROOT / "config" / "targets.yaml"


def load_targets(path: Optional[str] = None) -> Dict[str, Any]:
    p = Path(path) if path else REG_PATH
    if not p.exists():
        return {"targets": {}}
    data = yaml.safe_load(p.read_text()) or {}
    if "targets" not in data:
        data["targets"] = {}
    return data


def get_target(name: str, path: Optional[str] = None) -> Dict[str, Any]:
    return load_targets(path).get("targets", {}).get(name, {})


def fetch_uniprot_fasta(uniprot_id: str) -> str:
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.fasta"
    with urllib.request.urlopen(url) as r:
        text = r.read().decode("utf-8", errors="ignore")
    lines = [ln.strip() for ln in text.splitlines() if ln and not ln.startswith(">")]
    return "".join(lines)
