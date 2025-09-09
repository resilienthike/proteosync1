from __future__ import annotations
from pathlib import Path
from typing import Optional
from pydantic import BaseModel, Field
import yaml
from .paths import REPO_ROOT, ARTIFACTS_DIR, RAW_DIR, STD_DIR, POCKETS_DIR, DATA_DIR

class Project(BaseModel):
    name: str = "proteosync"
    version: str = "0.1.0"

class Paths(BaseModel):
    artifacts: Path = ARTIFACTS_DIR
    raw: Path = RAW_DIR
    std: Path = STD_DIR
    pockets: Path = POCKETS_DIR
    data: Path = DATA_DIR

class Run(BaseModel):
    log_level: str = "INFO"

class Settings(BaseModel):
    project: Project = Field(default_factory=Project)
    paths: Paths = Field(default_factory=Paths)
    run: Run = Field(default_factory=Run)

def load_settings(path: Optional[str] = None) -> Settings:
    cfg_path = Path(path) if path else (REPO_ROOT / "config" / "proteosync.yaml")
    data = {}
    if cfg_path.exists():
        import yaml as _yaml  # keep import local if file missing
        data = _yaml.safe_load(cfg_path.read_text()) or {}
    settings = Settings(**data)
    # ensure dirs
    for d in (settings.paths.artifacts, settings.paths.raw, settings.paths.std,
              settings.paths.pockets, settings.paths.data):
        Path(d).mkdir(parents=True, exist_ok=True)
    return settings
