from __future__ import annotations
from pathlib import Path
from typing import Optional
from pydantic import BaseModel, Field
import yaml
from .paths import REPO_ROOT, ARTIFACTS_DIR, RAW_DIR, STD_DIR, POCKETS_DIR


class Project(BaseModel):
    name: str = "proteosync"
    version: str = "0.1.0"


class Paths(BaseModel):
    artifacts: Path = ARTIFACTS_DIR
    raw: Path = RAW_DIR
    std: Path = STD_DIR
    pockets: Path = POCKETS_DIR


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
        data = yaml.safe_load(cfg_path.read_text()) or {}
    settings = Settings(**data)
    settings.paths.artifacts.mkdir(parents=True, exist_ok=True)
    settings.paths.raw.mkdir(parents=True, exist_ok=True)
    settings.paths.std.mkdir(parents=True, exist_ok=True)
    settings.paths.pockets.mkdir(parents=True, exist_ok=True)
    return settings
