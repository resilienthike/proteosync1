from __future__ import annotations
from abc import ABC, abstractmethod
from pathlib import Path

class BaseEngine(ABC):
    @abstractmethod
    def prepare(self, **kwargs) -> Path: ...
    @abstractmethod
    def run(self, **kwargs) -> Path: ...
    @abstractmethod
    def convert(self, **kwargs) -> Path: ...

class NullEngine(BaseEngine):
    def __init__(self, workdir: str | Path = "artifacts/md/null"):
        self.workdir = Path(workdir)
        self.workdir.mkdir(parents=True, exist_ok=True)
    def prepare(self, **kwargs) -> Path:
        p = self.workdir / "prepared.txt"
        p.write_text("prepared\n")
        return p
    def run(self, **kwargs) -> Path:
        p = self.workdir / "traj.dummy"
        p.write_text("traj\n")
        return p
    def convert(self, **kwargs) -> Path:
        p = self.workdir / "traj.xtc"
        p.write_text("xtc\n")
        return p
