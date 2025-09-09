from __future__ import annotations
from pathlib import Path

RAW_DIR = Path("artifacts/structures/raw")
STD_DIR = Path("artifacts/structures/std")


def download_or_load_structures(subtypes: list[str]) -> dict[str, Path]:
    raise NotImplementedError("Add raw structures later, then implement loader.")


def standardize_copy(subtype: str) -> Path:
    raise NotImplementedError("Implement after structures are added.")
