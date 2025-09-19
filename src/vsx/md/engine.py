"""MD Engine implementations."""

from pathlib import Path
from typing import Union


class NullEngine:
    """A null/mock MD engine that creates dummy files for testing."""

    def __init__(self, workdir: Union[str, Path]):
        """Initialize the NullEngine with a working directory.

        Args:
            workdir: Directory where files will be created
        """
        self.workdir = Path(workdir)
        self.workdir.mkdir(parents=True, exist_ok=True)

    def prepare(self) -> Path:
        """Prepare simulation - creates a dummy preparation file.

        Returns:
            Path to the preparation file
        """
        prep_file = self.workdir / "prepared.txt"
        prep_file.write_text("Simulation prepared")
        return prep_file

    def run(self) -> Path:
        """Run simulation - creates a dummy output file.

        Returns:
            Path to the simulation output file
        """
        run_file = self.workdir / "trajectory.dcd"
        run_file.write_text("Simulation trajectory data")
        return run_file

    def convert(self) -> Path:
        """Convert simulation output - creates a dummy converted file.

        Returns:
            Path to the converted file
        """
        convert_file = self.workdir / "converted.pdb"
        convert_file.write_text("Converted simulation data")
        return convert_file
