from vsx.utils.settings import load_settings
from typer.testing import CliRunner
from vsx.bench.cli import app


def test_settings_load():
    s = load_settings()
    assert s.project.name == "proteosync"
    assert s.paths.raw.exists()


def test_cli_paths():
    r = CliRunner().invoke(app, ["paths"])
    assert r.exit_code == 0
    assert "RAW_DIR=" in r.stdout
