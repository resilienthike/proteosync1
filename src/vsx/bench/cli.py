from __future__ import annotations
import json
import typer
from vsx import __version__
from vsx.utils.paths import REPO_ROOT, RAW_DIR, STD_DIR, POCKETS_DIR
from vsx.utils.settings import load_settings
from vsx.utils.logging import setup_logging
from vsx.feats.pocket_def import init_pocket_file, validate_pocket

app = typer.Typer(no_args_is_help=True, add_completion=False)

@app.callback()
def _root(version: bool = typer.Option(False, "--version", help="show version and exit")):
    setup_logging()
    if version:
        typer.echo(f"proteosync {__version__}")
        raise typer.Exit(0)

@app.command()
def status():
    raw = list(RAW_DIR.glob("*"))
    std = list(STD_DIR.glob("*"))
    typer.echo(f"raw_files={len(raw)}")
    typer.echo(f"std_files={len(std)}")

@app.command()
def paths():
    typer.echo(f"REPO_ROOT={REPO_ROOT}")
    typer.echo(f"RAW_DIR={RAW_DIR}")
    typer.echo(f"STD_DIR={STD_DIR}")
    typer.echo(f"POCKETS_DIR={POCKETS_DIR}")

@app.command("show-config")
def show_config():
    cfg = load_settings().model_dump()
    typer.echo(json.dumps(cfg, indent=2, default=str))

@app.command("init-pocket")
def init_pocket(spec: str):
    target, name = spec.split(":", 1)
    p = init_pocket_file(target, name, overwrite=False)
    typer.echo(str(p))

@app.command("validate-pocket")
def validate_pocket_cmd(spec: str):
    target, name = spec.split(":", 1)
    typer.echo(validate_pocket(target, name))

if __name__ == "__main__":
    app()
