from __future__ import annotations
import logging
from logging.handlers import RotatingFileHandler
from rich.logging import RichHandler
from .paths import ARTIFACTS_DIR
from .settings import load_settings


def setup_logging(level: str | None = None) -> logging.Logger:
    cfg = load_settings()
    lvl_name = (level or cfg.run.log_level).upper()
    lvl = getattr(logging, lvl_name, logging.INFO)

    root = logging.getLogger()
    root.handlers.clear()
    root.setLevel(lvl)

    log_dir = ARTIFACTS_DIR / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)
    file_handler = RotatingFileHandler(
        log_dir / "proteosync.log", maxBytes=1_000_000, backupCount=5
    )
    file_handler.setFormatter(
        logging.Formatter("%(asctime)s %(levelname)s %(name)s: %(message)s")
    )
    file_handler.setLevel(lvl)

    console = RichHandler(rich_tracebacks=True, markup=True)
    console.setLevel(lvl)

    root.addHandler(file_handler)
    root.addHandler(console)
    return logging.getLogger("proteosync")
