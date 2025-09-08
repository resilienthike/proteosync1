from __future__ import annotations
import argparse, json
from vsx.utils.paths import REPO_ROOT, RAW_DIR, STD_DIR
from vsx.utils.config import load_config

def main() -> None:
    p = argparse.ArgumentParser(description="Proteosync CLI")
    p.add_argument("--version", action="store_true")
    p.add_argument("--show-paths", dest="show_paths", action="store_true")
    p.add_argument("--show-config", dest="show_config", action="store_true")
    args = p.parse_args()
    if args.version:
        print("proteosync 0.1.0"); return
    if args.show_paths:
        print(f"REPO_ROOT={REPO_ROOT}")
        print(f"RAW_DIR={RAW_DIR}")
        print(f"STD_DIR={STD_DIR}")
        return
    if args.show_config:
        cfg = load_config()
        print(json.dumps(cfg, indent=2))
        return
    print("OK")

if __name__ == "__main__":
    main()
