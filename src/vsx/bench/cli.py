from __future__ import annotations
import argparse, json
from pathlib import Path
from vsx.utils.paths import REPO_ROOT, RAW_DIR, STD_DIR
from vsx.utils.config import load_config
from vsx.feats.pocket_def import init_pocket_file, validate_pocket

def main() -> None:
    p = argparse.ArgumentParser(description="Proteosync CLI")
    p.add_argument("--version", action="store_true")
    p.add_argument("--show-paths", dest="show_paths", action="store_true")
    p.add_argument("--show-config", dest="show_config", action="store_true")
    p.add_argument("--status", action="store_true")
    p.add_argument("--list-raw", action="store_true")
    p.add_argument("--list-std", action="store_true")
    p.add_argument("--init-pocket", metavar="TARGET:NAME")
    p.add_argument("--validate-pocket", metavar="TARGET:NAME")
    args = p.parse_args()

    if args.version:
        print("proteosync 0.1.0"); return
    if args.show_paths:
        print(f"REPO_ROOT={REPO_ROOT}"); print(f"RAW_DIR={RAW_DIR}"); print(f"STD_DIR={STD_DIR}"); return
    if args.show_config:
        print(json.dumps(load_config(), indent=2)); return
    if args.status:
        raw = list(Path(RAW_DIR).glob("*")); std = list(Path(STD_DIR).glob("*"))
        print(f"raw_files={len(raw)}"); print(f"std_files={len(std)}"); return
    if args.list_raw:
        for fp in sorted(Path(RAW_DIR).glob("*")): print(fp.name); return
    if args.list_std:
        for fp in sorted(Path(STD_DIR).glob("*")): print(fp.name); return
    if args.init_pocket:
        target, name = args.init_pocket.split(":", 1)
        print(init_pocket_file(target, name, overwrite=False)); return
    if args.validate_pocket:
        target, name = args.validate_pocket.split(":", 1)
        print(validate_pocket(target, name)); return
    print("OK")

if __name__ == "__main__":
    main()
