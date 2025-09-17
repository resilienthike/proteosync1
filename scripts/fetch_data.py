# scripts/fetch_data.py
import argparse
import subprocess
from pathlib import Path

# --- Configuration ---
# Use the same "Anyone" link you already have
TARGET_DATA_URLS = {
    "GLP1R": "https://netorg18511885-my.sharepoint.com/:u:/g/personal/partnerships_varosync_com/EZQou1R5skdAnJnSbvHYmREB8Xq7rybFSoLRnSe950Kjeg?e=nbQnUK"
}

# Define project paths
REPO_ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = REPO_ROOT / "artifacts" / "data"
# --------------------

def fetch_seed_structure(target_name: str):
    """
    Checks for a seed structure file and downloads it if not present using curl.
    """
    print(f"--- Ensuring seed structure exists for target: {target_name} ---")
    
    if target_name not in TARGET_DATA_URLS:
        raise ValueError(f"Target '{target_name}' not defined in TARGET_DATA_URLS.")

    url = TARGET_DATA_URLS[target_name]
    target_dir = DATA_DIR / target_name
    output_path = target_dir / "seed_structure.cif"

    if output_path.exists():
        print(f"✅ Seed file already exists: {output_path}")
        return

    print(f"--> Seed file not found. Downloading from your link using curl...")
    target_dir.mkdir(parents=True, exist_ok=True)
    
    # Use curl with -L to follow redirects and -o to specify the output file
    command = ["curl", "-L", "-o", str(output_path), url]
    
    try:
        subprocess.run(command, check=True, capture_output=True, text=True)
        print(f"✅ Download complete. File saved to: {output_path}")
    except subprocess.CalledProcessError as e:
        print(f"🔥 Error downloading file with curl. Return code: {e.returncode}")
        print(f"   Stderr: {e.stderr}")
        print("Please check the URL. If this still fails, the link may not be a direct download link.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Download seed structures for the project.")
    parser.add_argument("--target", "-t", required=True, help="Target name to fetch data for.")
    args = parser.parse_args()
    
    fetch_seed_structure(args.target)