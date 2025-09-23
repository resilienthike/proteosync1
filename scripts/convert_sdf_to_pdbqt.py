#!/usr/bin/env python3
"""
Convert SDF files to PDBQT format for AutoDock Vina.

This script uses OpenBabel to convert ligand structures from SDF format
to PDBQT format, which is required for AutoDock Vina docking.
"""

import argparse
import os
import sys
from pathlib import Path
from subprocess import run, PIPE
import logging

def setup_logging(verbose=False):
    """Set up logging configuration."""
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format='%(asctime)s - %(levelname)s - %(message)s',
        datefmt='%H:%M:%S'
    )

def convert_sdf_to_pdbqt(sdf_file, output_dir=None, verbose=False):
    """
    Convert SDF file to PDBQT format using OpenBabel.
    
    Args:
        sdf_file (str): Path to input SDF file
        output_dir (str): Output directory (default: same as input)
        verbose (bool): Enable verbose output
        
    Returns:
        str: Path to output PDBQT file
    """
    sdf_path = Path(sdf_file)
    
    if not sdf_path.exists():
        raise FileNotFoundError(f"SDF file not found: {sdf_file}")
    
    # Determine output directory and filename
    if output_dir:
        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)
        pdbqt_file = output_path / f"{sdf_path.stem}.pdbqt"
    else:
        pdbqt_file = sdf_path.with_suffix('.pdbqt')
    
    logging.info(f"Converting {sdf_file} to {pdbqt_file}")
    
    # Run OpenBabel conversion
    cmd = [
        'obabel',
        str(sdf_file),
        '-O', str(pdbqt_file),
        '--gen3d',  # Generate 3D coordinates if needed
        '-h'        # Add hydrogens
    ]
    
    if verbose:
        logging.debug(f"Running command: {' '.join(cmd)}")
    
    try:
        result = run(cmd, capture_output=True, text=True, check=True)
        
        if verbose and result.stderr:
            logging.debug(f"OpenBabel output: {result.stderr}")
            
        logging.info(f"Successfully converted to: {pdbqt_file}")
        return str(pdbqt_file)
        
    except Exception as e:
        logging.error(f"Conversion failed: {e}")
        raise

def convert_multiple_sdf(input_pattern, output_dir=None, verbose=False):
    """
    Convert multiple SDF files matching a pattern.
    
    Args:
        input_pattern (str): Glob pattern for SDF files
        output_dir (str): Output directory
        verbose (bool): Enable verbose output
        
    Returns:
        list: List of converted PDBQT files
    """
    from glob import glob
    
    sdf_files = glob(input_pattern)
    
    if not sdf_files:
        logging.warning(f"No SDF files found matching pattern: {input_pattern}")
        return []
    
    logging.info(f"Found {len(sdf_files)} SDF files to convert")
    
    converted_files = []
    failed_files = []
    
    for sdf_file in sdf_files:
        try:
            pdbqt_file = convert_sdf_to_pdbqt(sdf_file, output_dir, verbose)
            converted_files.append(pdbqt_file)
        except Exception as e:
            logging.error(f"Failed to convert {sdf_file}: {e}")
            failed_files.append(sdf_file)
    
    logging.info(f"Conversion complete: {len(converted_files)} successful, {len(failed_files)} failed")
    
    if failed_files:
        logging.warning(f"Failed files: {failed_files}")
    
    return converted_files

def main():
    parser = argparse.ArgumentParser(
        description="Convert SDF files to PDBQT format for AutoDock Vina",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Convert single SDF file
  python convert_sdf_to_pdbqt.py ligands.sdf
  
  # Convert multiple SDF files
  python convert_sdf_to_pdbqt.py "*.sdf" -o pdbqt_files/
  
  # Convert with verbose output
  python convert_sdf_to_pdbqt.py prepared_ligands.sdf --verbose
        """
    )
    
    parser.add_argument(
        'input',
        help='Input SDF file or glob pattern (e.g., "*.sdf")'
    )
    
    parser.add_argument(
        '-o', '--output-dir',
        help='Output directory for PDBQT files (default: same as input)'
    )
    
    parser.add_argument(
        '-v', '--verbose',
        action='store_true',
        help='Enable verbose output'
    )
    
    args = parser.parse_args()
    
    setup_logging(args.verbose)
    
    try:
        # Check if OpenBabel is available
        result = run(['obabel', '--version'], capture_output=True, text=True)
        logging.info(f"Using OpenBabel: {result.stdout.strip()}")
        
    except FileNotFoundError:
        logging.error("OpenBabel not found. Please install it:")
        logging.error("  conda install -c conda-forge openbabel")
        sys.exit(1)
    
    try:
        # Determine if input is a single file or pattern
        if '*' in args.input or '?' in args.input:
            # Multiple files with glob pattern
            converted_files = convert_multiple_sdf(args.input, args.output_dir, args.verbose)
            
            if converted_files:
                print(f"\nSuccessfully converted {len(converted_files)} files:")
                for f in converted_files:
                    print(f"  {f}")
            else:
                print("No files were converted.")
                sys.exit(1)
                
        else:
            # Single file
            if not Path(args.input).exists():
                logging.error(f"File not found: {args.input}")
                sys.exit(1)
                
            pdbqt_file = convert_sdf_to_pdbqt(args.input, args.output_dir, args.verbose)
            print(f"\nConversion successful: {pdbqt_file}")
            
    except Exception as e:
        logging.error(f"Conversion failed: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()