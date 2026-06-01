#!/usr/bin/env python3
import argparse
import subprocess
import shutil
import sys
from pathlib import Path


def clean_directory(dir_path: Path):
    """Creates the directory if it doesn't exist, or empties it if it does."""
    if not dir_path.exists():
        dir_path.mkdir(parents=True)
        return

    for item in dir_path.iterdir():
        if item.is_dir() and not item.is_symlink():
            shutil.rmtree(item)
        else:
            item.unlink()


def main():
    parser = argparse.ArgumentParser(
        description="Simulate a clonal matrix, usage matrix and read counts.")

    # Updated to include single-letter flags (and -r for seed)
    parser.add_argument('-m', '--mutations', required=True,
                        type=str, help='Number of mutations.')
    parser.add_argument('-s', '--samples', required=True,
                        type=str, help='Number of sequenced samples.')
    parser.add_argument('-c', '--coverage', required=True,
                        type=str, help='Expected sequencing coverage.')
    parser.add_argument('-l', '--loss', type=str, default='l2',
                        help='Loss function L_i(.) to use for optimization [nargs=0..1] [default: "l2"]')
    parser.add_argument('-r', '--seed', type=str,
                        help='Random seed.', default=None)

    args = parser.parse_args()

    # Define working paths
    sim_dir = Path('./sim')
    out_dir = Path('./out')

    # Prefix for simulate.py outputs
    sim_prefix = sim_dir / 'sim'

    # --- Step 1: Clean the sim directory ---
    clean_directory(sim_dir)

    # --- Step 2: Run simulate.py ---
    print("[*] Running simulate...")
    # We still pass the full --flag names to simulate.py as required by that specific script
    simulate_cmd = [
        "python", "scripts/simulate.py",
        "--mutations", args.mutations,
        "--samples", args.samples,
        "--coverage", args.coverage,
        "--output", str(sim_prefix)
    ]

    # Only append --seed if the user passed it in
    if args.seed is not None:
        simulate_cmd.extend(["--seed", args.seed])

    try:
        subprocess.run(simulate_cmd, check=True)
    except subprocess.CalledProcessError as e:
        print(f"[!] Error running simulate.py: {e}", file=sys.stderr)
        sys.exit(1)

    # --- Step 3: Clean/Create the output directory ---
    clean_directory(out_dir)

    # --- Step 4: Run fastppm ---
    print(f"[*] Running fastppm...")
    fastppm_cmd = [
        "./fastppm-cli",
        "-v", f"{sim_prefix}_variant_matrix.txt",
        "-d", f"{sim_prefix}_total_matrix.txt",
        "-t", f"{sim_prefix}_tree.txt",
        "-o", str(out_dir / "fastppm_results.json"),
        "-f", "verbose",
        "-l", args.loss
    ]

    try:
        subprocess.run(fastppm_cmd, check=True)
    except subprocess.CalledProcessError as e:
        print(f"[!] Error running fastppm: {e}", file=sys.stderr)
        sys.exit(1)

    # --- Step 5: Copy simulation files to the output directory ---

    files_to_copy = [
        "sim_variant_matrix.txt",
        "sim_total_matrix.txt",
        "sim_tree.txt",
        "sim_clonal_matrix.txt",
        "sim_usage_matrix.txt",
        "sim_frequency_matrix.txt"
    ]

    for file_name in files_to_copy:
        src_file = sim_dir / file_name
        if src_file.exists():
            shutil.copy2(src_file, out_dir / file_name)
        else:
            print(f"[!] Expected simulation file {
                  src_file} not found.", file=sys.stderr)
            sys.exit(1)

    # --- Step 6: Visualization ---
    print("[*] Running visualization...")
    try:
        subprocess.run(["python", "visual.py"], check=True)
    except subprocess.CalledProcessError as e:
        print(f"[!] Error running visual.py: {e}", file=sys.stderr)
        sys.exit(1)

    print("[*] Completed!")


if __name__ == "__main__":
    main()
