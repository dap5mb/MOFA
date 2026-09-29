#!/usr/bin/env python3
import os
import argparse
import sys
import pandas as pd
import h5py

import numpy as np

def get_final_elbo(model_path: str) -> float:
    """Extracts final non-NaN ELBO value from the saved MOFA HDF5 file."""
    with h5py.File(model_path, 'r') as f:
        if 'training_stats' in f and 'elbo' in f['training_stats']:
            elbo = f['training_stats']['elbo'][:]
            # Filter out NaNs before getting the final value
            valid_elbo = elbo[~np.isnan(elbo)]
            if len(valid_elbo) > 0:
                return float(valid_elbo[-1])
    return -float('inf')

def main():
    parser = argparse.ArgumentParser(description="Identify best MOFA model seed based on final ELBO.")
    parser.add_argument("--model_dir", required=True, help="Directory containing model_seed_*.hdf5 files.")
    args = parser.parse_args()
    
    if not os.path.exists(args.model_dir):
        print(f"Error: Directory '{args.model_dir}' does not exist.", file=sys.stderr)
        sys.exit(1)
        
    all_files = os.listdir(args.model_dir)
    print(f"Found {len(all_files)} total file(s) in '{args.model_dir}'")
    
    matching_files = [f for f in all_files if f.startswith("model_seed_") and f.endswith(".hdf5")]
    print(f"Found {len(matching_files)} matching 'model_seed_*.hdf5' file(s)")

    best_elbo = -float('inf')
    best_model = None
    
    for filename in matching_files:
        path = os.path.join(args.model_dir, filename)
        try:
            elbo = get_final_elbo(path)
            print(f"{filename}: Final ELBO = {elbo}")
            if elbo > best_elbo:
                best_elbo = elbo
                best_model = path
        except Exception as e:
            print(f"Could not read {filename}: {e}")
                
    if best_model:
        print(f"\nBest model identified: {best_model} (ELBO: {best_elbo})")
        with open(os.path.join(args.model_dir, "best_model_path.txt"), "w") as f:
            f.write(best_model)
    else:
        print("No valid models found with non-inf ELBO.")

if __name__ == "__main__":
    main()