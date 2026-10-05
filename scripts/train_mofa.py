from mofapy2.run.entry_point import entry_point
import pandas as pd
import numpy as np

import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import os

# Expectations saved to each HDF5. mofapy2's default is W, Z and Sigma only; Tau (per-feature noise
# precision) and AlphaZ (per-group factor activity) are needed for an exact projection of PDBP.
SAVED_EXPECTATIONS = ["W", "Z", "Sigma", "Tau", "AlphaZ"]

def train_mofa_with_seeds(mofa_input_file, output_dir, n_seeds=5, n_factors=10, n_iter=1000):
    """
    Trains MOFA with multiple random initializations to check stability.
    Saves ELBO diagnostic plots for each seed.
    """
    os.makedirs(output_dir, exist_ok=True)
    results = {}
    data = pd.read_csv(mofa_input_file)
    
    for seed in range(n_seeds):
        print(f"--- Training with Seed: {seed} ---")
        
        ent = entry_point()
        
        # Load from long format as specified in mofa_prep
        ent.set_data_df(data) 
        
        ent.set_data_options(scale_views=True)
        
        ent.set_model_options(factors=n_factors, 
                              ard_factors=True, 
                              ard_weights=True)
        ent.set_train_options(
            iter=n_iter, 
            seed=seed,
            convergence_mode="slow"
        )
        
        ent.build()
        
        ent.run()
        
        # Save Model
        model_path = os.path.join(output_dir, f"model_seed_{seed}.hdf5")
        ent.save(model_path, expectations=SAVED_EXPECTATIONS)
        
        # Capture and Plot ELBO
        elbo = ent.model.getTrainingStats()['elbo']
        results[seed] = elbo
        
        # Diagnostic Plotting
        clean_elbo = [e for e in elbo if not np.isnan(e)]  
        plt.figure()
        plt.plot(clean_elbo, marker='o')
        plt.gca().xaxis.set_major_locator(MaxNLocator(integer=True))
        plt.xlabel("Checkpoint", fontsize=11)
        plt.ylabel("ELBO (Evidence Lower Bound)", fontsize=11)
        plt.title(f"ELBO Convergence - Seed {seed}")
        plt.savefig(os.path.join(output_dir, f"elbo_plot_seed_{seed}.png"))
        plt.close()
        
        # mofapy2 stores the ELBO before training plus one entry per iteration, and truncates
        # the array on convergence, so a full-length array means training hit the cap
        n_run = len(elbo) - 1
        status = "stopped at the iteration cap; check convergence" if n_run >= n_iter else "converged"
        print(f"Final ELBO for seed {seed}: {clean_elbo[-1]} after {n_run} iterations ({status})")
        
    return results
