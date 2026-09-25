from mofapy2.run.entry_point import entry_point
import pandas as pd
import numpy as np

import matplotlib.pyplot as plt
import os

def train_mofa_with_seeds(mofa_input_file, output_dir, n_seeds=5, n_factors=10):
    """
    Trains MOFA with multiple random initializations to check stability.
    Saves ELBO diagnostic plots for each seed.
    """
    os.makedirs(output_dir, exist_ok=True)
    results = {}
    
    for seed in range(n_seeds):
        print(f"--- Training with Seed: {seed} ---")
        
        ent = entry_point()
        # Load from long format as specified in mofa_prep
        ent.load_data(mofa_input_file) 
        
        ent.set_data_options(
            view_names=["RNA_seq", "CSF_Proteomics", "Plasma_Proteomics"],
            scale_views=True
        )
        ent.set_model_options(factors=n_factors)
        ent.set_train_options(
            iter=20, 
            seed=seed,
            convergence_mode="slow"
        )
        
        ent.run()
        
        # Save Model
        model_path = os.path.join(output_dir, f"model_seed_{seed}.hdf5")
        ent.save(model_path)
        
        # Capture and Plot ELBO
        elbo = ent.model.nodes['ELBO'].value
        results[seed] = elbo
        
        # Diagnostic Plotting
        plt.figure()
        plt.plot(elbo)
        plt.title(f"ELBO Convergence - Seed {seed}")
        plt.savefig(os.path.join(output_dir, f"elbo_plot_seed_{seed}.png"))
        plt.close()
        
        print(f"Final ELBO for seed {seed}: {elbo[-1]}")
        
    return results
