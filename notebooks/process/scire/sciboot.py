from time import time
from scipy.stats import bootstrap 
import pandas as pd
import numpy as np
from sklearn.mixture import GaussianMixture
from concurrent.futures import ThreadPoolExecutor

# Define the main function
def main():
    # Define the g function, compute_bootstrap function, and stats list here

    # Define the g function
    def g(roi_data):
        gmm_model = GaussianMixture(n_components=2,
                                    init_params='random_from_data',
                                    covariance_type='full',
                                    n_init=100,
                                    random_state=42)
        gmm_model.fit(roi_data.reshape(-1, 1))
        means = gmm_model.means_.flatten()
        sds = np.sqrt(gmm_model.covariances_.flatten())
        weights = gmm_model.weights_.flatten()
        c1_idx = np.argmin(means)
        bic_score = gmm_model.bic(roi_data.reshape(-1, 1))
        return np.array([means[c1_idx], sds[c1_idx], weights[c1_idx], bic_score])

    # Define the function to compute bootstrap and store results
    def compute_bootstrap(nre, roi_data):
        bre = bootstrap((roi_data.flatten(),), g, random_state=42, vectorized=False, n_resamples=nre)
        distdf = pd.DataFrame(bre.bootstrap_distribution.T, columns=stats)
        sbre_df = pd.DataFrame(columns=['roi', 'stat', 'n_res', 'low', 'high', 'se'])
        for s in stats:
            i = stats.index(s)
            sbre_df = sbre_df.append({
                'roi': roi_name,
                'stat': [s],
                'n_res': nre,
                'low': bre.confidence_interval.low.tolist()[i],
                'high': bre.confidence_interval.high.tolist()[i],
                'se': bre.standard_error[i]
            }, ignore_index=True)
        return sbre_df, distdf

    # Define the list of statistics
    stats = ['mean', 'sd', 'weight', 'bic']

    # Load data
    # Define the path to your data file
    cpath = '/Users/linlin/Library/Mobile Documents/com~apple~CloudDocs/Desktop/2024/AD lab/dataexp/pySuStaIn/notebooks/dat/'
    lobes = ['L_MTL', 'R_MTL', 'L_temporal', 'R_temporal', 'L_frontal', 'R_frontal',
             'L_occipital', 'R_occipital', 'L_parietal', 'R_parietal']
    df = pd.read_csv(cpath + 'wide_data.csv')
    dat = df.loc[:, lobes].values

    # Define the ThreadPoolExecutor
    with ThreadPoolExecutor(max_workers=10) as executor:
        futures = []
        for roi_idx, roi_data in enumerate(dat.T):
            roi_name = lobes[roi_idx]
            for nre in [5000]:
                futures.append(executor.submit(compute_bootstrap, nre, roi_data))
        completed_tasks = 0
        # Store the results
        concatenated_distributions = []
        sbre_dfs = []
        for future in futures:
            sbre_df, distdf = future.result()
            sbre_dfs.append(sbre_df)
            concatenated_distributions.append(distdf)
            completed_tasks += 1
            print(f"Progress: {completed_tasks}/10 tasks completed.")

        # Concatenate and save results
        sbre_df = pd.concat(sbre_dfs)
        sbre_df.to_csv(cpath + 'scibootstrap_results2.csv', index=False)

        
        for i, distdf in enumerate(concatenated_distributions):
            distdf.to_csv(cpath + lobes[i] + '_dist2.csv', index=False)
    
        


if __name__ == "__main__":
    start_time = time()
    main()
    end_time = time()
    execution_time = end_time - start_time
    print(f"Execution Time: {execution_time} seconds")
