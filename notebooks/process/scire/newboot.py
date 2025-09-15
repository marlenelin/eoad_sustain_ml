import warnings
import numpy as np
import pandas as pd
from sklearn.mixture import GaussianMixture
from scipy import stats
from sklearn.exceptions import ConvergenceWarning
from concurrent.futures import ThreadPoolExecutor, as_completed
from time import time
from tqdm import tqdm

# Suppress convergence warnings
warnings.filterwarnings("ignore", category=ConvergenceWarning)

# Load data
cpath = '/Users/linlin/Library/Mobile Documents/com~apple~CloudDocs/Desktop/cap/dataexp/pySuStaIn/dat/'
lobes = ['L_MTL', 'R_MTL', 'L_temporal', 'R_temporal', 'L_frontal', 'R_frontal', 
         'L_occipital', 'R_occipital', 'L_parietal', 'R_parietal']  
df = pd.read_csv(cpath + 'wide_data.csv')   
dat = df.loc[:, lobes].values 
cn_idx = df.index[df['dx'] == 'CN']
eoad_idx = df.index[df['dx'] == 'EOAD'] 
eononad_idx = df.index[df['dx'] == 'EOnonAD'] 

# Perform bootstrap resampling
num_iterations = 5000
bootstrap_results = []

def fit_gmm(gmm_model, roi_data):
    gmm_model.fit(roi_data.reshape(-1, 1))
    means = gmm_model.means_.flatten()
    sds = np.sqrt(gmm_model.covariances_.flatten())
    weights = gmm_model.weights_.flatten()
    c1_idx = np.argmin(means)
    bic_score = gmm_model.bic(roi_data.reshape(-1, 1))
    z_scores = (roi_data - means[c1_idx]) / sds[c1_idx]
    slope, intercept, r_value, p_value, std_err = stats.linregress(z_scores, roi_data)
    return means[c1_idx], sds[c1_idx], weights[c1_idx], intercept, slope, bic_score

def resample_subject_level(data, seed):
    np.random.seed(seed)
    return data.sample(frac=1, replace=True, random_state=seed)

def bootstrap_iteration(i):
    seed = 42 + i  # Different seed for each iteration to ensure reproducibility

    # Resample datasets
    cn_data = df.loc[cn_idx]
    cn_resampled = resample_subject_level(cn_data, seed)
    data_no_eononad = df.loc[df.index.difference(eononad_idx)]
    data_with_eononad = df

    resampled_no_eononad = resample_subject_level(data_no_eononad, seed)
    resampled_with_eononad = resample_subject_level(data_with_eononad, seed)

    results = []
    for lobe in lobes:
        # Standardization for CN
        cn_mean = cn_resampled[lobe].mean()
        cn_sd = cn_resampled[lobe].std()

        cn_resampled['z_scores'] = (cn_resampled[lobe] - cn_mean) / cn_sd
        cn_slope, cn_intercept, cn_r_value, cn_p_value, cn_std_err = stats.linregress(cn_resampled['z_scores'], cn_resampled[lobe])

        # Fit GMM models and record values
        gmm_model_2_no_eononad = GaussianMixture(n_components=2,
                                                 init_params='random_from_data',
                                                 covariance_type='full',
                                                 n_init=100,
                                                 tol=0)
        gmm_model_2_with_eononad = GaussianMixture(n_components=2,
                                                   init_params='random_from_data',
                                                   covariance_type='full',
                                                   n_init=100,
                                                   tol=0)
        gmm_model_3_with_eononad = GaussianMixture(n_components=3,
                                                   init_params='random_from_data',
                                                   covariance_type='full',
                                                   n_init=100,
                                                   tol=0)

        gmm2_no_eononad_stats = fit_gmm(gmm_model_2_no_eononad, resampled_no_eononad[lobe].values)
        gmm2_with_eononad_stats = fit_gmm(gmm_model_2_with_eononad, resampled_with_eononad[lobe].values)
        gmm3_with_eononad_stats = fit_gmm(gmm_model_3_with_eononad, resampled_with_eononad[lobe].values)

        # Store results for this iteration
        results.append({
            'roi': lobe,
            'method': 'cn',
            'mean': cn_mean,
            'sd': cn_sd,
            'bic': None,
            'intercept': cn_intercept,
            'slope': cn_slope
        })
        results.append({
            'roi': lobe,
            'method': 'gmm2',
            'mean': gmm2_with_eononad_stats[0],
            'sd': gmm2_with_eononad_stats[1],
            'bic': gmm2_with_eononad_stats[5],
            'intercept': gmm2_with_eononad_stats[3],
            'slope': gmm2_with_eononad_stats[4]
        })
        results.append({
            'roi': lobe,
            'method': 'gmm2non',
            'mean': gmm2_no_eononad_stats[0],
            'sd': gmm2_no_eononad_stats[1],
            'bic': gmm2_no_eononad_stats[5],
            'intercept': gmm2_no_eononad_stats[3],
            'slope': gmm2_no_eononad_stats[4]
        })
        results.append({
            'roi': lobe,
            'method': 'gmm3non',
            'mean': gmm3_with_eononad_stats[0],
            'sd': gmm3_with_eononad_stats[1],
            'bic': gmm3_with_eononad_stats[5],
            'intercept': gmm3_with_eononad_stats[3],
            'slope': gmm3_with_eononad_stats[4]
        })
    return results

def main():
    with ThreadPoolExecutor() as executor:
        futures = [executor.submit(bootstrap_iteration, i) for i in range(num_iterations)]
        for i, future in tqdm(enumerate(as_completed(futures)), total=num_iterations):
            bootstrap_results.extend(future.result())
            print(f"Iteration {i + 1}/{num_iterations} completed.")

    # Create final dataframe
    results_df = pd.DataFrame(bootstrap_results)

    # Transform to long format and handle duplicates
    results_long = results_df.melt(id_vars=['roi', 'method'], var_name='stat', value_name='value')

    # Save the final dataframe
    results_long.to_csv(cpath + 'standardization_results.csv', index=False)

if __name__ == "__main__":
    start_time = time()
    main()
    end_time = time()
    execution_time = end_time - start_time
    print(f"Execution Time: {execution_time} seconds")
