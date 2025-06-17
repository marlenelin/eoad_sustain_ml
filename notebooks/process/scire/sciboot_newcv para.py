import numpy as np
import pandas as pd
from concurrent.futures import ThreadPoolExecutor
from sklearn.mixture import GaussianMixture
from datetime import datetime
from scipy.stats import bootstrap

n_re = 5

# Function to fit GMM and extract statistics for the first component given one dimensional data and gmm model
def fit_gmm(gmm_model, roi_data, idx): 
    roi_data = roi_data[idx]
    gmm_model.fit(roi_data.reshape(-1, 1))
    means = gmm_model.means_.flatten()
    sds = np.sqrt(gmm_model.covariances_.flatten())
    weights = gmm_model.weights_.flatten()
    c1_idx = np.argmin(means)
    bic_score = gmm_model.bic(roi_data.reshape(-1, 1))
    z =  (1 - means[c1_idx]) / sds[c1_idx]
    return means[c1_idx], sds[c1_idx], weights[c1_idx], bic_score, z 

# Function to perform bootstrapping and compute statistics
def bootstrap_gmm_stats(gmm_model,roi_data,n_re, rsl): #random_state list to rand sample.
    boot_means = []
    boot_sds = []
    boot_weights = []
    boot_bics = []
    boot_z = []
    for i in range(n_re):
        rs = rsl[i]
        dat_idx = np.random.RandomState(rs).randint(0,467, size=467)
        mean, sd, weight, bic, z = fit_gmm(gmm_model, roi_data, dat_idx)
        boot_means.append(mean)
        boot_sds.append(sd)
        boot_weights.append(weight)
        boot_bics.append(bic)
        boot_z.append(z)
    
    boot_means = np.array(boot_means)
    boot_sds = np.array(boot_sds)
    boot_weights = np.array(boot_weights)
    boot_bics = np.array(boot_bics)
    boot_z = np.array(boot_z)
    
    
 
    # Create dataframe for bootstrapped values
    boot_df = pd.DataFrame({ 
        'mean1': boot_means,
        'sd1': boot_sds,
        'weight1': boot_weights,
        'bic': boot_bics,
        'z': boot_z
    })

    return  boot_df
 

if __name__ == "__main__": 
    cpath = '/Users/linlin/Library/Mobile Documents/com~apple~CloudDocs/Desktop/2024/AD lab/dataexp/pySuStaIn/notebooks/dat/'
    lobes = ['L_MTL', 'R_MTL', 'L_temporal', 'R_temporal', 'L_frontal', 'R_frontal',
             'L_occipital', 'R_occipital', 'L_parietal', 'R_parietal']
    df = pd.read_csv(cpath + 'wide_data.csv')
    dat = df.loc[:, lobes].values

    bootdf = pd.DataFrame(columns=['roi', 'stat', 'value'])
    compk = {'L_MTL': 2, 'R_MTL': 2, 'L_temporal': 2, 'R_temporal': 2, 'L_frontal': 2,
             'R_frontal': 2, 'L_occipital': 3, 'R_occipital': 3, 'L_parietal': 2, 'R_parietal': 2}

    for roi_idx, roi_data in enumerate(dat.T):
        roi_name = lobes[roi_idx]  
        gmmm = GaussianMixture(n_components=compk[roi_name],
                                 init_params='random_from_data',
                                 covariance_type='full',
                                 n_init=100,
                                 random_state=42)
        roi_boot_df = bootstrap_gmm_stats(gmmm, roi_data, n_re, np.random.randint(0, 9999999, size=n_re))
        roi_df = pd.DataFrame(columns=['roi', 'stat', 'value'])  
        roi_df[['stat','value']] = pd.melt(roi_boot_df, var_name='stat', value_name='value')
        roi_df['roi'] = roi_name  
        bootdf = pd.concat([bootdf, roi_df], ignore_index=True) 

    bootdf.to_csv('bootstrapped_stats.csv', index=False)





