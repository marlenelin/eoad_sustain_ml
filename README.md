# EOAD SuStaIn Subtypes Project Repository

This repository contains codes and file specification for the project SuStaIn-based subtyping of early-onset Alzheimer’s disease (EOAD) using baseline Flortaucipir-PET data from the Longitudinal Early-Onset Alzheimer's Disease Study (LEADS). All components are cataloged in `file_description.xlsx`.

Manuscript is not published but you can read the e-poster for AAIC 2025 [here](https://github.com/marlenelin/eoad_sustain_ml/blob/main/MarleneLin_AAIC2025_e-poster.pdf).

---

## Code Specification

### `pySuStaIn/`
SuStaIn model source code and custom wrappers for model fitting and plotting (as of 202502)

### `notebooks/`
Python codes for data cleaning, running SuStaIn, and most of the analyses:
- **process/** — preprocessing/input preparation. run `suvr_calculate.ipynb` to generate Lobar ROI SUVRs for Flortaucipir PET, then run `bootstrap trial.ipynb` + `gmm_standardize.ipynb` to derive 2-Gaussian Mixture Model based z-scores as input to the SuStaIn model. 
- **data_merger/** — (not uploaded) notebooks for data-cleaning and merging
- **additional_sens_analysis/** — sensitivity analyses (e.g. alternative subtype numbers, subtype probability exclusion criteria, supplementary addition, participants who changed subtype at follow-ups)
- **`re_inter_c2mean_20.ipynb`** — running SuStaIn model and perform most of the subsequent (tabular) data analyses here

### `R_code/`
R codes for linear mixed-effect modeling (LME) of longitudinal cognitive data
- **`sustain_cogdata.Rmd`** — cognitive data modeling 
- **`table2.Rmd`** — generate table2 in the main text (slope and slope comparison from LME)


### `matlab_scripts/`

MATLAB scripts for image visualizations and voxelwise analyses

- **`baseline_subtype_ss_avg.mat`** — visualize average Flortaucipir PET SUVR images for subtypes and subtype and stage groups
- **baseline subtype-versus-rest comparison** — run `flexible_doublep_baseline_comparison.m` to perform baseline voxelwise comparison between one subtype and the rest of the subtypes combined for different covariate, significance threshold, and data modalities, then run `onevr_double_p_combo.m` to visualize the comparison t-value maps with double threshold (FWE corrected p < 0.05 and uncorrected p < 0.001)
- **voxelwise LME modeling** — change settings in `run_vs_script.m` to read in the correct filepaths for different data modalities, and then specify model equations in `batch_vs.m`, finally visualize the annual rate of change in the target variable in `day2bl_double_p_combo.m`
- **individual.m / poorly.m** — visualize individual participants or average Flortaucipir PET SUVR maps for those with poor model fit
- **brainnet config/** — configurations for brainnet viewer visualization, see `file_description.xlsx` >> `**`brainnet config` tab for details

---

## `file_description.xlsx`

Contains `code_spec`, `matlab_script`, `config_spec` (alternative SuStaIn model configuration and the corresponding out files + subtype and stage assignment + baseline subtype visualization storage location), `data_spec`, and `BrainNetViewer_colormap` tabs

---

## Project backup (self-note)

| Location            | Contents                                                                                                                                                                                                                                                                  |
| ------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| **Box**             | - `sustain_backup\` (model and result files backup)<br>- Old meeting notes                                                                                                                                                                                                 |
| **Network Drive (petcore personal)** | - Writing (AAIC, manuscript drafts, literature notes, presentations, figures/tables, capstone submission)<br>- File storage description<br>- Model and result files backup<br>- Python and R code backup<br>- MATLAB scripts<br>- Data backup |
| **Local work laptop**    | - R and Python code<br>- MATLAB script backup<br>- Final model files<br>- Data (**repository is based here**, but final model files and data are **not pushed to GitHub**)                                                                                                |
| **VM**              | - Nothing except raw/processed images stored under `/mnt/coredata/processing`                                                                                                                                                 
---

## Setup

```bash
conda env create -f sustain_eoad.yml
conda activate sustain_eoad
pip install -e .
