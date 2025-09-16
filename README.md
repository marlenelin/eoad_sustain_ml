# EOAD SuStaIn Subtypes Project Repository

This repository contains codes and file specification for the project <i>SuStaIn-based subtyping of early-onset Alzheimer’s disease (EOAD) using baseline Flortaucipir-PET data from the Longitudinal Early-Onset Alzheimer's Disease Study (LEADS)</i>. All components are cataloged in `file_description.xlsx`.

You can find the e-poster (for AAIC 2025) [here](https://github.com/marlenelin/eoad_sustain_ml/blob/main/MarleneLin_AAIC2025_e-poster.pdf).

---

## Code Specification

### `pySuStaIn/`
SuStaIn model source code and custom wrappers for model fitting and plotting (cloned from [ucl-pond/pySuStaIn](https://github.com/ucl-pond/pySuStaIn) prior to 202502)

### `notebooks/`
Python codes for data cleaning, running SuStaIn, and most of the analyses:
- **process/** — model input preparation: go through `suvr_calculate.ipynb` to generate Lobar ROI SUVRs for Flortaucipir PET, then run `bootstrap trial.ipynb` + `gmm_standardize.ipynb` to derive 2-Gaussian Mixture Model based z-scores as input to the SuStaIn model 
- **data_merger/** — (not uploaded) notebooks for data-cleaning and merging
- **additional_sens_analysis/** — notebooks for sensitivity analyses and analyses added during review rounds (e.g. alternative subtype numbers, subtype probability exclusion criteria, supplementary addition, characteristics of participants who changed subtype at follow-ups)
- **`re_inter_c2mean_20.ipynb`** — running SuStaIn model and performing most of the subsequent (tabular) data analyses + visualization here

### `R_code/`
R codes for linear modeling of cognitive data
- **`sustain_cogdata.Rmd`** — cognitive data modeling, includes both baseline (interaction between stage and subtype) as well as longitudinal (linear mixed-effect [LME] modeling with baseline subtype and time from baseline interactions). also includes missingness and correlation check among scores
- **`table2.Rmd`** — generate table2 in the main text (slope and slope comparison from LME, longitudinal participant count)


### `matlab_scripts/`

MATLAB scripts for image visualizations and voxelwise analyses

- **`baseline_subtype_ss_avg.mat`** — visualize average Flortaucipir PET SUVR images for subtypes and subtype and stage groups
- **baseline subtype-versus-rest comparison** — run `flexible_doublep_baseline_comparison.m` to perform baseline voxelwise comparison between one subtype and the rest of the subtypes combined for different covariate, significance threshold, and data modalities, then run `onevr_double_p_combo.m` to visualize the comparison t-value maps with double threshold (FWE corrected p < 0.05 and uncorrected p < 0.001)
- **voxelwise LME modeling** — change settings in `run_vs_script.m` to read in the correct filepaths for different data modalities, and then specify model equations in `batch_vs.m`, finally visualize the annual rate of change in the target variable in `day2bl_double_p_combo.m`
- **`individual.m` / `poorly.m`** — visualize individual participants or average Flortaucipir PET SUVR maps for those with poor model fit
- **brainnet config/** — configurations for brainnet viewer visualization, see `file_description.xlsx` >> `brainnet config` tab for details

---

## `file_description.xlsx`

Contains `code_spec` (specify what each notebook/R file is used for), `matlab_script` (specify what each matlab script and mat file is for), `config_spec` (alternative SuStaIn model configuration and the corresponding out files + subtype and stage assignment + baseline subtype visualization storage location), `data_spec` (specify what each csv/xlsx data file is for), and `BrainNetViewer_colormap` (configurations file for BrainNetViewer) tabs

---

## Project backup (self-note)

| Content Type               | Network Drive (petcore personal)                                                            | Work Laptop                                                               | VM                          |
| ------------------------- | -------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------- | ----------------------------- |
| Model and result files     | ✔️                                                                                           | final model files                                                          | —                            |
| Writing                   | ✔️                                                                                           | —                                                                         | —                            |
| File storage description  | ✔️                                                                                           | ✔️                                                                        | —                            |
| Python and R code         | ✔️ (backup)                                                                                  | ✔️                                                                       | —                            |
| MATLAB scripts            | ✔️                                                                                           | ✔️ (backup)                                                               | —                            |
| Data                      | ✔️ (backup)                                                                                  | ✔️                                                                       | —                            |
| Images                    | ✔️ (copied 202509)                                                                                            | -                                                                         | ✔️ (`/mnt/coredata`)         |
                                                                                                                                  
---

## Setup

```bash
conda env create -f sustain_eoad.yml
conda activate sustain_eoad
pip install -e .

# Follow pySuStaIn's instructions for installing the package-required dependencies