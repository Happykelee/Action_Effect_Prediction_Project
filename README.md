# Action-Effect Prediction Project

[![Preprint](https://img.shields.io/badge/Preprint-PsyArXiv-blue.svg)](https://psyarxiv.com/[Your-Link-Here])
[![Python](https://img.shields.io/badge/Python-3.11.8-blue.svg)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

This repository contains the custom code, statistical data, and intermediate datasets used for the behavioral, univariate CSD, and multivariate decoding analyses in our paper: **Neural Dynamics of Action-Effect Predictions within Effect-Related Regions Diverge Between Stimulus-Driven and Voluntary Actions"**.

## ⚠️ Important Note on Terminology and Variable Naming
> **Please read before diving into the code:** 
> The core data analysis pipeline for this project was developed and executed prior to the final drafting of the manuscript. Consequently, you may encounter slight discrepancies between the variable names used in the codebase and the finalized terminology adopted in the paper (e.g., specific abbreviations for experimental phases or conditions). 
> **These variables refer to the exact same concepts.** To ensure transparency and ease of understanding, we have included detailed comments and mapping dictionaries at the beginning of the relevant scripts (or as inline updates) to bridge the code-to-manuscript nomenclature.

---

## 📂 Repository Structure

The analytical pipeline is modularized and organized sequentially from raw data preprocessing to final figure generation.

### 📜 Main Scripts
| File Name | Description | Software |
| :--- | :--- | :--- |
| `AEP_01_EEG_preprocessing.m` | Core pipeline for raw EEG data preprocessing and cleaning. | MATLAB |
| `AEP_01b_extract_behavioral_rt.m` | Extraction of raw single-trial reaction times (RTs). | MATLAB |
| `AEP_01c_standardize_filenames.m` | Utility script for standardizing data filenames. | MATLAB |
| `AEP_01d_group_rt_cleaning.m` | Group-level behavioral RT cleaning and outlier rejection. | MATLAB |
| `AEP_02_Behavioral_Performance.ipynb` | Statistical analysis of behavioral data (RTs and accuracy). | Python/Jupyter |
| `AEP_03a_permutation_stats.py` | Worker script for cluster-based permutation statistical testing. | Python |
| `AEP_03_CSD_Univariate_Analysis.ipynb` | Univariate analyses of event-related Current Source Density (CSD). | Python/Jupyter |
| `AEP_04a_decoding_worker_pipeline.py` | Backend pipeline and models for multivariate pattern analysis (MVPA). | Python |
| `AEP_04_CSD_Multivariate_Decoding.ipynb` | Execution and evaluation of multivariate decoding across time. | Python/Jupyter |
| `AEP_05_Figures.ipynb` | Generates the exact main and supplementary figures used in the manuscript. | Python/Jupyter |

### 📁 Data Folders
* **`matdata/`**: Contains intermediate MATLAB structures necessary for behavioral and marker alignments.
  * `markers_des.mat`: Event marker descriptions.
  * `RT_infos.mat` & `RT_infos_update.mat`: Extracted and updated reaction time information.
* **`statdata/`**: Contains aggregated statistical outputs and behavioral files for Python analyses.
  * `Behavioral_subjects.xlsx`: Subject-level demographic and global behavioral data.
  * `Cued_Accuracy.csv`: Accuracy metrics specifically for the Cued condition.
  * `dict_permutation_outputs.npy`: Serialized dictionary containing the extensive computational outputs from the cluster-based permutation tests (loads directly into the Figures notebook).

---

## 🛠️ System Requirements & Dependencies

### MATLAB & EEGLAB Environment
The raw EEG data preprocessing (`AEP_01` series) was conducted using **MATLAB R2023b** and **EEGLAB v2024.0**. To ensure exact reproducibility of the artifact rejection and filtering pipelines, the following EEGLAB plugins are required:
* `PICARD` (v1.0) - *Pion-Tonachini et al., 2019* (faster independent component analysis)
* `ICLabel` (v1.6) - *Pion-Tonachini et al., 2019* (for detecting blink, muscle artifacts and other noises)
* `CSD Toolbox` (v1.1) - *Kayser & Tenke, 2006* (current source density transformation)


### Python Environment
The main statistical analyses, multivariate decoding, and visualizations were conducted using **Python 3.11.8** within the **Anaconda 24.3.0** platform. The following key packages are required:
* `numpy` (1.26.4)
* `pandas` (2.3.2)
* `scipy` (1.13.1) - *Virtanen et al., 2020*
* `pingouin` (0.5.5) - *Vallat, 2018* (for statistical testing)
* `scikit-learn` (1.4.1) - *Pedregosa et al., 2011* (for multivariate decoding)
* `mne` (1.10.0) - *Gramfort et al., 2013* (for EEG analysis and spatial representations)
* `matplotlib` (3.8.4) - for visualization

---

## 🚀 Usage Guide

1. **Pre-requisites:** Ensure that you have the intermediate data files placed in the correct `matdata/` and `statdata/` directories as cloned from this repository.
2. **Behavioral Data:** Run `AEP_02_Behavioral_Performance.ipynb` to reproduce the RT and accuracy statistics.
3. **EEG/CSD Results:** For univariate results, run the `AEP_03` series. Note that the permutation test (`AEP_03a`) may take considerable computational time depending on your CPU.
4. **Decoding:** Run the `AEP_04` series for the MVPA results. 
5. **Visualization:** Once all data frames and `.npy` outputs are generated, run `AEP_05_Figures.ipynb` to reproduce the exact topological maps, waveforms, and bar charts seen in the manuscript.

---

## 📖 Citation

If you use this code, the statistical pipeline, or the associated datasets in your research, please cite our preprint:

> [Your Name], [Co-authors]. (2026). [Your Paper Title]. *PsyArXiv*. [https://doi.org/10.31234/osf.io/XXXXX](https://psyarxiv.com/[Your-Link-Here])

For any questions, issues, or requests regarding the codebase or data, please feel free to open an issue in this repository or contact the corresponding author at caizhengting2858@gmail.com.
