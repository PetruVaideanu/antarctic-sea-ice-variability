# components-of-antarctic-sea-ice-variability


This repository contains the analysis code and figures used in the study:
**"Components of Antarctic sea ice variability"** (Vaideanu et al., 2025, under review).


## 📁 Repository Structure

- `CCA_ERAS_1950_2023/`: CCA analysis using ERA5 SIC and ERSSTv5 SST (Python)
- `CCA_NSIDC_obs_1980_2022/`: CCA using observed satellite NSIDC SIC (Python)
- `CCM_R_code/`: Causal analysis using Convergent Cross Mapping (R)
- `Model_Comparison/`: Comparison of nine ESM models and ERA5 patterns (Python)

One large .nc file (SIC_ERA5_reEOF19500223.nc) used in the CCA_ERAS_1950_2023/Antarctic_Peninsula, CCA_ERAS_1950_2023/Composite_maps exceed GitHub’s 100 MB file size limit and are therefore is not included in this repository. It can be obtained by running CCA_ERAS_1950_2023/Script_CCA_ERA5.ipynb or from Zenodo:https://doi.org/10.5281/zenodo.11565956


## 📊 Reproducing Results

Each folder contains scripts or notebooks to reproduce the paper’s main and supplementary figures.

## 🔓 License

This project is licensed under the MIT License.

## 📌 Citation

```
Vaideanu, P. , et al. (2025). Components of Antarctic sea ice variability. [DOI will be added when available]
```
