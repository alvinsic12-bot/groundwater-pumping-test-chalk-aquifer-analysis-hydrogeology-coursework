# Groundwater Pumping Test Analysis – Unconfined Chalk Aquifer (Trumpletts Farm)

This repository contains a fully reproducible Python-based workflow for the analysis of a pumping test conducted in an unconfined Chalk aquifer at Trumpletts Farm, UK.  
The project integrates groundwater-level processing, Jacob drawdown and recovery analysis, statistical uncertainty estimation, and potentiometric surface reconstruction.

This work was developed as part of the MSc Hydrology and Water Resources Management programme at Imperial College London.

---

## Project Overview

The objectives of this project are to:

- Process multi-well groundwater level time series during a controlled pumping test  
- Separate drawdown and recovery phases for multiple observation wells  
- Estimate transmissivity (T) and storativity (S) using Jacob straight-line methods  
- Quantify parameter uncertainty using 95% confidence and prediction bands  
- Reconstruct the potentiometric surface and regional hydraulic gradient  
- Interpret spatial heterogeneity and fracture-controlled flow in the Chalk aquifer  

The workflow is fully automated and reproducible using Python.

---

## Repository Structure

This repository currently contains the following main scripts:

### 1. `drawdown_recovery_processing.py`
Processes raw water-level data for six monitoring wells:
- Converts dip measurements to groundwater heads
- Identifies pumping onset automatically
- Separates drawdown and recovery phases
- Produces a 2 × 3 panel plot for all wells

Purpose:  
> Raw data pre-processing and quality-controlled visualisation of drawdown and recovery behaviour.

---

### 2. `hydraulic_gradient_analysis.py`
Reconstructs the potentiometric surface using four key wells:
- Plane fitting to groundwater heads
- Calculation of hydraulic gradient magnitude and flow direction
- 3D potentiometric surface
- Plan-view contour map
- East–west and north–south cross sections

Purpose:  
> Spatial interpretation of groundwater flow and hydraulic gradients in the unconfined Chalk aquifer.

---

### 3. `jacob_ci_and_prediction_bands.py`
Applies the Jacob straight-line method with statistical uncertainty:
- Drawdown and recovery regression in semi-log space
- Estimation of transmissivity from regression slopes
- 95% confidence intervals and prediction bands
- Joint visualisation of parameter uncertainty

Purpose:  
> Quantification of uncertainty associated with Jacob parameter estimates.

---

### 4. `simultaneous_jacob_fitting.py`
Performs simultaneous Jacob fitting of drawdown and recovery data:
- Non-linear curve fitting using analytical Theis–Jacob solutions
- Joint inversion of transmissivity (T) and storativity (S)
- Comparison of:
  - Drawdown-only fitting  
  - Recovery-only fitting  
  - Simultaneous joint fitting  
- Automatic result visualisation and parameter table output

Purpose:  
> Physically consistent joint inversion of aquifer parameters using the full pumping–recovery cycle.

---

## Methods Summary

- Groundwater heads are computed from casing elevations and dipper measurements  
- Pumping start time is detected automatically based on time proximity  
- Drawdown is analysed using the Jacob straight-line approximation of the Theis solution  
- Recovery is analysed using the Theis–Jacob residual drawdown formulation  
- Statistical uncertainty is quantified using regression-based confidence and prediction intervals  
- Potentiometric surface is reconstructed using least-squares plane fitting  
- Hydraulic gradients and flow directions are derived analytically from the fitted surface  

All analyses assume a laterally extensive, fracture-dominated unconfined Chalk aquifer.

---

## Data

The scripts are designed to read from a structured Excel file containing:

- Manual dipper measurements (PL10A, PL10B, PL10E, BBA)
- Logged groundwater level records (BBA_log, BBA_obs)
- Time stamps and reference elevations

The Excel file is not uploaded for data protection reasons but can be substituted with the same structure.

---

## Author

**Chengxu (Alvin) Shen**  
MSc Hydrology and Water Resources Management  
Imperial College London  

---

## Notes

- This repository is intended for academic demonstration and reproducibility.
- The scripts are modular and can be adapted to other pumping-test datasets with minimal modification.
- The results are directly linked to figures and interpretation presented in the associated MSc coursework submission.
