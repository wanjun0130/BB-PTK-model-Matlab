# BB-PTK-model-Matlab Code

## Overview

This repository provides the Matlab implementation of a burden-based pulmonary tissue kinetics (BB-PTK) model for predicting time-dependent lung dosimetry of organic chemicals following pulmonary deposition.

The model represents the lung as two regions simulated in parallel:

- Tracheobronchial region (TB)
- Alveolar region (AL)

Each region is described using an 8-compartment micro-anatomical transport framework:

- aEp: apical epithelial lining / surface layer  
- imEp: immune cells on epithelium  
- cEp: epithelial cytosol  
- int: interstitium  
- sm: smooth muscle  
- imInt: immune cells in interstitium  
- cEd: endothelial cytosol  
- p: plasma  

The BB-PTK model predicts pulmonary dosimetry based on minimal physicochemical inputs:

- logPN (neutral lipophilicity)
- pKa (acid dissociation constant)
- z (charge state at physiological pH)

---

## Methodological Improvements

This implementation extends the published pulmonary transport modeling framework of Yu and Rosania (2010) with several methodological and implementation-level improvements.

### 1. Stabilized Nernst–Planck formulation

The original Nernst–Planck factor,

φ(N) = N / (exp(N) − 1)

was reformulated using a numerically stable implementation:

- φ(N) ≈ 1 − N/2 + N²/12, for |N| < 1e−6
- φ(N) = N / (exp(N) − 1), otherwise

This modification avoids numerical instability caused by catastrophic cancellation when N is close to zero.

---

### 2. Explicit treatment of neutral molecules (z = 0)

Neutral species are handled explicitly by bypassing ionization partitioning:

```matlab
isNeutralMol = abs(z) < 1e-12;
if isNeutralMol
    fd = 0;
end
```

In the actual implementation, the neutral-only branch sets all ionized fractions to zero and bypasses Henderson–Hasselbalch splitting for z = 0 species. This avoids non-physical ionization-related artifacts in transport calculations.

### 3.Transformation into a high-throughput BB-PTK platform

The original compartmental transport framework was further extended into a high-throughput burden-based pulmonary tissue kinetics (BB-PTK) platform.

(a) Mass-based dosimetry representation

Model outputs are converted from concentration-based state variables to compartmental amounts using compartment volumes:

Z = Y × M_v

This enables direct calculation of burden-based pulmonary dosimetry metrics in molar units.

(b) Automated batch computation

The model was expanded from a single-compound computational framework to an automated batch-processing workflow for multiple compounds.

(c) Output of dosimetry metrics

The high-throughput implementation supports automated extraction of key dosimetry outputs, including:

Lung burden (mol)
Plasma burden (mol)
Tmax
AUP

Lung burden is defined as the sum of all lung-associated compartmental amounts:

Lung burden (mol) = Z(1) + Z(2) + Z(3) + Z(4) + Z(5) + Z(6) + Z(7)

These modifications improve numerical robustness, ensure mechanistic consistency for neutral compounds, and enable scalable dosimetry-oriented simulation across multiple chemicals.

---

## Repository Structure

BB-PTK-model-Matlab/

├─ README.md

├─ LICENSE

├─ src/

│ ├─ drugil_matlab_400nmol_final.m

│ ├─ LungModel_Rat_tb_D_15_dense_huge_free_final.m

│ ├─ LungModel_Rat_al_D_15_dense_huge_free_final.m

│ └─ LungModelOde.m

├─ input/

│ └─ input_template.xlsx

└─ output/

---

## Main Files

### 1. drugil_matlab_400nmol_final.m
Main batch execution script.

Functions:
- Reads compound information from Excel
- Performs input sanitation
- Runs TB and AL simulations
- Interpolates results to a unified time grid
- Combines TB + AL outputs
- Exports per-compound and summary results

---

### 2. LungModel_Rat_tb_D_15_dense_huge_free_final.m
BB-PTK model for tracheobronchial (TB) region.

---

### 3. LungModel_Rat_al_D_15_dense_huge_free_final.m
BB-PTK model for alveolar (AL) region.

---

### 4. LungModelOde.m

```matlab
function [dy] = LungModelOde(t,x,m,g)
dy = m*x + g;
end
```
---

## Model Inputs
Required columns in Excel:

  Column	       Description

  pKa            Acid dissociation constant
  
  logPN	         Neutral lipophilicity
  
  pKa	           Acid dissociation constant
  
  z	             Charge (-1 / 0 / +1)
  
  obs_T50	       Optional experimental T50

## Dose Setting
Total dose:dose_all = 400 nmol

Regional partition:
- TB: 280 nmol (70%)
- AL: 120 nmol (30%)

## Numerical Implementation
- ODE solver: ode15s
- Simulation duration: 7 days (604800 s)
- Tolerances:
  - RelTol = 1e-8
  - AbsTol = 1e-10
- Stability improvements:
  - expm1-based Nernst–Planck implementation
  - neutral molecule bypass (z = 0)
  - NaN/Inf pre-checks

## Model Outputs
Time-dependent outputs
- Lung burden (total & tissue)
- Plasma concentration
- Compartment-level profiles

Summary metrics
- T50
- Lung burden fractions
- Compartment contributions
- Peak metrics
- Regional metrics (TB & AL)

## Input File Format
Example:

  DrugName	 logPN	 pKa	   z	   obs_T50
  
  Chem1	     2.30	   8.10	   1	   35
  
  Chem2	     1.75	   6.50	   -1	   48

Notes:
- z must be -1 / 0 / +1
- pKa required for charged species

## How to Run
Step 1. Place files
Ensure all .m files are in /src

Step 2. Update paths
inputXlsx  = fullfile('input','input_template.xlsx');
folderPath = fullfile('output');

Step 3. Run in Matlab
drugil_matlab_400nmol_final

Step 4. Outputs
Generated in /output:
- summary_all.csv
- per-compound CSVs
- merged Excel results

## Matlab Version
Developed with:Matlab 2020a

## Notes
- Plasma acts as a one-way sink (no recirculation)
- Single-dose simulation
- Designed for high-throughput screening
- Users should verify input consistency before running

## Recommended Citation

If you use this model, please cite the associated thesis or manuscript.

(Citation details will be updated after publication)

## License and Contact

This code is provided for academic research use.

Users are strongly requested to contact the author by email before using this model.

Contact:

  Name: Wanjun Zhang
  
  Email: zhangwanjun0130@163.com
  
  See LICENSE for details.

---
