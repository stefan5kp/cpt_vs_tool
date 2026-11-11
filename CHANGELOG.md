# Changelog

## 2026-11-11
- The first release v0.1.0 beta

## 2026-11-10

- Introduced `topo.ipynb` and `groundwater.ipynb` to query USGS topo elevation and Ma et al. 2026 groundwater data in `dev/utilities`.

## 2026-11-05

- Introduced `geology.ipynb` to query Wills et al. (2015) surficial geology units in `dev/utilities`.

## 2026-11-01

- Updated the Zhang et al. model coefficients, Kc and m(Ic) functions, and geology-specific "a" values.

## 2026-10-27

- Consolidated all functions into a single scpt_tool.py at the repository root; each model folder now contains one example notebook and its sample data.
- Layer-Based Regression now uses the same regression algorithm and cost function as the all-in-one tool.

## 2026-08-25

- Updated the all-in-one-tool to the most recent version and added function for exporting the site specific Vs result.

## 2026-06-22

- Updated input parameters to be qc for the Zhang et al. model (without pore pressure correction).
- Updated the model offset "a" to be geology specific.

## 2026-06-01

- Updated input parameters for the ergodic model to be inverse-filtered qt and fs (Boulanger and DeJong, 2018).
- Added Kc plot and shading for Ic > 2.6 for ergodic model.
- Updated name to be Zhang et al. model.
- Updated the cost function for the layer-based regression tool.

## 2026-05-19

- Correction Made to the Robertson 2012 Calculation in the all-in-one tool.
- Ergodic model Tool is now available with the most recent version of the model and the corrected Robertson 2012 Calculation.
