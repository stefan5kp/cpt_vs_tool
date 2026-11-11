# cpt_vs_tool

This tool is designed to automate the Vs interpretation process from SCPT. `interpretation_tool` provides the all-in-one tool for interpreting Vs profiles from Zhang et al. (2026), Robertson (2012) and the layer-based regression approach. And `lbr` interpret Vs profiles from the conventional layer-based regression approach using only the travel time data. 

## Installation

```
pip install cpt-vs-tool
```

To run the example notebooks from a clone of the repository, install the package in editable mode with the notebook extras from the repository root:

```
pip install -e ".[notebooks]"
```

## Example

```python
import numpy as np
import pandas as pd
from cpt_vs_tool import compute_profiles, plot_profiles, compute_VsZ, a_coeff

cpt = pd.read_csv("avn_1.csv")
tt = pd.read_csv("avn_1_tt.csv")

p = compute_profiles(
    cpt["depth (m)"].to_numpy(),                # CPT depth
    cpt["qt (MPa)"].to_numpy(),                 # qt
    cpt["fs (MPa)"].to_numpy(),                 # fs
    dGWT=31.0,                                  # groundwater 
    TT_depth=tt["depth (m)"].to_numpy(),        # travel time depth
    Meas_TT=tt["traveltime (ms)"].to_numpy(),   # travel time measurement
    geology="Qal1",                             # surficial geology
)

plot_profiles(p)
print(compute_VsZ(p))
```

Pass `TT_depth=None, Meas_TT=None` when there is no travel time data; only Zhang et al. (geo) and Robertson (2012) are then computed. Use `LBR` for layer-based regression of travel time data alone.

## **User Guide**  
    
1. Purely interpret the traveltime data into stepped Vs profiles from SCPT, please use the layer-based regression `lbr` tool.  
2. Interpretation tool interpret Vs profile using Zhang et al. (2026) (geo) - only CPT data, Zhang et al. (2026) (site-specific) - CPT and travel time data, Robertson (2012) and lbr.  
3. Regression and data files for the mand utility tools are included in the `dev/` folder.  
    
## Repository Structure  

```
cpt_vs_tool/
├── pyproject.toml            # package metadata and dependencies (PyPI: cpt-vs-tool)
├── requirements.txt          # editable install of the package + notebook extras
├── LICENSE                   # BSD-3-Clause
├── src/cpt_vs_tool/          # the installable package
│   ├── __init__.py           # public API (from cpt_vs_tool import ...)
│   └── vs.py                 # all functions (CPT processing, Zhang et al. / Robertson models, layer-based regression)
├── tests/                    # pytest suite + sample data (avn_1)
├── .github/workflows/        # main.yml: build, test and publish to PyPI on a tag push
├── interpretation_tool/      # example: scpt_tool_v1.ipynb + sample data (avn_1)
├── lbr/                      # example: layerbasedregression.ipynb + sample data
└── dev/                      # model development (not needed to use the tools)
    ├── src_files/            # VSPDB query (vspdb_2_query.ipynb) and the CPT, travel time, groundwater and geology tables  
    ├── groundwater/          # groundwater table depths from Ma et al. (2026) and USGS topography, with validation  
    ├── data_observation/     # data_observation.ipynb + residual data (cpt_res_data_r3.csv, tt_res_data_r3.csv)  
    ├── regression/           # model regression: v6_mle.ipynb (maximum likelihood) and bayesian_pymc.ipynb  
    ├── results/              # results_comparison_gui.ipynb: compares the CPT-Vs models with the measured profiles  
    ├── sl_validation/        # sl_validation.ipynb: all-in-one tool applied to the Perris and UCD validation sites  
    └── utilities/            # geology.ipynb: Wills et al. (2015) surficial geology / Vs30 query (wills_shp/ shapefile)  
```

The tool notebooks import from the installed `cpt_vs_tool` package (see Installation), and all notebooks read their data through relative paths, so run them from inside their own folder. `pip install -r requirements.txt` is equivalent to `pip install -e ".[notebooks]"`; the `dev/` notebooks need additional packages (e.g. pymc, geopandas).  

## Update Log

See [CHANGELOG.md](https://github.com/stefan5kp/cpt_vs_tool/blob/main/CHANGELOG.md) for the full history of updates.

## Citations  
Zhang, ZS, SJ Brandenberg, JP Stewart (2026). Algorithm for interpreting velocity profiles from seismic Cone Penetration Test data, 13th National Conference on Earthquake Engineering, Portland, OR, EERI  

Zhang, ZS, SJ Brandenberg, JP Stewart (2027). CPT-Based Vs Estimation with or without Travel Time Measurements, ASCE2027, Philadelphia, PA, ASCE  