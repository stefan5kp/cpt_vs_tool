# cpt_vs_tool

This tool is designed to automate the Vs interpretation process from SCPT. `interpretation_tool` provides the all-in-one tool for interpreting Vs profiles from Zhang et al. (2026), Robertson (2012) and the layer-based regression approach. And `lbr` interpret Vs profiles from the conventional layer-based regression approach using only the travel time data. 

## **User Guide**  
    
1. Purely interpret the traveltime data into stepped Vs profiles from SCPT, please use the layer-based regression `lbr` tool.  
2. Interpretation tool interpret Vs profile using Zhang et al. (2026) (geo) - only CPT data, Zhang et al. (2026) (site-specific) - CPT and travel time data, Robertson (2012) and lbr.  
3. Regression and data files for the mand utility tools are included in the `dev/` folder.  
    
## Repository Structure  

```
cpt_vs_tool/
├── scpt_tool.py              # all shared functions (CPT processing, Zhang et al. / Robertson models, layer-based regression)
├── requirements.txt
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

Each tool notebook adds the repository root to its path and imports from `scpt_tool.py`, and the `dev/` notebooks read their data through relative paths, so run the notebooks from inside their own folder. Install dependencies with `pip install -r requirements.txt`.  

## Update Log 
    
2026-05-19  
	- Correction Made to the Robertson 2012 Calculation in the all-in-one tool.  
	- Ergodic model Tool is now available with the most recent version of the model and the corrected Robertson 2012 Calculation. 

2026-06-01  
    - Updated input parameters for the ergodic model to be inverse-filtered qt and fs (Boulanger and DeJong, 2018).  
    - Added Kc plot and shading for Ic > 2.6 for ergodic model.  
    - Updated name to be Zhang et al. model.  
    - Updated the cost function for the layer-based regression tool.  

2026-06-22  
	- Updated input parameters to be qc for the Zhang et al. model (without pore pressure correction).  
	- Updated the model offset "a" to be geology specifc.  

2026-08-25  
    - Updated the all-in-one-tool to the most recent version and added function for exporting the site specific Vs result.  

2026-10-27  
    - Consolidated all functions into a single scpt_tool.py at the repository root; each model folder now contains one example notebook and its sample data.  
    - Layer-Based Regression now uses the same regression algorithm and cost function as the all-in-one tool.  

2026-11-01  
    - Updated the Zhang et al. model coefficients, Kc and m(Ic) functions, and geology-specific "a" values.  

2026-11-05  
    - Introduced `geology.ipynb` to query Wills et al. (2015) surficial geology units in `dev/utilities`.  

2026-11-10  
    - Introduced `topo.ipynb` and `groundwater.ipynb` to query USGS topo elevation and Ma et al. 2026 groundwater data in `dev/utilities`. 

## Citations  
Zhang, ZS, SJ Brandenberg, JP Stewart (2026). Algorithm for interpreting velocity profiles from seismic Cone Penetration Test data, 13th National Conference on Earthquake Engineering, Portland, OR, EERI  

Zhang, ZS, SJ Brandenberg, JP Stewart (2027). CPT-Based Vs Estimation with or without Travel Time Measurements, ASCE2027, Philadelphia, PA, ASCE  