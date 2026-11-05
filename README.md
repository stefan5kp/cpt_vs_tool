# SCPT_tool

This tool is designed to automate the Vs interpretation process from SCPT. As of May 15, 2026, the layer-based regression tool is available in the folder with sample files. All-in-one tool is a comprehensive tool set for SCPT - Vs interpretation or Vs estimate from CPT data. The tool is currently at beta stage, more features will be added to the tool at its official release.  

## **User Guide**  
    
1. If you wanted to purely interpret the traveltime data into stepped Vs profiles from SCPT, please use the layer-based regression tool.  
2. If you wanted to estimate Vs profiles from purely cone penetration data, please use the ergodic model tool.  
3. Site-specific joint algorithm that leverage both the traveltime and cone penetration data is included in the all-in-one tool.  


## Repository Structure  

```
SCPT_tool/
├── SCPT_tool.py              # all shared functions (CPT processing, Zhang et al. / Robertson models, layer-based regression)
├── requirements.txt
├── All_in_one_tool/          # example: SCPT_Tool_V1.ipynb + sample data (AVN_1)
├── Ergodic_Model/            # example: Ergodic_CPT-Vs_Model.ipynb + sample data
└── Layer-Based_Regression/   # example: LayerBasedRegression.ipynb + sample data
```

Each example notebook adds the repository root to its path and imports from `SCPT_tool.py`, so run the notebooks from inside their own folder. Install dependencies with `pip install -r requirements.txt`.  

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
    - Consolidated all functions into a single SCPT_tool.py at the repository root; each model folder now contains one example notebook and its sample data.  
    - Layer-Based Regression now uses the same regression algorithm and cost function as the all-in-one tool.  

2026-11-01  
    - Updated the Zhang et al. model coefficients, Kc and m(Ic) functions, and geology-specific "a" values.  
    - Site-specific "a" is now fit to the measured travel time of each layer (each CPT point represents an equal share of its layer) instead of the interval velocity.  

## Citation  
Zhang, ZS, SJ Brandenberg, JP Stewart (2026). Algorithm for interpreting velocity profiles from seismic Cone Penetration Test data, 13th National Conference on Earthquake Engineering, Portland, OR, EERI