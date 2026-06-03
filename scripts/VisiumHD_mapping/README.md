# VisiumHD mapping

TACCO (Mages et al. 2023) was used for mapping the Human Female Reproductive System Cell Atlas v1 reference to VisiumHD data from different reproductive organs.

- 01_Prepare_VisiumHD.py -> Loads Space Ranger outputs, attaches spatial coordinates and sample metadata, applies minimal QC filtering, and writes one <donor>_raw.h5ad per donor. 
- 02_Run_Tacco.py -> Performs cell type deconvolution of VisiumHD spots.
- 03_Coocurrance_analysis.py -> Performs coocurrance analysis and plots the results.
