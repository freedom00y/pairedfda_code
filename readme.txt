## Install the R package “pairedfda”
1. Use command line tool to Navigate to the folder which “pairedfda_3.3.tar.gz”
2. Run “R CMD INSTALL pairedfda_3.3.tar.gz”

## Simulation
All files are listed in the “Simulation_ABS” folder
1. “onesimu.R” is the template for Scenario 1-3; “simu_template.R” can run the results of Scenario 1-3.
2. “onesimu_outlier.R” is the template for Scenario 4 (observation outliers); “simu_template_outlier.R” can run the results of Scenario 4.
3. “onesimu_curveoutlier.R” is the template for Scenario 5 (function outliers); “simu_template_curveoutlier.R” can run the results of Scenario 5.
4. “onesimu_kakb.R” is the template for selecting the number of PCs; “simu_kakb_template.R” can run the results of Table 3.

## Real Data Analysis
All files are listed in the “Real_Analysis_RI” folder
1. “real_data_latest.Rmd” analysis the type Ia supernova light curve data. The figures in the real data analysis can be reproduced by this file.