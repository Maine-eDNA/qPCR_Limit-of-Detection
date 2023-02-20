# qPCR_Limit-of-Detection
Calculates limit of detection and PCR efficiency from qPCR dilution series Ct values.

# Language 
R
# Required Packages
none
# Require Input
A .csv file with qPCR results for a dilution series; this file must include a column named "target_conc" which is the concentration of your target sequence/DNA (must be a number) and a column named "Ct" which is the cycle at which the sample crossed the threshold (if the sample didn't amplify, leave blank or put "NA"); see example input data
# Output
Standard curve +/- 95% confidence intervals, PCR efficiency, Limit of Blanks, Limit of Detection, and Limit of Quantification, and a plot of all this
# Optional Input/Output: 
You can upload Cts of unknown samples and estimate concentrations based on standard curve; see example unknown dataset
# References
Derived heavily from Klymus et al 2019 https://doi.org/10.1002/edn3.29 and with LOQ equation from Forootan et al. 2017 https://doi.org/10.1016%2Fj.bdq.2017.04.001
# Note
This script uses the simplest way to calculate LOD. See Klymus et al 2019 https://doi.org/10.1002/edn3.29 for more complicated model fitting options.
