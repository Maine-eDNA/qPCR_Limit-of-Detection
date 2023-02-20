# qPCR_Limit-of-Detection
Calculates limit of detection and PCR efficiency from qPCR dilution series Ct values.

# Language 
R
# Required Packages
none
# Require Input
A .csv file with dilution series qPCR results; this file must include a column named "target_conc" which is the concentration of your target sequence/DNA (must be a number) and a column named "Ct" which is the cycle at which the sample crossed the threshold (if the sample didn't amplify, leave blank or put "NA"); see example input data
# Output
Standarve curve +/- 95% confidence intervals, PCR efficiency, Limit of Blanks, Limit of Detection, and Limit of Quantitifaction, and a plot of all this
# Optional Input/Output: 
You can upload Cts of unknown samples and calculate target sequence detections and concentrations from them using the standard curve, LOD, and LOQ estimated above; see example unknown data; output will be a csv with estimated concentrations
