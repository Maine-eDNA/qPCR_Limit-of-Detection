###########################################
### qPCR Dilution Curve basic Analysis
### Erin Grey, 2022-11-28
### Required package: none
### Required files: input.csv with dilutions and ct values, option unknowns.csv for prediction
### Output: Standard curve, PCR efficiency, LOB, LOD, LOQ, option target concentration estimates for unknowns
##############################################

########### ADJUST THE STUFF BELOW FOR YOUR DATASET ##################
## 1. Housekeeping - clear memory, load packages, load input data  & set constants, etc.
rm(list = ls()) # clears the global memory
dat <- read.csv("example_input.csv", header=TRUE) #path to your input data, just using an example here
n_reps <- 3 # number of replicates of each dilution (minimum is 3), adjust as needed
max_cycle_n <- 40 # number of cycles in your thermocycler program, adjust as needed
conc_unit <- "copies/ul" # the units used for target concentrations - copies/ul, etc), adjust as needed
LOD_threshold <- 0.95 # set LOQ at 95% of detection rate, you can change based on your preference
LOQ_threshold <- 0.35 # set LOQ at 35% of CV, can change, you can change based on your preference
optional <- "yes" # if you want to upload unknown Cts and estimate target concentrations for them type "yes", otherwise this wont' happen
unknowns_file <- "example_unknown.csv" # if yes to optional_unknowns, adjust to path of your input data

########### AFTER THIS POINT YOU CAN JUST RUN THE WHOLE SCRIPT ##################
## 2. Format input data
dat$logCopy <- log10(dat$target_conc) #log the target concentration numbers
dat$Factor_conc <- factor(dat$target_conc) #make each concentration value a group for detection rate calculations
dat_noNTC <- subset(dat, target_conc >0) # subset of the data without NTCs, useful for later
dat_NTC <- subset(dat, target_conc == 0) # subset with only NTCs, useful for later

## 3. Estimate LOB, Standard Curve, & Efficiency
LOB <- min(dat[which(dat$target_conc==0),"Ct"], na.rm=TRUE) #Limit of Blank
mod_lm <- lm(data=subset(dat, target_conc >0) , Ct~logCopy) #fit a linear model to non-blank Ct values
yint <- mod_lm$coefficients[1]
slope <- mod_lm$coefficients[2]
r2 <- summary(mod_lm)$adj.r.squared
newx <- seq(min(dat_noNTC$logCopy, na.rm=TRUE), max(dat_noNTC$logCopy, na.rm=TRUE), by=0.05)
conf_intervals <- predict(mod_lm, newdata=data.frame(logCopy=newx), interval="confidence",level = 0.95)
E <- (10^(-1/coef(mod_lm)[["logCopy"]]))-1 #calculate efficiency

## 4. Estimate LOD and LOQ
detection_rates <- table(dat$Factor_conc[!is.na(dat$Ct)])/n_reps #calculate detection rate for each concentrations (=#detections/#replicates)
LOD <- min(as.numeric(names(which(detection_rates > LOD_threshold))))
sds_by_conc <- aggregate(dat$Ct, list(dat$Factor_conc), FUN=sd) 
sds_by_conc$cv <- sqrt((1+E)^((sds_by_conc$x^2)*log(1+E))-1) #CV function from Forootan et al. 2017 https://doi.org/10.1016%2Fj.bdq.2017.04.001 
LOQ <- min(as.numeric(as.character(sds_by_conc[which(sds_by_conc$x<LOQ_threshold),1])), na.rm=TRUE) #find concentration where CV < 0.35

## 5. Plot Standard Curve & stuff
plot(dat$Ct~dat$logCopy, pch=21, bg=rgb(red=0.5,blue=0.5,green=0.5,alpha=0.5),  cex=3, 
    ylim=c(min(dat$Ct, na.rm=TRUE), max_cycle_n),
    xlim=c(0, max(dat$logCopy,na.rm=TRUE)),
    xlab="log(Target Concentration)", ylab="Ct",
    main="Standard Curve")
abline(mod_lm, lwd=3, col="darkgrey") # add standard curve line
lines(newx, conf_intervals[,2], lty=2, lwd=2, col="darkgrey") #lower 95% CI line 
lines(newx, conf_intervals[,3], lty=2, lwd=2, col="darkgrey") #upper 95% CI line 
abline(v=log10(LOD), col="blue",lwd=2) # add Limit of Detection line
abline(h=yint+slope*log10(LOD), col="blue", lty=2, lwd=2) # add Limit of Detection predicted Cq line
abline(v=log10(LOQ), col="green", lwd=2) # add Limit of Quantification
abline(h=yint+slope*log10(LOQ), col="green", lty=2, lwd=2) # add Limit of Quantification predicted Cq line
# if blanks were detected, add them in red warning color
if (is.finite(LOB)) {
  points(rep(0, times=dim(dat_NTC)[1]), dat_NTC$Ct, pch=21, bg="red", cex=3)
  abline(h=LOB, col="red", lwd=3)
  abline(v=(LOB-yint)/slope, col="red", lwd=3, lty=2)
  LOB_format <- paste("LOB Ct=", signif(LOB, digits=2))
} else {
  LOB_format <- "NA"
}
legend(x="topright", pch=NULL, lty=NULL, cex=0.8,
       c(paste("Curve: y=", signif(yint,digits=2), "+", signif(slope,digits=2),"x"),
         paste("R-squared.adj =", signif(r2,digits=2)), 
         paste("Efficiency =", signif(E, digits=2)),
         LOB_format,
         paste("LOD =", signif(LOD, digits=2), conc_unit),
         paste("LOQ =", signif(LOQ, digits=2), conc_unit)))
legend(x="right", pch=NULL, lty=c(1,2,1,2,1,2,1,2), cex=0.8, col=c("darkgrey","darkgrey","red", "red","blue","blue","green", "green"),
       c("Standard Curve","Standard+/-95", "LOB", "LOB_predictx","LOD", "LOD_predicty", "LOQ", "LOQ_predicty"))

### 6. Optional
if(optional == "yes"){
  dat_unknowns = read.csv(unknowns_file, header=TRUE)
  dat_unknowns$log_target_conc_predict <- (dat_unknowns$Ct-yint)/slope
  dat_unknowns$target_conc_predict <- dat_unknowns$log_target_conc_predict^10
  dat_unknowns$better_than_blank <- ifelse(dat_unknowns$Ct>= LOB, "no", "yes")
  dat_unknowns$detectable <- ifelse(dat_unknowns$target_conc_predict>= LOD, "yes", "no")
  dat_unknowns$quantifiable <- ifelse(dat_unknowns$target_conc_predict>= LOQ, "yes", "no")
  write.table(dat_unknowns, "predicted_unknowns.csv", sep=",", row.names = FALSE)
}

### 7, Print Output
output <- c(paste("Standard Curve: Ct =", signif(yint,digits=3), "+", signif(slope,digits=3),"TargetConcentration"),
paste("R-squared.adj:", signif(r2,digits=3)),
paste("Efficiency =", signif(E, digits=3)),
paste("LOB Ct=", signif(LOB, digits=3)),
paste("LOD=", signif(LOD, digits=3), conc_unit),
paste("LOQ=", signif(LOQ, digits=3), conc_unit),
paste("unknowns predicted? ", optional))