
## Analysis of fsc results; compute mean values
##
## Based on "getCI.R" script, available at https://github.com/OB-lab/James_et_al._2021-MBE/blob/81cbade4c7260fc4c8c8b42c96f5937c239f9457/fastsimcoal2/getCI.R   (James_et_al._2021-MBE)
## 


####### load libraries
library(rcompanion)


##### read arguments
args <- commandArgs(TRUE)

## Load input file
fsc_results_file <- args[1]
fsc_results <- read.table(fsc_results_file, skip=1, header=FALSE, sep="\t", dec=".")

## Output file
OUTfile <- args[2]

## Output folder
OUTfolder <- args[3]



######### Compute mean value and lower and upper 95% confidence interval for each parameter
ci <- c()
for (i in seq(1,11)) {
   j1 <- (3*i)-2 
   j2 <- (3*i)-1
   j3 <- (3*i)
   aux1 <- paste("V",i, sep="")
   ci[j1] <- as.numeric(groupwiseMean(var =aux1, data=fsc_results, conf=0.95, digits=6)[3])		# Mean
   ci[j2] <- as.numeric(groupwiseMean(var =aux1, data=fsc_results, conf=0.95, digits=6)[5])		# Percentile.lower
   ci[j3] <- as.numeric(groupwiseMean(var =aux1, data=fsc_results, conf=0.95, digits=6)[6])		# Percentile.upper
}



##### Save results to file

# Parameters named in the same order that specified in the parameters file
ci <- t(as.data.frame(ci))
colnames(ci) <- c("PUBpopsize_Mean","PUBpopsize_min", "PUBpopsize_max", "PENDpopsize_Mean", "PENDpopsize_min", "PENDpopsize_max", "PLATYpopsize_Mean", "PLATYpopsize_min", "PLATYpopsize_max", "ANCpopsize1_Mean", "ANCpopsize1_min", "ANCpopsize1_max", "ANCpopsize2_Mean", "ANCpopsize2_min", "ANCpopsize2_max", 				  "TDIV1_Mean", "TDIV1_min", "TDIV1_max", "MIG0_pubpen_Mean", "MIG0_pubpen_min", "MIG0_pubpen_max", "TDIV2_Mean", "TDIV2_min", "TDIV2_max", "Nm_pubpen0_Mean", "Nm_pubpen0_min", "Nm_pubpen0_max", 
                  "MaxEstLhood_Mean", "MaxEstLhood_min", "MaxEstLhood_max", "MaxObsLhood_Mean", "MaxObsLhood_min", "MaxObsLhood_max")					  
					  

setwd(OUTfolder)					  
write.table(ci, OUTfile, sep="\t", quote=F, row.names=F)

