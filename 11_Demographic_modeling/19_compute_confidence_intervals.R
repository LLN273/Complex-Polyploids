
## Analysis of fsc results; estimation of confidence intervals
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



######### Compute lower and upper 95% confidence interval for each parameter
ci <- c()

for (i in seq(1,25)) {
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

colnames(ci) <- c("PUBpopsize_Mean","PUBpopsize_lower", "PUBpopsize_upper", "PENDpopsize_Mean", "PENDpopsize_lower", "PENDpopsize_upper", "PLATYpopsize_Mean", "PLATYpopsize_lower", "PLATYpopsize_upper", "ANCpopsize1_Mean", "ANCpopsize1_lower", "ANCpopsize1_upper", "ANCpopsize2_Mean", "ANCpopsize2_lower", "ANCpopsize2_upper", 				  "TDIV1_Mean", "TDIV1_lower", "TDIV1_upper", 
                  "MIG0_pubpen_Mean", "MIG0_pubpen_lower", "MIG0_pubpen_upper", 
                  "MIG0_pubpla_Mean", "MIG0_pubpla_lower", "MIG0_pubpla_upper", 
                  "MIG0_penpub_Mean", "MIG0_penpub_lower", "MIG0_penpub_upper", 
                  "MIG0_penpla_Mean", "MIG0_penpla_lower", "MIG0_penpla_upper", 
                  "MIG0_plapen_Mean", "MIG0_plapen_lower", "MIG0_plapen_upper", 
                  "MIG0_plapub_Mean", "MIG0_plapub_lower", "MIG0_plapub_upper", 
                  "MIG1_pubpen_Mean", "MIG1_pubpen_lower", "MIG1_pubpen_upper", 
                  "MIG1_penpub_Mean", "MIG1_penpub_lower", "MIG1_penpub_upper", 
                  "TDIV2_Mean", "TDIV2_lower", "TDIV2_upper", 
                  "Nm_pubpen0_Mean", "Nm_pubpen0_lower", "Nm_pubpen0_upper", 
                  "Nm_pubpla0_Mean", "Nm_pubpla0_lower", "Nm_pubpla0_upper", 
                  "Nm_penpub0_Mean", "Nm_penpub0_lower", "Nm_penpub0_upper", 
                  "Nm_penpla0_Mean", "Nm_penpla0_lower", "Nm_penpla0_upper", 	
                  "Nm_plapen0_Mean", "Nm_plapen0_lower", "Nm_plapen0_upper", 	
                  "Nm_plapub0_Mean", "Nm_plapub0_lower", "Nm_plapub0_upper", 
                  "Nm_pubpen1_Mean", "Nm_pubpen1_lower", "Nm_pubpen1_upper", 
                  "Nm_penpub1_Mean", "Nm_penpub1_lower", "Nm_penpub1_upper", 				  				  				  			  				  
                  "MaxEstLhood_Mean", "MaxEstLhood_lower", "MaxEstLhood_upper", "MaxObsLhood_Mean", "MaxObsLhood_lower", "MaxObsLhood_upper")				  
					  
setwd(OUTfolder)					  
write.table(ci, OUTfile, sep="\t", quote=F, row.names=F)

