
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
for (i in seq(1,25)) {
   #j1 <- (3*i)-2 
   #j2 <- (3*i)-1
   #j3 <- (3*i)
   aux1 <- paste("V",i, sep="")
   ci[i] <- as.numeric(groupwiseMean(var =aux1, data=fsc_results, conf=0.95, digits=6)[3])		# Mean
   #ci[j2] <- as.numeric(groupwiseMean(var =aux1, data=fsc_results, conf=0.95, digits=6)[5])		# Percentile.lower
   #ci[j3] <- as.numeric(groupwiseMean(var =aux1, data=fsc_results, conf=0.95, digits=6)[6])		# Percentile.upper
}



##### Save results to file

# Parameters named in the same order that specified in the parameters file
ci <- t(as.data.frame(ci))
		
colnames(ci) <- c("PUBpopsize_Mean", "PENDpopsize_Mean", "PLATYpopsize_Mean", "ANCpopsize1_Mean", "ANCpopsize2_Mean", "TDIV1_Mean", 
                  "MIG0_pubpen_Mean",  
                  "MIG0_pubpla_Mean", 
                  "MIG0_penpub_Mean", 
                  "MIG0_penpla_Mean", 
                  "MIG0_plapen_Mean", 
                  "MIG0_plapub_Mean",  
                  "MIG1_pubpen_Mean",  
                  "MIG1_penpub_Mean",  
                  "TDIV2_Mean",  
                  "Nm_pubpen0_Mean", 
                  "Nm_pubpla0_Mean", 
                  "Nm_penpub0_Mean",  
                  "Nm_penpla0_Mean",  	
                  "Nm_plapen0_Mean",  	
                  "Nm_plapub0_Mean",  
                  "Nm_pubpen1_Mean",  
                  "Nm_penpub1_Mean",  				  				  				  			  				  
                  "MaxEstLhood_Mean", "MaxObsLhood_Mean")				  

setwd(OUTfolder)					  
write.table(ci, OUTfile, sep="\t", quote=F, row.names=F)

