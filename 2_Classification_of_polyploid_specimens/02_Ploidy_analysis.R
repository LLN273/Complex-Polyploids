
## Analysis of allelic counts for each sample; used to determine ploidy level
## Luis Leal 2021


##### Clear all states
if(!is.null(dev.list())) dev.off()
par(mfrow=c(1,1))
rm(list=ls(all=TRUE))



#library(BiocManager)
#install.packages("PerformanceAnalytics")
library(PerformanceAnalytics)


#Verify versions of R and other loaded software:
#sessionInfo()






## read arguments
args <- commandArgs(TRUE)


## Import Allelic Counts file (for biallelic SNPs)
AC_file <- args[1]
AllelicCounts <- read.table(AC_file, header=TRUE, sep="\t", dec=".")

# output folder
OUTFolder <- args[2]

# output file (root)
OUTFile <- args[3]







###### Process Allelic Counts file

# remove SNPs if # read < 10 for any of the alleles
AllelicCounts_clean <- AllelicCounts[AllelicCounts$REF_COUNT > 9 & AllelicCounts$ALT_COUNT > 9, ]

# Reference read-frequency at biallelic SNPs
AllelicCounts_clean$REF_freq <- AllelicCounts_clean$REF_COUNT/(AllelicCounts_clean$REF_COUNT + AllelicCounts_clean$ALT_COUNT)

# remove Reference read-frequency > 0.8 or < 0.2
AllelicCounts_clean <- AllelicCounts_clean[AllelicCounts_clean$REF_freq <= 0.8 &  AllelicCounts_clean$REF_freq >= 0.2, ]

##### histogram
myhyst <- hist(AllelicCounts_clean$REF_freq, breaks = 50)

# position of histogram peak [0.5: diploid; 0.75: tetraploid; 0.66: triploid]
hist_max <- which(myhyst$counts == max(myhyst$counts))
max_allelicFreq <- myhyst$breaks[hist_max]			

# save histogram
setwd(OUTFolder)
OutFile_plot <- paste(OUTFile, "jpeg", sep = ".")
jpeg(file=OutFile_plot)
hist(AllelicCounts_clean$REF_freq, breaks = 50)
dev.off()



# Calculate Kurtosis
kurt <- kurtosis(AllelicCounts_clean$REF_freq, method = "excess")
#kurtosis(AllelicCounts_clean$REF_freq, method = "excess")


# "method"	
# A character string which specifies the method of computation. 
# These are either "moment", "fisher", or "excess". If "excess" is selected, then the value of the kurtosis is 
# computed by the "moment" method and a value of 3 will be subtracted. The "moment" method is based on the 
# definitions of kurtosis for distributions; these forms should be used when resampling (bootstrap or jackknife). 
# The "fisher" method correspond to the usual "unbiased" definition of sample variance, although in the case of 
# kurtosis exact unbiasedness is not possible. The "sample" method gives the sample kurtosis of the distribution.







############ Save output to file


outData_df <- data.frame(max_allelicFreq=max_allelicFreq, kurtosis=kurt)

setwd(OUTFolder)

outFile_data <- paste(OUTFile, "txt", sep = ".")

write.table(outData_df, 
            file = outFile_data, sep = "\t", quote = FALSE, row.names=FALSE)


