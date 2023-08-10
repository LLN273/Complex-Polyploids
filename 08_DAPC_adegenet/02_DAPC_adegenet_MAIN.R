
## Discriminant Analysis of Principal Components with adegenet
## Luis Leal 2021


##### Clear all states
#if(!is.null(dev.list())) dev.off()
par(mfrow=c(1,1))
rm(list=ls(all=TRUE))

#library(BiocManager)
#install.packages("adegenet")
#install.packages("pegas")
library(adegenet)
library(pegas)	# required to import vcf files
library(vcfR)

## Reference adegenet package
## https://github.com/thibautjombart/adegenet/wiki
## https://bmcgenomdata.biomedcentral.com/articles/10.1186/1471-2156-11-94#Sec8

#Verify versions of R and other loaded software:
sessionInfo()


##### Input folder
IN_FOLDER <- "/home/snufkin/Dropbox/PlantEco/03_Paper_Origin_Bpubescens/01_DAPC_covariant_markers"

##### Input VCF
InVCF <- "MASTER_SingleRuns_ALLsamples_GOODonly_NANA_SILVER_WHITE_HUMILIS_PLATYPHYLLA-filtered_tranche_2_PASS-BIALLELIC_FINAL_NoPrivateAlleles_NONCODING_LD_CLEAN_100GENOT_REDUX.vcf.gz"

##### file with ploidy for each sample
PLD <- "00_SamplePloidy_nana.pendula.pubescens.platy_FORMAS-GeneTREE.txt"


########################### Load VCF file

# Notes:
# The input VCF file can be compressed (*.gz).
# from, to:		the loci to read; by default, the first 10,000.

setwd(IN_FOLDER)
my.loci <- read.vcf(InVCF, to = 50000)   

# sample names
sample_list <- row.names(my.loci)

# number of samples
N_SAMPLES <- length(sample_list)



########################### Read ploidy file; create vector with ploidy values

SPLOIDY <- read.table(PLD, header=T)

ploidy_df <- data.frame(matrix(ncol = 3, nrow = N_SAMPLES))
colnames(ploidy_df) <- c('sample', 'ploidy','pop')
ploidy_df$sample <- sample_list
ploidy_df$ploidy <- "0"
ploidy_df$pop <- "0"

for(i in 1:N_SAMPLES) {
  aux_pl <- which(SPLOIDY$Sample==ploidy_df[i,1])
  ploidy_df[i,2] <- SPLOIDY[aux_pl,2]   
  ploidy_df[i,3] <- SPLOIDY[aux_pl,3] 
}

# check whether there are missing values
sum(is.na(ploidy_df))

ploidy_df$ploidy <- as.factor(ploidy_df$ploidy)
ploidy_df$pop <- as.factor(ploidy_df$pop)
str(ploidy_df)

# vector with ploidy for each sample
ploidy_vec <- as.vector(ploidy_df$ploidy)

# vector with population info
pop_vec <- as.vector(ploidy_df$pop)



########################### Convert VCf to genind format

genind_data <- loci2genind(my.loci, ploidy = ploidy_vec, na.alleles = c("."), unphase = TRUE)

# Add population info
genind_data$pop <- as.factor(pop_vec)



#### Perform Discriminant Analysis of Principal Component (DAPC)
## DAPC is implemented by the function dapc, which first transforms the data using PCA,
## and then performs a Discriminant Analysis on the retained principal components.

dapc1 <- dapc(genind_data, genind_data$pop, n.pca = 11, n.da = nPop(genind_data) - 1)

dapc1

### plot

# rotate plot axis (if required)
dapc1$ind.coord <- -dapc1$ind.coord    

myCol <- c("darkorange1","purple","plum1","azure3","green","darkgreen","peachpuff3","black","red","firebrick4","cyan","blue",
           "gold","cyan4","darkcyan","dimgray","chocolate4","antiquewhite3","lightpink","red","deeppink",
           "seagreen1","sandybrown","dimgray")
scatter(dapc1, xax=1, yax=2, scree.da=FALSE, bg="white", pch=20, cell=0, cstar=0, col=myCol, solid=0.7, cex=3,clab=0, leg=TRUE, cleg=1)




####### Find optimal number of clusters
####### https://raw.githubusercontent.com/thibautjombart/adegenet/master/tutorials/tutorial-dapc.pdf

## find optimal number of clusters, based on BIC criteria
grp <- find.clusters(genind_data, max.n.clust=20) 

grp$size   # number of elements in each group

