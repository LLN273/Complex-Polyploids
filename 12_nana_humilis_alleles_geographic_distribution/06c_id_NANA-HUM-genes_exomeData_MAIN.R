
## Script used to identify B. pubescens alleles of B. nana or B. humilis origin
## Luis Leal 2023


##### Clear all states
if(!is.null(dev.list())) dev.off()
par(mfrow=c(1,1))
rm(list=ls(all=TRUE))

#library(BiocManager)
#BiocManager::install('ape')
#BiocManager::install('ips')
library("ape")
library("ips")

## Reference APE package
## Paradis E. et al.  (2004 ) APE: analyses of phylogenetics and evolution in R language . Bioinformatics , 20 , 289 –290.


#Verify versions of R and other loaded software:
#sessionInfo()


## read arguments
args <- commandArgs(TRUE)

## Import pre-ASTRAL file containing gene trees estimated using IQTREE2
gene_Trees_file <- args[1]
gene_Trees <- read.table(gene_Trees_file, header=FALSE, sep="", dec=".")

# Polyploid
POLY_ID <- args[2]

# Outgroup
outgroup_ID <- args[3]

# output folder
OUTFolder <- args[4]

# output file
OUTFile <- args[5]

#gene list
GL_file <- args[6]
GL <- read.table(GL_file, header=FALSE, sep="", dec=".")



## for each  gene tree, identify sister species

nana_list <- c()
hum_list <- c()
pend_list <- c()
pendplaty_list <- c()

for (i in c(1:nrow(gene_Trees))){

   mytr <- read.tree(text = gene_Trees[i,1])
   mytr <- root(mytr , outgroup_ID)

   #### Find sister species to true parental species
   sister_ANC1 <- sister(mytr, POLY_ID, type = "terminal", label = TRUE)
   sister_ANC1 <- paste(sort(sister_ANC1), collapse = '_')
   #print(sister_ANC1)


   if (sister_ANC1 == "NANA_FINLAND_ENONTEKIO.A002_2") {
      nana_list <- append(nana_list, GL[i,1] )
   }

   if (sister_ANC1 == "HUM02") {
      hum_list <- append(hum_list, GL[i,1] )
   }
   
   if (sister_ANC1 == "PENDULA_FINLAND_LOIMAA.LOIMAA") {
      pend_list <- append(pend_list, GL[i,1] )
   }
   
   if (sister_ANC1 == "PENDULA_FINLAND_LOIMAA.LOIMAA_PLATYPHYLLA_RUSSIA.A001_1") {
      pendplaty_list <- append(pendplaty_list, GL[i,1] )
   }

}



## convert list to df

nana_list_df <- as.data.frame(table(nana_list))
nana_list_df <- nana_list_df[,-c(2)]

hum_list_df <- as.data.frame(table(hum_list))
hum_list_df <- hum_list_df[,-c(2)]

pend_list_df <- as.data.frame(table(pend_list))
pend_list_df <- pend_list_df[,-c(2)]

pendplaty_list_df <- as.data.frame(table(pendplaty_list))
pendplaty_list_df <- pendplaty_list_df[,-c(2)]


############ Save output to file

setwd(OUTFolder)

OUTF_1 <- paste(OUTFile, "-nana.txt", sep = "")
OUTF_2 <- paste(OUTFile, "-humilis.txt", sep = "")
OUTF_3 <- paste(OUTFile, "-pendula.txt", sep = "")
OUTF_4 <- paste(OUTFile, "-pendplaty.txt", sep = "")

write.table(nana_list_df, 
            file = OUTF_1, sep = "\t", quote = FALSE, row.names=FALSE, col.names=FALSE)

write.table(hum_list_df, 
            file = OUTF_2, sep = "\t", quote = FALSE, row.names=FALSE, col.names=FALSE)

write.table(pend_list_df, 
            file = OUTF_3, sep = "\t", quote = FALSE, row.names=FALSE, col.names=FALSE)

write.table(pendplaty_list_df, 
            file = OUTF_4, sep = "\t", quote = FALSE, row.names=FALSE, col.names=FALSE)


