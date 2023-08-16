# Cluster analysis of B. pubescens alleles of B. nana or B. humilis origin (based on k-means clustering analysis carried out using MATLAB)
# Luis Leal, 2021

#install.packages("matrixStats")

library(gplots)
library(pheatmap)
library(circlize)
library(grid)
library(gridExtra)
library(readr)
library(matrixStats)
library(dplyr)
library(cowplot)

#Verify versions of R and other loaded software:
sessionInfo()

########## Clear all states and remove all plots
#if(!is.null(dev.list())) dev.off()
rm(list=ls(all=TRUE))


###### Input folder 
IN_FOLDER <- "/Users/luisleal/Dropbox/PlantEco/03_Paper_Origin_Bpubescens/04_PHYLOGENIES"   

###### in files
infile_nana_ALT1 <- "ID_genes-nana_summary_ALT1.txt"      # pendula polarization
infile_hum_ALT1 <- "ID_genes-humilis_summary_ALT1.txt"
infile_pend_ALT1 <- "ID_genes-pendula_summary_ALT1.txt"
infile_pendplaty_ALT1 <- "ID_genes-pendplaty_summary_ALT1.txt"

infile_nana_ALT2 <- "ID_genes-nana_summary_ALT2.txt"      # nana polarization
infile_hum_ALT2 <- "ID_genes-humilis_summary_ALT2.txt"
infile_pend_ALT2 <- "ID_genes-pendula_summary_ALT2.txt"
infile_pendplaty_ALT2 <- "ID_genes-pendplaty_summary_ALT2.txt" 

infile_nana_ALT3 <- "ID_genes-nana_summary_ALT3.txt"      # platyphylla polarization
infile_hum_ALT3 <- "ID_genes-humilis_summary_ALT3.txt"
infile_pend_ALT3 <- "ID_genes-pendula_summary_ALT3.txt"
infile_pendplaty_ALT3 <- "ID_genes-pendplaty_summary_ALT3.txt"

infile_nana_ALT4 <- "ID_genes-nana_summary_ALT4.txt"      # humilis polarization
infile_hum_ALT4 <- "ID_genes-humilis_summary_ALT4.txt"
infile_pend_ALT4 <- "ID_genes-pendula_summary_ALT4.txt"
infile_pendplaty_ALT4 <- "ID_genes-pendplaty_summary_ALT4.txt"




#### load data

no_col_nana <- max(count.fields(file.path(IN_FOLDER, infile_nana_ALT1), sep = ","))   # finds the maximum number of columns (each row has a different number of columns)
data_nana_ALT1 <- read.table(file.path(IN_FOLDER, infile_nana_ALT1), header = FALSE, sep = ",", fill=TRUE, col.names=c("Sample",1:(no_col_nana-1)))

no_col_hum <- max(count.fields(file.path(IN_FOLDER, infile_hum_ALT1), sep = ","))   # finds the maximum number of columns (each row has a different number of columns)
data_hum_ALT1 <- read.table(file.path(IN_FOLDER, infile_hum_ALT1), header = FALSE, sep = ",", fill=TRUE, col.names=c("Sample",1:(no_col_hum-1)))

no_col_pend <- max(count.fields(file.path(IN_FOLDER, infile_pend_ALT1), sep = ","))   # finds the maximum number of columns (each row has a different number of columns)
data_pend_ALT1 <- read.table(file.path(IN_FOLDER, infile_pend_ALT1), header = FALSE, sep = ",", fill=TRUE, col.names=c("Sample",1:(no_col_pend-1)))

no_col_pendplaty <- max(count.fields(file.path(IN_FOLDER, infile_pendplaty_ALT1), sep = ","))   # finds the maximum number of columns (each row has a different number of columns)
data_pendplaty_ALT1 <- read.table(file.path(IN_FOLDER, infile_pendplaty_ALT1), header = FALSE, sep = ",", fill=TRUE, col.names=c("Sample",1:(no_col_pendplaty-1)))


no_col_nana <- max(count.fields(file.path(IN_FOLDER, infile_nana_ALT2), sep = ","))   # finds the maximum number of columns (each row has a different number of columns)
data_nana_ALT2 <- read.table(file.path(IN_FOLDER, infile_nana_ALT2), header = FALSE, sep = ",", fill=TRUE, col.names=c("Sample",1:(no_col_nana-1)))

no_col_hum <- max(count.fields(file.path(IN_FOLDER, infile_hum_ALT2), sep = ","))   # finds the maximum number of columns (each row has a different number of columns)
data_hum_ALT2 <- read.table(file.path(IN_FOLDER, infile_hum_ALT2), header = FALSE, sep = ",", fill=TRUE, col.names=c("Sample",1:(no_col_hum-1)))

no_col_pend <- max(count.fields(file.path(IN_FOLDER, infile_pend_ALT2), sep = ","))   # finds the maximum number of columns (each row has a different number of columns)
data_pend_ALT2 <- read.table(file.path(IN_FOLDER, infile_pend_ALT2), header = FALSE, sep = ",", fill=TRUE, col.names=c("Sample",1:(no_col_pend-1)))

no_col_pendplaty <- max(count.fields(file.path(IN_FOLDER, infile_pendplaty_ALT2), sep = ","))   # finds the maximum number of columns (each row has a different number of columns)
data_pendplaty_ALT2 <- read.table(file.path(IN_FOLDER, infile_pendplaty_ALT2), header = FALSE, sep = ",", fill=TRUE, col.names=c("Sample",1:(no_col_pendplaty-1)))


no_col_nana <- max(count.fields(file.path(IN_FOLDER, infile_nana_ALT3), sep = ","))   # finds the maximum number of columns (each row has a different number of columns)
data_nana_ALT3 <- read.table(file.path(IN_FOLDER, infile_nana_ALT3), header = FALSE, sep = ",", fill=TRUE, col.names=c("Sample",1:(no_col_nana-1)))

no_col_hum <- max(count.fields(file.path(IN_FOLDER, infile_hum_ALT3), sep = ","))   # finds the maximum number of columns (each row has a different number of columns)
data_hum_ALT3 <- read.table(file.path(IN_FOLDER, infile_hum_ALT3), header = FALSE, sep = ",", fill=TRUE, col.names=c("Sample",1:(no_col_hum-1)))

no_col_pend <- max(count.fields(file.path(IN_FOLDER, infile_pend_ALT3), sep = ","))   # finds the maximum number of columns (each row has a different number of columns)
data_pend_ALT3 <- read.table(file.path(IN_FOLDER, infile_pend_ALT3), header = FALSE, sep = ",", fill=TRUE, col.names=c("Sample",1:(no_col_pend-1)))

no_col_pendplaty <- max(count.fields(file.path(IN_FOLDER, infile_pendplaty_ALT3), sep = ","))   # finds the maximum number of columns (each row has a different number of columns)
data_pendplaty_ALT3 <- read.table(file.path(IN_FOLDER, infile_pendplaty_ALT3), header = FALSE, sep = ",", fill=TRUE, col.names=c("Sample",1:(no_col_pendplaty-1)))


no_col_nana <- max(count.fields(file.path(IN_FOLDER, infile_nana_ALT4), sep = ","))   # finds the maximum number of columns (each row has a different number of columns)
data_nana_ALT4 <- read.table(file.path(IN_FOLDER, infile_nana_ALT4), header = FALSE, sep = ",", fill=TRUE, col.names=c("Sample",1:(no_col_nana-1)))

no_col_hum <- max(count.fields(file.path(IN_FOLDER, infile_hum_ALT4), sep = ","))   # finds the maximum number of columns (each row has a different number of columns)
data_hum_ALT4 <- read.table(file.path(IN_FOLDER, infile_hum_ALT4), header = FALSE, sep = ",", fill=TRUE, col.names=c("Sample",1:(no_col_hum-1)))

no_col_pend <- max(count.fields(file.path(IN_FOLDER, infile_pend_ALT4), sep = ","))   # finds the maximum number of columns (each row has a different number of columns)
data_pend_ALT4 <- read.table(file.path(IN_FOLDER, infile_pend_ALT4), header = FALSE, sep = ",", fill=TRUE, col.names=c("Sample",1:(no_col_pend-1)))

no_col_pendplaty <- max(count.fields(file.path(IN_FOLDER, infile_pendplaty_ALT4), sep = ","))   # finds the maximum number of columns (each row has a different number of columns)
data_pendplaty_ALT4 <- read.table(file.path(IN_FOLDER, infile_pendplaty_ALT4), header = FALSE, sep = ",", fill=TRUE, col.names=c("Sample",1:(no_col_pendplaty-1)))


row.names(data_nana_ALT1) <- data_nana_ALT1[,1]
row.names(data_hum_ALT1) <- data_hum_ALT1[,1]
row.names(data_pend_ALT1) <- data_pend_ALT1[,1]
row.names(data_pendplaty_ALT1) <- data_pendplaty_ALT1[,1]

row.names(data_nana_ALT2) <- data_nana_ALT2[,1]
row.names(data_hum_ALT2) <- data_hum_ALT2[,1]
row.names(data_pend_ALT2) <- data_pend_ALT2[,1]
row.names(data_pendplaty_ALT2) <- data_pendplaty_ALT2[,1]

row.names(data_nana_ALT3) <- data_nana_ALT3[,1]
row.names(data_hum_ALT3) <- data_hum_ALT3[,1]
row.names(data_pend_ALT3) <- data_pend_ALT3[,1]
row.names(data_pendplaty_ALT3) <- data_pendplaty_ALT3[,1]

row.names(data_nana_ALT4) <- data_nana_ALT4[,1]
row.names(data_hum_ALT4) <- data_hum_ALT4[,1]
row.names(data_pend_ALT4) <- data_pend_ALT4[,1]
row.names(data_pendplaty_ALT4) <- data_pendplaty_ALT4[,1]


data_nana_ALT1 <- data_nana_ALT1[,-c(1)]
data_hum_ALT1 <- data_hum_ALT1[,-c(1)]
data_pend_ALT1 <- data_pend_ALT1[,-c(1)]
data_pendplaty_ALT1 <- data_pendplaty_ALT1[,-c(1)]

data_nana_ALT2 <- data_nana_ALT2[,-c(1)]
data_hum_ALT2 <- data_hum_ALT2[,-c(1)]
data_pend_ALT2 <- data_pend_ALT2[,-c(1)]
data_pendplaty_ALT2 <- data_pendplaty_ALT2[,-c(1)]

data_nana_ALT3 <- data_nana_ALT3[,-c(1)]
data_hum_ALT3 <- data_hum_ALT3[,-c(1)]
data_pend_ALT3 <- data_pend_ALT3[,-c(1)]
data_pendplaty_ALT3 <- data_pendplaty_ALT3[,-c(1)]

data_nana_ALT4 <- data_nana_ALT4[,-c(1)]
data_hum_ALT4 <- data_hum_ALT4[,-c(1)]
data_pend_ALT4 <- data_pend_ALT4[,-c(1)]
data_pendplaty_ALT4 <- data_pendplaty_ALT4[,-c(1)]


## transpose
data_nana_t_ALT1 <- as.data.frame(t(data_nana_ALT1))
data_hum_t_ALT1 <- as.data.frame(t(data_hum_ALT1))
data_pend_t_ALT1 <- as.data.frame(t(data_pend_ALT1))
data_pendplaty_t_ALT1 <- as.data.frame(t(data_pendplaty_ALT1))

data_nana_t_ALT2 <- as.data.frame(t(data_nana_ALT2))
data_hum_t_ALT2 <- as.data.frame(t(data_hum_ALT2))
data_pend_t_ALT2 <- as.data.frame(t(data_pend_ALT2))
data_pendplaty_t_ALT2 <- as.data.frame(t(data_pendplaty_ALT2))

data_nana_t_ALT3 <- as.data.frame(t(data_nana_ALT3))
data_hum_t_ALT3 <- as.data.frame(t(data_hum_ALT3))
data_pend_t_ALT3 <- as.data.frame(t(data_pend_ALT3))
data_pendplaty_t_ALT3 <- as.data.frame(t(data_pendplaty_ALT3))

data_nana_t_ALT4 <- as.data.frame(t(data_nana_ALT4))
data_hum_t_ALT4 <- as.data.frame(t(data_hum_ALT4))
data_pend_t_ALT4 <- as.data.frame(t(data_pend_ALT4))
data_pendplaty_t_ALT4 <- as.data.frame(t(data_pendplaty_ALT4))


### get list of genes (remove '' entry in list of genes)
GL_nana_ALT1 <- unique(unlist(data_nana_t_ALT1))
GL_nana_ALT2 <- unique(unlist(data_nana_t_ALT2))
GL_nana_ALT3 <- unique(unlist(data_nana_t_ALT3))
GL_nana_ALT4 <- unique(unlist(data_nana_t_ALT4))
GL_nana <- unique(c(GL_nana_ALT1, GL_nana_ALT2, GL_nana_ALT3, GL_nana_ALT4))

GL_hum_ALT1 <- unique(unlist(data_hum_t_ALT1))
GL_hum_ALT2 <- unique(unlist(data_hum_t_ALT2))
GL_hum_ALT3 <- unique(unlist(data_hum_t_ALT3))
GL_hum_ALT4 <- unique(unlist(data_hum_t_ALT4))
GL_hum <- unique(c(GL_hum_ALT1, GL_hum_ALT2, GL_hum_ALT3, GL_hum_ALT4))

GL_pend_ALT1 <- unique(unlist(data_pend_t_ALT1))
GL_pend_ALT2 <- unique(unlist(data_pend_t_ALT2))
GL_pend_ALT3 <- unique(unlist(data_pend_t_ALT3))
GL_pend_ALT4 <- unique(unlist(data_pend_t_ALT4))
GL_pend <- unique(c(GL_pend_ALT1, GL_pend_ALT2, GL_pend_ALT3, GL_pend_ALT4))

GL_pendplaty_ALT1 <- unique(unlist(data_pendplaty_t_ALT1))
GL_pendplaty_ALT2 <- unique(unlist(data_pendplaty_t_ALT2))
GL_pendplaty_ALT3 <- unique(unlist(data_pendplaty_t_ALT3))
GL_pendplaty_ALT4 <- unique(unlist(data_pendplaty_t_ALT4))
GL_pendplaty <- unique(c(GL_pendplaty_ALT1, GL_pendplaty_ALT2, GL_pendplaty_ALT3, GL_pendplaty_ALT4))

GL_nana <- GL_nana[! GL_nana %in% c('')]    # list of genes containing an allele of nana origin
GL_hum <- GL_hum[! GL_hum %in% c('')]       # list of genes containing an allele of humilis origin
GL_pend <- GL_pend[! GL_pend %in% c('')]
GL_pendplaty <- GL_pendplaty[! GL_pendplaty %in% c('')]



### Number of unique genes (when nana, hum, pend, pendplaty datasets are concatenated)
GL_total <- c(GL_nana, GL_hum, GL_pend, GL_pendplaty)
N_gene_total <- length(unique(GL_total))



### For each sample, evaluate whether gene has been detected in nana, humilis, or pendula dataset
aux_nrows <- N_gene_total 
aux_ncol <- ncol(data_nana_t_ALT1)

# create all zero dataframes
db_nana <-as.data.frame(matrix(0, nrow = aux_nrows, ncol = aux_ncol))
db_hum <- as.data.frame(matrix(0, nrow = aux_nrows, ncol = aux_ncol))
db_pend <- as.data.frame(matrix(0, nrow = aux_nrows, ncol = aux_ncol))
db_pendplaty <-as.data.frame(matrix(0, nrow = aux_nrows, ncol = aux_ncol))

row.names(db_nana) <- sort(unique(GL_total))
row.names(db_hum) <- sort(unique(GL_total))
row.names(db_pend) <- sort(unique(GL_total))
row.names(db_pendplaty) <- sort(unique(GL_total))

colnames(db_nana) <- colnames(data_nana_t_ALT1)
colnames(db_hum) <- colnames(data_hum_t_ALT1)
colnames(db_pend) <- colnames(data_hum_t_ALT1)
colnames(db_pendplaty) <- colnames(data_pendplaty_t_ALT1)


# scan ALT1 data
for (i in row.names(db_nana)) {        # for each gene i
   for (j in colnames(db_nana)) {      # check if detected in sample j
      if (i %in% data_nana_t_ALT1[,j]) {
         db_nana[i,j] <- 1  
      } 
   }
}

for (i in row.names(db_hum)) {        # for each gene i
  for (j in colnames(db_hum)) {      # check if detected in sample j
    if (i %in% data_hum_t_ALT1[,j]) {
       db_hum[i,j] <- 1  
    } 
  }
}

for (i in row.names(db_pend)) {        # for each gene i
  for (j in colnames(db_pend)) {      # check if detected in sample j
    if (i %in% data_pend_t_ALT1[,j]) {
      db_pend[i,j] <- 1  
    } 
  }
}

for (i in row.names(db_pendplaty)) {        # for each gene i
  for (j in colnames(db_pendplaty)) {      # check if detected in sample j
    if (i %in% data_pendplaty_t_ALT1[,j]) {
      db_pendplaty[i,j] <- 1  
    } 
  }
}

sum(db_nana == "1")
sum(db_hum == "1")
sum(db_pend == "1")
sum(db_pendplaty == "1")


# scan ALT2 data
for (i in row.names(db_nana)) {        # for each gene i
  for (j in colnames(db_nana)) {      # check if detected in sample j
    if (i %in% data_nana_t_ALT2[,j]) {
      db_nana[i,j] <- 1  
    } 
  }
}

for (i in row.names(db_hum)) {        # for each gene i
  for (j in colnames(db_hum)) {      # check if detected in sample j
    if (i %in% data_hum_t_ALT2[,j]) {
      db_hum[i,j] <- 1  
    } 
  }
}

for (i in row.names(db_pend)) {        # for each gene i
  for (j in colnames(db_pend)) {      # check if detected in sample j
    if (i %in% data_pend_t_ALT2[,j]) {
      db_pend[i,j] <- 1  
    } 
  }
}

for (i in row.names(db_pendplaty)) {        # for each gene i
  for (j in colnames(db_pendplaty)) {      # check if detected in sample j
    if (i %in% data_pendplaty_t_ALT2[,j]) {
      db_pendplaty[i,j] <- 1  
    } 
  }
}

sum(db_nana == "1")
sum(db_hum == "1")
sum(db_pend == "1")
sum(db_pendplaty == "1")


# scan ALT3 data
for (i in row.names(db_nana)) {        # for each gene i
  for (j in colnames(db_nana)) {      # check if detected in sample j
    if (i %in% data_nana_t_ALT3[,j]) {
      db_nana[i,j] <- 1  
    } 
  }
}

for (i in row.names(db_hum)) {        # for each gene i
  for (j in colnames(db_hum)) {      # check if detected in sample j
    if (i %in% data_hum_t_ALT3[,j]) {
      db_hum[i,j] <- 1  
    } 
  }
}

for (i in row.names(db_pend)) {        # for each gene i
  for (j in colnames(db_pend)) {      # check if detected in sample j
    if (i %in% data_pend_t_ALT3[,j]) {
      db_pend[i,j] <- 1  
    } 
  }
}

for (i in row.names(db_pendplaty)) {        # for each gene i
  for (j in colnames(db_pendplaty)) {      # check if detected in sample j
    if (i %in% data_pendplaty_t_ALT3[,j]) {
      db_pendplaty[i,j] <- 1  
    } 
  }
}

sum(db_nana == "1")
sum(db_hum == "1")
sum(db_pend == "1")
sum(db_pendplaty == "1")


# scan ALT4 data
for (i in row.names(db_nana)) {        # for each gene i
  for (j in colnames(db_nana)) {      # check if detected in sample j
    if (i %in% data_nana_t_ALT4[,j]) {
      db_nana[i,j] <- 1  
    } 
  }
}

for (i in row.names(db_hum)) {        # for each gene i
  for (j in colnames(db_hum)) {      # check if detected in sample j
    if (i %in% data_hum_t_ALT4[,j]) {
      db_hum[i,j] <- 1  
    } 
  }
}

for (i in row.names(db_pend)) {        # for each gene i
  for (j in colnames(db_pend)) {      # check if detected in sample j
    if (i %in% data_pend_t_ALT4[,j]) {
      db_pend[i,j] <- 1  
    } 
  }
}

for (i in row.names(db_pendplaty)) {        # for each gene i
  for (j in colnames(db_pendplaty)) {      # check if detected in sample j
    if (i %in% data_pendplaty_t_ALT4[,j]) {
      db_pendplaty[i,j] <- 1  
    } 
  }
}

sum(db_nana == "1")
sum(db_hum == "1")
sum(db_pend == "1")
sum(db_pendplaty == "1")





### For each population, compute fraction of individuals for which a certain gene is detected

# Populations
pops <- unique(substr(colnames(db_nana), 0, 3))
#pops <- str_replace(pops, c("-"), c(""))

# remove populations SE pop
pops <- pops[! pops %in% c('SE-')]
#pops <- append(pops,"FI-BP-19")
#aux_order <- c("SPN","ES-","NO-","DJU","SKA","VAL","MOK","LIN","GNA","DOK","BRA","BUR","PIT","JOK","NED","KAU","NOR","FI-BP-19","LT-","YT-","URL","NSV")

# reorder populations
aux_order <- c("SPN","POL","YT-","URL","LT-","NSV","NO-","FI-","DJU","SKA","VAL","ASP","MOK","LIN","GNA","DOK","BRA","BUR","PIT","JOK","NED","KAU","NOR")
pops <- pops[order(match(pops, aux_order))]

# create population df
db_nana_pop <- as.data.frame(matrix(0, nrow = aux_nrows, ncol = length(pops)))
row.names(db_nana_pop) <- row.names(db_nana)
colnames(db_nana_pop) <- pops

db_hum_pop <- as.data.frame(matrix(0, nrow = aux_nrows, ncol = length(pops)))
row.names(db_hum_pop) <- row.names(db_hum)
colnames(db_hum_pop) <- pops


# function used to create table where, for each gene, it shows the percentage of individuals on each population where said gene has been detected.

pop_presence_nana <- function(pops) {
  db_pop <- db_nana[ , grepl( paste(pops,collapse="|") , names( db_nana ) ) ]
  if (class(db_pop) == "data.frame") {
     db_pop$norm <- rowSums(db_pop)/ncol(db_pop)
  } else {
    db_pop$norm <- db_pop
  }
  return(db_pop$norm)
}

pop_presence_hum <- function(pops) {
  db_pop <- db_hum[ , grepl( paste(pops,collapse="|") , names( db_hum ) ) ]
  if (class(db_pop) == "data.frame") {
    db_pop$norm <- rowSums(db_pop)/ncol(db_pop)
  } else {
    db_pop$norm <- db_pop
  }
  return(db_pop$norm)
}

pop_presence_pend <- function(pops) {
  db_pop <- db_pend[ , grepl( paste(pops,collapse="|") , names( db_pend ) ) ]
  if (class(db_pop) == "data.frame") {
    db_pop$norm <- rowSums(db_pop)/ncol(db_pop)
  } else {
    db_pop$norm <- db_pop
  }
  return(db_pop$norm)
}

pop_presence_pendplaty <- function(pops) {
  db_pop <- db_pendplaty[ , grepl( paste(pops,collapse="|") , names( db_pendplaty ) ) ]
  if (class(db_pop) == "data.frame") {
    db_pop$norm <- rowSums(db_pop)/ncol(db_pop)
  } else {
    db_pop$norm <- db_pop
  }
  return(db_pop$norm)
}




############ compute averages for different population groups

pops2 <- c("ES", "Central_Europe", "Central_Asia","SVsouth", "SVcenter", "SVnorth", "Arctic")

db_nana_pop_REDUX <-as.data.frame(matrix(0, nrow = nrow(db_nana) , ncol = length(pops2)))
colnames(db_nana_pop_REDUX) <- pops2
row.names(db_nana_pop_REDUX) <- row.names(db_nana)

db_hum_pop_REDUX <-as.data.frame(matrix(0, nrow = nrow(db_hum) , ncol = length(pops2)))
colnames(db_hum_pop_REDUX) <- pops2
row.names(db_hum_pop_REDUX) <- row.names(db_hum)

db_pend_pop_REDUX <-as.data.frame(matrix(0, nrow = nrow(db_pend) , ncol = length(pops2)))
colnames(db_pend_pop_REDUX) <- pops2
row.names(db_pend_pop_REDUX) <- row.names(db_pend)

db_pendplaty_pop_REDUX <-as.data.frame(matrix(0, nrow = nrow(db_pendplaty) , ncol = length(pops2)))
colnames(db_pendplaty_pop_REDUX) <- pops2
row.names(db_pendplaty_pop_REDUX) <- row.names(db_pendplaty)


db_nana_pop_REDUX$ES <- pop_presence_nana(c('SPN'))
db_nana_pop_REDUX$Central_Europe <- pop_presence_nana(c("POL","YT-","LT-"))
db_nana_pop_REDUX$Central_Asia <- pop_presence_nana(c("URL","NSV"))
db_nana_pop_REDUX$SVsouth<- pop_presence_nana(c("DJU","SKA","VAL"))
db_nana_pop_REDUX$SVcenter<- pop_presence_nana(c("ASP","MOK","LIN","GNA","DOK"))
db_nana_pop_REDUX$SVnorth<- pop_presence_nana(c("BRA","BUR","PIT","JOK"))
db_nana_pop_REDUX$Arctic<- pop_presence_nana(c("NED","KAU","NOR"))

db_hum_pop_REDUX$ES <- pop_presence_hum(c('SPN'))
db_hum_pop_REDUX$Central_Europe <- pop_presence_hum(c("POL","YT-","LT-"))
db_hum_pop_REDUX$Central_Asia <- pop_presence_hum(c("URL","NSV"))
db_hum_pop_REDUX$SVsouth<- pop_presence_hum(c("DJU","SKA","VAL"))
db_hum_pop_REDUX$SVcenter<- pop_presence_hum(c("ASP","MOK","LIN","GNA","DOK"))
db_hum_pop_REDUX$SVnorth<- pop_presence_hum(c("BRA","BUR","PIT","JOK"))
db_hum_pop_REDUX$Arctic<- pop_presence_hum(c("NED","KAU","NOR"))

db_pend_pop_REDUX$ES <- pop_presence_pend(c('SPN'))
db_pend_pop_REDUX$Central_Europe <- pop_presence_pend(c("POL","YT-","LT-"))
db_pend_pop_REDUX$Central_Asia <- pop_presence_pend(c("URL","NSV"))
db_pend_pop_REDUX$SVsouth<- pop_presence_pend(c("DJU","SKA","VAL"))
db_pend_pop_REDUX$SVcenter<- pop_presence_pend(c("ASP","MOK","LIN","GNA","DOK"))
db_pend_pop_REDUX$SVnorth<- pop_presence_pend(c("BRA","BUR","PIT","JOK"))
db_pend_pop_REDUX$Arctic<- pop_presence_pend(c("NED","KAU","NOR"))

db_pendplaty_pop_REDUX$ES <- pop_presence_pendplaty(c('SPN'))
db_pendplaty_pop_REDUX$Central_Europe <- pop_presence_pendplaty(c("POL","YT-","LT-"))
db_pendplaty_pop_REDUX$Central_Asia <- pop_presence_pendplaty(c("URL","NSV"))
db_pendplaty_pop_REDUX$SVsouth<- pop_presence_pendplaty(c("DJU","SKA","VAL"))
db_pendplaty_pop_REDUX$SVcenter<- pop_presence_pendplaty(c("ASP","MOK","LIN","GNA","DOK"))
db_pendplaty_pop_REDUX$SVnorth<- pop_presence_pendplaty(c("BRA","BUR","PIT","JOK"))
db_pendplaty_pop_REDUX$Arctic<- pop_presence_pendplaty(c("NED","KAU","NOR"))

db_nana_pop_REDUX <- round(db_nana_pop_REDUX[,], digits = 3)
db_hum_pop_REDUX <- round(db_hum_pop_REDUX[,], digits = 3)
db_pend_pop_REDUX <- round(db_pend_pop_REDUX[,], digits = 3)
db_pendplaty_pop_REDUX <- round(db_pendplaty_pop_REDUX[,], digits = 3)


##### FILTERING STEP I
#####
##### remove gene if it is not detected in any of the samples

db_nana_pop_REDUX$avg <- rowMeans(db_nana_pop_REDUX[,], na.rm=TRUE)   # filter
db_nana_pop_REDUX <- db_nana_pop_REDUX[db_nana_pop_REDUX$avg > 0,]

db_hum_pop_REDUX$avg <- rowMeans(db_hum_pop_REDUX[,], na.rm=TRUE) 
db_hum_pop_REDUX <- db_hum_pop_REDUX[db_hum_pop_REDUX$avg > 0,]

db_pend_pop_REDUX$avg <- rowMeans(db_pend_pop_REDUX[,], na.rm=TRUE) 
db_pend_pop_REDUX <- db_pend_pop_REDUX[db_pend_pop_REDUX$avg > 0,]

db_pendplaty_pop_REDUX$avg <- rowMeans(db_pendplaty_pop_REDUX[,], na.rm=TRUE)   
db_pendplaty_pop_REDUX <- db_pendplaty_pop_REDUX[db_pendplaty_pop_REDUX$avg > 0,]

db_nana_pop_REDUX_store <- db_nana_pop_REDUX  # store variable
db_hum_pop_REDUX_store <- db_hum_pop_REDUX
db_pend_pop_REDUX_store <- db_pend_pop_REDUX
db_pendplaty_pop_REDUX_store <- db_pendplaty_pop_REDUX



##### FILTERING STEP II
#####
##### delete gene if gene is not observed in at least XXXX% of all individuals in at least ONE population   

#ref_value <- 0.65
ref_value <- 0 
 
db_nana_pop_REDUX$max <- rowMaxs(as.matrix(db_nana_pop_REDUX[,c(1:(ncol(db_nana_pop_REDUX)-1))]))
db_nana_pop_REDUX <- db_nana_pop_REDUX[db_nana_pop_REDUX$max > ref_value,]

db_hum_pop_REDUX$max <- rowMaxs(as.matrix(db_hum_pop_REDUX[,c(1:(ncol(db_hum_pop_REDUX)-1))]))
db_hum_pop_REDUX <- db_hum_pop_REDUX[db_hum_pop_REDUX$max > ref_value,]

db_pend_pop_REDUX$max <- rowMaxs(as.matrix(db_pend_pop_REDUX[,c(1:(ncol(db_pend_pop_REDUX)-1))]))
db_pend_pop_REDUX <- db_pend_pop_REDUX[db_pend_pop_REDUX$max > ref_value,]

db_pendplaty_pop_REDUX$max <- rowMaxs(as.matrix(db_pendplaty_pop_REDUX[,c(1:(ncol(db_pendplaty_pop_REDUX)-1))]))
db_pendplaty_pop_REDUX <- db_pendplaty_pop_REDUX[db_pendplaty_pop_REDUX$max > ref_value,]





######
###### sava dataset and carry-out k-mean cluster analysis using MATLAB
######

data_out <- db_nana_pop_REDUX
data_out <- data_out[,-c(ncol(data_out))]
data_out <- data_out[,-c(ncol(data_out))]
data_out$gene <- rownames(data_out)
data_out2 <- data_out %>% select(gene, everything())
setwd(IN_FOLDER)
write.table(data_out2, 
            file = "db_nana_pop_REDUX.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

data_out3 <- db_hum_pop_REDUX
data_out3 <- data_out3[,-c(ncol(data_out3))]
data_out3 <- data_out3[,-c(ncol(data_out3))]
data_out3$gene <- rownames(data_out3)
data_out4 <- data_out3 %>% select(gene, everything())
setwd(IN_FOLDER)
write.table(data_out4, 
            file = "db_hum_pop_REDUX.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)







######
###### save data (all genes, including cases where allele frequencies remain low across all populations)
######
data_out <- db_nana_pop_REDUX_store
data_out <- data_out[,-c(ncol(data_out))]
data_out$gene <- rownames(data_out)
data_out2 <- data_out %>% select(gene, everything())
setwd(IN_FOLDER)
write.table(data_out2, 
            file = "db_nana_pop_ALL.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

data_out <- db_hum_pop_REDUX_store
data_out <- data_out[,-c(ncol(data_out))]
data_out$gene <- rownames(data_out)
data_out2 <- data_out %>% select(gene, everything())
setwd(IN_FOLDER)
write.table(data_out2, 
            file = "db_hum_pop_ALL.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

data_out <- db_pend_pop_REDUX_store
data_out <- data_out[,-c(ncol(data_out))]
data_out$gene <- rownames(data_out)
data_out2 <- data_out %>% select(gene, everything())
setwd(IN_FOLDER)
write.table(data_out2, 
            file = "db_pend_pop_ALL.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

data_out <- db_pendplaty_pop_REDUX_store
data_out <- data_out[,-c(ncol(data_out))]
data_out$gene <- rownames(data_out)
data_out2 <- data_out %>% select(gene, everything())
setwd(IN_FOLDER)
write.table(data_out2, 
            file = "db_pendplaty_pop_ALL.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)





#### Import MATLAB results (after k-means analysis)
## primer on k-means clustering analysis with R
# https://statsandr.com/blog/clustering-analysis-k-means-and-hierarchical-clustering-by-hand-and-in-r/

setwd(IN_FOLDER)
db_nana_kmeans <- read.csv("presence_matrix_after_kmeans_nana_ALLgenes.txt", row.names="gene")    
db_hum_kmeans <- read.csv("presence_matrix_after_kmeans_hum_ALLgenes.txt", row.names="gene")      


#Confirm Number of elements per cluster 
table(db_nana_kmeans$clustering_new)
table(db_hum_kmeans$clustering_new)


# reorder factor order, if needed
db_nana_kmeans$clustering_new <- as.factor(db_nana_kmeans$clustering_new)
db_hum_kmeans$clustering_new <- as.factor(db_hum_kmeans$clustering_new)

my_oder <- c("1","2","3","4","5","6","7")

# nana
ge_aux = db_nana_kmeans[order(factor(db_nana_kmeans$clustering_new,levels=c(my_oder ))),]
db_nana_kmeans_matrix <- as.matrix(ge_aux[,1:7])

nana_topGo_df <- db_nana_kmeans  # to be used when saving gene lists for topGO analysis
nana_topGo_df$cluster <- NA

blockA = db_nana_kmeans[db_nana_kmeans$clustering_new == "3",]
blockA = blockA[order(-blockA$Arctic),]
nana_topGo_df[nana_topGo_df$clustering_new == "3", "cluster"] <- "N1"

blockB = db_nana_kmeans[db_nana_kmeans$clustering_new == "4",]
blockB = blockB[order(-blockB$Arctic),]
nana_topGo_df[nana_topGo_df$clustering_new == "4", "cluster"] <- "N2"

blockC = db_nana_kmeans[db_nana_kmeans$clustering_new == "1",]
blockC = blockC[order(-blockC$Arctic),]
nana_topGo_df[nana_topGo_df$clustering_new == "1", "cluster"] <- "N3"

blockD = db_nana_kmeans[db_nana_kmeans$clustering_new == "2",]
blockD = blockD[order(-blockD$Arctic),]
nana_topGo_df[nana_topGo_df$clustering_new == "2", "cluster"] <- "N4"

blockE = db_nana_kmeans[db_nana_kmeans$clustering_new == "6",]
blockE = blockE[order(-blockE$Arctic),]
nana_topGo_df[nana_topGo_df$clustering_new == "6", "cluster"] <- "N5"

blockF = db_nana_kmeans[db_nana_kmeans$clustering_new == "5",]
blockF = blockF[order(-blockF$Arctic),]
nana_topGo_df[nana_topGo_df$clustering_new == "5", "cluster"] <- "N6"

db_nana_kmeans_sort <- rbind(blockA, blockB, blockC, blockD, blockE, blockF)
db_nana_kmeans_matrix <- as.matrix(db_nana_kmeans_sort [,1:7])


# hum
ge_aux = db_hum_kmeans[order(factor(db_hum_kmeans$clustering_new,levels=c(my_oder ))),]
db_hum_kmeans_matrix <- as.matrix(ge_aux[,1:7])

hum_topGo_df <- db_hum_kmeans  # to be used when saving gene lists for topGO analysis
hum_topGo_df$cluster <- NA

blockA2 = db_hum_kmeans[db_hum_kmeans$clustering_new == "6",]
blockA2 = blockA2[order(-blockA2$Central_Asia),]
hum_topGo_df[hum_topGo_df$clustering_new == "6", "cluster"] <- "H1"

blockB2 = db_hum_kmeans[db_hum_kmeans$clustering_new == "3",]
blockB2 = blockB2[order(-blockB2$Central_Asia),]
hum_topGo_df[hum_topGo_df$clustering_new == "3", "cluster"] <- "H2"

blockC2 = db_hum_kmeans[db_hum_kmeans$clustering_new == "5",]
blockC2 = blockC2[order(-blockC2$Central_Asia),]
hum_topGo_df[hum_topGo_df$clustering_new == "5", "cluster"] <- "H3"

blockD2 = db_hum_kmeans[db_hum_kmeans$clustering_new == "2",]
blockD2 = blockD2[order(-blockD2$Central_Asia),]
hum_topGo_df[hum_topGo_df$clustering_new == "2", "cluster"] <- "H4"

blockE2 = db_hum_kmeans[db_hum_kmeans$clustering_new == "4",]
blockE2 = blockE2[order(-blockE2$Central_Asia),]
hum_topGo_df[hum_topGo_df$clustering_new == "4", "cluster"] <- "H5"

blockF2 = db_hum_kmeans[db_hum_kmeans$clustering_new == "1",]
blockF2 = blockF2[order(-blockF2$Central_Asia),]
hum_topGo_df[hum_topGo_df$clustering_new == "1", "cluster"] <- "H6"

db_hum_kmeans_sort <- rbind(blockA2, blockB2, blockC2, blockD2, blockE2, blockF2)
db_hum_kmeans_matrix <- as.matrix(db_hum_kmeans_sort [,1:7])




########### heatmaps

### font size
fontsize = 12            # font size for word 'Treatment' (top-right corner)
fontsize_col = 12        # Treatment name (T0, L1, ...) font size

# heatmap color range
my_colors <- c("mediumblue", "yellow", "red")
color_3 <- colorRampPalette(my_colors, space = "Lab")
color_range <- color_3(100)



### nana

# color list (treatments)        <<<<< ADJUST FUNCTION NAMES AND COLORS AS REQUIRED
my_colour = list(
  Cluster= c(N1 = "deepskyblue", N2 = "violet", N3 = "mistyrose1", N4 = "springgreen4", N5 = "green", N6 = "lightslategray")
)

## cell width and height
mycellwidth = 15
mycellheight = 1.5

# cluster ID
my_cluster_ID <- c(rep("N1", nrow(blockA)),rep("N2", nrow(blockB)),rep("N3", nrow(blockC)),
                   rep("N4", nrow(blockD)),rep("N5", nrow(blockE)),rep("N6", nrow(blockF)))
my_sample_ID_df <- data.frame(Cluster = my_cluster_ID)
row.names(my_sample_ID_df) <- row.names(db_nana_kmeans_sort)
colnames(db_nana_kmeans_matrix) <- c("ES","C. Europe","C. Asia","SV south","SV center","SV north","Arctic")

pushViewport(viewport(width = 0.9, height = 0.9))
my_heatmap <- pheatmap(db_nana_kmeans_matrix, 
                       col = color_range,
                       cluster_cols=FALSE, 
                       cluster_rows=FALSE,
                       show_rownames=FALSE, 
                       annotation_colors = my_colour,
                       annotation_row = my_sample_ID_df,
                       fontsize = fontsize, 
                       fontsize_col = fontsize_col, 
                       treeheight_row = 0, 
                       treeheight_col = 0,
                       cellwidth = mycellwidth, 
                       cellheight = mycellheight,
                       angle_col = "90",
                       legend = TRUE,
                       annotation_legend = TRUE,
                       clustering_method = "average",
                       main = substitute(paste("Allele origin:", italic(' B. nana'), sep=" ")))






### humilis

# color list (treatments)        <<<<< ADJUST FUNCTION NAMES AND COLORS AS REQUIRED
my_colour2 = list(
  Cluster= c(H1 = "deepskyblue", H2 = "violet", H3 = "mistyrose1", H4 = "springgreen4", H5 = "green", H6 = "lightslategray")
)

## cell width and height
mycellwidth = 15
mycellheight = 1.5

# cluster ID
my_cluster_ID2 <- c(rep("H1", nrow(blockA2)),rep("H2", nrow(blockB2)),rep("H3", nrow(blockC2)),
                    rep("H4", nrow(blockD2)),rep("H5", nrow(blockE2)),rep("H6", nrow(blockF2)))
my_sample_ID_df2 <- data.frame(Cluster = my_cluster_ID2)
row.names(my_sample_ID_df2) <- row.names(db_hum_kmeans_sort)
colnames(db_hum_kmeans_matrix) <- c("ES","C. Europe","C. Asia","SV south","SV center","SV north","Arctic")

pushViewport(viewport(width = 0.9, height = 0.9))
my_heatmap2 <- pheatmap(db_hum_kmeans_matrix, 
                       col = color_range,
                       cluster_cols=FALSE, 
                       cluster_rows=FALSE,
                       show_rownames=FALSE, 
                       annotation_colors = my_colour2,
                       annotation_row = my_sample_ID_df2,
                       fontsize = fontsize, 
                       fontsize_col = fontsize_col, 
                       treeheight_row = 0, 
                       treeheight_col = 0,
                       cellwidth = mycellwidth, 
                       cellheight = mycellheight,
                       angle_col = "90",
                       legend = TRUE,
                       annotation_legend = TRUE,
                       clustering_method = "average",
                       main = substitute(paste("Allele origin:", italic(' B. humilis'), sep=" ")))



plot_grid(my_heatmap[[4]], my_heatmap2[[4]], nrow = 1, align = 'v')




