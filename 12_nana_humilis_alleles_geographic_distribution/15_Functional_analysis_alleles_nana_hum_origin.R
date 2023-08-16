## Functional analysis of B. pubescens alleles of B. nana or B. humilis origin
## Luis Leal, 2022


#Libraries
library(plotly)
library(gplots)
library(ggplot2)
library(cowplot)
library(colorspace)
library(gridExtra)
library(ggpubr)
library(grid)
library(Metrics)
library(stringr)
library(scales)
library(matrixStats)


#Verify versions of R and other loaded software:
sessionInfo()

########## Clear all states and remove all plots
#if(!is.null(dev.list())) dev.off()
rm(list=ls(all=TRUE))

###### Input folder 
IN_FOLDER <- "/Users/luisleal/Dropbox/PlantEco/03_Paper_Origin_Bpubescens/04_PHYLOGENIES"     # iMac

###### in files
infile_nana <- "presence_matrix_after_kmeans_nana_annotated_0.65_NEW2.txt"   # allele observed in >=65% of all individuals in at least one population; k-means annotated
infile_hum <- "presence_matrix_after_kmeans_hum_annotated_0.65_NEW2.txt"

infile_nana_ALL <- "db_nana_pop_ALL.txt"     # list of alleles observed in at least one individual 
infile_hum_ALL <- "db_hum_pop_ALL.txt"
infile_pend_ALL <- "db_pend_pop_ALL.txt"

annot_file <- "presence_matrix_after_kmeans_annotated_genes.txt"   # annotated birch genes (for gene set of interest only)


#### load data

data_nana <- read.delim(file.path(IN_FOLDER, infile_nana), header = TRUE, sep = "\t")
data_hum <- read.delim(file.path(IN_FOLDER, infile_hum), header = TRUE, sep = "\t")

data_nana_ALL <- read.delim(file.path(IN_FOLDER, infile_nana_ALL), header = TRUE, sep = "\t")
data_hum_ALL <- read.delim(file.path(IN_FOLDER, infile_hum_ALL), header = TRUE, sep = "\t")
data_pend_ALL <- read.delim(file.path(IN_FOLDER, infile_pend_ALL), header = TRUE, sep = "\t")

anott_data <- read.delim(file.path(IN_FOLDER, annot_file), header = TRUE, sep = "\t")



### Clean tables
anott_data <-anott_data[,-c(2:13,15:18,49)]

data_nana_clean <- data_nana[, -c(9)]  
data_nana_clean <- merge(data_nana_clean,anott_data,by="gene",all.x = TRUE)

data_hum_clean <- data_hum[, -c(9)]  
data_hum_clean <- merge(data_hum_clean,anott_data,by="gene",all.x = TRUE)



### Flag genes without GO info
data_nana_clean$sum <- rowSums(data_nana_clean[ ,c(11:ncol(data_nana_clean))], na.rm = TRUE)
data_hum_clean$sum <- rowSums(data_hum_clean[ ,c(11:ncol(data_hum_clean))], na.rm = TRUE)

#data_nana_clean <- data_nana_clean[data_nana_clean$sum > 0,]
#data_hum_clean <- data_hum_clean[data_hum_clean$sum > 0,]



# recode GO terms (remove GO codes)
colnames(data_nana_clean) <- gsub("..GO.[0-9]*.", "", colnames(data_nana_clean))
colnames(data_hum_clean) <- gsub("..GO.[0-9]*.", "", colnames(data_hum_clean))



###################### Auxiliary funtions

dataProcessFun <- function(data_nana_N1p, mycolname1, mycolname3)
{
  aux_1 <- dim(data_nana_N1p)[1] * dim(data_nana_N1p)[2]
  data_plot_1 <- data.frame(X = rep(NA,aux_1), Gene = rep(NA,aux_1), Y = rep(NA,aux_1))
  names(data_plot_1)[names(data_plot_1)=="X"] <- mycolname1
  names(data_plot_1)[names(data_plot_1)=="Y"] <- mycolname3
  
  n <- 1
  for (i in colnames(data_nana_N1p)) {
    for (j in rownames(data_nana_N1p)) {
      data_plot_1[n,1] <- i
      data_plot_1[n,2] <- j
      data_plot_1[n,3] <- data_nana_N1p[j,i]
      n <- n + 1
    }
  }
  
  return(data_plot_1)
}


FractionPlot <- function(data_plot_1)
{
  level_order <- c("ES","C. Europe","C. Asia","SV south","SV center","SV north","Arctic")
  
  plotN1 <-ggplot(data_plot_1 ,
                  aes(x = factor(Location, level = level_order), y = Fraction, color = Gene, group = Gene)) + 
    xlab("Location") +
    ylab("Fraction of individuals with introgressed allele") +
    theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
    theme(panel.background = element_rect(fill = 'white', colour = 'black')) +
    geom_line() +
    scale_y_continuous(limits = c(0, 1))
}


dataGo <- function(data_nana_N1, data_nana_N1p)
{
  data_nana_N1go <- data_nana_N1[data_nana_N1$gene_name_1 %in% rownames(data_nana_N1p),]
  row.names(data_nana_N1go) <- data_nana_N1go$gene_name_1
  data_nana_N1go_p <- data_nana_N1go[, c(11:ncol(data_nana_N1go))]
  data_nana_N1go_p <- data_nana_N1go_p[,-c(ncol(data_nana_N1go_p))]
  data_nana_N1go_p[is.na(data_nana_N1go_p)] <- 0
  
  # clean column names
  colnames(data_nana_N1go_p) <- gsub("\\.\\.", " (", colnames(data_nana_N1go_p))
  colnames(data_nana_N1go_p) <- gsub("\\.$", ")", colnames(data_nana_N1go_p))
  colnames(data_nana_N1go_p) <- gsub("\\.", " ", colnames(data_nana_N1go_p))
  colnames(data_nana_N1go_p) <- gsub("\\(GO ", "(GO\\:", colnames(data_nana_N1go_p))
  
  return(data_nana_N1go_p)
}


ColorGo <- function(data_plot_Ngo)
{
  data_plot_Ngo$mycolor <- "black"
  data_plot_Ngo$mycolor[grepl("gametophyte development", data_plot_Ngo$GO_term)] <- "chartreuse"
  data_plot_Ngo$mycolor[grepl("reproductive system development", data_plot_Ngo$GO_term)] <- "chartreuse"
  data_plot_Ngo$mycolor[grepl("recognition of pollen", data_plot_Ngo$GO_term)] <- "chartreuse"
  data_plot_Ngo$mycolor[grepl("seed germination", data_plot_Ngo$GO_term)] <- "chartreuse4"
  data_plot_Ngo$mycolor[grepl("seed development", data_plot_Ngo$GO_term)] <- "chartreuse"
  data_plot_Ngo$mycolor[grepl("root system development", data_plot_Ngo$GO_term)] <- "chartreuse4"
  data_plot_Ngo$mycolor[grepl("shoot system development", data_plot_Ngo$GO_term)] <- "chartreuse4"
  data_plot_Ngo$mycolor[grepl("response to cold", data_plot_Ngo$GO_term)] <- "orange"
  data_plot_Ngo$mycolor[grepl("response to heat", data_plot_Ngo$GO_term)] <- "orange"
  data_plot_Ngo$mycolor[grepl("response to water deprivation", data_plot_Ngo$GO_term)] <- "orange"
  data_plot_Ngo$mycolor[grepl("response to salt stress", data_plot_Ngo$GO_term)] <- "orange"
  data_plot_Ngo$mycolor[grepl("response to light stimulus", data_plot_Ngo$GO_term)] <- "violet"
  data_plot_Ngo$mycolor[grepl("circadian rhythm", data_plot_Ngo$GO_term)] <- "violet"
  data_plot_Ngo$mycolor[grepl("response to fungus", data_plot_Ngo$GO_term)] <- "red"
  data_plot_Ngo$mycolor[grepl("response to bacterium", data_plot_Ngo$GO_term)] <- "red"
  data_plot_Ngo$mycolor[grepl("response to virus", data_plot_Ngo$GO_term)] <- "red"
  data_plot_Ngo$mycolor[grepl("defense response to insect", data_plot_Ngo$GO_term)] <- "red"
  data_plot_Ngo$mycolor[grepl("mitotic DNA replication checkpoint signaling", data_plot_Ngo$GO_term)] <- "deepskyblue"
  data_plot_Ngo$mycolor[grepl("DNA repair", data_plot_Ngo$GO_term)] <- "deepskyblue"
  data_plot_Ngo$mycolor[grepl("DNA recombination", data_plot_Ngo$GO_term)] <- "deepskyblue"
  data_plot_Ngo$mycolor[grepl("chromosome organization", data_plot_Ngo$GO_term)] <- "deepskyblue"
  data_plot_Ngo$mycolor[grepl("leaf senescence", data_plot_Ngo$GO_term)] <- "seagreen1"
  data_plot_Ngo$mycolor[grepl("abscisic acid pathways", data_plot_Ngo$GO_term)] <- "lightcyan3"
  data_plot_Ngo$mycolor[grepl("auxin pathways", data_plot_Ngo$GO_term)] <- "lightcyan3"
  data_plot_Ngo$mycolor[grepl("jasmonic acid pathways", data_plot_Ngo$GO_term)] <- "lightcyan3"
  data_plot_Ngo$mycolor[grepl("ethylene pathways", data_plot_Ngo$GO_term)] <- "lightcyan3"
  data_plot_Ngo$mycolor[grepl("gibberellic acid pathways", data_plot_Ngo$GO_term)] <- "lightcyan3"
  data_plot_Ngo$mycolor[grepl("regulation of gene expression", data_plot_Ngo$GO_term)] <- "lightpink"
  data_plot_Ngo$mycolor[grepl("histone acetylation", data_plot_Ngo$GO_term)] <- "lightpink"
  data_plot_Ngo$mycolor[grepl("histone deacetylation", data_plot_Ngo$GO_term)] <- "lightpink"
  data_plot_Ngo$mycolor[grepl("histone methylation", data_plot_Ngo$GO_term)] <- "lightpink"
  
  return(data_plot_Ngo)
}


GOplotFUN <- function(data_plot_Ngo, level_order_go, level_order_gene, myPlotMargins, plot_title)
{
  plotN1go <- ggplot(data_plot_Ngo, aes(x=factor(GO_term, level = level_order_go), y=factor(Gene, level = level_order_gene), label=Gene)) + 
    geom_point(stat='identity', aes(col=Flag), size=data_plot_Ngo$mydotsize, colour = data_plot_Ngo$mycolor)  +
    ggtitle(plot_title) +
    theme(axis.title.x = element_blank()) +
    theme(axis.title.y = element_blank()) +
    theme(axis.text.y = element_text(size=10, angle=0)) +
    theme(axis.text.x = element_text(size=10, angle=90)) +
    theme(panel.background = element_rect(fill = 'white', colour = 'black')) +
    theme(panel.grid.major = element_line(size = 0.3, linetype = 'solid', colour = "bisque")) +
    theme(axis.text.x = element_text(vjust = 0.5, hjust = 0.95)) +
    theme(plot.margin=unit(myPlotMargins, 'cm')) +
    labs(title = plot_title) +
    theme(plot.title = element_text(vjust = 0.5, hjust = 0.01, size = 16, face = "plain")) +
    coord_flip()

   return(plotN1go)
}


plotDots <- function(data_plot_N1_dots, level_order, level_order_gene, myPlotMargins, myDotSize)
{
  
  #my_colors <- c("mediumblue", "yellow", "red")
  my_colors <- c("red", "gold", "green3")
  color_3 <- colorRampPalette(my_colors, space = "Lab")
  color_range <- color_3(100)
  
  data_plot_N1_dots$mycolors <- "black" 
  data_plot_N1_dots$aux <- round((100/myDotSize)*data_plot_N1_dots$mydotsize)
  data_plot_N1_dots$aux[data_plot_N1_dots$aux == 0] <- 1
  data_plot_N1_dots$mycolors <- color_range[data_plot_N1_dots$aux]
  
  plotN1_dots <- ggplot(data_plot_N1_dots, aes(x=factor(Location, level = level_order), y=factor(Gene, level = level_order_gene), label=Gene)) + 
    geom_point(stat='identity', aes(col=Flag), size=data_plot_N1_dots$mydotsize, colour = data_plot_N1_dots$mycolors)  +
    theme(axis.title.y = element_text(size = rel(1.8), angle = 90)) +
    theme(axis.title.x = element_blank()) +
    theme(axis.title.y = element_blank()) +
    theme(axis.ticks.x = element_blank()) +
    theme(axis.text.y = element_text(size=10, angle=0)) +
    theme(axis.text.x = element_text(size=6, angle=40)) +
    theme(panel.background = element_rect(fill = 'white', colour = 'white')) +
    theme(panel.grid.major = element_line(size = 0.3, linetype = 'solid', colour = "bisque")) +
    theme(axis.text.x = element_text(vjust = 0.95, hjust = 0.95)) +
    theme(plot.margin=unit(myPlotMargins, 'cm')) +
    coord_flip()
  
  return(plotN1_dots)
}






#################### 1. All genes with an allele of nana origin and containing a GO annotation

data_nana_all <- data_nana_clean[data_nana_clean$sum > 0,]

# select genes for which the nana allele has been detected in at least 65% of all individuals in a particular population
data_nana_all$maxf <- rowMaxs(as.matrix(data_nana_all[,c(2:8)]))
data_nana_all <- data_nana_all[data_nana_all$maxf >= 0.65, ]
data_nana_all <- data_nana_all[,-c(ncol(data_nana_all))]

## covert dataset to ggplot-friendly format
data_nana_allp <- data_nana_all
row.names(data_nana_allp) <- data_nana_allp$gene_name_1
data_nana_allp <- data_nana_allp[, c(2:8)]

colnames(data_nana_allp) <- c("ES","C. Europe","C. Asia","SV south","SV center","SV north","Arctic")
data_plot_all_nana <- dataProcessFun(data_nana_allp,"Location","Fraction")

# get dataset with GO terms
data_nana_allgo_p <- dataGo(data_nana_all, data_nana_allp)

# covert dataset to ggplot-friendly format
data_plot_Nallgo <- dataProcessFun(data_nana_allgo_p, "GO_term", "Flag")

# color GO terms
data_plot_Nallgo <- ColorGo(data_plot_Nallgo)

# Select only observed GO/gene pairs 
data_plot_Nallgo <- data_plot_Nallgo[data_plot_Nallgo$Flag == "1",]

# dot size
myDotSize <- 4
data_plot_Nallgo$mydotsize <- myDotSize*(data_plot_Nallgo$Flag) 

# plot margins
myPlotMargins <- c(4.5,1,1,1)

# plot title
plot_title <- expression(paste("Allele origin:", italic(' nana'), sep=" "))

level_order_go <- c("gametophyte development","reproductive system development","recognition of pollen", "seed development",
                    "seed germination", "root system development", "shoot system development", 
                    "response to light stimulus", "circadian rhythm",
                    "leaf senescence",
                    "response to water deprivation","response to cold", "response to heat","response to salt stress",
                    "response to fungus","response to bacterium","response to virus","defense response to insect",
                    "mitotic DNA replication checkpoint signaling","DNA repair","DNA recombination","chromosome organization",
                    "regulation of gene expression","histone acetylation","histone deacetylation","histone methylation",
                    "abscisic acid pathways","auxin pathways","jasmonic acid pathways","ethylene pathways","gibberellic acid pathways")
level_order_go <- rev(level_order_go)

# gene order
level_order_gene <- c("CYP86A8","PWR","MIP2","ERL1","RST1","FAB1B","GRV2","APRR7","MAX2","RGI2","TPR2","ARR12","CAMTA2","DIN10","At3g21300 ","GWD1","DCL4","NHX7","RCF3",
                      "4CLL7","KCS11 ","PGDH1","SIR","RF298","GFS12",
                      "EXO70-G1","At4g28260",
                      "DRM1","GH3.6","PCMP-E47","MCM4",
                      "SPT6","HULK3","VIP6","CPL1","PAT1","GLR3.3 ","PLDDELTA",
                      "LD","VIP1","PHYB","TMK3","SNL2","OBE3","ARF5","EAF1A","LHW","GLU1","ALA4",
                      "ERECTA","CDC5","CuAOzeta","EDR1","SHD","FPA","GLR3.4","AFB2","FLS2","ERS1","ACC1","MBR2",
                      "MED12","DDM1","ZEP","CHR12","AHK4","LOX5","RBOHD",
                      "CPL3","NPR1","DCL1","SRFR1","SNL4",
                      "CALS11","PAL1",
                      "TRN1","ABCG37","SPA2","SWEETIE","SIZ1",
                      "CNGC15","MBD9","BRL2",
                      "RBOHB",
                      "PUB13",
                      "POL2A",
                      "LHY","XRN4","PHOT2","BHLH155","SBT6.1","RCA","EIN3","AT4G03260","DCL3","RPK2","VIL2","AGO7","MDN1")

# plot
(plotNallgo <- GOplotFUN(data_plot_Nallgo, level_order_go, level_order_gene, myPlotMargins, plot_title))

## dot plots data
data_plot_Nall_dots <- data_plot_all_nana

# GO and dots plots
myPlotMargins <- c(0.5,0.5,0,0.1)
plot_title <- ""
plotNallgo <- GOplotFUN(data_plot_Nallgo, level_order_go, level_order_gene, myPlotMargins, plot_title)

# Dots
myPlotMargins <- c(0.2,0.5,7.5,0.1)

myDotSize <- 4   ## dot size
data_plot_Nall_dots$mydotsize <- myDotSize*(data_plot_Nall_dots$Fraction) 
level_order <- c("ES","C. Europe","C. Asia","SV south","SV center","SV north","Arctic")

data_plot_Nall_dots <- data_plot_Nall_dots[data_plot_Nall_dots$Gene %in% data_plot_Nallgo$Gene,]
plotNall_dots <- plotDots(data_plot_Nall_dots, level_order, level_order_gene, myPlotMargins, myDotSize)

plot_grid(plotNallgo, plotNall_dots, ncol = 1, align = 'v')






#################### 2. All genes with an allele of humilis origin and containing a GO annotation

data_hum_all <- data_hum_clean[data_hum_clean$sum > 0,]

# select genes for which the humilis allele has been detected in at least 65% of all individuals in a particular population
data_hum_all$maxf <- rowMaxs(as.matrix(data_hum_all[,c(2:8)]))
data_hum_all <- data_hum_all[data_hum_all$maxf >= 0.65, ]
data_hum_all <- data_hum_all[,-c(ncol(data_hum_all))]

## covert dataset to ggplot-friendly format
data_hum_allp <- data_hum_all
row.names(data_hum_allp) <- data_hum_allp$gene_name_1
data_hum_allp <- data_hum_allp[, c(2:8)]

colnames(data_hum_allp) <- c("ES","C. Europe","C. Asia","SV south","SV center","SV north","Arctic")
data_plot_all_hum <- dataProcessFun(data_hum_allp,"Location","Fraction")



###### GO analysis plot (for same gene set in plotN1)

# get dataset with GO terms
data_hum_allgo_p <- dataGo(data_hum_all, data_hum_allp)

# covert dataset to ggplot-friendly format
data_plot_Hallgo <- dataProcessFun(data_hum_allgo_p, "GO_term", "Flag")

# color GO terms
data_plot_Hallgo <- ColorGo(data_plot_Hallgo)

# Select only observed GO/gene pairs 
data_plot_Hallgo <- data_plot_Hallgo[data_plot_Hallgo$Flag == "1",]



# dot size
myDotSize <- 4
data_plot_Hallgo$mydotsize <- myDotSize*(data_plot_Hallgo$Flag) 

# plot margins
myPlotMargins <- c(4.5,1,1,1)

# plot title
plot_title <- expression(paste("Allele origin:", italic(' humilis'), sep=" "))

# Go term order
level_order_go <- c("gametophyte development","reproductive system development","recognition of pollen", "seed development",
                    "seed germination", "root system development", "shoot system development", 
                    "response to light stimulus", "circadian rhythm",
                    "leaf senescence",
                    "response to water deprivation","response to cold", "response to heat","response to salt stress",
                    "response to fungus","response to bacterium","response to virus","defense response to insect",
                    "mitotic DNA replication checkpoint signaling","DNA repair","DNA recombination","chromosome organization",
                    "regulation of gene expression","histone acetylation","histone deacetylation","histone methylation",
                    "abscisic acid pathways","auxin pathways","jasmonic acid pathways","ethylene pathways","gibberellic acid pathways")
level_order_go <- rev(level_order_go)

# gene order
level_order_gene <- c("XRN4","CALS12","LHY","RAD17","RECQL4A","UVH3","DOT2","DRD1","MNB8.8","AP22.19","At1g74770","RDR6",
                      "DME","RBOHD","TPST","PLDDELTA","SAB","ERL1","LOX5","DTX48","AVP1","CMT2","AGO7","NMAT4","DEK1","SUVR5","JAR1","BHLH155",
                      "HSR4 ","BAM3","BRM","SPY","SIZ1","PIF3",
                      "CHR12","GN","EIN2","CLPB4","AAO2","HOS1","GLR3.3","ATXR7","SNL4",
                      "CIPK23","PLDBETA1","LOX3", "BAG6","CAMTA3","BIR2",
                      "ARR12","MTB","NPR1","AHK2",
                      "CRY1","MKP1","NPF6.3",
                      "IRE1B","PWR","CNGC15")


# plot
(plotHallgo <- GOplotFUN(data_plot_Hallgo, level_order_go, level_order_gene, myPlotMargins, plot_title))

## dot plots data
data_plot_Hall_dots <- data_plot_all_hum

# GO and dots plots
myPlotMargins <- c(0.5,0.5,0,0.1)
plot_title <- ""
plotHallgo <- GOplotFUN(data_plot_Hallgo, level_order_go, level_order_gene, myPlotMargins, plot_title)

# Dots
myPlotMargins <- c(0.2,0.5,7.5,0.1)

myDotSize <- 4   ## dot size
data_plot_Hall_dots$mydotsize <- myDotSize*(data_plot_Hall_dots$Fraction) 
level_order <- c("ES","C. Europe","C. Asia","SV south","SV center","SV north","Arctic")

data_plot_Hall_dots <- data_plot_Hall_dots[data_plot_Hall_dots$Gene %in% data_plot_Hallgo$Gene,]
plotHall_dots <- plotDots(data_plot_Hall_dots, level_order, level_order_gene, myPlotMargins, myDotSize)

plot_grid(plotHallgo, plotHall_dots, ncol = 1, align = 'v')







############## Scatter plots for individual genes


AlleleFractionPlot <- function(data_plot_1, gname, legend_position)
{
  level_order <- c("ES", "C. Europe", "C. Asia", "SV south", "SV center", "SV north", "Arctic" )
  
  plotN1 <-ggplot(data_plot_1 ,
                  aes(x = factor(Location, level = level_order), y = Fraction, color = Origin, group = Origin)) + 
    geom_point(aes(color = Origin, fill = Origin), size=5.5, alpha = 0.5, shape = 21) + 
    xlab("") +
    ylab("Fraction of individuals \n with introgressed allele") +
    scale_colour_manual(values = darken(cbbPalette, 0.3)) +
    scale_fill_manual(values = cbbPalette) +
    theme(axis.text=element_text(size=12), axis.title=element_text(size=15,face="plain")) +
    theme(axis.text.x = element_text(angle = 60, hjust=1, vjust=1))  +
    theme(panel.background = element_rect(fill = "white", colour = "white", size = 0.75)) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    theme(axis.line.x.bottom=element_line(color="black")) +
    theme(axis.line.y.left=element_line(color="black")) +
    geom_line() +
    scale_y_continuous(limits = c(0, 1)) +
    labs(title = gname) +
    theme(plot.title = element_text(vjust = 0.5, hjust = 0.01, size = 16, face = "plain")) +
    theme(legend.text=element_text(size=10),
                  legend.position = legend_position,
                  legend.direction="vertical",
                  legend.background = element_rect(fill = "white"),
                  legend.key = element_rect(fill = "white"))
}

cbbPalette <- c("bisque","violet","grey87","green2","darkorange","chartreuse4","turquoise3")





####### PWR (PWR) (POWERDRESS)
gname <- "PWR"
gname_Bp <- "Bpev01.c0043.g0028"

allele_freq_nana <- data_nana_ALL[data_nana_ALL$gene == gname_Bp, c(2:8)]
rownames(allele_freq_nana) <- "nana"

allele_freq_hum <- data_hum_ALL[data_hum_ALL$gene == gname_Bp, c(2:8)]
rownames(allele_freq_hum) <- "humilis"

allele_freq_pend <- data_pend_ALL[data_pend_ALL$gene == gname_Bp, c(2:8)]
rownames(allele_freq_pend) <- "pendula"

#allele_freq_data <- rbind(allele_freq_nana,allele_freq_hum,allele_freq_pend)
allele_freq_data <- rbind(allele_freq_nana,allele_freq_hum)
colnames(allele_freq_data) <- c("ES", "C. Europe", "C. Asia", "SV south", "SV center", "SV north", "Arctic")

## covert dataset to ggplot-friendly format
data_plot_G1 <- dataProcessFun(allele_freq_data,"Location","Fraction")
colnames(data_plot_G1)[2] <- "Origin"

# plot
legend_position = c(0.7, 0.2)
(plotG1 <- AlleleFractionPlot(data_plot_G1, gname, legend_position))





####### CALS12 (CALLOSE SYNTHASE 12)
gname <- "CALS12"
gname_Bp <- "Bpev01.c0845.g0002"

allele_freq_nana <- data_nana_ALL[data_nana_ALL$gene == gname_Bp, c(2:8)]
allele_freq_nana[1,] <- c(0,0,0,0,0,0,0)   # allele was not detected in any of the samples
rownames(allele_freq_nana) <- "nana"

allele_freq_hum <- data_hum_ALL[data_hum_ALL$gene == gname_Bp, c(2:8)]
rownames(allele_freq_hum) <- "humilis"

allele_freq_pend <- data_pend_ALL[data_pend_ALL$gene == gname_Bp, c(2:8)]
allele_freq_pend[1,] <- c(0,0,0,0,0,0,0)   # allele was not detected in any of the samples
rownames(allele_freq_pend) <- "pendula"

#allele_freq_data <- rbind(allele_freq_nana,allele_freq_hum,allele_freq_pend)
allele_freq_data <- rbind(allele_freq_nana,allele_freq_hum)
colnames(allele_freq_data) <- c("ES", "C. Europe", "C. Asia", "SV south", "SV center", "SV north", "Arctic")

## covert dataset to ggplot-friendly format
data_plot_G2 <- dataProcessFun(allele_freq_data,"Location","Fraction")
colnames(data_plot_G2)[2] <- "Origin"

# plot
legend_position = c(0.25, 0.5)
(plotG2 <- AlleleFractionPlot(data_plot_G2, gname, legend_position))





####### PHYB (Phytochrome B)
gname <- "PHYB"
gname_Bp <- "Bpev01.c2035.g0006"

allele_freq_nana <- data_nana_clean[data_nana_clean$gene_name_1 == gname, c(2:8)]
rownames(allele_freq_nana) <- "nana"

allele_freq_hum <- data_hum_ALL[data_hum_ALL$gene == gname_Bp, c(2:8)]
rownames(allele_freq_hum) <- "humilis"

allele_freq_pend <- data_pend_ALL[data_pend_ALL$gene == gname_Bp, c(2:8)]
rownames(allele_freq_pend) <- "pendula"

#allele_freq_data <- rbind(allele_freq_nana,allele_freq_hum,allele_freq_pend)
allele_freq_data <- rbind(allele_freq_nana,allele_freq_hum)
colnames(allele_freq_data) <- c("ES", "C. Europe", "C. Asia", "SV south", "SV center", "SV north", "Arctic")

## covert dataset to ggplot-friendly format
data_plot_G3 <- dataProcessFun(allele_freq_data,"Location","Fraction")
colnames(data_plot_G3)[2] <- "Origin"

# plot
legend_position = c(0.25, 0.8)
(plotG3 <- AlleleFractionPlot(data_plot_G3, gname, legend_position))




####### NPR1 (Regulatory protein NPR1)
gname <- "NPR1"
gname_Bp <- "Bpev01.c0051.g0060"

allele_freq_nana <- data_nana_ALL[data_nana_ALL$gene == gname_Bp, c(2:8)]
rownames(allele_freq_nana) <- "nana"

allele_freq_hum <- data_hum_clean[data_hum_clean$gene_name_1 == gname, c(2:8)]
rownames(allele_freq_hum) <- "humilis"

allele_freq_pend <- data_pend_ALL[data_pend_ALL$gene == gname_Bp, c(2:8)]
rownames(allele_freq_pend) <- "pendula"

#allele_freq_data <- rbind(allele_freq_nana,allele_freq_hum,allele_freq_pend)
allele_freq_data <- rbind(allele_freq_nana,allele_freq_hum)
colnames(allele_freq_data) <- c("ES", "C. Europe", "C. Asia", "SV south", "SV center", "SV north", "Arctic")

## covert dataset to ggplot-friendly format
data_plot_G4 <- dataProcessFun(allele_freq_data,"Location","Fraction")
colnames(data_plot_G4)[2] <- "Origin"

# plot
legend_position = c(0.7, 0.8)
(plotG4 <- AlleleFractionPlot(data_plot_G4, gname, legend_position))




(plot_grid_set1 <- plot_grid(plotG1, plotG2, plotG3, plotG4, nrow = 1, ncol = 4))





