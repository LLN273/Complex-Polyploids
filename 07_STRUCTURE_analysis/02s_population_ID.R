# Associate sample to population
# Luis Leal, 2021


##### Clear all states
if(!is.null(dev.list())) dev.off()
par(mfrow=c(1,1))
rm(list=ls(all=TRUE))

library("stringr")

#Verify versions of R and other loaded software:
#sessionInfo()

## read arguments
args <- commandArgs(TRUE)

## Ploidy for each sample
Ploidy_file <- args[1]
SPLOIDY <- read.table(Ploidy_file, header = F, sep =" ", na.strings = "NA")

## Sample list
SampleList_file <- args[2]
SampleList <- read.table(SampleList_file, header = F, sep =" ", na.strings = "NA")

# output folder
OUTFolder <- args[3]

# output file, population
OUTFile_pop <- args[4]

# output file, ploidy
OUTFile_ploidy <- args[5]


### Ploidy

SampleList$ploidy <- "0"

SampleList[, 1] <- sapply(SampleList[, 1], as.character)

for(i in 1:nrow(SampleList)) {
  (aux_pl <-which(SPLOIDY$V1==SampleList[i,1]))
  SampleList[i,"ploidy"] <- SPLOIDY[aux_pl,2]   
}

SampleList$pop <- substr(SampleList$V1, 0, 3)
SampleList$pop2 <- substr(SampleList$V1, 0, 8)

SampleList$pop_ploidy <- paste(SampleList$ploidy, SampleList$pop, sep = "")
SampleList$pop_ploidy2 <- paste(SampleList$ploidy, SampleList$pop2, sep = "")       # used for FI population

SampleList$pop_order <- "0"


# Population ID

SampleList$pop_order[grepl("Pendula-Germany.Kyynel", SampleList$V1)] <- "1"
SampleList$pop_order[grepl("Pendula-Finland-Loimaa.Loimaa", SampleList$V1)] <- "1"
SampleList$pop_order[grepl("Pendula-Finland-Helsinki.Ritva", SampleList$V1)] <- "1"
SampleList$pop_order[grepl("Pendula-Finland-Kuopio.Kanttarelli", SampleList$V1)] <- "1"
SampleList$pop_order[grepl("Pendula-Finland-Luuta.Luuta", SampleList$V1)] <- "1"
SampleList$pop_order[grepl("Pendula-Finland-Oulu.Poyta", SampleList$V1)] <- "1"
SampleList <- transform(SampleList, pop_order = ifelse(pop_ploidy == "2SPN", "2", pop_order))
SampleList <- transform(SampleList, pop_order = ifelse(pop_ploidy == "2ES-", "2", pop_order))
SampleList <- transform(SampleList, pop_order = ifelse(pop_ploidy2 == "2FR-BP-04", "3", pop_order)) 	# FR, inland
SampleList <- transform(SampleList, pop_order = ifelse(pop_ploidy2 == "2FR-BP-03", "4", pop_order)) 	# FR, Alps
SampleList <- transform(SampleList, pop_order = ifelse(pop_ploidy2 == "2FR-BP-21", "4", pop_order)) 	# FR, Alps
SampleList <- transform(SampleList, pop_order = ifelse(pop_ploidy == "2CH-", "5", pop_order))
SampleList <- transform(SampleList, pop_order = ifelse(pop_ploidy == "2IT-", "6", pop_order))
SampleList <- transform(SampleList, pop_order = ifelse(pop_ploidy == "2GB-", "7", pop_order)) 
SampleList <- transform(SampleList, pop_order = ifelse(pop_ploidy == "2DE-", "8", pop_order)) 
SampleList <- transform(SampleList, pop_order = ifelse(pop_ploidy == "2NO-", "9", pop_order))
SampleList$pop_order[grepl("Pendula-Norway.Drobak", SampleList$V1)] <- "9O"
SampleList <- transform(SampleList, pop_order = ifelse(pop_ploidy == "2DJU", "10", pop_order))
SampleList <- transform(SampleList, pop_order = ifelse(pop_ploidy == "2SKA", "10", pop_order))
SampleList <- transform(SampleList, pop_order = ifelse(pop_ploidy == "2VAL", "10", pop_order))
SampleList <- transform(SampleList, pop_order = ifelse(pop_ploidy == "2ASP", "10", pop_order))
SampleList <- transform(SampleList, pop_order = ifelse(pop_ploidy == "2MOK", "11", pop_order))
SampleList <- transform(SampleList, pop_order = ifelse(pop_ploidy == "2LIN", "11", pop_order))
SampleList <- transform(SampleList, pop_order = ifelse(pop_ploidy == "2GNA", "11", pop_order))
SampleList <- transform(SampleList, pop_order = ifelse(pop_ploidy == "2DOK", "11", pop_order))
SampleList <- transform(SampleList, pop_order = ifelse(pop_ploidy == "2BRA", "12", pop_order))
SampleList <- transform(SampleList, pop_order = ifelse(pop_ploidy == "2BUR", "12", pop_order))
SampleList <- transform(SampleList, pop_order = ifelse(pop_ploidy == "2PIT", "12", pop_order))
SampleList <- transform(SampleList, pop_order = ifelse(pop_ploidy == "2JOK", "12", pop_order))
SampleList <- transform(SampleList, pop_order = ifelse(pop_ploidy == "2SE-", "12", pop_order)) 
SampleList <- transform(SampleList, pop_order = ifelse(pop_ploidy2 == "2FI-BP-19", "13", pop_order)) 
SampleList$pop_order[grepl("Pendula-Finland-Posio", SampleList$V1)] <- "13"
SampleList$pop_order[grepl("Pendula-Finland-Rovaniemi", SampleList$V1)] <- "13"
SampleList$pop_order[grepl("Pendula-Finland-Kittila", SampleList$V1)] <- "13"
SampleList <- transform(SampleList, pop_order = ifelse(pop_ploidy2 == "2FI-BP-20", "14", pop_order)) 
SampleList$pop_order[grepl("Pendula-Finland-Punkaharju", SampleList$V1)] <- "14"
SampleList$pop_order[grepl("Pendula-Finland-Vehmersalmi", SampleList$V1)] <- "14"
SampleList$pop_order[grepl("Pendula-Finland-Loppi", SampleList$V1)] <- "15"
SampleList <- transform(SampleList, pop_order = ifelse(pop_ploidy == "2LT-", "16", pop_order))
SampleList <- transform(SampleList, pop_order = ifelse(pop_ploidy == "2POL", "17", pop_order))
SampleList$pop_order[grepl("2YT-23", SampleList$pop_ploidy2)] <- "18"
SampleList$pop_order[grepl("2YT-24", SampleList$pop_ploidy2)] <- "18"
SampleList$pop_order[grepl("2YT-27", SampleList$pop_ploidy2)] <- "19"
SampleList$pop_order[grepl("2YT-28", SampleList$pop_ploidy2)] <- "19"
SampleList$pop_order[grepl("Pendula-Russia-Voronezh", SampleList$V1)] <- "20"
SampleList <- transform(SampleList, pop_order = ifelse(pop_ploidy == "2URL", "22", pop_order))
SampleList$pop_order[grepl("Pendula-Russia-Krasnoyarsk", SampleList$V1)] <- "21"
SampleList <- transform(SampleList, pop_order = ifelse(pop_ploidy == "2NSV", "22", pop_order))
SampleList$pop_order[grepl("2YT-86", SampleList$pop_ploidy2)] <- "22"
SampleList$pop_order[grepl("2YT-87", SampleList$pop_ploidy2)] <- "22"
SampleList$pop_order[grepl("Pendula-Russia-Yekaterinburg", SampleList$V1)] <- "22"
SampleList$pop_order[grepl("Pendula-Russia-Yakutsk", SampleList$V1)] <- "23"

SampleList <- transform(SampleList, pop_order = ifelse(pop_ploidy == "2OBL", "24", pop_order))
SampleList$pop_order[grepl("Platyphylla-Russia", SampleList$V1)] <- "25"
SampleList <- transform(SampleList, pop_order = ifelse(pop_ploidy == "2MAN", "25", pop_order))

SampleList$pop_order[grepl("Pubescens-Finland-Oulu", SampleList$V1)] <- "26"
SampleList <- transform(SampleList, pop_order = ifelse(pop_ploidy == "4SPN", "27", pop_order))
SampleList <- transform(SampleList, pop_order = ifelse(pop_ploidy == "4ES-", "27", pop_order))
SampleList <- transform(SampleList, pop_order = ifelse(pop_ploidy == "4FR-", "28", pop_order))
SampleList <- transform(SampleList, pop_order = ifelse(pop_ploidy == "4CH-", "29", pop_order))
SampleList <- transform(SampleList, pop_order = ifelse(pop_ploidy == "4IT-", "30", pop_order))
SampleList$pop_order[grepl("Pendula-Ireland.", SampleList$V1)] <- "31"
SampleList <- transform(SampleList, pop_order = ifelse(pop_ploidy == "4GB-", "32", pop_order))
SampleList <- transform(SampleList, pop_order = ifelse(pop_ploidy == "4DE-", "33", pop_order))
SampleList <- transform(SampleList, pop_order = ifelse(pop_ploidy == "4NO-", "34", pop_order))
SampleList <- transform(SampleList, pop_order = ifelse(pop_ploidy == "4DJU", "35", pop_order))
SampleList <- transform(SampleList, pop_order = ifelse(pop_ploidy == "4SKA", "35", pop_order))
SampleList <- transform(SampleList, pop_order = ifelse(pop_ploidy == "4VAL", "35", pop_order))
SampleList <- transform(SampleList, pop_order = ifelse(pop_ploidy == "4ASP", "35", pop_order))
SampleList <- transform(SampleList, pop_order = ifelse(pop_ploidy == "4MOK", "35", pop_order))
SampleList <- transform(SampleList, pop_order = ifelse(pop_ploidy == "4LIN", "36", pop_order))
SampleList <- transform(SampleList, pop_order = ifelse(pop_ploidy == "4GNA", "36", pop_order))
SampleList <- transform(SampleList, pop_order = ifelse(pop_ploidy == "4DOK", "36", pop_order))
SampleList <- transform(SampleList, pop_order = ifelse(pop_ploidy == "4BRA", "37", pop_order))
SampleList <- transform(SampleList, pop_order = ifelse(pop_ploidy == "4BUR", "37", pop_order))
SampleList <- transform(SampleList, pop_order = ifelse(pop_ploidy == "4SE-", "37", pop_order))
SampleList <- transform(SampleList, pop_order = ifelse(pop_ploidy == "4PIT", "38", pop_order))
SampleList <- transform(SampleList, pop_order = ifelse(pop_ploidy == "4JOK", "38", pop_order))
SampleList <- transform(SampleList, pop_order = ifelse(pop_ploidy == "4NED", "39", pop_order)) 
SampleList <- transform(SampleList, pop_order = ifelse(pop_ploidy == "4KAU", "39", pop_order)) 
SampleList <- transform(SampleList, pop_order = ifelse(pop_ploidy == "4NOR", "39", pop_order)) 
SampleList <- transform(SampleList, pop_order = ifelse(pop_ploidy2 == "4FI-BP-19", "40", pop_order))
SampleList <- transform(SampleList, pop_order = ifelse(pop_ploidy2 == "4FI-BP-20", "41", pop_order))
SampleList <- transform(SampleList, pop_order = ifelse(pop_ploidy == "4LT-", "42", pop_order))
SampleList <- transform(SampleList, pop_order = ifelse(pop_ploidy == "4POL", "43", pop_order))
SampleList$pop_order[grepl("4YT-23", SampleList$pop_ploidy2)] <- "44"
SampleList$pop_order[grepl("4YT-24", SampleList$pop_ploidy2)] <- "44"
SampleList$pop_order[grepl("4YT-27", SampleList$pop_ploidy2)] <- "45"
SampleList$pop_order[grepl("4YT-28", SampleList$pop_ploidy2)] <- "45"
SampleList <- transform(SampleList, pop_order = ifelse(pop_ploidy == "4URL", "46", pop_order))
SampleList <- transform(SampleList, pop_order = ifelse(pop_ploidy == "4NSV", "47", pop_order))
SampleList$pop_order[grepl("4YT-86", SampleList$pop_ploidy2)] <- "47"

SampleList$pop_order[grepl("Nana-Finland-Enontekio", SampleList$V1)] <- "48"
SampleList <- transform(SampleList, pop_order = ifelse(pop_ploidy == "2NAN", "48", pop_order))

SampleList <- transform(SampleList, pop_order = ifelse(pop_ploidy == "2HUM", "49", pop_order))

SampleList[, "pop_order"] <- sapply(SampleList[, "pop_order"], as.character)

## ploidy
SamplePloidy <- SampleList$ploidy

## Population
SampleList <- SampleList[,-c(2:6)]

### Write output files
setwd(OUTFolder)
write.table(SamplePloidy,file = OUTFile_ploidy, sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(SampleList,file = OUTFile_pop, sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
