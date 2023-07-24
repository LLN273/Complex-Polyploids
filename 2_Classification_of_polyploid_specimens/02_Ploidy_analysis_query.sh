#!/bin/bash
#
#SBATCH -J ploidy
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 00:10:00
##SBATCH -A snic2021-22-727
#SBATCH -A snic2022-22-589
#SBATCH -M snowy
#SBATCH --mail-user luis.leal@ebc.uu.se
#SBATCH --mail-type=FAIL
ulimit -c unlimited
#set -eo pipefail



### Script used to perform ploidy analysis for each sample



module load bioinfo-tools
module load R/4.0.0




############### input data files and paths

#remember initial path
SRCDIR_INI=$(pwd)

# sample name
SampleName=${1:?msg}	

# allelic Counts file
ACfile=${2:?msg}	

# input folder 
AA=${3:?msg}

# output folder
RR=${4:?msg}



## output file (r0ot name)
OUTFile1="${SampleName}_ploidy_analysis"


echo
echo "Sample:" $SampleName
echo "Allelic counts file:" $ACfile
echo "Input folder:" $AA
echo "Output folder:" $RR
echo "Outfile (root):" $OUTFile1
echo






## Perform ploidy analysis

## Input file (remove headers and masked sites)
cat $AA/$ACfile | grep -v '^@'| grep -wv N > $SNIC_TMP/__ACfile.txt



Rscript --no-save /crex1/proj/snic2017-7-149/private/Luis/P06_birch_exome_IGA-Sweden_2019/15_ploidy_analysis/02_Ploidy_analysis.R $SNIC_TMP/__ACfile.txt \
                                                                                                                                          $RR \
                                                                                                                                          $OUTFile1




