#!/bin/bash
#
#SBATCH -J maxg1
#SBATCH -p devcore
#SBATCH -n 1
#SBATCH -t 01:00:00
#SBATCH -A snic2021-22-727
#SBATCH -M snowy
#SBATCH --mail-user luis.leal@ebc.uu.se
#SBATCH --mail-type=FAIL
ulimit -c unlimited
#set -eo pipefail


## Generate list of loci that pass all tests (determined using 02_get_number_informative_genes2_exomeData_MAIN.sh)


#### paths and folder names

#remember initial path
SRCDIR_INI=$(pwd)                                           	 

###################################### input files 

declare -A AA
			
AA[1]=${SRCDIR_INI}/INFORMATIVE_SITES_exomeData_silver_TETRA_100_SLIM
AA[2]=${SRCDIR_INI}/INFORMATIVE_SITES_exomeData_white_TETRA_100_SLIM

# input subfolders
SUBFOLDER[1]=SINGLE_TRANSCRIPTS
SUBFOLDER[2]=SINGLE_INTRONS

# MAIN sample list (REDUX: excludes diploid samples used in phylogeny backbone, which be included in a different list, below)
SPL[1]=/crex1/proj/snic2017-7-149/private/birch_phylogeny_scripts/P03_MS_alignment/00_silver_ploidyAnalysis_REDUX.txt			
SPL[2]=/crex1/proj/snic2017-7-149/private/birch_phylogeny_scripts/P03_MS_alignment/00_white_ploidyAnalysis_REDUX.txt			

# reference sequence ID [can be any sequence in the MSA]
REFGEN="PENDULA_FINLAND_LOIMAA.LOIMAA"




######## read sample list (MAIN)

declare -A SAMPLE_LIST

i=1

while read -r LINE                                                                      
do
   SAMPLE_LIST[1,$i]=$LINE
   SAMPLE_LIST_1[$i]=$LINE
   let "i+=1"
done < ${SPL[1]}


### Number of samples
N_SAMPLES[1]=${#SAMPLE_LIST_1[@]}



i=1

while read -r LINE                                                                      
do
   SAMPLE_LIST[2,$i]=$LINE
   SAMPLE_LIST_2[$i]=$LINE
   let "i+=1"
done < ${SPL[2]}


### Number of samples
N_SAMPLES[2]=${#SAMPLE_LIST_2[@]}




#### FEATURE: GENE EXON INTRON
FT[1]=exons
FT[2]=introns

# model name
MODEL_outfilename=MFP_ModelFinder



for n in `seq 1 1 1`; do		# genomic feature (exon intron)

   GO_FEATURE=${FT[$n]}

   for j in `seq 1 1 2`; do  					         # sample set  

      for m in `seq 1 1 ${N_SAMPLES[$j]}`; do 			# polyploid ID

         GO_SAMPLE=${SAMPLE_LIST[$j,$m]}
         old='-'
         new='_'
         GO_POLYP="${GO_SAMPLE//${old}/${new}}"	# Polyploid ID ('-' replaced by '_')
         GO_POLYP=${GO_POLYP^^}			# uppercase
         GO_REF=${REFGEN}

         GO_AA=${AA[$j]}/REFGEN_${GO_REF}
         GO_RR=$GO_AA

         echo $GO_POLYP $GO_FEATURE

         # NO SLIM
         #AUX_prefix_out=""

         # No MAN02, No NAN01
         AUX_prefix_out="_SLIM-NoMAN02NoNAN01"

         #input file
         Infile=INFORMATIVE_SITES_03.${GO_FEATURE}_UFBoot_${MODEL_outfilename}_${GO_POLYP}${AUX_prefix_out}.txt

         #outfile
         OUT_FILE=00_Birch_gene_list_${GO_FEATURE}_${MODEL_outfilename}_${GO_POLYP}${AUX_prefix_out}_TRIM_MAX-00000.txt
         
         cat $GO_AA/${Infile} | grep -v FAIL | cut -f1 -d' ' > $GO_RR/$OUT_FILE

      done
   done
done


