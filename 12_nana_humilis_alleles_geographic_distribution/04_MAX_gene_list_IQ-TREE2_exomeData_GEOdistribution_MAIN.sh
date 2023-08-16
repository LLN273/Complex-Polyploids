#!/bin/bash
#
#SBATCH -J maxg
#SBATCH -p devcore
#SBATCH -n 1
#SBATCH -t 01:00:00
#SBATCH -A snic2022-22-909
#SBATCH -M snowy
#SBATCH --mail-user luis.leal@ebc.uu.se
#SBATCH --mail-type=FAIL
ulimit -c unlimited


## Script used to generate list of loci that pass all tests (determined using 03_get_number_informative_genes2_exomeData_GEOdistribution_MAIN.sh)


##################################### paths and folder names

#remember initial path
SRCDIR_INI=$(pwd)                                           	 

# input files 
AA[1]=${SRCDIR_INI}/INFORMATIVE_SITES_exomeData_white_TETRA_100_SLIM_ALT1
AA[2]=${SRCDIR_INI}/INFORMATIVE_SITES_exomeData_white_TETRA_100_SLIM_ALT2
AA[3]=${SRCDIR_INI}/INFORMATIVE_SITES_exomeData_white_TETRA_100_SLIM_ALT3
AA[4]=${SRCDIR_INI}/INFORMATIVE_SITES_exomeData_white_TETRA_100_SLIM_ALT4

# input subfolders
SUBFOLDER[1]=SINGLE_TRANSCRIPTS

# MAIN sample list
SPL=/crex1/proj/snic2017-7-149/private/birch_phylogeny_scripts/P03_MS_alignment/00_white_ploidyAnalysis_REDUX.txt			

# reference sequence ID [can be any sequence in the MSA]
REFGEN[1]="PENDULA_FINLAND_LOIMAA.LOIMAA"
REFGEN[2]="NANA_FINLAND_ENONTEKIO.A002_2"
REFGEN[3]="PLATYPHYLLA_RUSSIA.A001_1"
REFGEN[4]="HUM02"

#### FEATURE: GENE EXON INTRON
FT[1]=exons

# model name
MODEL_outfilename=MFP_ModelFinder



######## read sample list (MAIN)

declare -A SAMPLE_LIST

i=1

while read -r LINE                                                                      
do
   SAMPLE_LIST[$i]=$LINE
   SAMPLE_LIST_1[$i]=$LINE
   let "i+=1"
done < ${SPL}


### Number of samples
N_SAMPLES=${#SAMPLE_LIST_1[@]}





for n in `seq 1 1 1`; do		# genomic feature (exons)

   GO_FEATURE=${FT[$n]}

   for j in `seq 1 1 4`; do  					         # polarization  

      for m in `seq 1 1 ${N_SAMPLES}`; do 			# polyploid ID

         GO_SAMPLE=${SAMPLE_LIST[$m]}
         old='-'
         new='_'
         GO_POLYP="${GO_SAMPLE//${old}/${new}}"	# Polyploid ID ('-' replaced by '_')
         GO_POLYP=${GO_POLYP^^}			# uppercase
         GO_REF=${REFGEN[$j]}

         GO_AA=${AA[$j]}/REFGEN_${GO_REF}
         GO_RR=$GO_AA

         echo $GO_POLYP $GO_FEATURE

         AUX_prefix_out="_SLIM-NoMAN02NoNAN01"

         #input file
         Infile=INFORMATIVE_SITES_03.${GO_FEATURE}_UFBoot_${MODEL_outfilename}_${GO_POLYP}${AUX_prefix_out}.txt

         #outfile
         OUT_FILE=00_Birch_gene_list_${GO_FEATURE}_${MODEL_outfilename}_${GO_POLYP}${AUX_prefix_out}_TRIM_MAX-00000.txt
         

         cat $GO_AA/${Infile} | grep -v FAIL | cut -f1 -d' ' > $GO_RR/$OUT_FILE

      done
   done
done


echo
echo "Done!"









