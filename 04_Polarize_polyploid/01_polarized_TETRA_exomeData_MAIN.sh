#!/bin/bash
#

## Script used to polarize tetraploid sequence in the MSA.
## Polarization done separately for each locus.


#### paths and folder names

# remember initial path
SRCDIR_INI=$(pwd)      

# Path to folder containing input data
AA=/crex1/proj/snic2020-6-184/private/Luis/P06_birch_exome_IGA-Sweden_2019/42_polarizationProtocol_MS_alignment_IUPAC_JG/BWA_unmasked

# input/output subfolder 1
RRsub[1]=exome_WGS

# output folder
RR=/crex1/proj/snic2020-6-184/private/Luis/P06_birch_exome_IGA-Sweden_2019/43a_polarizationProtocol_Polarize_IUPAC_ALT1/BWA_unmasked

# output subfolders
SUBFOLDER[1]=SINGLE_TRANSCRIPTS
SUBFOLDER[2]=SINGLE_INTRONS

# MAIN sample list
SPL[1]=/crex1/proj/snic2017-7-149/private/birch_phylogeny_scripts/P03_MS_alignment/00_white_ploidyAnalysis_REDUX.txt





# reference sequence ID [can be any sequence in the MSA]
REFGEN="PENDULA_FINLAND_LOIMAA.LOIMAA"

# Number of species (including root species)
NSPEC=14	

#file containing gene list (or transcript list)
GL=/crex1/proj/snic2017-7-149/private/birch_phylogeny_scripts/P03_MS_alignment/00_gene_list.txt	



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




# genomic feature
GF[1]=exons
GF[2]=introns



for n in `seq 1 1 1`; do		# genomic feature (exon or intron)

   GO_FEATURE=${GF[$n]}

   for k in `seq 1 1 1`; do 			# sample set

      for j in `seq 1 1 ${N_SAMPLES[$k]}`; do 			# polyploid ID

         GO_SAMPLE=${SAMPLE_LIST[$k,$j]}
         old='-'
         new='_'
         GO_POLYP="${GO_SAMPLE//${old}/${new}}"	# Polyploid ID ('-' replaced by '_')
         GO_POLYP=${GO_POLYP^^}			# uppercase
         GO_REF=${REFGEN}
         GO_NSPEC=${NSPEC}

         GO_AA=${AA}/${RRsub[$k]}/${GO_SAMPLE}
         GO_RR=${RR}/REFGEN_${GO_REF}/${SUBFOLDER[$n]}/${GO_SAMPLE}
         mkdir -p $GO_RR

         cd $GO_RR
         rsync -ah $SRCDIR_INI/01s_polarizeTETRA.py .
   
         echo $GO_NSPEC $GO_FEATURE
         echo $GO_POLYP
         echo $GO_REF
         echo $GL
         echo $GO_AA
         echo $GO_RR

         cd $SRCDIR_INI

         sbatch ./01_polarized_TETRA_exomeData_query.sh $GO_AA \
                                                        $GO_RR \
                                                        $GO_NSPEC \
                                                        $GO_POLYP \
                                                        $GO_REF \
                                                        $GL \
                                                        $GO_FEATURE

         echo

      done
   done     
done





