#!/bin/bash
#
#SBATCH -J sfreq1
#SBATCH -p devcore
#SBATCH -n 1
#SBATCH -t 01:00:00
#SBATCH -A snic2021-22-727
#SBATCH -M snowy
#SBATCH --mail-user luis.leal@ebc.uu.se
#SBATCH --mail-type=FAIL
ulimit -c unlimited
 

####### Script used to identify the sister species of the polyploid


#### load modules
module load bioinfo-tools
module load R/4.0.0


#### paths and folder names

#remember initial path
SRCDIR_INI=$(pwd)                                           	 

# input folder 

declare -A AA
declare -A RR

AA=/crex1/proj/snic2020-6-184/private/Luis/P06_birch_exome_IGA-Sweden_2019/43c_polarizationProtocol_ASTRAL_TETRA_ALT1_SLIM/BWA_unmasked

# input subfolders
SUBFOLDER[1]=SINGLE_TRANSCRIPTS
SUBFOLDER[2]=SINGLE_INTRONS


# Sample list 
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





#### model name
MODEL_outfilename=MFP_ModelFinder

#### FEATURE: GENE EXON INTRON
FT[1]=exons
FT[2]=introns



for n in `seq 1 1 1`; do		# genomic feature (exon intron)

   GO_FEATURE=${FT[$n]}

   for j in `seq 2 1 2`; do  					         # sample sets 

      for m in `seq 1 1 ${N_SAMPLES[$j]}`; do 			# polyploid ID

         GO_SAMPLE=${SAMPLE_LIST[$j,$m]}
         old='-'
         new='_'
         GO_POLYP="${GO_SAMPLE//${old}/${new}}"	# Polyploid ID ('-' replaced by '_')
         GO_POLYP=${GO_POLYP^^}			# uppercase
         GO_REF=${REFGEN}

         GO_AA=${AA}/REFGEN_${GO_REF}/${SUBFOLDER[$n]}/${GO_SAMPLE}
         GO_RR=$GO_AA
      
         GO_POLYPT=${GO_POLYP}

         # NO SLIM
         #AUX_prefix_out=""

         # No MAN02, No NAN01
         AUX_prefix_out="_SLIM-NoMAN02NoNAN01"

         ## file with gene trees
         gene_Trees_file=$GO_AA/GENE_TREES.${GO_FEATURE}_trimN_SLIM-${GO_POLYP}-ALT.fasta_UFBoot_${MODEL_outfilename}${AUX_prefix_out}.newick

         ## tetraploid ID
         TETRAPLOID=${GO_POLYPT}

         ## outgroup
         outgroup_ID="ALNUS_GLUTINOSA.V70901"

         ## output file
         OUTFile="sister_ID_analysis_GENE_TREES.${GO_FEATURE}_trimN_SLIM-${GO_POLYP}-ALT.fasta_UFBoot_${MODEL_outfilename}${AUX_prefix_out}.txt"

         echo $GO_POLYP $GO_FEATURE $GO_AA

         Rscript --no-save /crex1/proj/snic2017-7-149/private/Luis/P10_SIMULATED_FASTA_PHYLOGENY/05_phylogeny_ASTRAL/06_IQTREE_gene_tree_DISTANCE_ALT.R $gene_Trees_file \
                                                                  $TETRAPLOID \
                                                                  $outgroup_ID \
                                                                  $GO_RR \
                                                                  $OUTFile


      done
   done
done

