#!/bin/bash
#
#SBATCH -J preASTR1
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 24:00:00
#SBATCH -A snic2021-22-727
#SBATCH -M snowy
#SBATCH --mail-user luis.leal@ebc.uu.se
#SBATCH --mail-type=FAIL
ulimit -c unlimited
 

## Script used to create ASTRAL input file.
## It collects the consensus trees obtained for each locus, using IQ-TEST2, that passed the filtering criteria (previous step).



## load modules
module load bioinfo-tools


############################################################################## paths and folder names

####################################### remember initial path
SRCDIR_INI=$(pwd)  

####################################### input folder 

declare -A AA
declare -A RR

AA=/crex1/proj/snic2020-6-184/private/Luis/P06_birch_exome_IGA-Sweden_2019/43b_polarizationProtocol_IQ-TREE2_IUPAC_ALT1_SLIM/BWA_unmasked

# input subfolders
SUBFOLDER[1]=SINGLE_TRANSCRIPTS
SUBFOLDER[2]=SINGLE_INTRONS


###################################### output folder
RR=/crex1/proj/snic2020-6-184/private/Luis/P06_birch_exome_IGA-Sweden_2019/43c_polarizationProtocol_ASTRAL_TETRA_ALT1_SLIM/BWA_unmasked


##################################### Gene lists folder
GLF[1]=/crex1/proj/snic2017-7-149/private/birch_phylogeny_scripts/P05_phylogeny_IQ-TREE2/INFORMATIVE_SITES_exomeData_silver_TETRA_100_SLIM
GLF[2]=/crex1/proj/snic2017-7-149/private/birch_phylogeny_scripts/P05_phylogeny_IQ-TREE2/INFORMATIVE_SITES_exomeData_white_TETRA_100_SLIM


##################################### Sample list
SPL[1]=/crex1/proj/snic2017-7-149/private/birch_phylogeny_scripts/P03_MS_alignment/00_silver_ploidyAnalysis_REDUX.txt			
SPL[2]=/crex1/proj/snic2017-7-149/private/birch_phylogeny_scripts/P03_MS_alignment/00_white_ploidyAnalysis_REDUX.txt			


##################################### reference sequence ID [can be any sequence in the MSA]
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
         GO_RR=${RR}/REFGEN_${GO_REF}/${SUBFOLDER[$n]}/${GO_SAMPLE}
         mkdir -p $GO_RR

         # NO SLIM
         #AUX_prefix_out=""

         # No MAN02, No NAN01
         AUX_prefix_out="_SLIM-NoMAN02NoNAN01"

         # Gene list
         GL=${GLF[$j]}/REFGEN_${GO_REF}/00_Birch_gene_list_${GO_FEATURE}_${MODEL_outfilename}_${GO_POLYP}${AUX_prefix_out}_TRIM_MAX-00000.txt


         ######## read gene list

         unset GENE_LIST
         x=1

         while read -r LINE                                                                      
         do
            GENE_LIST[$x]=$LINE
            let "x+=1"
         done < $GL

         N_GENES=${#GENE_LIST[@]}

        
         # input file suffix
         InFileSUF=${GO_FEATURE}_trimN-ALT.fasta_UFBoot_${MODEL_outfilename}${AUX_prefix_out}.iqtree

         #output file
         OUTfile=GENE_TREES.${GO_FEATURE}_trimN_SLIM-${GO_POLYP}-ALT.fasta_UFBoot_${MODEL_outfilename}${AUX_prefix_out}.newick

         rm -f ${GO_RR}/$OUTfile

         echo $InFileSUF
         echo $OUTfile
         echo $GL
         echo $GO_AA
         echo $GO_RR
         echo "Total number of genes:" $N_GENES
         echo

         cd $GO_AA

         for i in `seq 1 1 $N_GENES`; do                 # cycle through genes

            GO_GENE=${GENE_LIST[$i]}
         
            # get IQ-TEST2 consensus tree for each gene           
            cat ${GO_GENE}.${InFileSUF} | grep -A 2 'Consensus tree in newick format:' | tail -1 >> ${GO_RR}/${OUTfile}
      
         done
      done
   done
done


echo "Done!"


