#!/bin/bash
#
#SBATCH -J preASTR
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 24:00:00
#SBATCH -A snic2022-22-589
#SBATCH -M snowy
#SBATCH --mail-user luis.leal@ebc.uu.se
#SBATCH --mail-type=FAIL
ulimit -c unlimited

 
## Script used to create file with IQ-TEST2 consensus trees obtained for each gene & sample | Use gene list that passed all filtering steps

echo
echo "Starting Uppmax jobs ..."
echo

# load modules
module load bioinfo-tools


############################################################################## paths and folder names

####################################### remember initial path
SRCDIR_INI=$(pwd)  
                                        
####################################### input folder 
AA[1]=/crex1/proj/snic2020-6-184/private/Luis/P06_birch_exome_IGA-Sweden_2019/43b_polarizationProtocol_IQ-TREE2_IUPAC_ALT1_SLIM/BWA_unmasked
AA[2]=/crex1/proj/snic2020-6-184/private/Luis/P06_birch_exome_IGA-Sweden_2019/44b_polarizationProtocol_IQ-TREE2_IUPAC_ALT2_SLIM/BWA_unmasked
AA[3]=/crex1/proj/snic2020-6-184/private/Luis/P06_birch_exome_IGA-Sweden_2019/46b_polarizationProtocol_IQ-TREE2_IUPAC_ALT3_SLIM/BWA_unmasked
AA[4]=/crex1/proj/snic2020-6-184/private/Luis/P06_birch_exome_IGA-Sweden_2019/45b_polarizationProtocol_IQ-TREE2_IUPAC_ALT4_SLIM/BWA_unmasked

# input subfolders
SUBFOLDER[1]=SINGLE_TRANSCRIPTS

##################################### Gene lists folder
GLF[1]=${SRCDIR_INI}/INFORMATIVE_SITES_exomeData_white_TETRA_100_SLIM_ALT1
GLF[2]=${SRCDIR_INI}/INFORMATIVE_SITES_exomeData_white_TETRA_100_SLIM_ALT2
GLF[3]=${SRCDIR_INI}/INFORMATIVE_SITES_exomeData_white_TETRA_100_SLIM_ALT3
GLF[4]=${SRCDIR_INI}/INFORMATIVE_SITES_exomeData_white_TETRA_100_SLIM_ALT4

###################################### output folder
RR[1]=/crex1/proj/snic2020-6-184/private/Luis/P06_birch_exome_IGA-Sweden_2019/48_phylogeny_NAN_HUM_PEND_genes_geo_distribution/ALT1
RR[2]=/crex1/proj/snic2020-6-184/private/Luis/P06_birch_exome_IGA-Sweden_2019/48_phylogeny_NAN_HUM_PEND_genes_geo_distribution/ALT2
RR[3]=/crex1/proj/snic2020-6-184/private/Luis/P06_birch_exome_IGA-Sweden_2019/48_phylogeny_NAN_HUM_PEND_genes_geo_distribution/ALT3
RR[4]=/crex1/proj/snic2020-6-184/private/Luis/P06_birch_exome_IGA-Sweden_2019/48_phylogeny_NAN_HUM_PEND_genes_geo_distribution/ALT4

###################################### MAIN sample list
SPL=/crex1/proj/snic2017-7-149/private/birch_phylogeny_scripts/P03_MS_alignment/00_white_ploidyAnalysis_REDUX.txt			

###################################### reference sequence ID [can be any sequence in the MSA]
REFGEN[1]="PENDULA_FINLAND_LOIMAA.LOIMAA"
REFGEN[2]="NANA_FINLAND_ENONTEKIO.A002_2"
REFGEN[3]="PLATYPHYLLA_RUSSIA.A001_1"
REFGEN[4]="HUM02"

# model name
MODEL_outfilename=MFP_ModelFinder

#### FEATURE: GENE EXON INTRON
FT[1]=exons



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

   for j in `seq 4 1 4`; do  					         # polarization  

      for m in `seq 1 1 ${N_SAMPLES}`; do 			# polyploid ID

         GO_SAMPLE=${SAMPLE_LIST[$m]}
         old='-'
         new='_'
         GO_POLYP="${GO_SAMPLE//${old}/${new}}"	# Polyploid ID ('-' replaced by '_')
         GO_POLYP=${GO_POLYP^^}			# uppercase
         GO_REF=${REFGEN[$j]}

         GO_AA=${AA[$j]}/REFGEN_${GO_REF}/${SUBFOLDER[$n]}/${GO_SAMPLE}
         GO_RR=${RR[$j]}/REFGEN_${GO_REF}/${SUBFOLDER[$n]}/${GO_SAMPLE}
         mkdir -p $GO_RR

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




