#!/bin/bash
#
#SBATCH -J getAlleles
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 02:00:00
#SBATCH -A snic2022-22-909
#SBATCH -M snowy
#SBATCH --mail-user luis.leal@ebc.uu.se
#SBATCH --mail-type=FAIL
ulimit -c unlimited
 

### Script used to identify B. pubescens alleles of B. nana or B. humilis origin 


echo
echo "Starting Uppmax jobs ..."
echo

# load modules
module load bioinfo-tools
module load R/4.0.0

#### paths and folder names

#remember initial path
SRCDIR_INI=$(pwd)                                           	 


####################################### input folder 
AA[1]=/crex1/proj/snic2020-6-184/private/Luis/P06_birch_exome_IGA-Sweden_2019/48_phylogeny_NAN_HUM_PEND_genes_geo_distribution/ALT1
AA[2]=/crex1/proj/snic2020-6-184/private/Luis/P06_birch_exome_IGA-Sweden_2019/48_phylogeny_NAN_HUM_PEND_genes_geo_distribution/ALT2
AA[3]=/crex1/proj/snic2020-6-184/private/Luis/P06_birch_exome_IGA-Sweden_2019/48_phylogeny_NAN_HUM_PEND_genes_geo_distribution/ALT3
AA[4]=/crex1/proj/snic2020-6-184/private/Luis/P06_birch_exome_IGA-Sweden_2019/48_phylogeny_NAN_HUM_PEND_genes_geo_distribution/ALT4

# input subfolders
SUBFOLDER[1]=SINGLE_TRANSCRIPTS

##################################### Gene lists folder
GLF[1]=${SRCDIR_INI}/INFORMATIVE_SITES_exomeData_white_TETRA_100_SLIM_ALT1
GLF[2]=${SRCDIR_INI}/INFORMATIVE_SITES_exomeData_white_TETRA_100_SLIM_ALT2
GLF[3]=${SRCDIR_INI}/INFORMATIVE_SITES_exomeData_white_TETRA_100_SLIM_ALT3
GLF[4]=${SRCDIR_INI}/INFORMATIVE_SITES_exomeData_white_TETRA_100_SLIM_ALT4

# MAIN sample list
SPL=/crex1/proj/snic2017-7-149/private/birch_phylogeny_scripts/P03_MS_alignment/00_white_ploidyAnalysis_REDUX.txt

# reference sequence ID
REFGEN[1]="PENDULA_FINLAND_LOIMAA.LOIMAA"
REFGEN[2]="NANA_FINLAND_ENONTEKIO.A002_2"
REFGEN[3]="PLATYPHYLLA_RUSSIA.A001_1"
REFGEN[4]="HUM02"

# Outfile (root)
OUTFroot[1]="ALT1"
OUTFroot[2]="ALT2"
OUTFroot[3]="ALT3"
OUTFroot[4]="ALT4"

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

   for j in `seq 1 1 4`; do  					  # polarization  

      for m in `seq 1 1 ${N_SAMPLES}`; do 			# polyploid ID
	  
            GO_SAMPLE=${SAMPLE_LIST[$m]}
            old='-'
            new='_'
            GO_POLYP="${GO_SAMPLE//${old}/${new}}"	# Polyploid ID ('-' replaced by '_')
            GO_POLYP=${GO_POLYP^^}			# uppercase
            GO_REF=${REFGEN[$j]}

            GO_AA=${AA[$j]}/REFGEN_${GO_REF}/${SUBFOLDER[$n]}/${GO_SAMPLE}
            GO_RR=$GO_AA
      
            GO_POLYPT=${GO_POLYP}

            AUX_prefix_out="_SLIM-NoMAN02NoNAN01"

            # Gene list
            GL=${GLF[$j]}/REFGEN_${GO_REF}/00_Birch_gene_list_${GO_FEATURE}_${MODEL_outfilename}_${GO_POLYP}${AUX_prefix_out}_TRIM_MAX-00000.txt

            ## file with gene trees
            gene_Trees_file=$GO_AA/GENE_TREES.${GO_FEATURE}_trimN_SLIM-${GO_POLYP}-ALT.fasta_UFBoot_${MODEL_outfilename}${AUX_prefix_out}.newick

            ## tetraploid ID
            TETRAPLOID=${GO_POLYPT}

            ## outgroup
            outgroup_ID="ALNUS_GLUTINOSA.V70901"

            ## output file (root)
            OUTFile_aux="ID_genes"
            OUTFile=${OUTFile_aux}_${OUTFroot[$j]}

            echo $GO_POLYP $GO_FEATURE $GO_AA

            Rscript --no-save $SRCDIR_INI/06c_id_NANA-HUM-genes_exomeData_MAIN.R $gene_Trees_file \
                                                                                 $TETRAPLOID \
                                                                                 $outgroup_ID \
                                                                                 $GO_RR \
                                                                                 $OUTFile \
                                                                                 $GL


      done
   done
done


echo
echo "Done!"







