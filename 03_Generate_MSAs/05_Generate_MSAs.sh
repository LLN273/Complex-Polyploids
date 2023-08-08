#!/bin/bash
#
#SBATCH -J MSA
#SBATCH -p core 
#SBATCH -n 1
#SBATCH -t 12:00:00
#SBATCH -A snic2022-22-589
#SBATCH -M snowy
#SBATCH --mail-user luis.leal@ebc.uu.se
#SBATCH --mail-type=FAIL
ulimit -c unlimited


# Script used to generate a MSA for each locus (based on exonic sequences).


# load modules
module load bioinfo-tools


#### paths and folder names

#remember initial path
SRCDIR_INI=$(pwd)                                           	 

########### folder containing input fasta files
AA=/crex1/proj/snic2020-6-184/private/Luis/P06_birch_exome_IGA-Sweden_2019/41_polarizationProtocol_fetch_gene_sequences_IUPAC_JG/BWA_unmasked

########### input/output subfolder 1
RRsub[1]=exome_WGS
RRsub[2]=exome_WGS


###########  output folder
RR=/crex1/proj/snic2020-6-184/private/Luis/P06_birch_exome_IGA-Sweden_2019/42_polarizationProtocol_MS_alignment_IUPAC_JG/BWA_unmasked


########### MAIN sample list 
SPL[1]=/crex1/proj/snic2017-7-149/private/birch_phylogeny_scripts/P03_MS_alignment/00_silver_ploidyAnalysis_REDUX.txt			
SPL[2]=/crex1/proj/snic2017-7-149/private/birch_phylogeny_scripts/P03_MS_alignment/00_white_ploidyAnalysis_REDUX.txt			






######## Phylogeny: tree backbone

## Input folder
BB[1]=/crex1/proj/snic2020-6-184/private/Luis/P09_birch_phylogeny/P02a_fetch_gene_sequences_IUPAC_SG/BWA_unmasked
BB[2]=/crex1/proj/snic2020-6-184/private/Luis/P06_birch_exome_IGA-Sweden_2019/41_polarizationProtocol_fetch_gene_sequences_IUPAC_JG/BWA_unmasked/exome_WGS

## Sample lists
SPL_BB[1]=/crex1/proj/snic2017-7-149/private/birch_phylogeny_scripts/P03_MS_alignment/00_sampleList_BASE_diploid-species.txt			
SPL_BB[2]=/crex1/proj/snic2017-7-149/private/birch_phylogeny_scripts/P03_MS_alignment/00_sampleList_BASE_diploid-species_extra.txt	

## file containing gene list (or transcript list)
GL[1]=/crex1/proj/snic2017-7-149/private/birch_phylogeny_scripts/P03_MS_alignment/00_gene_list.txt		

			               			

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




######## read sample list (backbone)

declare -A SAMPLE_LIST_BB

i=1

while read -r LINE                                                                      
do
   SAMPLE_LIST_BB[1,$i]=$LINE
   SAMPLE_LIST_BB_1[$i]=$LINE
   let "i+=1"
done < ${SPL_BB[1]}


### Number of samples
N_SAMPLES_BB[1]=${#SAMPLE_LIST_BB_1[@]}


i=1

while read -r LINE                                                                      
do
   SAMPLE_LIST_BB[2,$i]=$LINE
   SAMPLE_LIST_BB_2[$i]=$LINE
   let "i+=1"
done < ${SPL_BB[2]}


### Number of samples
N_SAMPLES_BB[2]=${#SAMPLE_LIST_BB_2[@]}





######## read gene list

declare -A GENE_LIST

i=1

while read -r LINE                                                                      
do
   GENE_LIST[$i]=$LINE
   GENE_LIST_1[$i]=$LINE
   let "i+=1"
done < ${GL[1]}


### Number of samples
N_GENES=${#GENE_LIST_1[@]}






##### reads each filename at a time, stores names in array (MAIN)
declare -A READS  

for k in `seq 1 1 2`; do 
   for i in `seq 1 1 ${N_SAMPLES[$k]}`; do       #loop through samples     
      READS[$k,$i]=${SAMPLE_LIST[$k,$i]}.exons.concatenated.fa
   done
done



##### reads each filename at a time, stores names in array (Backbone)
declare -A READS_BB  

for k in `seq 1 1 2`; do 
   for i in `seq 1 1 ${N_SAMPLES_BB[$k]}`; do       #loop through samples     
      READS_BB[$k,$i]=${SAMPLE_LIST_BB[$k,$i]}.exons.concatenated.fa
   done
done





########### Create fasta file containing sequence for a specific gene for all species

echo 'Create fasta file for each gene'
for k in `seq 1 1 2`; do 				# sample set (MAIN)

   GO_AA=${AA}/${RRsub[$k]} 
   
   for n in `seq 1 1 $N_GENES`; do 			# loop through gene list

      GO_GENE=${GENE_LIST[$n]}

      for i in `seq 1 1 ${N_SAMPLES[$k]}`; do  				# loop through samples (MAIN)
   
         GO_SAMPLE=${SAMPLE_LIST[$k,$i]}
         GO_RR=${RR}/${RRsub[$k]}/${GO_SAMPLE}   
         mkdir -p $GO_RR 

         GO_OUTFILE=${GO_GENE}.exons.fa
         rm -f $GO_RR/$GO_OUTFILE

	       echo $GO_GENE $GO_SAMPLE

         #### Add sample being tested (from MAIN list)
         DNAFILE_exons=${READS[$k,$i]} 
         #paste sequences associated to a gene into a single pasta file
         header_AUX=${GO_SAMPLE}
         body_AUX="$(grep -A 1 $GO_GENE $GO_AA/$DNAFILE_exons | grep -v "^>")"
	       echo ">"${header_AUX} >> $GO_RR/$GO_OUTFILE 
	       echo ${body_AUX} >> $GO_RR/$GO_OUTFILE

      
         #### Add backbone samples, set #1
         for j in `seq 1 1 ${N_SAMPLES_BB[1]}`; do  			# loop through samples
   
            GO_BB=${BB[1]}
	          DNAFILE_exons_BB=${READS_BB[1,$j]} 
            GO_SAMPLE_BB=${SAMPLE_LIST_BB[1,$j]}
         
            # paste sequences associated to a gene into a single pasta file
            header_AUX_BB=${GO_SAMPLE_BB}
	          body_AUX_BB="$(grep -A 1 $GO_GENE $GO_BB/$DNAFILE_exons_BB | grep -v "^>")"
	  
	          echo ">"${header_AUX_BB} >> $GO_RR/$GO_OUTFILE 
	          echo ${body_AUX_BB} >> $GO_RR/$GO_OUTFILE

         done

         #### Add backbone samples, set #2
         for j in `seq 1 1 ${N_SAMPLES_BB[2]}`; do  			# loop through samples
   
            GO_BB=${BB[2]}
	          DNAFILE_exons_BB=${READS_BB[2,$j]} 
            GO_SAMPLE_BB=${SAMPLE_LIST_BB[2,$j]}
         
            # paste sequences associated to a gene into a single pasta file
            header_AUX_BB=${GO_SAMPLE_BB}
	          body_AUX_BB="$(grep -A 1 $GO_GENE $GO_BB/$DNAFILE_exons_BB | grep -v "^>")"
	  
	          echo ">"${header_AUX_BB} >> $GO_RR/$GO_OUTFILE 
	          echo ${body_AUX_BB} >> $GO_RR/$GO_OUTFILE
            echo

         done
      done
   done
done




