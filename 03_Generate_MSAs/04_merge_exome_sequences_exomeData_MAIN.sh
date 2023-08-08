#!/bin/bash
#
#SBATCH -J mergeExome
#SBATCH -p core 
#SBATCH -n 1
#SBATCH -t 12:00:00
#SBATCH -A snic2022-22-589
#SBATCH -M snowy
#SBATCH --mail-user luis.leal@ebc.uu.se
#SBATCH --mail-type=FAIL
ulimit -c unlimited



## Scrip used to concatenate exonic sequences associated to the same gene.
## Performed independently for each sample.



# load modules
module load bioinfo-tools




##### Paths and folders

#remember initial path
SRCDIR_INI=$(pwd) 

# Path to folder containing input data
AA=/crex1/proj/snic2020-6-184/private/Luis/P06_birch_exome_IGA-Sweden_2019/41_polarizationProtocol_fetch_gene_sequences_IUPAC_JG/BWA_unmasked

# input/output subfolder 1
RRsub[1]=exome_WGS
RRsub[2]=exome_WGS
RRsub[3]=exome_WGS

# Sample list
SPL[1]=/crex1/proj/snic2017-7-149/private/Luis/P06_birch_exome_IGA-Sweden_2019/15_ploidy_analysis/00_silver_ploidyAnalysis.txt			   
SPL[2]=/crex1/proj/snic2017-7-149/private/Luis/P06_birch_exome_IGA-Sweden_2019/15_ploidy_analysis/00_white_ploidyAnalysis.txt			     
SPL[3]=/crex1/proj/snic2017-7-149/private/Luis/P06_birch_exome_IGA-Sweden_2019/15_ploidy_analysis/00_platyphylla_ploidyAnalysis.txt		




######## read sample list

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





i=1

while read -r LINE                                                                      
do
   SAMPLE_LIST[3,$i]=$LINE
   SAMPLE_LIST_3[$i]=$LINE
   let "i+=1"
done < ${SPL[3]}



### Number of samples
N_SAMPLES[3]=${#SAMPLE_LIST_3[@]}


# Genomic feature
Gfeature[1]=exons
Gfeature[2]=introns



cd $SRCDIR_INI


for n in `seq 1 1 1`; do 						# loop through genomic features (exons, introns)

   rm -f $SRCDIR_INI/slurm_merge_${Gfeature[$n]}_sequences_exomeData.txt 
  
   for k in `seq 1 1 3`; do 						# cycle through datasets
      for i in `seq 1 1 ${N_SAMPLES[$k]}`; do                 		# cycle through samples

         GO_AA=${AA}/${RRsub[$k]}
         GO_RR=$GO_AA     
         DNAFILE_A=${SAMPLE_LIST[$k,$i]}.${Gfeature[$n]}.fa
	       GO_SAMPLE=${SAMPLE_LIST[$k,$i]}  
	     
         echo $GO_SAMPLE   
         echo $DNAFILE_A   
         echo $GO_AA
         echo $GO_RR
                                           
	       # Order fasta headers alphabetically 
         perl -pe 's/[\r\n]+/;/g; s/>/\n>/g' $GO_AA/$DNAFILE_A | sort -t"[" -k2,2V | sed 's/;/\n/g' | sed '/^$/d' > $GO_RR/${DNAFILE_A%.fa}_aux1.fa

         # reformat fasta headers (keep only gene name)
         cat $GO_RR/${DNAFILE_A%.fa}_aux1.fa | awk -F".m" '/>/{$0=">"$1}1' > $GO_RR/${DNAFILE_A%.fa}_aux2.fa

		     #### concatenate exon sequences associated to the same gene (for each individual transcript)
		     awk '/^>/ {if(prev!=$0) {prev=$0;printf("\n%s\n",$0);} next;} {printf("%s",$0);} END {printf("\n");}' $GO_RR/${DNAFILE_A%.fa}_aux2.fa > $GO_RR/_AUX_${GO_SAMPLE}.${Gfeature[$n]}.concatenated.fa
		     cat $GO_RR/_AUX_${GO_SAMPLE}.${Gfeature[$n]}.concatenated.fa | tail -n +2 > $GO_RR/${GO_SAMPLE}.${Gfeature[$n]}.concatenated.fa

         rm -f $GO_RR/_AUX_${GO_SAMPLE}.${Gfeature[$n]}.concatenated.fa
         rm -f $GO_RR/${DNAFILE_A%.fa}_aux*.fa

         # Check number of nucleotides and number of headers in fasta file
         aux_1="$(cat $GO_AA/$DNAFILE_A | grep -v '^>' | grep -io [A-T] | wc -l  )"
         aux_1_concat="$(cat $GO_RR/${GO_SAMPLE}.${Gfeature[$n]}.concatenated.fa | grep -v '^>' | grep -io [A-T] | wc -l  )"
         aux_headers="$(cat $GO_RR/${GO_SAMPLE}.${Gfeature[$n]}.concatenated.fa | grep '^>' | wc -l  )"
         echo $GO_RR/${GO_SAMPLE}.${Gfeature[$n]}.concatenated.fa $aux_1 $aux_1_concat $aux_headers >> $SRCDIR_INI/slurm_merge_${Gfeature[$n]}_sequences_exomeData.txt 
                                 
         echo

      done
   done
done


echo done!
echo



