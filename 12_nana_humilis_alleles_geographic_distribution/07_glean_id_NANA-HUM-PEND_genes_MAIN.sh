#!/bin/bash
#
#SBATCH -J geneID
#SBATCH -p devcore
#SBATCH -n 1
#SBATCH -t 1:00:00
#SBATCH -A snic2022-22-909
#SBATCH -M snowy
#SBATCH --mail-user luis.leal@ebc.uu.se
#SBATCH --mail-type=FAIL
ulimit -c unlimited

 

## Generates summary file listing B. humilis genes of B. nana or B. humilis origin for each sample

echo
echo "Starting Uppmax jobs ..."
echo


# load modules
module load bioinfo-tools


#### paths and folder names

#remember initial path
SRCDIR_INI=$(pwd)                                           	 

# input folder 
AA[1]=/crex1/proj/snic2020-6-184/private/Luis/P06_birch_exome_IGA-Sweden_2019/48_phylogeny_NAN_HUM_PEND_genes_geo_distribution/ALT1
AA[2]=/crex1/proj/snic2020-6-184/private/Luis/P06_birch_exome_IGA-Sweden_2019/48_phylogeny_NAN_HUM_PEND_genes_geo_distribution/ALT2
AA[3]=/crex1/proj/snic2020-6-184/private/Luis/P06_birch_exome_IGA-Sweden_2019/48_phylogeny_NAN_HUM_PEND_genes_geo_distribution/ALT3
AA[4]=/crex1/proj/snic2020-6-184/private/Luis/P06_birch_exome_IGA-Sweden_2019/48_phylogeny_NAN_HUM_PEND_genes_geo_distribution/ALT4

# input subfolders
SUBFOLDER[1]=SINGLE_TRANSCRIPTS

# MAIN sample list (triploids and 4X-Bpendula removed)
SPL=/crex1/proj/snic2017-7-149/private/birch_phylogeny_scripts/P11_phylogeny_NAN_HUM_PEND_genes_geo_distribution/00_white_ploidyAnalysis_REDUX_CLEAN.txt

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




echo

for n in `seq 1 1 1`; do		# genomic feature (exon intron)

   GO_FEATURE=${FT[$n]}
   
   GO_RR=$SRCDIR_INI

   for j in `seq 1 1 4`; do  					         # polarization  


      ## output file
      outFile_1="ID_genes-nana_summary_"${OUTFroot[$j]}".txt"
      outFile_2="ID_genes-humilis_summary_"${OUTFroot[$j]}".txt"
      outFile_3="ID_genes-pendula_summary_"${OUTFroot[$j]}".txt"
      outFile_4="ID_genes-pendplaty_summary_"${OUTFroot[$j]}".txt"
      rm -f $GO_RR/$outFile_1
      rm -f $GO_RR/$outFile_2
      rm -f $GO_RR/$outFile_3
      rm -f $GO_RR/$outFile_4

      touch $GO_RR/$outFile_1
      touch $GO_RR/$outFile_2
      touch $GO_RR/$outFile_3
      touch $GO_RR/$outFile_4
      
      for m in `seq 1 1 ${N_SAMPLES}`; do 			# polyploid ID

         GO_SAMPLE=${SAMPLE_LIST[$m]}
         old='-'
         new='_'
         GO_POLYP="${GO_SAMPLE//${old}/${new}}"	# Polyploid ID ('-' replaced by '_')
         GO_POLYP=${GO_POLYP^^}			# uppercase
         GO_REF=${REFGEN[$j]}

         GO_AA=${AA[$j]}/REFGEN_${GO_REF}/${SUBFOLDER[$n]}/${GO_SAMPLE}

         GO_POLYPT=${GO_POLYP}

         ## input file
         In_file_1="ID_genes_"${OUTFroot[$j]}"-nana.txt"
         In_file_2="ID_genes_"${OUTFroot[$j]}"-humilis.txt"
         In_file_3="ID_genes_"${OUTFroot[$j]}"-pendula.txt"
         In_file_4="ID_genes_"${OUTFroot[$j]}"-pendplaty.txt"

         echo $GO_POLYP $GO_FEATURE

         # concatenate files, one at a time

         echo $GO_SAMPLE > $SNIC_TMP/_inFile_1.txt
         cat $GO_AA/$In_file_1 >> $SNIC_TMP/_inFile_1.txt
         awk 'BEGIN { ORS = " " } { print }' $SNIC_TMP/_inFile_1.txt > $SNIC_TMP/_inFile_1B.txt
         echo "" >> $SNIC_TMP/_inFile_1B.txt 

         echo $GO_SAMPLE > $SNIC_TMP/_inFile_2.txt
         cat $GO_AA/$In_file_2 >> $SNIC_TMP/_inFile_2.txt
         awk 'BEGIN { ORS = " " } { print }' $SNIC_TMP/_inFile_2.txt > $SNIC_TMP/_inFile_2B.txt
         echo "" >> $SNIC_TMP/_inFile_2B.txt 
		 
         echo $GO_SAMPLE > $SNIC_TMP/_inFile_3.txt
         cat $GO_AA/$In_file_3 >> $SNIC_TMP/_inFile_3.txt
         awk 'BEGIN { ORS = " " } { print }' $SNIC_TMP/_inFile_3.txt > $SNIC_TMP/_inFile_3B.txt
         echo "" >> $SNIC_TMP/_inFile_3B.txt 
		 
         echo $GO_SAMPLE > $SNIC_TMP/_inFile_4.txt
         cat $GO_AA/$In_file_4 >> $SNIC_TMP/_inFile_4.txt
         awk 'BEGIN { ORS = " " } { print }' $SNIC_TMP/_inFile_4.txt > $SNIC_TMP/_inFile_4B.txt
         echo "" >> $SNIC_TMP/_inFile_4B.txt 

         cat $GO_RR/$outFile_1 $SNIC_TMP/_inFile_1B.txt >  $SNIC_TMP/_outFile_1_aux.txt
         cat $GO_RR/$outFile_2 $SNIC_TMP/_inFile_2B.txt >  $SNIC_TMP/_outFile_2_aux.txt
         cat $GO_RR/$outFile_3 $SNIC_TMP/_inFile_3B.txt >  $SNIC_TMP/_outFile_3_aux.txt
         cat $GO_RR/$outFile_4 $SNIC_TMP/_inFile_4B.txt >  $SNIC_TMP/_outFile_4_aux.txt

         mv $SNIC_TMP/_outFile_1_aux.txt $GO_RR/$outFile_1
         mv $SNIC_TMP/_outFile_2_aux.txt $GO_RR/$outFile_2
         mv $SNIC_TMP/_outFile_3_aux.txt $GO_RR/$outFile_3
         mv $SNIC_TMP/_outFile_4_aux.txt $GO_RR/$outFile_4

      done
	  
	  # replace space by comma (first remove space at end of each line)
	  cat $GO_RR/$outFile_1 | sed 's/ *$//' | sed -e 's/ /,/g' > $SNIC_TMP/_aux1
	  cat $GO_RR/$outFile_2 | sed 's/ *$//' | sed -e 's/ /,/g' > $SNIC_TMP/_aux2
	  cat $GO_RR/$outFile_3 | sed 's/ *$//' | sed -e 's/ /,/g' > $SNIC_TMP/_aux3
	  cat $GO_RR/$outFile_4 | sed 's/ *$//' | sed -e 's/ /,/g' > $SNIC_TMP/_aux4

	  mv $SNIC_TMP/_aux1 $GO_RR/$outFile_1
	  mv $SNIC_TMP/_aux2 $GO_RR/$outFile_2
	  mv $SNIC_TMP/_aux3 $GO_RR/$outFile_3
	  mv $SNIC_TMP/_aux4 $GO_RR/$outFile_4
	  
	  # clean auxiliary files
	  rm -f $SNIC_TMP/_*
	  
   done
done


echo
echo "Done!"







