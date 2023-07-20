#!/bin/bash -l
#SBATCH -J checkLIBS
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 02:00:00
#SBATCH -A snic2021-22-291
#SBATCH -M snowy
#SBATCH --mail-user luis.leal@ebc.uu.se
#SBATCH --mail-type=FAIL
ulimit -c unlimited
set -eo pipefail




##### Script used to get detailed information about each Illumina library (instrument name, run number, flow-cell ID, etc)



# load modules
module load bioinfo-tools



######### PATHS and filenames

#initial folder
SRCDIR=$(pwd)                                                   

#read name of files containing paired reads (with path)
READ1=${1:?msg}		
READ2=${2:?msg}

#read name of files containing paired reads (without path)
rawReads1=${3:?msg}     
rawReads2=${4:?msg}

echo
echo "Libraries:"
echo $READ1
echo $READ2


# output folder (root name)
RR=${5:?msg}

#sample name
SNAME=${6:?msg}

echo
echo "Sample name:"
echo $SNAME
echo


# output file
OUTF1=${SNAME}_check_read-names_OUT.txt
OUTF2=${SNAME}_check_index_seq_OUT.txt

rm -f $RR/$OUTF1
rm -f $RR/$OUTF2






######### Check libraries

DNAFILE_A=$READ1                 
DNAFILE_B=$READ2
   

#### get instrument name, run number, and flow-cell ID
zcat $AA/$DNAFILE_A | sed -n 1~4p | cut -d ' ' -f1 | cut -d ':' -f1-3 | sort | uniq -c | sort -k1 -nr > $SNIC_TMP/runs_list_A.txt
zcat $AA/$DNAFILE_B | sed -n 1~4p | cut -d ' ' -f1 | cut -d ':' -f1-3 | sort | uniq -c | sort -k1 -nr > $SNIC_TMP/runs_list_B.txt

#### get number of runs (based on sequencing instrument name and run id)
runs_list_aux1=$(cat $SNIC_TMP/runs_list_A.txt | wc -l )   
runs_list_aux2=$(cat $SNIC_TMP/runs_list_B.txt | wc -l )   
   
#### get number of reads in library
reads_aux1=$(zcat $AA/$DNAFILE_A | sed -n 1~4p | wc -l)  
reads_aux2=$(zcat $AA/$DNAFILE_B | sed -n 1~4p | wc -l)
   
#### get number of unique headers
reads_uniq1=$(zcat $AA/$DNAFILE_A | sed -n 1~4p | sort | uniq | wc -l)  
reads_uniq2=$(zcat $AA/$DNAFILE_B | sed -n 1~4p | sort | uniq | wc -l)  
   
#### difference
aux1=2
diff_aux1="$(awk '{print $1-$2}' <<<"$reads_aux1 $reads_uniq1")"
diff_aux2="$(awk '{print $1-$2}' <<<"$reads_aux2 $reads_uniq2")"   
   
#### get list of index sequences 
zcat $AA/$DNAFILE_A | grep '^@' | cut -d ' ' -f2 | cut -d ':' -f4 | sort | uniq -c | sort -k1 -nr > $SNIC_TMP/index_seq_1.txt
zcat $AA/$DNAFILE_B | grep '^@' | cut -d ' ' -f2 | cut -d ':' -f4 | sort | uniq -c | sort -k1 -nr > $SNIC_TMP/index_seq_2.txt
   
#### get total number of index sequences and index with most hits
index_top_hits1=$(cat $SNIC_TMP/index_seq_1.txt | head -1 | cut -d ' ' -f2)  
nb_index1=$(cat $SNIC_TMP/index_seq_1.txt | wc -l)  

index_top_hits2=$(cat $SNIC_TMP/index_seq_2.txt | head -1 | cut -d ' ' -f2)  
nb_index2=$(cat $SNIC_TMP/index_seq_2.txt | wc -l) 
   
 


#### save results to file
echo "Sample_name" "Library_name" "No_runs_R1" "No_runs_R2" "Reads_R1" "Reads_R2" "Unique_read_names_R1" "Unique_read_names_R2" "difference_R1" "difference_R2"
echo "Sample_name" "Library_name" "No_runs_R1" "No_runs_R2" "Reads_R1" "Reads_R2" "Unique_read_names_R1" "Unique_read_names_R2" "difference_R1" "difference_R2" > $RR/$OUTF1

echo ${SNAME} ${rawReads1%_R1_001.trim.fastq.gz} $runs_list_aux1 $runs_list_aux2 $reads_aux1 $reads_aux2 $reads_uniq1 $reads_uniq2 $diff_aux1 $diff_aux2
echo ${SNAME} ${rawReads1%_R1_001.trim.fastq.gz} $runs_list_aux1 $runs_list_aux2 $reads_aux1 $reads_aux2 $reads_uniq1 $reads_uniq2 $diff_aux1 $diff_aux2 >> $RR/$OUTF1


echo
echo "Sample_name" "Library_name" "Number_of_indexes_R1" "Number_of_indexes_R2" "Index_top_hits_R1" "Index_top_hits_R2"
echo "Sample_name" "Library_name" "Number_of_indexes_R1" "Number_of_indexes_R2" "Index_top_hits_R1" "Index_top_hits_R2" > $RR/$OUTF2

echo ${SNAME} ${rawReads1%_R1_001.trim.fastq.gz} $nb_index1 $nb_index2 $index_top_hits1 $index_top_hits2
echo ${SNAME} ${rawReads1%_R1_001.trim.fastq.gz} $nb_index1 $nb_index2 $index_top_hits1 $index_top_hits2 >> $RR/$OUTF2

echo
echo "Run list (R1):"
cat $SNIC_TMP/runs_list_A.txt

echo >> $RR/$OUTF1
echo "Run list (R1):" >> $RR/$OUTF1
cat $SNIC_TMP/runs_list_A.txt >> $RR/$OUTF1

echo
echo "Run list (R2):"
cat $SNIC_TMP/runs_list_B.txt

echo >> $RR/$OUTF1
echo "Run list (R2):" >> $RR/$OUTF1
cat $SNIC_TMP/runs_list_B.txt >> $RR/$OUTF1

echo
echo "Index list (R1):"
cat $SNIC_TMP/index_seq_1.txt

echo >> $RR/$OUTF2
echo "Index list (R1):" >> $RR/$OUTF2
cat $SNIC_TMP/index_seq_1.txt >> $RR/$OUTF2

echo
echo "Index list (R2):"
cat $SNIC_TMP/index_seq_2.txt

echo >> $RR/$OUTF2
echo "Index list (R2):" >> $RR/$OUTF2
cat $SNIC_TMP/index_seq_2.txt >> $RR/$OUTF2




## clean scratch disk
rm -f $SNIC_TMP/*.txt






