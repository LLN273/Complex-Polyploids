#!/bin/bash -l
#SBATCH -J SplitInst
#SBATCH -p core
#SBATCH -n 6
#SBATCH -t 06:00:00
#SBATCH -A snic2021-22-291
#SBATCH -M snowy
#SBATCH --mail-user luis.leal@ebc.uu.se
#SBATCH --mail-type=FAIL
ulimit -c unlimited
set -eo pipefail


## Takes a fastq file as input and produces a separate fastq for each sequencing-machine.


# load modules
module load bioinfo-tools
#module load FastQC/0.11.5
#module load MultiQC/1.2




######## PATHS and filenames

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


# Folder containing results from check-libs
ISEQ=${7:?msg}








######## create list of machine/runs/flowcells
cat $ISEQ/${SNAME}_check_read-names_OUT.txt | grep '@' | cut -f2 -d' ' | sort | uniq > $SNIC_TMP/_AUX_runs.txt

echo
echo "list of sequencing machine/runs/flowcells:"
cat $SNIC_TMP/_AUX_runs.txt
echo





######## read list of sequencing machine/runs/flowcells

i=1

while read -r LINE                                                                      
do
   SEQ_INSTR[$i]=$LINE
   let "i+=1"
done < $SNIC_TMP/_AUX_runs.txt


#Number of sequencing machine/runs/flowcells
N_INSTR=${#SEQ_INSTR[@]}

#echo
#echo $N_INSTR




########  Split library by sequencing machine/run/flowcell

for i in `seq 1 1 $N_INSTR`; do   

   ##### grep reads on interest
   echo 
   echo "Searching for run"${i} ":" ${SEQ_INSTR[$i]}
   
   zcat $READ1 | grep -A 3 "${SEQ_INSTR[$i]}" > $RR/${SNAME}-RUN${i}_R1.fastq
   zcat $READ2 | grep -A 3 "${SEQ_INSTR[$i]}" > $RR/${SNAME}-RUN${i}_R2.fastq

   gzip -c $RR/${SNAME}-RUN${i}_R1.fastq > $RR/${SNAME}-RUN${i}_R1.fastq.gz
   gzip -c $RR/${SNAME}-RUN${i}_R2.fastq > $RR/${SNAME}-RUN${i}_R2.fastq.gz


done







######## clean scratch disk
cd $SNIC_TMP
rm -f $SNIC_TMP/*gz
rm -f $SNIC_TMP/_AUX_runs.txt



