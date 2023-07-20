#!/bin/bash -l
#SBATCH -J splitbylane
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 06:00:00
#SBATCH -A snic2021-22-291
#SBATCH -M snowy
#SBATCH --mail-user luis.leal@ebc.uu.se
#SBATCH --mail-type=FAIL
ulimit -c unlimited
set -eo pipefail


## Script takes a fastq file as input and produces a separate fastq for each sequencing lane.
## It is assumed that input fastq files have already been split by sequencing-machine/run/flow-cell, if required.




# load modules
module load bioinfo-tools


# input folder
AA=${1:?msg}		

echo
echo "Input folder:"
echo $AA

# output folder
RR=${2:?msg}			

#sample name
SNAME=${3:?msg}

echo
echo "Sample:"
echo $SNAME

# Folder containing results from check-libs
ISEQ=${4:?msg}

# current path
SRCDIR=$(pwd)                                                                    




###### create list of machine/runs/flowcells
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




###  Split Reads For Different Flowcell Lanes
###  Read name assumed to be formatted as follows
### 	@NS500198:189:H3VJ5BGXY:1:11101:1168:18287
###  with flow cell lane indicated in the 4th field (1, in the example above)


for j in `seq 1 1 $N_INSTR`; do 

  echo 
  echo "Processing run "${j} "of" $N_INSTR

  gunzip -c $AA/${SNAME}-RUN${j}_R1.fastq.gz > $SNIC_TMP/rawReads1
  gunzip -c $AA/${SNAME}-RUN${j}_R2.fastq.gz > $SNIC_TMP/rawReads2
	
	cd $SNIC_TMP
	
	for i in `seq 1 1 2`; do
	
	  rm -f $SNIC_TMP/lane.*.fastq
	
    # cerate fastq file for each lane
    awk 'BEGIN {FS = ":"} {lane=$4 ; print > "lane."lane".fastq" ; for (i = 1; i <= 3; i++) {getline ; print > "lane."lane".fastq"}}' < rawReads${i}
		
	  # get list of new fastq files (one for each lane) 
	  find lane.*.fastq -maxdepth 0 -type f > __aux1.txt
		
		# copy results to output folder
		while read -r LINE                                                                      
		do
		   aux_n1=${LINE%.fastq}
		   gzip -c $SNIC_TMP/$LINE > $RR/${SNAME}-RUN${j}-${aux_n1}_R${i}.paired.fastq.gz
		done < $SNIC_TMP/__aux1.txt

    done
	
	rm -f $SNIC_TMP/rawReads*
	rm -f $SNIC_TMP/*fastq*
	rm -f $SNIC_TMP/__aux1.txt

done


rm -f $SNIC_TMP/_AUX_runs.txt







