#!/bin/bash
#
#SBATCH -J MarkDuplicates
#SBATCH -p core 
#SBATCH -n 1
#SBATCH -t 02:00:00
#SBATCH -A snic2021-22-291
#SBATCH -M snowy
#SBATCH --mail-user luis.leal@ebc.uu.se
#SBATCH --mail-type=FAIL
ulimit -c unlimited
#set -eo pipefail



## Script used to identify read pairs that are likely to have originated from duplicates of the same original DNA fragments through some artifactual processes. 
## This step eliminates PCR duplicates (arising from library preparation) across all lanes in addition to optical duplicates (which are by definition only per-lane).


# Load software
module load bioinfo-tools
module load picard/2.10.3
module load samtools/1.12




########### input data


#read input folder name
AA=${1:?msg}                            

# output folder
RR=${2:?msg}			

#sample name
SNAME=${3:?msg}

echo
echo "Sample:" $SNAME
echo

echo "Output folder:" $RR
echo

# Run info
RUNinfo=${4:?msg}

# current path
SRCDIR_INI=$(pwd)   




### Determine number of runs & lanes for this sample
N_FILES="$(cat $RUNinfo | grep -w ${SNAME} | wc -l)"

### get info for each run & lane
cat $RUNinfo | grep -w ${SNAME} > $SNIC_TMP/__INFO_aux.txt

echo
echo "Runs and lanes info"
cat $SNIC_TMP/__INFO_aux.txt
echo





### Get list of BAM files to be merged
bamlist=$(for f in $AA/${SNAME}*-MergeBamAlignment.bam; do echo -n "I=$f " ; done)

## List of bam files being processed
echo $bamlist
echo



########## Merge SAM files associated to the same sample/run

echo
java -Xmx5G -jar $PICARD_HOME/picard.jar MergeSamFiles \
	                                       $bamlist \
										                     SORT_ORDER=coordinate \
	                                       O=$SNIC_TMP/${SNAME}.bam




########## Mark duplicates
echo
cd $SNIC_TMP
mkdir __RD_aux1

java -Xmx5G -jar $PICARD_HOME/picard.jar MarkDuplicates \
                                         INPUT=$SNIC_TMP/${SNAME}.bam \
                                         OUTPUT=$SNIC_TMP/${SNAME}_marked_duplicates.bam \
                                         M=$SNIC_TMP/${SNAME}_marked_dup_metrics.txt \
                                         TMP_DIR=$SNIC_TMP/__RD_aux1 
   



########## index BAM file
echo
java -Xmx5G -jar $PICARD_HOME/picard.jar BuildBamIndex \
                                         INPUT=$SNIC_TMP/${SNAME}_marked_duplicates.bam
		
		
		
		
		
								 
							 
########## copy results to output folder
echo
rsync -ah $SNIC_TMP/${SNAME}_marked_duplicates.bam $RR/${SNAME}_marked_duplicates.bam
rsync -ah $SNIC_TMP/${SNAME}_marked_duplicates.bai $RR/${SNAME}_marked_duplicates.bai
rsync -ah $SNIC_TMP/${SNAME}_marked_dup_metrics.txt $RR/${SNAME}_marked_dup_metrics.txt




   
########### glean statistics about BAM file
echo
samtools flagstat $RR/${SNAME}_marked_duplicates.bam  > $RR/${SNAME}_marked_duplicates.flagstat



## clean scratch disk
cd $SNIC_TMP
rm -rf $SNIC_TMP/__RD_aux1
rm -f $SNIC_TMP/*.bam
rm -f $SNIC_TMP/*.bai
rm -f $SNIC_TMP/*.txt


