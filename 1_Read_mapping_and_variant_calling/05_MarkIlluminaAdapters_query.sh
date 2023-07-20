#!/bin/bash -l
#SBATCH -J MarkIlluminaAdapters
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 02:00:00
#SBATCH -A snic2021-22-291
#SBATCH -M snowy
#SBATCH --mail-user luis.leal@ebc.uu.se
#SBATCH --mail-type=FAIL
ulimit -c unlimited
#set -eo pipefail



## Script used to identify and mark position of Illumina adapters


# load modules
module load bioinfo-tools
module load picard/2.10.3


#read input folder name
AA=${1:?msg}                            

# output folder
RR=${2:?msg}			

#sample name
SNAME=${3:?msg}

echo
echo "Sample:"
echo $SNAME
echo

# Run info
RUNinfo=${4:?msg}

# current path
SRCDIR=$(pwd)                                                                    




### Determine number of runs & lanes for this sample
N_FILES="$(cat $RUNinfo | grep -w ${SNAME} | wc -l)"

### get info for each run & lane
cat $RUNinfo | grep -w ${SNAME} > $SNIC_TMP/__INFO_aux.txt

echo
echo "Runs and lanes info"
cat $SNIC_TMP/__INFO_aux.txt
echo





########## Mark adapter sequences for each run and lane
####
#### MORE INFO:
#### https://gatkforums.broadinstitute.org/gatk/discussion/6483
#### https://gatkforums.broadinstitute.org/gatk/discussion/6484#latest#top
## Extra notes:
## picard's MarkIlluminaAdapters produces two files. 
## (1) The metrics file, XXXXX_metrics.txt bins the number of tagged adapter bases versus the number of reads. 
## (2) The XXXXX_markilluminaadapters.bam file is identical to the input BAM, XXXXX_revertsam.bam, except reads with adapter sequences will be marked with a tag in XT:i:# format, where # denotes the 5' starting position of the adapter sequence. 
##     At least six bases are required to mark a sequence. Reads without adapter sequence remain untagged.


while read -r LINE                                                                      
do 

   cd $SNIC_TMP
   
   run_number=$(echo $LINE | cut -f2 -d' ')												                # run number
   flowcell_lane=$(echo $LINE | cut -f4 -d' ' | cut -f4 -d':')				            # flowcell lane	

   ########## file name
   BAM1=${SNAME}-${run_number}-lane.${flowcell_lane}-fastq2ubam.bam
   
   echo
   echo "Input file:" $BAM1
   echo "run number:" $run_number
   echo "flowcell_lane:" $flowcell_lane
   echo 



   ########## Mark adapter sequences using MarkIlluminaAdapters
   java -Xmx5G -Djava.io.tmpdir=$SNIC_TMP -jar $PICARD_HOME/picard.jar MarkIlluminaAdapters \
                                               INPUT=$AA/$BAM1 \
                                               OUTPUT=$SNIC_TMP/${SNAME}-${run_number}-lane.${flowcell_lane}-markilluminaadapters.bam \
                                               METRICS=$RR/${SNAME}-${run_number}-lane.${flowcell_lane}-markilluminaadapters_metrics.txt

   
   
   ########## copy uBAM to results folder
   rsync -ah $SNIC_TMP/${SNAME}-${run_number}-lane.${flowcell_lane}-markilluminaadapters.bam $RR/${SNAME}-${run_number}-lane.${flowcell_lane}-markilluminaadapters.bam
   
   
   ########## clean temporary folder
   #rm -rf $SNIC_TMP/__MIA_aux
   
   echo

done < $SNIC_TMP/__INFO_aux.txt                    
   
   
   
   
   
# notes (picard.jar FastqToSam): 
# INPUT (File)				Required.
# OUTPUT (File)				If output is not specified, just the metrics are generated Default value: null.
# METRICS (File)			Histogram showing counts of bases_clipped in how many reads Required.
# TMP_DIR					     Working folder used to process large files



## clean scratch disk
cd $SNIC_TMP/
rm -f $SNIC_TMP/*.bam
rm -f $SNIC_TMP/__INFO_aux.txt








