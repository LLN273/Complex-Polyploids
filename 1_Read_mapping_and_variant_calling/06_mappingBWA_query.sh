#!/bin/bash
#
#SBATCH -J map_bwa
#SBATCH -p core 
#SBATCH -n 10
#SBATCH -t 02:00:00
#SBATCH -A snic2021-22-291
#SBATCH -M snowy
#SBATCH --mail-user luis.leal@ebc.uu.se
#SBATCH --mail-type=FAIL
ulimit -c unlimited
#set -eo pipefail



## Script used to convert bam to fastq, then pipe into BWA. Finally, merge original uBAM and BAM produced by BWA to generate a clean BAM


# Load software
module load bioinfo-tools
module load picard/2.10.3
module load bwa/0.7.17
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

# reference genome
refGenome=${4:?msg}			

# reference genome folder
refGenomefolder=${5:?msg}

echo "Reference genome:" $refGenome
echo "Reference genome folder:" $refGenomefolder
echo
echo "Output folder:" $RR

# Run info
RUNinfo=${6:?msg}

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



########### copy input files to scratch disk

cd $SNIC_TMP
rsync -ah $refGenomefolder/${refGenome%.fa*}* .







########### convert bam to fastq, then pipe into BWA. Finally, merge original uBAM and BAM produced by BWA to generate a clean BAM



while read -r LINE                                                                      
do 

   cd $SNIC_TMP

   run_number=$(echo $LINE | cut -f2 -d' ')												                # run number

   flowcell_lane=$(echo $LINE | cut -f4 -d' ' | cut -f4 -d':')				            # flowcell lane	

   ########## file name
   BAM1=${SNAME}-${run_number}-lane.${flowcell_lane}-markilluminaadapters.bam

   echo
   echo "Input file:" $BAM1
   echo "Input folder:" $AA
   echo "run number:" $run_number
   echo "flowcell_lane:" $flowcell_lane
   echo 

   
   cd $SNIC_TMP
   mkdir __BWA_aux1
   mkdir __BWA_aux3
   

   
   java -Xmx40G -jar $PICARD_HOME/picard.jar SamToFastq \
                                             INPUT=$AA/$BAM1 \
                                             FASTQ=/dev/stdout \
                                             CLIPPING_ATTRIBUTE=XT \
   											                     CLIPPING_ACTION=2 \
											                       INTERLEAVE=true \
											                       NON_PF=true \
                                             TMP_DIR=$SNIC_TMP/__BWA_aux1 | \
   bwa mem -t 10 \
   	       -M \
   		     -p \
   		     $SNIC_TMP/$refGenome \
           /dev/stdin | \
   java -Xmx40G -jar $PICARD_HOME/picard.jar  MergeBamAlignment \
                                              ALIGNED_BAM=/dev/stdin \
                                              UNMAPPED_BAM=$AA/$BAM1 \
                                              OUTPUT=$SNIC_TMP/${SNAME}-${run_number}-lane.${flowcell_lane}-MergeBamAlignment.bam \
                                              REFERENCE_SEQUENCE=$SNIC_TMP/$refGenome \
    										                      CREATE_INDEX=true \
											                        ADD_MATE_CIGAR=true \
                                              CLIP_ADAPTERS=false \
 											                        CLIP_OVERLAPPING_READS=true \
                                              INCLUDE_SECONDARY_ALIGNMENTS=true \
 											                        MAX_INSERTIONS_OR_DELETIONS=-1 \
                                              PRIMARY_ALIGNMENT_STRATEGY=MostDistant \
 											                        ATTRIBUTES_TO_RETAIN=XS \
                                              TMP_DIR=$SNIC_TMP/__BWA_aux3




   # copy results to output folder
   rsync -ah $SNIC_TMP/${SNAME}-${run_number}-lane.${flowcell_lane}-MergeBamAlignment.bam $RR/${SNAME}-${run_number}-lane.${flowcell_lane}-MergeBamAlignment.bam
   rsync -ah $SNIC_TMP/${SNAME}-${run_number}-lane.${flowcell_lane}-MergeBamAlignment.bai $RR/${SNAME}-${run_number}-lane.${flowcell_lane}-MergeBamAlignment.bai





   ########### glean statistics about BAM file
   samtools flagstat $RR/${SNAME}-${run_number}-lane.${flowcell_lane}-MergeBamAlignment.bam > $RR/${SNAME}-${run_number}-lane.${flowcell_lane}-MergeBamAlignment.flagstat


   ########## clean temporary folder
   rm -rf $SNIC_TMP/__BWA_aux1
   rm -rf $SNIC_TMP/__BWA_aux3

   
   echo


done < $SNIC_TMP/__INFO_aux.txt






####### clean scratch disk
cd $SNIC_TMP/
rm -f $SNIC_TMP/*.bam
rm -f $SNIC_TMP/*.sam
rm -f $SNIC_TMP/*.fa*
rm -f $SNIC_TMP/__INFO_aux.txt
rm -f $SNIC_TMP/${refGenome%.fa*}*











