#!/bin/bash -l
#SBATCH -J FqTOuBAM
#SBATCH -p core
#SBATCH -n 2
#SBATCH -t 02:00:00
#SBATCH -A snic2021-22-291
#SBATCH -M snowy
#SBATCH --mail-user luis.leal@ebc.uu.se
#SBATCH --mail-type=FAIL
ulimit -c unlimited
#set -eo pipefail



# Script converts fastq files to unmapped BAM (uBAM) format.
# Adds read-group tags.


# load modules
module load bioinfo-tools
module load picard/2.10.3



#read input folder name
AA=${1:?msg}                            

# output folder
RR=${2:?msg}			

#sample name
SNAME=${3:?msg}

# Run info
RUNinfo=${4:?msg}

echo
echo "Sample:"
echo $SNAME

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



##################################################################################### Create uBAM for each run and lane
####
#### MORE INFO:
#### https://gatkforums.broadinstitute.org/gatk/discussion/6483
#### https://gatkforums.broadinstitute.org/gatk/discussion/6484#latest#top
#### https://gatkforums.broadinstitute.org/gatk/discussion/6472/read-groups
#### https://broadinstitute.github.io/picard/command-line-overview.html#FastqToSam


while read -r LINE                                                                      
do 

   run_number=$(echo $LINE | cut -f2 -d' ')

   ### get library info
   ### example of line format: @NS500198:23:H0VUMAGXX:3:11401:11643:1016 1:N:0:5

   flowcell_ID=$(echo $LINE | cut -f4 -d' ' | cut -f3 -d':')				              # flowcell id 						(H0VUMAGXX)
   flowcell_lane=$(echo $LINE | cut -f4 -d' ' | cut -f4 -d':')				            # flowcell lane						(3)
   RG_ID=${flowcell_ID}.${SNAME}.${flowcell_lane}										              # Because some barcodes are troublesome characters (eg '<' or ':'), we will use the sample name as the barcode
   SM=$SNAME																			                                # sample name
   LB=${SNAME}_LIB1																		                            # Library identifier (during DNA preparation) >> use a different lib number (eg LIB2) if two DNA libraries are prepared for the same sample 
   PU=${flowcell_ID}.${SNAME}.${flowcell_lane}											              # Because some barcodes are troublesome characters (eg '<' or ':'), we will use the sample name as the barcode
   PL=ILLUMINA																			                              # platform used to generate the data
   SC=IGA_ITALY																		                                # sequencing center (all samples sequenced at IGA Tech, ITALY; 
																						                                      # see for example https://www.ebi.ac.uk/ena/data/view/ERS1792007)

                                                            
   ######## file names  
   READ1=${SNAME}-${run_number}-lane.${flowcell_lane}_R1.paired.fastq.gz
   READ2=${SNAME}-${run_number}-lane.${flowcell_lane}_R2.paired.fastq.gz
   
   echo
   echo "Input file 1:" $READ1       
   echo "Input file 2:" $READ2
   echo "run number:" $run_number
   echo "flowcell id:" $flowcell_ID
   echo "flowcell_lane:" $flowcell_lane
   echo "read group identifier:" $RG_ID
   echo "sample name:" $SM
   echo "Library identifier:" $LB
   echo "platform unit:" $PU
   echo 

   ######## Convert FASTQ to unmapped BAM (uBAM) and add read group information

   java -Xmx10G -Djava.io.tmpdir=$SNIC_TMP -jar $PICARD_HOME/picard.jar FastqToSam \
                                                FASTQ=$AA/$READ1 \
                                                FASTQ2=$AA/$READ2 \
                                                OUTPUT=$SNIC_TMP/${SNAME}-${run_number}-lane.${flowcell_lane}-fastq2ubam.bam \
                                                READ_GROUP_NAME=${RG_ID} \
                                                SAMPLE_NAME=${SM} \
                                                LIBRARY_NAME=${LB} \
                                                PLATFORM_UNIT=${PU} \
                                                PLATFORM=${PL} \
                                                SEQUENCING_CENTER=${SC} \
											                          SORT_ORDER=queryname
   
   

   ########## copy uBAM to results folder
   rsync -ah $SNIC_TMP/${SNAME}-${run_number}-lane.${flowcell_lane}-fastq2ubam.bam $RR/${SNAME}-${run_number}-lane.${flowcell_lane}-fastq2ubam.bam
   

done < $SNIC_TMP/__INFO_aux.txt



# notes (picard.jar FastqToSam): 
# FASTQ (File)					Input fastq file (optionally gzipped) for single end data, or first read in paired end data. Required.
# FASTQ2 (File)					Input fastq file (optionally gzipped) for the second read of paired end data. Default value: null.
# OUTPUT (File)					Output SAM/BAM file. Required.
# READ_GROUP_NAME (String)		Read group name Default value: A. This option can be set to 'null' to clear the default value.
# SAMPLE_NAME (String)			Sample name to insert into the read group header Required.
# LIBRARY_NAME (String)			The library name to place into the LB attribute in the read group header Default value: null.
# PLATFORM_UNIT (String)		The platform unit (often run_barcode.lane) to insert into the read group header Default value: null.
# PLATFORM (String)				The platform type (e.g. illumina, solid) to insert into the read group header Default value: null.
# SEQUENCING_CENTER (String)	The sequencing center from which the data originated Default value: null.
# SORT_ORDER (SortOrder)		The sort order for the output sam/bam file. Default value: queryname. 
#								This option can be set to 'null' to clear the default value. Possible values: {unsorted, queryname, coordinate, duplicate, unknown}





## clean scratch disk
cd $SNIC_TMP/
rm -f $SNIC_TMP/*.bam
rm -f $SNIC_TMP/__INFO_aux.txt









