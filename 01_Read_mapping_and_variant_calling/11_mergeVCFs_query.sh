#!/bin/bash
#
#SBATCH -J MergeVcfs
#SBATCH -p core
#SBATCH -n 2
#SBATCH -t 06:00:00
#SBATCH -A snic2021-22-727
#SBATCH -M snowy
#SBATCH --mail-user luis.leal@ebc.uu.se
#SBATCH --mail-type=FAIL
ulimit -c unlimited



### Script used to combine all VCF files in one single file


module load bioinfo-tools
#module load GATK/4.2.0.0
module load picard/2.10.3



############### input data files and paths

#remember initial path
SRCDIR_INI=$(pwd)   

# input folder
AA=${1:?msg}
	
#output file
RR=${2:?msg}

# reference genome
refGenome=${3:?msg}			

# reference genome folder
refGenomefolder=${4:?msg}

# contig list
SPL=${5:?msg}

# Output file name
outFile_root=${6:?msg}




echo
echo "Sample list (BAM):" $SPL
echo "Input folder (BAM):" $AA
echo "Output folder:" $RR
echo "Output file (root):" db_${outFile_root}

echo "Reference genome:" $refGenome
echo "Reference genome folder:" $refGenomefolder
echo



########### copy reference genome to scratch disk

cd $SNIC_TMP
rsync -ah $refGenomefolder/${refGenome%.fa*}* .




########### Merge VCF files

java -Xmx8G -jar $PICARD_HOME/picard.jar MergeVcfs \
                                         REFERENCE_SEQUENCE=$SNIC_TMP/$refGenome \
                                         I=$SPL \
                                         O=$RR/${outFile_root}.vcf.gz




########### clean scratch disk
rm -r $SNIC_TMP/*.fa*





