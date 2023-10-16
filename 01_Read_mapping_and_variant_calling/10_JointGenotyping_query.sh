#!/bin/bash
#
#SBATCH -J GATK_JG
#SBATCH -p core 
#SBATCH -n 3
#SBATCH -t 12:00:00
#SBATCH -A snic2022-22-589
#SBATCH -M snowy
#SBATCH --mail-user luis.leal@ebc.uu.se
#SBATCH --mail-type=FAIL
ulimit -c unlimited


## Script used to perform variant calling using GATK; second step: GenotypeGVCFs


module load bioinfo-tools
module load GATK/4.2.0.0
#module load samtools/1.6


############### input data files and paths

#remember initial path
SRCDIR_INI=$(pwd)   

# path to GenomicsDB workspace
AA=${1:?msg}

# GenomicsDBImport folders
IN_FILE=${2:?msg}
	
#output folder
RR=${3:?msg}

# Output file name
outFile_root=${4:?msg}

# reference genome
refGenome=${5:?msg}			

# reference genome folder
refGenomefolder=${6:?msg}

# Intervals list (BED format) exome data 
ANNF=${7:?msg}

# Contig list
GO_CONTIGlist=${8:?msg}


echo
echo "Input folder:" $AA
echo "Output folder:" $RR
echo "Input file (merged gVCF):" $IN_FILE
echo "Output filename:" ${outFile_root}.vcf.gz
echo "Annotation file:" $ANNF
echo "Contig list:" $GO_CONTIGlist
echo "Reference genome:" $refGenome
echo "Reference genome folder:" $refGenomefolder
echo




########### copy reference genome to scratch disk

cd $SNIC_TMP
rsync -ah $refGenomefolder/${refGenome%.fa*}* .





################# Create intervals list (target genes)
################# File must have ".intervals" extension (gatk crashes otherwise) [see: https://gatk.broadinstitute.org/hc/en-us/articles/360035531852?id=11009]

### Convert position 0 to position 1 (GATK does not accept position 0)
echo
awk 'BEGIN {OFS="\t"} $2 == "0" {$2 = "1"}; 1'  $ANNF > $SNIC_TMP/__aux0.txt

### Convert intervals list to chr1:100-200 format
echo
cat $SNIC_TMP/__aux0.txt | cut -f1 > $SNIC_TMP/__aux1.txt
cat $SNIC_TMP/__aux0.txt | cut -f2 > $SNIC_TMP/__aux2.txt
cat $SNIC_TMP/__aux0.txt | cut -f3 > $SNIC_TMP/__aux3.txt
paste -d':' $SNIC_TMP/__aux1.txt $SNIC_TMP/__aux2.txt > $SNIC_TMP/__auxA.txt
paste -d'-' $SNIC_TMP/__auxA.txt $SNIC_TMP/__aux3.txt > $SNIC_TMP/__targets_aux.interval_list


# select subset of contigs
cat $GO_CONTIGlist | sort > $SNIC_TMP/__GO_CONTIGlist_sort.txt
grep -f $SNIC_TMP/__GO_CONTIGlist_sort.txt $SNIC_TMP/__targets_aux.interval_list > $SNIC_TMP/__targets.intervals







################# Perform joint genotyping on GenomicsDB workspace created with GenomicsDBImport

mkdir -p $SNIC_TMP/_TEMPgatk

gatk --java-options "-Xmx15g" \
    GenotypeGVCFs \
   --reference $SNIC_TMP/$refGenome \
   --variant gendb://${AA}/db_${IN_FILE}  \
   --intervals $SNIC_TMP/__targets.intervals \
   --include-non-variant-sites \
   --output $RR/${outFile_root}.vcf.gz \
   --tmp-dir $SNIC_TMP/_TEMPgatk \
   --verbosity ERROR
 







#### clean scratch disk
rm -r $SNIC_TMP/*.fa*
rm -r $SNIC_TMP/_*
rm -rf $SNIC_TMP/_TEMPgatk
rm -f $SNIC_TMP/__GO_CONTIGlist_sort.txt







