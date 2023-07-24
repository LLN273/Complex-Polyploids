#!/bin/bash
#
#SBATCH -J AllelicC
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 1:00:00
#SBATCH -A snic2022-22-589
#SBATCH -M snowy
#SBATCH --mail-user luis.leal@ebc.uu.se
#SBATCH --mail-type=FAIL
ulimit -c unlimited



### Script used to collect reference and alternate allele counts at biallelic sites




module load bioinfo-tools
module load GATK/4.2.0.0





############### input data files and paths

#remember initial path
SRCDIR_INI=$(pwd)

#sample name
SampleName=${1:?msg}	

#BAM file 
InBAM=${2:?msg}	

#VCF file 
InVCF=${3:?msg}	

# input folder (BAM files)
BB=${4:?msg}

# input folder (VCF files)
AA=${5:?msg}

# output folder
RR=${6:?msg}

# reference genome
refGenome=${7:?msg}			

# reference genome folder
refGenomefolder=${8:?msg}




echo
echo "Sample:" $SampleName
echo "BAM file:" $InBAM
echo "VCF file:" $InVCF
echo "Input folder (BAM):" $BB
echo "Input folder (VCF):" $AA
echo "Output folder:" $RR
echo
echo



##### Move data files to scratch disk 
cd $SNIC_TMP
rsync -ah $refGenomefolder/${refGenome%fa*}* .
zcat $AA/$InVCF > $SNIC_TMP/_biallelic.vcf







##### CollectAllelicCounts

gatk --java-options "-Xmx5g" \
          CollectAllelicCounts \
          -I $BB/$InBAM \
          -R $SNIC_TMP/${refGenome} \
          -L $SNIC_TMP/_biallelic.vcf \
          --sites-only-vcf-output true \
          --verbosity ERROR \
          -O $RR/${SampleName}.allelicCounts.tsv


# Notes:
# Overview
#Collects reference and alternate allele counts at specified sites. The alt count is defined as the total count minus the ref count, and the alt nucleotide is defined as the non-ref base with the highest count, with ties broken by the order of the bases in AllelicCountCollector#BASES. Only reads that pass the specified read filters and bases that exceed the specified minimum-base-quality will be counted.

#Input
#    SAM format read data
#    Reference FASTA file
#    Sites at which allelic counts will be collected

#Output
#    Allelic-counts file. This is a tab-separated values (TSV) file with a SAM-style header containing a read group sample name, a sequence dictionary, a row specifying the column headers contained in AllelicCountCollection.AllelicCountTableColumn, and the corresponding entry rows.








# Clean scratch disk
rm -f $SNIC_TMP/${refGenome%fa}*
rm -f $SNIC_TMP/_*


