#!/bin/bash
#
#SBATCH -J GATK_HC
#SBATCH -p core 
#SBATCH -n 3
#SBATCH -t 24:00:00
#SBATCH -A snic2022-22-589
#SBATCH -M snowy
#SBATCH --mail-user luis.leal@ebc.uu.se
#SBATCH --mail-type=FAIL
ulimit -c unlimited



## Script used to perform variant calling using GATK; first step: HaplotypeCaller


module load bioinfo-tools
module load GATK/4.2.0.0
module load samtools/1.12



############### input data files and paths

#remember initial path
SRCDIR_INI=$(pwd)   

#read BAM file (with path)
READ1=${1:?msg}

#read BAM file (without path)
BAM1=${2:?msg}		

#output file (root name)
RR=${3:?msg}

#sample name
SNAME=${4:?msg}

# reference genome
refGenome=${5:?msg}			

# reference genome folder
refGenomefolder=${6:?msg}

# ploidy
SPLOIDY=${7:?msg}

# Intervals list (BED format) exome data 
ANNF=${8:?msg}





echo
echo "Input samples (BAM):" $READ1
echo "Sample name:" $SNAME
echo "Reference genome:" $refGenome
echo "Reference genome folder:" $refGenomefolder
echo "Output folder:" $RR
echo "Sample ploidy:" $SPLOIDY
echo "Intervals list:" $ANNF
echo




############### copy files to scratch disk

cd $SNIC_TMP
rsync -ah $refGenomefolder/${refGenome%.fa*}* .




############### Create intervals list (target genes)

### Convert position 0 to position 1 (GATK does not accept position 0)
echo
awk 'BEGIN {OFS="\t"} $2 == "0" {$2 = "1"}; 1'  $ANNF > $SNIC_TMP/__aux0.txt

### Convert intervals list to chr1:100-200 format
echo
cat $SNIC_TMP/__aux0.txt | cut -f1 > $SNIC_TMP/__aux1.txt
cat $SNIC_TMP/__aux0.txt | cut -f2 > $SNIC_TMP/__aux2.txt
cat $SNIC_TMP/__aux0.txt | cut -f3 > $SNIC_TMP/__aux3.txt
paste -d':' $SNIC_TMP/__aux1.txt $SNIC_TMP/__aux2.txt > $SNIC_TMP/__auxA.txt
paste -d'-' $SNIC_TMP/__auxA.txt $SNIC_TMP/__aux3.txt > $SNIC_TMP/__targets_intervals.list






##### Run HaplotypeCaller variant calling in GVCF mode



gatk --java-options "-Xmx16G" \
	   HaplotypeCaller \
	   --emit-ref-confidence GVCF \
     --sample-ploidy $SPLOIDY \
	   --minimum-mapping-quality 10 \
	   --create-output-variant-index true \
     --reference $SNIC_TMP/${refGenome} \
     --intervals $SNIC_TMP/__targets_intervals.list \
     --input $READ1 \
     --output $SNIC_TMP/${SNAME}-individualHaplotypeCaller.g.vcf






##### Copy output files back to Uppmax folder 
cd $RR
rsync -ah $SNIC_TMP/${SNAME}-individualHaplotypeCaller.g.vcf .
rsync -ah $SNIC_TMP/${SNAME}-individualHaplotypeCaller.g.vcf.idx .



#### clean scratch disk
cd $SNIC_TMP
rm -r $SNIC_TMP/${SNAME}*
rm -f $SNIC_TMP/${refGenome%.fa*}*



