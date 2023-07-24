#!/bin/bash
#
#SBATCH -J GATKextract
#SBATCH -p core
#SBATCH -n 3
#SBATCH -t 06:00:00
#SBATCH -A snic2022-22-589
#SBATCH -M snowy
#SBATCH --mail-user luis.leal@ebc.uu.se
#SBATCH --mail-type=FAIL
ulimit -c unlimited
set -eo pipefail



### Script used to create individual VCF files for each sample

module load bioinfo-tools
module load GATK/4.2.0.0





#remember initial path
SRCDIR_INI=$(pwd)                                           	


#input folder
InFile1_path=${1:?msg}		

#input vcf file (root name)
InFile1=${2:?msg}		    

#sample name
SampleName=${3:?msg}		

#output folder
RR=${4:?msg}	            

# reference genome
refGenome=${5:?msg}			

# reference genome folder
refGenomefolder=${6:?msg}

# interval list (BED format) exome data [target genes]
ANNF=${7:?msg}


#output file (root name)
outFile_root=${SampleName}-JG


echo
echo "Sample name:" $SampleName
echo "Input file (root):" $InFile1
echo "Input folder:" $InFile1_path
echo "Output folder:" $RR
echo "Interval list:" $ANNF
echo "Reference genome:" $refGenomefolder/$refGenome
echo




##### Move data files to scratch disk 
cd $SNIC_TMP
rsync -ah $refGenomefolder/${refGenome%fa*}* .
rsync -ah ${InFile1_path}/${InFile1}_snps-filtered_tranche_2_PHYLO.vcf.gz $SNIC_TMP/_VCF_snp.vcf.gz
rsync -ah ${InFile1_path}/${InFile1}_indels-filtered_tranche_2_PHYLO.vcf.gz $SNIC_TMP/_VCF_indels.vcf.gz
sleep 1
rsync -ah ${InFile1_path}/${InFile1}_snps-filtered_tranche_2_PHYLO.vcf.gz.tbi $SNIC_TMP/_VCF_snp.vcf.gz.tbi
rsync -ah ${InFile1_path}/${InFile1}_indels-filtered_tranche_2_PHYLO.vcf.gz.tbi $SNIC_TMP/_VCF_indels.vcf.gz.tbi



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
paste -d'-' $SNIC_TMP/__auxA.txt $SNIC_TMP/__aux3.txt > $SNIC_TMP/__targets.intervals








##### Run SelectVariants to extract individual samples from full VCF file


## snps
echo
gatk --java-options "-Xmx5g" \
     SelectVariants \
     -R $SNIC_TMP/${refGenome} \
     -V $SNIC_TMP/_VCF_snp.vcf.gz \
     --sample-name $SampleName \
     --intervals $SNIC_TMP/__targets.intervals \
     --remove-unused-alternates true \
     --verbosity ERROR \
     -O $RR/${outFile_root}_snps.vcf.gz




## indels
echo
gatk --java-options "-Xmx5g" \
     SelectVariants \
     -R $SNIC_TMP/${refGenome} \
     -V $SNIC_TMP/_VCF_indels.vcf.gz \
     --sample-name $SampleName \
     --intervals $SNIC_TMP/__targets.intervals \
     --select-type-to-include INDEL \
     --remove-unused-alternates true \
     --exclude-non-variants true \
     --exclude-filtered true \
     --verbosity ERROR \
     -O $RR/${outFile_root}_indels.vcf.gz




    

##### Cleaning
rm -f $SNIC_TMP/${refGenome%fa}*
rm -f $SNIC_TMP/_*



