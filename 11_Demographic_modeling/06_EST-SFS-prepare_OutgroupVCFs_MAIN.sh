#!/bin/bash
#
#SBATCH -J preEST
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 03:00:00
#SBATCH -A snic2022-22-909
#SBATCH -M snowy
#SBATCH --mail-user luis.leal@ebc.uu.se
#SBATCH --mail-type=FAIL
ulimit -c unlimited


# Script used to prepare data required to run EST-SFS (step 1)
# Cleans VCFs associated to the three outgroup species (B. nana, B. populifolia, and B. occidentalis)



##### load modules
module load bioinfo-tools
module load GATK/4.2.0.0
module load htslib/1.9
module load tabix/0.2.6
module load picard/2.10.3



####################################### paths

#remember current path
SRCDIR=$(pwd)  

### Folder containing VCF files associated to outgroup species
BB=/crex1/proj/snic2020-6-184/nobackup/private/Luis/P1_DemSim/15d_Individual_VCFs_SG/BWA_unmasked

### input/output subfolder 1
RRsub=exome_WGS

## VCF files associated to outgroup species
VCF_OUTG[1]=Nana-Finland-Enontekio.A002-2_snps-CLEAN.vcf.gz				# B. nana
VCF_OUTG[2]=Populifolia-Canada.A005-11_snps-CLEAN.vcf.gz				# B. populifolia
VCF_OUTG[3]=Occidentalis-Canada-Alberta.A009-17_snps-CLEAN.vcf.gz		# B. occidentalis

### output folder
RR=/crex1/proj/snic2020-6-184/private/Luis/P06_birch_exome_IGA-Sweden_2019/32b_fastsimcoal2_ancestral_state/BWA_unmasked





############## 1. Clean VCF files (outgroup species only)
echo
echo "1A. Cleaning VCF files (outgroup species only)."

GO_RR=$RR/$RRsub
mkdir -p $GO_RR

for k in `seq 1 1 3`; do 			# dataset (root)
   echo
   echo ${k}...
   vcf_aux=${VCF_OUTG[$k]}
   zcat $BB/$vcf_aux | grep '^#' > $SNIC_TMP/${vcf_aux%.vcf.gz}_READY.vcf
   zcat $BB/$vcf_aux | grep -v '^#' > $SNIC_TMP/__aux_pass.txt
   awk '($4 == "A" || $4 == "C" || $4 == "G" || $4 == "T") { print $0 }' $SNIC_TMP/__aux_pass.txt > $SNIC_TMP/__aux_pass1.txt
   awk '$5 != "*" { print $0 }' $SNIC_TMP/__aux_pass1.txt > $SNIC_TMP/__aux_pass2.txt
   cat $SNIC_TMP/__aux_pass2.txt | grep PASS >> $SNIC_TMP/${vcf_aux%.vcf.gz}_READY.vcf
done

################# Compress output file
echo
echo "1B. Compress VCF files"

for k in `seq 1 1 3`; do 			# dataset (root)
echo
   echo ${k}...
   vcf_aux=${VCF_OUTG[$k]}
   bgzip -c $SNIC_TMP/${vcf_aux%.vcf.gz}_READY.vcf > $GO_RR/${vcf_aux%.vcf.gz}_READY.vcf.gz
done


################# Index vcf file
echo
echo "1C. Index VCF file"

for k in `seq 1 1 3`; do 			# dataset (root)
echo
   echo ${k}...
   vcf_aux=${VCF_OUTG[$k]}
   gatk --java-options "-Xmx5G" \
	            IndexFeatureFile \
	            -I $GO_RR/${vcf_aux%.vcf.gz}_READY.vcf.gz
done

 
## clean auxiliary files
rm -f $SNIC_TMP/*

echo
echo 'Done!'






