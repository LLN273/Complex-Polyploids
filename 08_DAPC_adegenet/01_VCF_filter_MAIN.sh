#!/bin/bash
#
#SBATCH -J VCFf
#SBATCH -p devcore
#SBATCH -n 1
#SBATCH -t 01:00:00
#SBATCH -A snic2022-22-909
#SBATCH --mail-user luis.leal@ebc.uu.se
#SBATCH --mail-type=FAIL
ulimit -c unlimited

 
#### Prepare VCF file for DAPC/adegenet analysis


#### load modules
module load bioinfo-tools
module load GATK/4.2.0.0


#### paths and folder names

#remember initial path
SRCDIR_INI=$(pwd)                                           	 

#input folder
AA=/crex1/proj/snic2020-6-184/private/Luis/P06_birch_exome_IGA-Sweden_2019/24_Data_for_PCA_and_Admixture/BWA_unmasked

#output folder
RR=/crex1/proj/snic2020-6-184/private/Luis/P06_birch_exome_IGA-Sweden_2019/31_DAPC_adegenet/BWA_unmasked

# input/output subfolder 1
RRsub=exome_WGS

# samples to exclude
STE=/crex1/proj/snic2017-7-149/private/Luis/P06_birch_exome_IGA-Sweden_2019/31_DAPC_adegenet/00_samples_exclude_REDUX.args	

# VCF list
VCFin[1]=MASTER_SingleRuns_ALLsamples_GOODonly_NANA_SILVER_WHITE_HUMILIS_PLATYPHYLLA-filtered_tranche_2_PASS-BIALLELIC_FINAL_NoPrivateAlleles_NONCODING_LD_CLEAN.vcf.gz                  	 



for k in `seq 1 1 1`; do 			# dataset (root)	

   GO_AA=${AA}/${RRsub}
   GO_RR=${RR}/${RRsub}
   mkdir -p $GO_RR

   # VCF file
   GO_VCF=${VCFin[$k]}

   ## output file (root)
   OUTFile=${GO_VCF%.vcf.gz}_100GENOT_REDUX.vcf.gz

   echo $GO_VCF
   echo $GO_AA
   echo 

   ####### include only variants if genotyping has been done in all samples.   
   echo
   echo
   echo "Filter per nocall-fraction"
   echo

   gatk --java-options "-Xmx5g" \
     SelectVariants \
     -V $GO_AA/$GO_VCF \
     --max-nocall-fraction 0.0 \
     --exclude-sample-name $STE \
     --remove-unused-alternates true \
     --exclude-non-variants true \
     --exclude-filtered true \
     --verbosity ERROR \
     -O $GO_RR/$OUTFile

   # --max-nocall-fraction 	Maximum fraction of samples with no-call genotypes (default: 1.0)
   # --remove-unused-alternates	Remove alternate alleles not present in any genotypes (default: false)
   # --exclude-non-variants	Don't include non-variant sites (default: false)
   # --exclude-filtered		Don't include filtered sites. If this flag is enabled, sites that have been marked as filtered (i.e. have anything other than `.` or `PASS` in the FILTER field) will be excluded from the output. (default: false)


done

