#!/bin/bash
#
#SBATCH -J fVCF
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 03:00:00
#SBATCH -A snic2022-22-909
#SBATCH -M snowy
#SBATCH --mail-user luis.leal@ebc.uu.se
#SBATCH --mail-type=FAIL
ulimit -c unlimited


### Scrip used to prepare VCF file for STRUCTURE analysis
### Include only sites if genotyped in 100% of all samples
### Select 10000 random SNPS


##### load modules
module load bioinfo-tools
module load GATK/4.2.0.0
module load htslib/1.9
module load tabix/0.2.6


####################################### paths

#remember current path
SRCDIR=$(pwd)  

#input folder
AA=/crex1/proj/snic2020-6-184/private/Luis/P06_birch_exome_IGA-Sweden_2019/24_Data_for_PCA_and_Admixture/BWA_unmasked

# input/output subfolder 1
RRsub[1]=exome_WGS

# VCF list
VCF[1]=MASTER_SingleRuns_ALLsamples_GOODonly_NANA_SILVER_WHITE_HUMILIS_PLATYPHYLLA-filtered_tranche_2_PASS-BIALLELIC_FINAL_NoPrivateAlleles_NONCODING

# input extension
OEXT[1]="_ReducedPEND_ReducedSW_LD_CLEAN.vcf.gz"

# output folder
RR=/crex1/proj/snic2020-6-184/private/Luis/P06_birch_exome_IGA-Sweden_2019/30_STRUCTURE/BWA_unmasked

# Sample list (balanced datasets)
SL[1]=/crex1/proj/snic2017-7-149/private/Luis/P06_birch_exome_IGA-Sweden_2019/30_STRUCTURE/00_Sample_list_NANA_SILVER_WHITE_HUMILIS_PLATYPHYLLA_ReducedPEND_ReducedSW_SHORT.args    # 126 samples

# FLAGS reference genome (unmasked:1; masked:0)
FLAG_unmasked=1					


########### Reference genome (unmasked)
if [[ ${FLAG_unmasked} = "1" ]] ; then
   refGenome=Koivu_gapfilled_scaffolds.fasta
   refGenomefolder=/crex1/proj/snic2017-7-149/private/Luis/birch_reference_2017
fi



for k in `seq 1 1 1`; do 			# dataset (root)	

   GO_AA=${AA}/${RRsub[$k]}
   GO_RR=${RR}/${RRsub[$k]}
   GO_SL=${SL[$k]}
   GO_FILE=${VCF[$k]}${OEXT[$k]}
   mkdir -p $GO_RR

   echo
   echo $GO_AA
   echo $GO_RR
   echo $GO_FILE
   echo 

   cd $GO_RR

   # Output file
   OUTfile=${GO_FILE%CLEAN.vcf.gz}

   
   # Select samples of interest
   # Select sites if genotyped in 100% of all samples
   echo
   echo "Selecting samples. Selecting sites genotyped in 100% of all samples..."
   gatk --java-options "-Xmx5g" \
        SelectVariants \
        -V ${GO_AA}/$GO_FILE \
        --restrict-alleles-to BIALLELIC \
        --sample-name $GO_SL \
        --max-nocall-fraction 0.0 \
        --remove-unused-alternates true \
        --exclude-non-variants true \
        --exclude-filtered true \
        --verbosity ERROR \
        -O $SNIC_TMP/_myVCF_1.vcf


   ## Update vcf fields for each variant after removing samples
   echo
   echo 'Update VCF fields...'

   module unload GATK/4.2.0.0
   module load GATK/3.8-0

    #### Annotate variant calls with context information (updates INFO filed based on current set of samples)
    #### Note: the following fields are automatically updated (the list below is not exhaustive...):
    ####       AF: Allele Frequency, for each ALT allele, in the same order as listed
    ####       AC: Allele count in genotypes, for each ALT allele, in the same order as listed


   java -Xmx5G -jar $GATK_HOME/GenomeAnalysisTK.jar \
                 -T VariantAnnotator \
                 -R $refGenomefolder/${refGenome} \
                 -V $SNIC_TMP/_myVCF_1.vcf \
                 -A Coverage \
                 -o $SNIC_TMP/_myVCF_1A.vcf


   module unload GATK/3.8-0
   module load GATK/4.2.0.0

   ## Identify and remove position of private alleles (singletons only) in VCF file associated to tetraploid samples
   ## There is currently no tool that can do this on datasets that include tetraploid samples. Instead, we will remove sites with AC=1.
   ## This will remove sites where alternative allele, in biallelic sites, is observed in just one sample. 
   ## This method will not remove  private tetratons (1/1/1/1).
   ## We also remove sites where all samples have the same alternate allele (AF=1.0)

   echo
   echo 'Flag private alleles...'

   gatk --java-options "-Xmx5g" \
      VariantFiltration \
      -V $SNIC_TMP/_myVCF_1A.vcf \
      --filter-name "AC_filter" \
      --filter-expression "AC < 2.0" \
      --filter-name "AF_filter" \
      --filter-expression "AF == 1.0" \
      --verbosity ERROR \
      -O $SNIC_TMP/_myVCF_1F.vcf


   # Check filtering; count number of variants
   cat $SNIC_TMP/_myVCF_1F.vcf | grep '^#' > $SNIC_TMP/_myVCF_2.vcf
   cat $SNIC_TMP/_myVCF_1F.vcf | grep -v '^#' | grep -v ';AF=1.00;' | grep -v ';AF=0.00;' |  grep -Pv '\tAC=1;' | grep PASS >> $SNIC_TMP/_myVCF_2.vcf
   NVAR_SNP="$(cat $SNIC_TMP/_myVCF_2.vcf | grep -v '^#' | wc -l)"
   echo
   echo "Number of variants (genotyped in 100% of all samples):" $NVAR_SNP
   echo


   #Select 10000 random SNPS
   echo
   echo "Selecting 10k SNPs..."
   if [[ ${NVAR_SNP} -gt 10000 ]] ; then

      aux_10k=11000
      RF="$(awk '{print $2/$1}' <<<"$NVAR_SNP $aux_10k")"

      gatk --java-options "-Xmx5g" \
           SelectVariants \
           -V $SNIC_TMP/_myVCF_2.vcf \
           --select-random-fraction ${RF} \
           --remove-unused-alternates true \
           --exclude-non-variants true \
           --exclude-filtered true \
           --verbosity ERROR \
           -O $SNIC_TMP/_myVCF_3.vcf

   cat $SNIC_TMP/_myVCF_3.vcf | grep '^#' > $SNIC_TMP/_myVCF_4.vcf
   cat $SNIC_TMP/_myVCF_3.vcf | grep -v '^#' | head -n 10000 >> $SNIC_TMP/_myVCF_4.vcf   
   bgzip -c $SNIC_TMP/_myVCF_4.vcf > $GO_RR/${OUTfile}10kSNP.vcf.gz   

   else
      bgzip -c $SNIC_TMP/_myVCF_2.vcf > $GO_RR/${OUTfile}10kSNP.vcf.gz
   fi

   NVAR_SNP2="$(zcat $GO_RR/${OUTfile}10kSNP.vcf.gz | grep -v '^#' | wc -l)"
   echo
   echo "Number of variants (final):" $NVAR_SNP2
   echo


      # --max-nocall-fraction 	Maximum fraction of samples with no-call genotypes (default: 1.0)
      # --remove-unused-alternates	Remove alternate alleles not present in any genotypes (default: false)
      # --exclude-non-variants	Don't include non-variant sites (default: false)
      # --exclude-filtered		Don't include filtered sites. If this flag is enabled, sites that have been marked as filtered (i.e. have anything other than `.` or `PASS` in the FILTER field) will be excluded from the output. (default: false)
      # --select-random-fraction / -fraction		Select a fraction of variants at random from the input. The value of this argument should be a number between 0 and 1 specifying the fraction of total variants to be randomly selected from the input callset. Note that this is done using a probabilistic function, so the final result is not guaranteed to carry the exact fraction requested. Can be used for large fractions. 

   rm -f $SNIC_TMP/_*

done


echo
echo "Done!"




