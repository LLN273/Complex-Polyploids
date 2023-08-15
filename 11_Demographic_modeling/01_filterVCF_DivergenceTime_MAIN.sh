#!/bin/bash
#
#SBATCH -J fVCF
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 24:00:00
#SBATCH -A snic2022-22-909
#SBATCH -M snowy
#SBATCH --mail-user luis.leal@ebc.uu.se
#SBATCH --mail-type=FAIL
ulimit -c unlimited



### Scrip used to find variant and invariant sites used in fastsimcoal2 analysis
### This is done in three steps:
###  1. find loci of interest in intronic regions (assumed to be neutral). Include also invariant sites.
###  2. find loci of interest in exonic regions (4-fold degenerate sites). Include also invariant sites.
###  3. join the two datasets.


##### load modules
module load bioinfo-tools
module load GATK/4.2.0.0
module load htslib/1.9
module load tabix/0.2.6
module load picard/2.10.3


## SnpSift (SnpSift version 5.0)
## https://pcingola.github.io/SnpEff/ss_filter/
## Citation: "Cingolani P., Platts A., Wang L.L., Coon M., Nguyen T., Wang L., Land S.J., Lu X., Ruden D.M. 2012. A program for annotating and predicting the effects of single nucleotide polymorphisms, SnpEff. Fly. 6:80–92."
SnpSift=/crex1/proj/snic2017-7-149/private/Luis/z_APPS/snpEff/SnpSift.jar
SnpEff=/crex1/proj/snic2017-7-149/private/Luis/z_APPS/snpEff/snpEff.jar



####################################### paths

# remember current path
SRCDIR=$(pwd)  

# input folder
AA=/crex1/proj/snic2020-6-184/private/Luis/P06_birch_exome_IGA-Sweden_2019/18_VCFgatk_c_HardFiltering/BWA_unmasked

# input/output subfolder 1
RRsub=exome_WGS

# VCF list
VCF[1]=MASTER_SingleRuns_ALLsamples_GOODonly_SILVER_WHITE_PLATYPHYLLA
VCF[2]=${VCF[1]}
VCF[3]=${VCF[1]}
VCF[4]=${VCF[1]}

# input extension
OEXT[1]=".vcf.gz"
OEXT[2]=${OEXT[1]}
OEXT[3]=${OEXT[1]}
OEXT[4]=${OEXT[1]}

# output folder
RR=/crex1/proj/snic2020-6-184/private/Luis/P06_birch_exome_IGA-Sweden_2019/32a_fastsimcoal2_Bpubescens_DivergenceTime_VCFfiles/BWA_unmasked

# output extension
OUTEXT[1]="PubSpain"
OUTEXT[2]="PubCentralEurope"
OUTEXT[3]="PubSVsouth"
OUTEXT[4]="PubCentralAsia"

# Depth filtering >> adjust also lines 220 and 454
DPfilter=8

# Intervals list (BED format) exome data, target genes
ANNF_intronic=/crex1/proj/snic2017-7-149/private/Luis/P06_birch_exome_IGA-Sweden_2019/07_Probe_set_MASKS/Bpendula.annotation-targetGenes_introns_CLEAN.bed		# TARGET INTRONS ONLY
ANNF_exonic=/crex1/proj/snic2017-7-149/private/Luis/P06_birch_exome_IGA-Sweden_2019/07_Probe_set_MASKS/Bpendula.annotation-targetGenes_exons_CLEAN.bed		# TARGET EXONS ONLY

# Sample set
sample_set="4_8"

# Sample list 
SL[1]=/crex1/proj/snic2017-7-149/private/Luis/P06_birch_exome_IGA-Sweden_2019/32_fastsimcoal2_Bpubescens_DivergenceTime/00_samples_SILVER_WHITE_PLATYPHYLLA_divergenceTime_${OUTEXT[2]}_${sample_set}.args	# Bpubescens from Spain
SL[2]=/crex1/proj/snic2017-7-149/private/Luis/P06_birch_exome_IGA-Sweden_2019/32_fastsimcoal2_Bpubescens_DivergenceTime/00_samples_SILVER_WHITE_PLATYPHYLLA_divergenceTime_${OUTEXT[3]}_${sample_set}.args	# Bpubescens from Central Europe
SL[3]=/crex1/proj/snic2017-7-149/private/Luis/P06_birch_exome_IGA-Sweden_2019/32_fastsimcoal2_Bpubescens_DivergenceTime/00_samples_SILVER_WHITE_PLATYPHYLLA_divergenceTime_${OUTEXT[4]}_${sample_set}.args	# Bpubescens from Southern Sweden
SL[4]=/crex1/proj/snic2017-7-149/private/Luis/P06_birch_exome_IGA-Sweden_2019/32_fastsimcoal2_Bpubescens_DivergenceTime/00_samples_SILVER_WHITE_PLATYPHYLLA_divergenceTime_${OUTEXT[5]}_${sample_set}.args	# Bpubescens from Central Asia



#### List of loci of B. pendula or B. platyphylla ancestry 

if [[ ${sample_set} = "4_8" ]] ; then
   ##(4/8 dataset)
   Bpa[1]=/crex1/proj/snic2017-7-149/private/Luis/P06_birch_exome_IGA-Sweden_2019/32_fastsimcoal2_Bpubescens_DivergenceTime/presence_matrix_annotated_genes_pendula_${OUTEXT[2]}-Samples_4.txt
   Bpa[2]=/crex1/proj/snic2017-7-149/private/Luis/P06_birch_exome_IGA-Sweden_2019/32_fastsimcoal2_Bpubescens_DivergenceTime/presence_matrix_annotated_genes_pendula_${OUTEXT[3]}-Samples_4.txt
   Bpa[3]=/crex1/proj/snic2017-7-149/private/Luis/P06_birch_exome_IGA-Sweden_2019/32_fastsimcoal2_Bpubescens_DivergenceTime/presence_matrix_annotated_genes_pendula_${OUTEXT[4]}-Samples_4.txt
   Bpa[4]=/crex1/proj/snic2017-7-149/private/Luis/P06_birch_exome_IGA-Sweden_2019/32_fastsimcoal2_Bpubescens_DivergenceTime/presence_matrix_annotated_genes_pendula_${OUTEXT[5]}-Samples_4.txt
fi


# FLAGS reference genome (unmasked:1; masked:0)
FLAG_unmasked=1					

########### Reference genome (unmasked)
if [[ ${FLAG_unmasked} = "1" ]] ; then
   refGenome=Koivu_gapfilled_scaffolds.fasta
   refGenomefolder=/crex1/proj/snic2017-7-149/private/Luis/birch_reference_2017
fi




for k in `seq 1 1 4`; do 			# dataset (root)	

   GO_AA=${AA}/${RRsub}
   GO_RR=${RR}/${RRsub}/${DPfilter}X/${sample_set}
   GO_SL=${SL[$k]}
   GO_FILE=${VCF[$k]}${OEXT[$k]}
   GO_TargetList=${Bpa[$k]}
   mkdir -p $GO_RR

   echo
   echo $GO_AA
   echo $GO_RR
   echo $GO_FILE
   echo $GO_SL
   echo 

   # Output file
   OUTfile=${GO_FILE%.vcf.gz}_${OUTEXT[$k]}

   echo 
   echo $OUTfile > $SRCDIR/Glean_number_variants_01_filterVCF_DivergenceTime_${DPfilter}X_${OUTEXT[$k]}_${sample_set}.txt




   #############################
   ############################# 1. Intronic regions
   #############################
   echo
   echo "I. INTRONIC REGIONS"



   ################# Select intronic sites of B. pendula ancestry
   grep -w -f $GO_TargetList $ANNF_intronic > $GO_RR/Bpendula.annotation-targetGenes_introns_CLEAN_${OUTEXT[$k]}.bed


   ################# Create intervals list (target genes)
   ################# File must have ".intervals" extension (gatk crashes otherwise) [see: https://gatk.broadinstitute.org/hc/en-us/articles/360035531852?id=11009]

   rm -f $SNIC_TMP/__targets.intervals

   ### Convert position 0 to position 1 (GATK does not accept position 0)
   echo
   awk 'BEGIN {OFS="\t"} $2 == "0" {$2 = "1"}; 1'  $GO_RR/Bpendula.annotation-targetGenes_introns_CLEAN_${OUTEXT[$k]}.bed > $SNIC_TMP/__aux0.txt

   ### Convert intervals list to chr1:100-200 format
   echo
   cat $SNIC_TMP/__aux0.txt | cut -f1 > $SNIC_TMP/__aux1.txt
   cat $SNIC_TMP/__aux0.txt | cut -f2 > $SNIC_TMP/__aux2.txt
   cat $SNIC_TMP/__aux0.txt | cut -f3 > $SNIC_TMP/__aux3.txt
   paste -d':' $SNIC_TMP/__aux1.txt $SNIC_TMP/__aux2.txt > $SNIC_TMP/__auxA.txt
   paste -d'-' $SNIC_TMP/__auxA.txt $SNIC_TMP/__aux3.txt > $SNIC_TMP/__targets.intervals

   cd $GO_RR   
   
   
   ################# Select samples of interest
   ################# Select intronic regions of B. pendula ancestry only

   echo
   echo "1. Selecting samples. Selecting intronic regions of B. pendula ancestry."
   
   gatk --java-options "-Xmx5g" \
        SelectVariants \
        -V ${GO_AA}/$GO_FILE \
        --intervals $SNIC_TMP/__targets.intervals \
        --select-type-to-include SNP \
        --select-type-to-include NO_VARIATION \
        --sample-name $GO_SL \
        --remove-unused-alternates true \
        --exclude-filtered true \
        --verbosity ERROR \
        -O $SNIC_TMP/_myVCF_1.vcf
   
   N_records="$(cat $SNIC_TMP/_myVCF_1.vcf | grep -v '^#' | wc -l)"
   echo $N_records


   ################# Update vcf fields for each variant after removing samples
   echo
   echo '2. Updating VCF fields...'

   module unload GATK/4.2.0.0
   module load GATK/3.8-0

    #### Annotate variant calls with context information (updates INFO filed based on current set of samples)

   java -Xmx5G -jar $GATK_HOME/GenomeAnalysisTK.jar \
                 -T VariantAnnotator \
                 -R $refGenomefolder/${refGenome} \
                 -V $SNIC_TMP/_myVCF_1.vcf \
                 -A Coverage \
                 -o $SNIC_TMP/_myVCF_2.vcf
   
   N_records="$(cat $SNIC_TMP/_myVCF_2.vcf | grep -v '^#' | wc -l)"
   echo $N_records


   module unload GATK/3.8-0
   module load GATK/4.2.0.0
   
   ################# Filter VCF file based on genotype fields (failed genotypes set to ./.):   		
   #################       allele min depth (8)              

   echo
   echo "3A. DP filtering"

   cat $SNIC_TMP/_myVCF_2.vcf | java -Xmx5G -jar $SnpSift \
                                                 filter "( GEN[*].DP < 8 )" \
                                                 --addFilter FAIL \
                                                 > $SNIC_TMP/${OUTfile}_DPfilter.vcf
   
   N_records="$(cat $SNIC_TMP/${OUTfile}_DPfilter.vcf | grep -v '^#' | wc -l)"
   echo $N_records
   
   # Notes (SnpSift)
   # GEN[*] means we are filetring all samples based on genotype fields (DP, GQ, etc) 
   # -a|--addFilter   : Add a string to FILTER VCF field if 'expression' is true. Default: '' (none)
   
   
   ################# Compress file
   echo
   echo "3B. Compressing VCF file"
   bgzip -c $SNIC_TMP/${OUTfile}_DPfilter.vcf > $SNIC_TMP/${OUTfile}_DPfilter.vcf.gz
   
   ################# Index vcf file
   echo
   echo "3C. Index VCF file"
   gatk --java-options "-Xmx5G" \
                       IndexFeatureFile \
                       -I $SNIC_TMP/${OUTfile}_DPfilter.vcf.gz

    
   
   ################# Create separate VCF files for variant and invariant sites
   ################# Select sites if genotyped in 100% of all samples
   ################# Variant sites: Biallelic sites only
   echo
   echo '4A. Creating separate VCF files for variant and invariant sites. Select sites genotyped in 100% of all samples.'
   
   #  Variant sites: select biallelic sites genotyped in 100% of all samples...
   gatk --java-options "-Xmx5g" \
        SelectVariants \
        -R $refGenomefolder/${refGenome} \
        -V $SNIC_TMP/${OUTfile}_DPfilter.vcf.gz \
        --intervals $SNIC_TMP/__targets.intervals \
        --select-type-to-include SNP \
        --restrict-alleles-to BIALLELIC \
        --max-nocall-fraction 0.0 \
        --remove-unused-alternates true \
        --exclude-filtered true \
        --verbosity ERROR \
        -O $SNIC_TMP/${OUTfile}_intronic_SNP_aux.vcf.gz
   
   N_records="$(zcat $SNIC_TMP/${OUTfile}_intronic_SNP_aux.vcf.gz | grep -v '^#' | wc -l)"
   echo $N_records
   
   #  Invariant sites: select sites genotyped in 100% of all samples...
   echo
   gatk --java-options "-Xmx5g" \
        SelectVariants \
        -R $refGenomefolder/${refGenome} \
        -V $SNIC_TMP/${OUTfile}_DPfilter.vcf.gz \
        --intervals $SNIC_TMP/__targets.intervals \
        --select-type-to-include NO_VARIATION \
        --max-nocall-fraction 0.0 \
        --remove-unused-alternates true \
        --exclude-filtered true \
        --verbosity ERROR \
        -O $SNIC_TMP/${OUTfile}_intronic_nonVariant_aux.vcf.gz

   N_records="$(zcat $SNIC_TMP/${OUTfile}_intronic_nonVariant_aux.vcf.gz | grep -v '^#' | wc -l)"
   echo $N_records
   
   

   ################# Clean VCF files
   ################# 
   echo
   echo '5A. Cleaning VCF files'
   
   # variant sites
   zcat $SNIC_TMP/${OUTfile}_intronic_SNP_aux.vcf.gz | grep '^#' > $SNIC_TMP/${OUTfile}_intronic_SNP.vcf
   zcat $SNIC_TMP/${OUTfile}_intronic_SNP_aux.vcf.gz | grep -v '^#' > $SNIC_TMP/__aux_pass.txt
   awk '$5 != "." { print $0 }' $SNIC_TMP/__aux_pass.txt > $SNIC_TMP/__aux_pass1.txt
   awk '$5 != "*" { print $0 }' $SNIC_TMP/__aux_pass1.txt > $SNIC_TMP/__aux_pass2.txt
   cat $SNIC_TMP/__aux_pass2.txt | grep -v ';AF=1.00;' | grep -v ';AF=0.00;' | grep PASS >> $SNIC_TMP/${OUTfile}_intronic_SNP.vcf
   
   
   # invariant sites
   zcat $SNIC_TMP/${OUTfile}_intronic_nonVariant_aux.vcf.gz | grep '^#' > $SNIC_TMP/${OUTfile}_intronic_nonVariant.vcf
   zcat $SNIC_TMP/${OUTfile}_intronic_nonVariant_aux.vcf.gz | grep -v '^#' > $SNIC_TMP/__aux_pass.txt
   awk '$5 == "." { print $0 }' $SNIC_TMP/__aux_pass.txt > $SNIC_TMP/__aux_pass1.txt
   awk '($4 == "A" || $4 == "C" || $4 == "G" || $4 == "T") { print $0 }' $SNIC_TMP/__aux_pass1.txt > $SNIC_TMP/__aux_pass2.txt
   cat $SNIC_TMP/__aux_pass2.txt | grep PASS >> $SNIC_TMP/${OUTfile}_intronic_nonVariant.vcf
   

   ################# Compress output file
   echo
   echo "5B. Compress VCF files"
   bgzip -c $SNIC_TMP/${OUTfile}_intronic_SNP.vcf > $GO_RR/${OUTfile}_intronic_SNP.vcf.gz
   bgzip -c $SNIC_TMP/${OUTfile}_intronic_nonVariant.vcf > $GO_RR/${OUTfile}_intronic_nonVariant.vcf.gz
   
   
   ################# Index vcf file
   echo
   echo "5C. Index VCF file"
   gatk --java-options "-Xmx5G" \
                       IndexFeatureFile \
                       -I $GO_RR/${OUTfile}_intronic_SNP.vcf.gz
   
   echo
   gatk --java-options "-Xmx5G" \
                       IndexFeatureFile \
                       -I $GO_RR/${OUTfile}_intronic_nonVariant.vcf.gz




   
   ################# Check filtering; count number of variant sites
   zcat $GO_RR/${OUTfile}_intronic_SNP.vcf.gz | grep -v '^#' > $SNIC_TMP/__aux_pass.txt
   awk '$5 != "." { print $0 }' $SNIC_TMP/__aux_pass.txt > $SNIC_TMP/__aux_pass1.txt
   awk '$5 != "*" { print $0 }' $SNIC_TMP/__aux_pass1.txt > $SNIC_TMP/__aux_pass2.txt
   NVAR_SNP_introns="$(cat $SNIC_TMP/__aux_pass2.txt | grep -v ';AF=1.00;' | grep -v ';AF=0.00;' | grep PASS | wc -l)"
   echo
   echo "Number of variants, intronic:" $NVAR_SNP_introns
   echo
   echo "Number of variants, intronic:" $NVAR_SNP_introns >> $SRCDIR/Glean_number_variants_01_filterVCF_DivergenceTime_${DPfilter}X_${OUTEXT[$k]}_${sample_set}.txt


   ################# Check filtering; count number of non-variant sites
   zcat $GO_RR/${OUTfile}_intronic_nonVariant.vcf.gz | grep -v '^#' > $SNIC_TMP/__aux_pass.txt
   #awk '$5 != "." { print $0 }' $SNIC_TMP/__aux_pass.txt > $SNIC_TMP/__aux_pass1.txt
   awk '$5 != "*" { print $0 }' $SNIC_TMP/__aux_pass.txt > $SNIC_TMP/__aux_pass2.txt
   NVAR_nonVar_introns="$(cat $SNIC_TMP/__aux_pass2.txt | grep PASS | wc -l)"
   echo
   echo "Number of non-variants, intronic:" $NVAR_nonVar_introns
   echo
   echo "Number of non-variants, intronic:" $NVAR_nonVar_introns >> $SRCDIR/Glean_number_variants_01_filterVCF_DivergenceTime_${DPfilter}X_${OUTEXT[$k]}_${sample_set}.txt








   #############################
   ############################# 2. Exonic regions (4FD only)
   #############################
   echo
   echo
   echo "II. EXONIC REGIONS"

   ################# Select exonic sites of B. pendula ancestry
   ################# Note: use -w option to get only exact matches
   grep -w -f $GO_TargetList $ANNF_exonic > $GO_RR/Bpendula.annotation-targetGenes_exons_CLEAN_${OUTEXT[$k]}.bed


   ################# Create intervals list (target genes)
   ################# File must have ".intervals" extension (gatk crashes otherwise) [see: https://gatk.broadinstitute.org/hc/en-us/articles/360035531852?id=11009]

   rm -f $SNIC_TMP/__targets.intervals
   rm -f $SNIC_TMP/__intervals.bed

   ### Convert position 0 to position 1 (GATK does not accept position 0)
   echo
   awk 'BEGIN {OFS="\t"} $2 == "0" {$2 = "1"}; 1'  $GO_RR/Bpendula.annotation-targetGenes_exons_CLEAN_${OUTEXT[$k]}.bed > $SNIC_TMP/__aux0.txt

   ### Convert intervals list to chr1:100-200 format
   echo
   cat $SNIC_TMP/__aux0.txt | cut -f1 > $SNIC_TMP/__aux1.txt
   cat $SNIC_TMP/__aux0.txt | cut -f2 > $SNIC_TMP/__aux2.txt
   cat $SNIC_TMP/__aux0.txt | cut -f3 > $SNIC_TMP/__aux3.txt
   paste -d':' $SNIC_TMP/__aux1.txt $SNIC_TMP/__aux2.txt > $SNIC_TMP/__auxA.txt
   paste -d'-' $SNIC_TMP/__auxA.txt $SNIC_TMP/__aux3.txt > $SNIC_TMP/__targets.intervals
   
   # annotation file when running snpEff
   paste $SNIC_TMP/__aux1.txt $SNIC_TMP/__aux2.txt > $SNIC_TMP/__auxA2.txt
   paste $SNIC_TMP/__auxA2.txt $SNIC_TMP/__aux3.txt > $SNIC_TMP/__intervals.bed




   ################# Select samples of interest
   ################# Select exonic regions of B. pendula ancestry only

   echo
   echo "6. Selecting samples. Selecting exonic sites of B. pendula ancestry."
   gatk --java-options "-Xmx5g" \
        SelectVariants \
        -V ${GO_AA}/$GO_FILE \
        --intervals $SNIC_TMP/__targets.intervals \
        --select-type-to-include SNP \
        --select-type-to-include NO_VARIATION \
        --sample-name $GO_SL \
        --remove-unused-alternates true \
        --exclude-filtered true \
        --verbosity ERROR \
        -O $SNIC_TMP/_myVCF_E1.vcf
   
   N_records="$(cat $SNIC_TMP/_myVCF_E1.vcf | grep -v '^#' | wc -l)"
   echo $N_records
   
   
   
   ################# Update vcf fields for each variant after removing samples
   echo
   echo '7. Updating VCF fields...'

   module unload GATK/4.2.0.0
   module load GATK/3.8-0

    #### Annotate variant calls with context information (updates INFO filed based on current set of samples)
    #### Note: the following fields are automatically updates (the list below is not exhaustive...):
    ####       AF: Allele Frequency, for each ALT allele, in the same order as listed
    ####       AC: Allele count in genotypes, for each ALT allele, in the same order as listed


   java -Xmx5G -jar $GATK_HOME/GenomeAnalysisTK.jar \
                    -T VariantAnnotator \
                    -R $refGenomefolder/${refGenome} \
                    -V $SNIC_TMP/_myVCF_E1.vcf \
                    -A Coverage \
                    -o $SNIC_TMP/_myVCF_E2.vcf
   
   N_records="$(cat $SNIC_TMP/_myVCF_E2.vcf | grep -v '^#' | wc -l)"
   echo $N_records


   module unload GATK/3.8-0
   module load GATK/4.2.0.0

   
   ################# Filter VCF file based on genotype fields (failed genotypes set to ./.):   		
   #################       allele min depth (6 or 8)              

   echo
   echo "8A. DP filtering"

   cat $SNIC_TMP/_myVCF_E2.vcf | java -Xmx5G -jar $SnpSift \
                                                  filter "( GEN[*].DP < 8 )" \
                                                  --addFilter FAIL \
                                                  > $SNIC_TMP/${OUTfile}_DPfilter_exonic.vcf
   
   N_records="$(cat $SNIC_TMP/${OUTfile}_DPfilter_exonic.vcf | grep -v '^#' | wc -l)"
   echo $N_records
   
   # Notes (SnpSift)
   # GEN[*] means we are filetring all samples based on genotype fields (DP, GQ, etc) 
   # -a|--addFilter   : Add a string to FILTER VCF field if 'expression' is true. Default: '' (none)
   
   
   ################# Compress file
   echo
   echo "8B. Compressing VCF file"
   bgzip -c $SNIC_TMP/${OUTfile}_DPfilter_exonic.vcf > $SNIC_TMP/${OUTfile}_DPfilter_exonic.vcf.gz
   
   ################# Index vcf file
   echo
   echo "8C. Index VCF file"
   gatk --java-options "-Xmx5G" \
                       IndexFeatureFile \
                       -I $SNIC_TMP/${OUTfile}_DPfilter_exonic.vcf.gz




   ################# Select sites if genotyped in 100% of all samples
   ################# Exonic variants: Biallelic sites only
   echo
   echo "9A. Selecting biallelic exonic sites genotyped in 100% of all samples..."
   gatk --java-options "-Xmx5g" \
        SelectVariants \
        -V $SNIC_TMP/${OUTfile}_DPfilter_exonic.vcf.gz \
        --intervals $SNIC_TMP/__targets.intervals \
        --select-type-to-include SNP \
        --restrict-alleles-to BIALLELIC \
        --max-nocall-fraction 0.0 \
        --remove-unused-alternates true \
        --exclude-filtered true \
        --verbosity ERROR \
        -O $SNIC_TMP/_myVCF_E3.vcf
   
   N_records="$(cat $SNIC_TMP/_myVCF_E3.vcf | grep -v '^#' | wc -l)"
   echo $N_records


   ################# Exonic invariant sites
   echo
   echo "9B. Selecting exonic invariant sites genotyped in 100% of all samples..."
   gatk --java-options "-Xmx5g" \
        SelectVariants \
        -V $SNIC_TMP/${OUTfile}_DPfilter_exonic.vcf.gz \
        --intervals $SNIC_TMP/__targets.intervals \
        --select-type-to-include NO_VARIATION \
        --max-nocall-fraction 0.0 \
        --remove-unused-alternates true \
        --exclude-filtered true \
        --verbosity ERROR \
        -O $SNIC_TMP/${OUTfile}_exonic_nonVariant_aux.vcf.gz
   
   N_records="$(zcat $SNIC_TMP/${OUTfile}_exonic_nonVariant_aux.vcf.gz | grep -v '^#' | wc -l)"
   echo $N_records
   
   
   ################# Clean VCF files (non variant only)
   ################# 
   echo
   echo '10A. Cleaning VCF files'
   
   # invariant sites
   zcat $SNIC_TMP/${OUTfile}_exonic_nonVariant_aux.vcf.gz | grep '^#' > $SNIC_TMP/${OUTfile}_exonic_nonVariant.vcf
   zcat $SNIC_TMP/${OUTfile}_exonic_nonVariant_aux.vcf.gz | grep -v '^#' > $SNIC_TMP/__aux_pass.txt
   awk '$5 == "." { print $0 }' $SNIC_TMP/__aux_pass.txt > $SNIC_TMP/__aux_pass1.txt
   awk '($4 == "A" || $4 == "C" || $4 == "G" || $4 == "T") { print $0 }' $SNIC_TMP/__aux_pass1.txt > $SNIC_TMP/__aux_pass2.txt
   cat $SNIC_TMP/__aux_pass2.txt | grep PASS >> $SNIC_TMP/${OUTfile}_exonic_nonVariant.vcf
   

   ################# Compress output file
   echo
   echo "10B. Compress VCF files"
   bgzip -c $SNIC_TMP/${OUTfile}_exonic_nonVariant.vcf > $GO_RR/${OUTfile}_exonic_nonVariant.vcf.gz
   
   
   ################# Index vcf file
   echo
   echo "10C. Index VCF file"
   
   gatk --java-options "-Xmx5G" \
        IndexFeatureFile \
        -I $GO_RR/${OUTfile}_exonic_nonVariant.vcf.gz
   
   
   
   
   ################# Check filtering; count total number of exonic variants
   ## exclude non-variant sites and sites that FAILED
   cat $SNIC_TMP/_myVCF_E3.vcf | grep -v '^#' > $SNIC_TMP/__aux_pass.txt
   awk '$5 != "." { print $0 }' $SNIC_TMP/__aux_pass.txt > $SNIC_TMP/__aux_pass1.txt
   awk '$5 != "*" { print $0 }' $SNIC_TMP/__aux_pass1.txt > $SNIC_TMP/__aux_pass2.txt
   NVAR_SNP_exonic_RAW="$(cat $SNIC_TMP/__aux_pass2.txt | grep -v ';AF=1.00;' | grep -v ';AF=0.00;' | grep PASS | wc -l)"
   echo
   echo "Number of variants, exonic (RAW):" $NVAR_SNP_exonic_RAW
   echo
   echo "Number of variants, exonic (RAW):" $NVAR_SNP_exonic_RAW >> $SRCDIR/Glean_number_variants_01_filterVCF_DivergenceTime_${DPfilter}X_${OUTEXT[$k]}_${sample_set}.txt


   ################# Check filtering; count total number of exonic invariants
   ## exclude sites that FAILED
   zcat $GO_RR/${OUTfile}_exonic_nonVariant.vcf.gz | grep -v '^#' > $SNIC_TMP/__aux_pass.txt
   #awk '$5 != "." { print $0 }' $SNIC_TMP/__aux_pass.txt > $SNIC_TMP/__aux_pass1.txt
   awk '$5 != "*" { print $0 }' $SNIC_TMP/__aux_pass.txt > $SNIC_TMP/__aux_pass2.txt
   NVAR_invar_exonic="$(cat $SNIC_TMP/__aux_pass2.txt | grep PASS | wc -l)"
   echo
   echo "Number of invariant sites, exonic:" $NVAR_invar_exonic
   echo
   echo "Number of invariant sites, exonic:" $NVAR_invar_exonic >> $SRCDIR/Glean_number_variants_01_filterVCF_DivergenceTime_${DPfilter}X_${OUTEXT[$k]}_${sample_set}.txt




   ################# Annotate VCF file using snpEff
   ################# Note: Bpendula snpEff database was created beforehand (located in /crex1/proj/snic2017-7-149/private/Luis/z_APPS/snpEff/data/Bpendula)
   echo
   echo "11. Annotate VCF file using snpEff"
   
   java -Xmx5g -jar $SnpEff \
                    -fi $SNIC_TMP/__intervals.bed \
			                    -stats snpEff_summary_${OUTEXT[$k]}.html \
                    -noLog \
                    -classic \
                    Bpendula \
                    $SNIC_TMP/_myVCF_E3.vcf > $SNIC_TMP/${OUTfile}_exonic_DPfilter_snpEff.vcf
   
   N_records="$(cat $SNIC_TMP/${OUTfile}_exonic_DPfilter_snpEff.vcf | grep -v '^#' | wc -l)"
   echo $N_records
   
   
   # snpEff options:
   # -fi , -filterInterval  <file>   : Only analyze changes that intersect with the intervals specified in this file (you may use this option many times)
   # -s , -stats                     : Name of stats file (summary). Default is 'snpEff_summary.html'
   # -csvStats                       : Create CSV summary file instead of HTML
   # -noLog							 : SnpEff will try to log usage statistics to our "log server". This is useful for us to understand user's needs and have some statistics on what users are doing with the program
   # 								   (e.g. decide whether a command or option is useful or not). Logging can be deactivated by using the -noLog command line option.
   # -no-downstream                  : Do not show DOWNSTREAM changes
   # -no-upstream                    : Do not show UPSTREAM changes
   # -classic                        : Use old style annotations instead of Sequence Ontology and Hgvs.






   ################# Keep only 4-fold degenerate sites
   ################# Scan VCF file and grep lines the contain SYNONYMOUS_CODING
   ################# From manual: https://pcingola.github.io/SnpEff/se_inputoutput/  >> SYNONYMOUS_CODING: Variant causes a codon that produces the same amino acid (e.g.: Ttg/Ctg, L/L)
   
   echo
   echo "12. Keep only 4-fold degenerate sites"

   cat $SNIC_TMP/${OUTfile}_exonic_DPfilter_snpEff.vcf | grep '^#' >  $SNIC_TMP/${OUTfile}_exonic_DPfilter_snpEff_4FD_aux.vcf
   cat $SNIC_TMP/${OUTfile}_exonic_DPfilter_snpEff.vcf | grep -v '^#' | grep -v 'NON_SYNONYMOUS_CODING' | grep 'SYNONYMOUS_CODING(LOW|SILENT' >> $SNIC_TMP/${OUTfile}_exonic_DPfilter_snpEff_4FD_aux.vcf

   N_records="$(cat $SNIC_TMP/${OUTfile}_exonic_DPfilter_snpEff_4FD_aux.vcf | grep -v '^#' | wc -l)"
   echo $N_records


   ################# Clean VCF files
   ################# 
   
   echo
   echo "13A. Reformat INFO field (remove comments inserted by snpEff)"
   cat $SNIC_TMP/${OUTfile}_exonic_DPfilter_snpEff_4FD_aux.vcf | grep '^#' > $SNIC_TMP/${OUTfile}_exonic_DPfilter_snpEff_4FD_aux2.vcf
   cat $SNIC_TMP/${OUTfile}_exonic_DPfilter_snpEff_4FD_aux.vcf | grep -v '^#' > $SNIC_TMP/__aux_pass.txt
   awk -v OFS="\t" '{sub(/\;EFF.*$/,"",$8); print}' $SNIC_TMP/__aux_pass.txt >> $SNIC_TMP/${OUTfile}_exonic_DPfilter_snpEff_4FD_aux2.vcf
   
   N_records="$(cat $SNIC_TMP/${OUTfile}_exonic_DPfilter_snpEff_4FD_aux2.vcf | grep -v '^#' | wc -l)"
   echo $N_records
   
   
   echo
   echo '13B. Cleaning VCF file'
   
   # variant sites
   cat $SNIC_TMP/${OUTfile}_exonic_DPfilter_snpEff_4FD_aux2.vcf | grep '^#' > $SNIC_TMP/${OUTfile}_exonic_DPfilter_snpEff_4FD.vcf
   cat $SNIC_TMP/${OUTfile}_exonic_DPfilter_snpEff_4FD_aux2.vcf | grep -v '^#' > $SNIC_TMP/__aux_pass.txt
   awk '$5 != "." { print $0 }' $SNIC_TMP/__aux_pass.txt > $SNIC_TMP/__aux_pass1.txt
   awk '$5 != "*" { print $0 }' $SNIC_TMP/__aux_pass1.txt > $SNIC_TMP/__aux_pass2.txt
   cat $SNIC_TMP/__aux_pass2.txt | grep -v ';AF=1.00;' | grep -v ';AF=0.00;' | grep PASS >> $SNIC_TMP/${OUTfile}_exonic_DPfilter_snpEff_4FD.vcf


   ################# Compress output file
   echo
   echo "13C. Compressing VCF files"
   bgzip -c $SNIC_TMP/${OUTfile}_exonic_DPfilter_snpEff_4FD.vcf > $GO_RR/${OUTfile}_exonic_snpEff_4FD.vcf.gz


   ################# Index vcf files
   echo
   echo "13D. Indexing VCF file"
   gatk --java-options "-Xmx5G" \
   	            IndexFeatureFile \
   	            -I $GO_RR/${OUTfile}_exonic_snpEff_4FD.vcf.gz



   ################# Check filtering; number of 4FD variants
   ## exclude non-variant sites and sites that FAILED
   cat $SNIC_TMP/${OUTfile}_exonic_DPfilter_snpEff_4FD.vcf | grep -v '^#' > $SNIC_TMP/__aux_pass.txt
   awk '$5 != "." { print $0 }' $SNIC_TMP/__aux_pass.txt > $SNIC_TMP/__aux_pass1.txt
   awk '$5 != "*" { print $0 }' $SNIC_TMP/__aux_pass1.txt > $SNIC_TMP/__aux_pass2.txt
   NVAR_SNP_exonic_4FD="$(cat $SNIC_TMP/__aux_pass2.txt | grep -v ';AF=1.00;' | grep -v ';AF=0.00;' | grep PASS | wc -l)"
   echo
   echo "Number of variants, exonic 4FD:" $NVAR_SNP_exonic_4FD
   echo
   echo "Number of variants, exonic 4FD:" $NVAR_SNP_exonic_4FD >> $SRCDIR/Glean_number_variants_01_filterVCF_DivergenceTime_${DPfilter}X_${OUTEXT[$k]}_${sample_set}.txt



   ################# Merge intronic and exonic datasets


   # VCf file list: variant sites
   echo
   echo "14A. Merging intronic and 4-fold degenerate variant datasets"
   
   rm -f $SNIC_TMP/_VCF.list
   echo $GO_RR/${OUTfile}_intronic_SNP.vcf.gz > $SNIC_TMP/_VCF.list
   echo $GO_RR/${OUTfile}_exonic_snpEff_4FD.vcf.gz >> $SNIC_TMP/_VCF.list
      
   java -Xmx8G -jar $PICARD_HOME/picard.jar MergeVcfs \
                                            REFERENCE_SEQUENCE=$refGenomefolder/${refGenome} \
                                            I=$SNIC_TMP/_VCF.list \
                                            O=$GO_RR/${OUTfile}_intronic_and_4FD_DPfilter.vcf.gz

		  
   # VCf file list: non-variant sites
   echo
   echo "14B. Merging intronic and exonic non-variant datasets"
   
   rm -f $SNIC_TMP/_VCF.list
   echo $GO_RR/${OUTfile}_intronic_nonVariant.vcf.gz > $SNIC_TMP/_VCF.list
   echo $GO_RR/${OUTfile}_exonic_nonVariant.vcf.gz >> $SNIC_TMP/_VCF.list
   
   java -Xmx8G -jar $PICARD_HOME/picard.jar MergeVcfs \
                                            REFERENCE_SEQUENCE=$refGenomefolder/${refGenome} \
                                            I=$SNIC_TMP/_VCF.list \
                                            O=$GO_RR/${OUTfile}_intronic_exonic_nonVariant.vcf.gz

	  

   
   ################# Check filtering; total number of variants (intronic + 4FD)
   ## exclude non-variant sites and sites that FAILED
   zcat $GO_RR/${OUTfile}_intronic_and_4FD_DPfilter.vcf.gz | grep -v '^#' > $SNIC_TMP/__aux_pass.txt
   awk '$5 != "." { print $0 }' $SNIC_TMP/__aux_pass.txt > $SNIC_TMP/__aux_pass1.txt
   awk '$5 != "*" { print $0 }' $SNIC_TMP/__aux_pass1.txt > $SNIC_TMP/__aux_pass2.txt
   NVAR_variants="$(cat $SNIC_TMP/__aux_pass2.txt | grep -v ';AF=1.00;' | grep -v ';AF=0.00;' | grep PASS | wc -l)"
   echo
   echo "Number of variants, intronic + exonic 4FD:" $NVAR_variants
   echo
   echo "Number of variants, intronic + exonic 4FD:" $NVAR_variants >> $SRCDIR/Glean_number_variants_01_filterVCF_DivergenceTime_${DPfilter}X_${OUTEXT[$k]}_${sample_set}.txt
   
   
   ################# Check filtering; total number of non-variant sites
   ## exclude non-variant sites and sites that FAILED
   zcat $GO_RR/${OUTfile}_intronic_exonic_nonVariant.vcf.gz | grep -v '^#' > $SNIC_TMP/__aux_pass.txt
   awk '$5 != "*" { print $0 }' $SNIC_TMP/__aux_pass.txt > $SNIC_TMP/__aux_pass2.txt
   NVAR_nonVariants="$(cat $SNIC_TMP/__aux_pass2.txt | grep PASS | wc -l)"
   echo
   echo "Number of non-variants, total:" $NVAR_nonVariants
   echo
   echo "Number of non-variants, total:" $NVAR_nonVariants >> $SRCDIR/Glean_number_variants_01_filterVCF_DivergenceTime_${DPfilter}X_${OUTEXT[$k]}_${sample_set}.txt


   rm -f $SNIC_TMP/_*
   rm -f $SNIC_TMP/${OUTfile}*

done


echo
echo "Done!"




