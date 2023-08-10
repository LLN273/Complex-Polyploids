## setwd() # make sure wd contains shell scripts invoked here!

require(tidyverse)
require(data.table)

source("ld_easy.R")


# read contig list
contig_df <- read.table("00_contig_list-WHITE_mixploidy.txt", header=F)
contig_list <- contig_df$V1


prepVCF(
    vcfdir="/crex1/proj/snic2020-6-184/private/Luis/P06_birch_exome_IGA-Sweden_2019/24_Data_for_PCA_and_Admixture/BWA_unmasked/exome_WGS/MASTER_SingleRuns_ALLsamples_GOODonly_NANA_SILVER_WHITE_HUMILIS_PLATYPHYLLA-filtered_tranche_2_PASS-BIALLELIC_FINAL_NoPrivateAlleles_NONCODING_ReducedPEND_ReducedSW.vcf.gz",
    bcftools.loadpath="module load gcc/4.8.2 perl/5.18.4 samtools/1.9",
    chromosomes=contig_list,
    temp.dir="/crex1/proj/snic2020-6-184/private/Luis/P06_birch_exome_IGA-Sweden_2019/24_Data_for_PCA_and_Admixture/BWA_unmasked/exome_WGS/_NANA_SILVER_WHITE_HUMILIS_PLATYPHYLLA"
)

pruneChromosomes(
    temp.dir="/crex1/proj/snic2020-6-184/private/Luis/P06_birch_exome_IGA-Sweden_2019/24_Data_for_PCA_and_Admixture/BWA_unmasked/exome_WGS/_NANA_SILVER_WHITE_HUMILIS_PLATYPHYLLA",
    windowsize=1000,
    rlimit=0.1
)

postVCF(
    vcfdir="/crex1/proj/snic2020-6-184/private/Luis/P06_birch_exome_IGA-Sweden_2019/24_Data_for_PCA_and_Admixture/BWA_unmasked/exome_WGS/MASTER_SingleRuns_ALLsamples_GOODonly_NANA_SILVER_WHITE_HUMILIS_PLATYPHYLLA-filtered_tranche_2_PASS-BIALLELIC_FINAL_NoPrivateAlleles_NONCODING_ReducedPEND_ReducedSW.vcf.gz",
    bcftools.loadpath="module load gcc/4.8.2 perl/5.18.4 samtools/1.9",
    temp.dir="/crex1/proj/snic2020-6-184/private/Luis/P06_birch_exome_IGA-Sweden_2019/24_Data_for_PCA_and_Admixture/BWA_unmasked/exome_WGS/_NANA_SILVER_WHITE_HUMILIS_PLATYPHYLLA"
)


# bcftools.loadpath="module load gcc/4.8.2 gdc perl/5.18.4 samtools/1.9",
