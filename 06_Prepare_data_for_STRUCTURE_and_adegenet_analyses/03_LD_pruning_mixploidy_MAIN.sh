#!/bin/bash
#
#SBATCH -J LDtetra
#SBATCH -p core 
#SBATCH -n 3
#SBATCH -t 10:00:00
#SBATCH -A snic2021-22-727
#SBATCH -M snowy
#SBATCH --mail-user luis.leal@ebc.uu.se
#SBATCH --mail-type=FAIL
ulimit -c unlimited
set -eo pipefail



# Script used to id SNPs in linkage disequilibrium (LD)
# r2 in this case is a measure of the association (correlation) between two genotypes associated to two different loci located within a certain window

# Using a modifies version of scripts developed for LD analysis of tetraploids by Weitz et al.
# https://github.com/LZeitler/tetrarenosa

## Citation
# Male meiotic recombination rate varies with seasonal temperature fluctuations in wild populations of autotetraploid Arabidopsis arenosa
# Andrew P. Weitz, Marinela Dukic, Leo Zeitler, Kirsten Bomblies
# First published: 17 July 2021
# https://doi.org/10.1111/mec.16084

# IMPORTANT: 
# 'ld_easy_prepvcf.sh' script modified in line 12:
# no restrictions on MAF or missing data, besides those already in use (no private alleles; missing < 50%) The problem here is that some populations are represented by one or two samples only


module load bioinfo-tools
module load bcftools/1.8
module load samtools/1.9
module load htslib/1.9
module load tabix/0.2.6
module load R/3.6.0			# R scripts used here were developed for R 3.6.0


#remember initial path
SRCDIR_INI=$(pwd)    

#input folder
AA=/crex1/proj/snic2020-6-184/private/Luis/P06_birch_exome_IGA-Sweden_2019/24_Data_for_PCA_and_Admixture/BWA_unmasked

# input/output subfolder 1
RRsub=exome_WGS

###### VCF list and extensions
VCF[1]=MASTER_SingleRuns_ALLsamples_GOODonly_NANA_SILVER_WHITE_HUMILIS_PLATYPHYLLA-filtered_tranche_2_PASS-BIALLELIC_FINAL_NoPrivateAlleles_NONCODING
Vext[1]="_ReducedPEND_ReducedSW.vcf"



for k in `seq 1 1 1`; do 			# dataset 

   GO_AA=${AA}/${RRsub}
   GO_RR=$GO_AA

   GO_SET=${VCF[$k]}${Vext[$k]}
   echo
   echo $GO_SET

   # compress and index vcf file
   bgzip -c $GO_AA/${GO_SET} > $GO_AA/${GO_SET}.gz
   tabix $GO_AA/${GO_SET}.gz

   # >>>>> IMPORTANT <<<<<<<
   #CONTIGS REMOVED (not enough variant sites) [LD script crashes]:
   cat $GO_AA/${GO_SET} | grep -v '^#' | cut -f1 | sort -t g -k 2 -g | uniq | grep -vx Contig141 | grep -vx Contig992 | grep -vx Contig340 | grep -vx Contig1586 | grep -vx Contig4423 | grep -vx Contig230 | grep -vx Contig795 | grep -vx Contig1593 | grep -vx Contig3916 | grep -vx Contig318 | grep -vx Contig130 | grep -vx Contig500 | grep -vx Contig2165 | grep -vx Contig426 | grep -vx Contig445 | grep -vx Contig1320 | grep -vx Contig1323 | grep -vx Contig1328 | grep -vx Contig1329 | grep -vx Contig1334 | grep -vx Contig1335 | grep -vx Contig1337 > $SRCDIR_INI/00_contig_list-NANA_SILVER_WHITE_HUMILIS_PLATYPHYLLA_mixploidy.txt

   # >>>>>>>>>>>>>>>>>>>>>>>>>>>  NOTE: 03s_ld_LARGEdatasets.R scrip must be adjusted depending on dataset <<<<<<<<<<<<<<<<<<<<<<<
   Rscript --no-save /crex1/proj/snic2017-7-149/private/Luis/P06_birch_exome_IGA-Sweden_2019/24_LD_prunning/03s_ld_LARGEdatasets.R

   #copy results to output folder
   cp $GO_AA/_NANA_SILVER_WHITE_HUMILIS_PLATYPHYLLA/ld_easy_temp/final.vcf.gz $GO_RR/${GO_SET%.vcf}_LD.vcf.gz

done


# clean auxiliary files
rm -f $SNIC_TMP/_*

echo
echo "Done!"



