#!/bin/bash
#
#SBATCH -J nonCoding
#SBATCH -p core 
#SBATCH -n 1
#SBATCH -t 5:00:00
#SBATCH -A snic2021-22-727
#SBATCH -M snowy
#SBATCH --mail-user luis.leal@ebc.uu.se
#SBATCH --mail-type=FAIL
ulimit -c unlimited
set -eo pipefail



# Script used to remove coding regions from a VCF file

module load bioinfo-tools
module load GATK/4.2.0.0


#remember initial path
SRCDIR_INI=$(pwd)    


#input folder
AA=/crex1/proj/snic2020-6-184/private/Luis/P06_birch_exome_IGA-Sweden_2019/18_VCFgatk_c_HardFiltering/BWA_unmasked

# input/output subfolder 1
RRsub[1]=exome_WGS

# VCF file (root)
VCF[1]=MASTER_SingleRuns_ALLsamples_GOODonly_NANA_SILVER_WHITE_HUMILIS_PLATYPHYLLA

# VCF extensions
#Vext[1]=filtered_tranche_2_PHYLO.vcf.gz
#Vext[2]=filtered_tranche_2_PASS_FINAL.vcf.gz
Vext[3]=filtered_tranche_2_PASS-BIALLELIC_FINAL_NoPrivateAlleles.vcf.gz

### Output folder
RR=/crex1/proj/snic2020-6-184/private/Luis/P06_birch_exome_IGA-Sweden_2019/24_Data_for_PCA_and_Admixture/BWA_unmasked

##### Sample list
SPL[1]=/crex1/proj/snic2017-7-149/private/Luis/P06_birch_exome_IGA-Sweden_2019/24_LD_prunning/00_nana_silver_white_humilis_platy_ploidyAnalysis.args		

# Input annotation files (exons only)
ANNF_exons=/crex1/proj/snic2017-7-149/private/Luis/P06_birch_exome_IGA-Sweden_2019/07_Probe_set_MASKS/Bpendula.annotation-targetGenes_exons.bed           

# FLAGS reference genome (unmasked:1; masked:0)
FLAG_unmasked=1					


########### Reference genome (unmasked)
if [[ ${FLAG_unmasked} = "1" ]] ; then
   refGenome=Koivu_gapfilled_scaffolds.fasta
   refGenomefolder=/crex1/proj/snic2017-7-149/private/Luis/birch_reference_2017
fi



############################################################## Remove coding regions from VCF files


for k in `seq 1 1 1`; do 			# dataset (root)	

   GO_AA=${AA}/${RRsub[$k]} 
   GO_RR=${RR}/${RRsub[$k]}

   ################# Create intervals list (target genes)
   ################# File must have ".intervals" extension (gatk crashes otherwise) [see: https://gatk.broadinstitute.org/hc/en-us/articles/360035531852?id=11009]

   ### Convert position 0 to position 1 (GATK does not accept position 0)
   echo
   awk 'BEGIN {OFS="\t"} $2 == "0" {$2 = "1"}; 1'  ${ANNF_exons} > $SNIC_TMP/__aux0.txt

   ### Convert intervals list to chr1:100-200 format
   echo
   cat $SNIC_TMP/__aux0.txt | cut -f1 > $SNIC_TMP/__aux1.txt
   cat $SNIC_TMP/__aux0.txt | cut -f2 > $SNIC_TMP/__aux2.txt
   cat $SNIC_TMP/__aux0.txt | cut -f3 > $SNIC_TMP/__aux3.txt
   paste -d':' $SNIC_TMP/__aux1.txt $SNIC_TMP/__aux2.txt > $SNIC_TMP/__auxA.txt
   paste -d'-' $SNIC_TMP/__auxA.txt $SNIC_TMP/__aux3.txt > $SNIC_TMP/__targets.intervals

   mkdir -p $GO_RR

   for j in `seq 3 1 3`; do 			# dataset (extension)

      GO_SET[1]=${VCF[$k]}-${Vext[$j]}

      echo

      GO_VCF1=${GO_SET[1]}

      # Index vcf file
      echo
      echo "Index VCF file"
      gatk --java-options "-Xmx5G" \
	            IndexFeatureFile \
	            -I $GO_AA/${GO_VCF1}

      echo
      echo "remove coding regions"
      gatk --java-options "-Xmx5g" \
           VariantFiltration \
           -R $refGenomefolder/${refGenome} \
           -V $GO_AA/${GO_VCF1} \
           --exclude-intervals $SNIC_TMP/__targets.intervals \
           --verbosity ERROR \
           -O $SNIC_TMP/_NONCODING.vcf

      echo
      echo "include only samples of interest"
      gatk --java-options "-Xmx5g" \
        SelectVariants \
        -R $refGenomefolder/${refGenome} \
        -V $SNIC_TMP/_NONCODING.vcf \
        --sample-name ${SPL[$k]} \
        --exclude-intervals $SNIC_TMP/__targets.intervals \
        --remove-unused-alternates true \
        --exclude-non-variants true \
        --verbosity ERROR \
        -O $GO_RR/${GO_VCF1%.vcf.gz}_NONCODING.vcf

      rm -f $SNIC_TMP/_NONCODING.vcf


   done
done






