#!/bin/bash
#
#SBATCH -J AlleleC
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 12:00:00
#SBATCH -A snic2022-22-589
#SBATCH -M snowy
#SBATCH --mail-user luis.leal@ebc.uu.se
#SBATCH --mail-type=FAIL
ulimit -c unlimited


## Identify and remove private alleles (singletons only).



##### load modules
module load bioinfo-tools
module load GATK/4.2.0.0


##### paths

#remember current path
SRCDIR=$(pwd)                                               				

#input folder
AA=/crex1/proj/snic2020-6-184/private/Luis/P06_birch_exome_IGA-Sweden_2019/18_VCFgatk_c_HardFiltering/BWA_unmasked

# input/output subfolder 1
RRsub[1]=exome_WGS

# VCF list
VCF[1]=MASTER_SingleRuns_ALLsamples_GOODonly_NANA_SILVER_WHITE_HUMILIS_PLATYPHYLLA

# VCF extensions
#Vext[1]=filtered_tranche_2_PHYLO.vcf.gz
#Vext[2]=filtered_tranche_2_PASS_FINAL.vcf.gz
Vext[3]=filtered_tranche_2_PASS-BIALLELIC_FINAL.vcf.gz

# Intervals list (BED format) exome data 
ANNF=/crex1/proj/snic2017-7-149/private/Luis/P06_birch_exome_IGA-Sweden_2019/07_Probe_set_MASKS/Bpendula.annotation-targetGenes_Gene-AND-2Kupstream_genes.bed	

# Samples to be removed
SRMV[1]=/crex1/proj/snic2017-7-149/private/Luis/P06_birch_exome_IGA-Sweden_2019/18_clean_VCF_VariantSites_datasets/00_RemoveSamples_E_MASTER_SingleRuns_ALLsamples_GOODonly_NANA_SILVER_WHITE_HUMILIS_PLATYPHYLLA.args

# FLAGS reference genome (unmasked:1; masked:0)
FLAG_unmasked=1					

########### Reference genome (unmasked)
if [[ ${FLAG_unmasked} = "1" ]] ; then
   refGenome=Koivu_gapfilled_scaffolds.fasta
   refGenomefolder=/crex1/proj/snic2017-7-149/private/Luis/birch_reference_2017
fi



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



############################## Remove private alleles 

for k in `seq 1 1 1`; do 			# dataset (root)	

   GO_AA=${AA}/${RRsub[$k]}
   GO_RR=${GO_AA}

   for j in `seq 3 1 3`; do 			# dataset (extension)

      GO_SET=${VCF[$k]}-${Vext[$j]}
      echo
      echo 'Dataset:' $GO_SET

      # remove WGS samples (apart from standards)
      echo
      echo "Removing WGS samples..."
      gatk --java-options "-Xmx5g" \
           SelectVariants \
           -R $refGenomefolder/${refGenome} \
           -V $GO_AA/${GO_SET} \
           --intervals $SNIC_TMP/__targets.intervals \
           --exclude-sample-name ${SRMV[$k]} \
           --select-type-to-include SNP \
           --exclude-non-variants true \
           --exclude-filtered true \
           --remove-unused-alternates true \
           --verbosity ERROR \
           -O $SNIC_TMP/_SFfilter.vcf.gz

      echo
      echo 'Removing private alleles...'

      gatk --java-options "-Xmx5g" \
      VariantFiltration \
      -R $refGenomefolder/${refGenome} \
      -V $SNIC_TMP/_SFfilter.vcf.gz \
      --intervals $SNIC_TMP/__targets.intervals \
      --filter-name "AC_filter" \
      --filter-expression "AC < 2.0" \
      --filter-name "AF_filter" \
      --filter-expression "AF == 1.0" \
      --verbosity ERROR \
      -O $SNIC_TMP/_NoPrivateAlleles.vcf.gz

      echo

      gatk --java-options "-Xmx5g" \
      SelectVariants \
      -R $refGenomefolder/${refGenome} \
      -V $SNIC_TMP/_NoPrivateAlleles.vcf.gz \
      --intervals $SNIC_TMP/__targets.intervals \
      --remove-unused-alternates true \
      --exclude-non-variants true \
      --exclude-filtered true \
      -O $GO_RR/${GO_SET%.vcf.gz}_NoPrivateAlleles.vcf.gz

      rm -f $SNIC_TMP/_SFfilter.vcf.gz
      rm -f $SNIC_TMP/_NoPrivateAlleles.vcf.gz
			    
   done
done


# remove auxiliary files									    
rm -f $SNIC_TMP/_*

echo
echo 'Done!'



