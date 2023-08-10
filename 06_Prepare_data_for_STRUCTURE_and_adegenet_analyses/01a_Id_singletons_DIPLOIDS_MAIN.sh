#!/bin/bash
#
#SBATCH -J singletons
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 12:00:00
#SBATCH -A snic2021-22-727
#SBATCH -M snowy
#SBATCH --mail-user luis.leal@ebc.uu.se
#SBATCH --mail-type=FAIL
ulimit -c unlimited


## Identify and remove position of private alleles (singletons and private doubletons) in VCF file >> can be applied to diploid samples only
## Remove also sites with AF=1 (indicates all samples have only the alternative allele)


##### load modules
module load bioinfo-tools
module load vcftools/0.1.16
module load GATK/4.2.0.0


########### paths

#remember current path
SRCDIR=$(pwd)                                               				

#input folder
AA=/crex1/proj/snic2020-6-184/private/Luis/P06_birch_exome_IGA-Sweden_2019/18_VCFgatk_c_HardFiltering/BWA_unmasked

# input/output subfolder 1
RRsub[1]=exome_WGS
RRsub[2]=exome_WGS
RRsub[3]=exome_WGS
RRsub[4]=exome_WGS
RRsub[5]=exome_WGS


# VCF list
VCF[1]=MASTER_SingleRuns_ALLsamples_GOODonly_SILVER_snps
VCF[2]=MASTER_SingleRuns_ALLsamples_GOODonly_SILVER_PLATYPHYLLA
VCF[3]=MASTER_SingleRuns_ALLsamples_GOODonly_NANA_snps
VCF[4]=MASTER_SingleRuns_ALLsamples_GOODonly_PLATYPHYLLA_snps
VCF[5]=MASTER_SingleRuns_ALLsamples_GOODonly_HUMILIS_snps


# VCF extensions
#Vext[1]=filtered_tranche_2_PHYLO.vcf.gz
Vext[2]=filtered_tranche_2_PASS_FINAL.vcf.gz
Vext[3]=filtered_tranche_2_PASS-BIALLELIC_FINAL.vcf.gz


###### Intervals list (BED format) exome data 
ANNF[1]=/crex1/proj/snic2017-7-149/private/Luis/P06_birch_exome_IGA-Sweden_2019/07_Probe_set_MASKS/Bpendula.annotation-targetGenes_Gene-AND-2Kupstream_genes.bed
ANNF[2]=${ANNF[1]}
ANNF[3]=${ANNF[1]}
ANNF[4]=${ANNF[1]}
ANNF[5]=${ANNF[1]}


###### Samples to be removed
SRMV[1]=/crex1/proj/snic2017-7-149/private/Luis/P06_birch_exome_IGA-Sweden_2019/18_clean_VCF_VariantSites_datasets/00_RemoveSamples_1_MASTER_SingleRuns_ALLsamples_GOODonly_SILVER_snps.args
SRMV[2]=/crex1/proj/snic2017-7-149/private/Luis/P06_birch_exome_IGA-Sweden_2019/18_clean_VCF_VariantSites_datasets/00_RemoveSamples_2_MASTER_SingleRuns_ALLsamples_GOODonly_SILVER_PLATYPHYLLA.args


###### FLAGS reference genome (unmasked:1; masked:0)
FLAG_unmasked=1					

###### Reference genome (unmasked)
if [[ ${FLAG_unmasked} = "1" ]] ; then
   refGenome=Koivu_gapfilled_scaffolds.fasta
   refGenomefolder=/crex1/proj/snic2017-7-149/private/Luis/birch_reference_2017
fi




################## Remove private alleles 


for k in `seq 1 1 5`; do 			# dataset (root)	

   GO_AA=${AA}/${RRsub[$k]}
   GO_RR=${GO_AA}

   ################# Create intervals list (target genes)
   ################# File must have ".intervals" extension (gatk crashes otherwise) [see: https://gatk.broadinstitute.org/hc/en-us/articles/360035531852?id=11009]

   ### Convert position 0 to position 1 (GATK does not accept position 0)
   echo
   awk 'BEGIN {OFS="\t"} $2 == "0" {$2 = "1"}; 1'  ${ANNF[$k]} > $SNIC_TMP/__aux0.txt

   ### Convert intervals list to chr1:100-200 format
   echo
   cat $SNIC_TMP/__aux0.txt | cut -f1 > $SNIC_TMP/__aux1.txt
   cat $SNIC_TMP/__aux0.txt | cut -f2 > $SNIC_TMP/__aux2.txt
   cat $SNIC_TMP/__aux0.txt | cut -f3 > $SNIC_TMP/__aux3.txt
   paste -d':' $SNIC_TMP/__aux1.txt $SNIC_TMP/__aux2.txt > $SNIC_TMP/__auxA.txt
   paste -d'-' $SNIC_TMP/__auxA.txt $SNIC_TMP/__aux3.txt > $SNIC_TMP/__targets.intervals


   for j in `seq 2 1 3`; do 			# dataset (extension)

      GO_SET=${VCF[$k]}-${Vext[$j]}
      echo
      echo 'Dataset:' $GO_SET

      # for SILVER and SILVER_PLATYPHYLLA exome datasets, remove WGS samples
      if [[ $k -le "2" ]] ; then
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
      else
         rsync -ahv $GO_AA/${GO_SET} $SNIC_TMP/_SFfilter.vcf.gz
      fi

      echo
      echo 'Removing private alleles...'

      # remove sites with AF=1
      gatk --java-options "-Xmx5g" \
                       VariantFiltration \
                       -R $refGenomefolder/${refGenome} \
                       -V $SNIC_TMP/_SFfilter.vcf.gz \
                       --intervals $SNIC_TMP/__targets.intervals \
                       --filter-name "AF_filter" \
                       --filter-expression "AF == 1.0" \
                       --verbosity ERROR \
                       -O $SNIC_TMP/_AFfilter.vcf.gz

      echo

      gatk --java-options "-Xmx5g" \
                          SelectVariants \
                          -R $refGenomefolder/${refGenome} \
                          -V $SNIC_TMP/_AFfilter.vcf.gz \
                          --intervals $SNIC_TMP/__targets.intervals \
                          --remove-unused-alternates true \
                          --exclude-non-variants true \
                          --exclude-filtered true \
                          -O $SNIC_TMP/_NoInvariant.vcf.gz
 
      echo

      # id position private alleles
      vcftools --gzvcf \
               $SNIC_TMP/_NoInvariant.vcf.gz \
               --singletons \
               --out $GO_RR/${GO_SET%.vcf.gz}

      cat $GO_RR/${GO_SET%.vcf.gz}.singletons | tail -n +2 | cut -f1-2 > $SNIC_TMP/_singletons.txt

      # remove private alleles
      vcftools --gzvcf \
               $SNIC_TMP/_NoInvariant.vcf.gz \
               --exclude-positions $SNIC_TMP/_singletons.txt \
               --recode \
               --recode-INFO-all \
              --out $GO_RR/${GO_SET%.vcf.gz}_NoPrivateAlleles

      rm -f $SNIC_TMP/_SFfilter.vcf.gz
      rm -f $SNIC_TMP/_AFfilter.vcf.gz
      rm -f $SNIC_TMP/_NoInvariant.vcf.gz
      rm -f $SNIC_TMP/_singletons.txt

      echo
      echo
			    
   done
done


# Notes
# --singletons		This option will generate a file detailing the location of singletons, and the individual they occur in. The file reports both true singletons, and private doubletons (i.e. SNPs where the minor allele only occurs in a single individual and that 
#			individual is homozygotic for that allele). The output file has the suffix ".singletons".

# --exclude-positions <filename>  Exclude a set of sites on the basis of a list of positions in a file. Each line of the input file should contain a (tab-separated) chromosome and position. The file can have comment lines that start with a "#", they will be ignored.
# --recode		This option is used to generate a new file in VCF from the input VCF file after applying the filtering options specified by the user. The output file has the suffix ".recode.vcf". By default, the INFO fields are removed from the output file, 
#			as the INFO values may be invalidated by the recoding (e.g. the total depth may need to be recalculated if individuals are removed). This behavior may be overriden by the following options. By default, BCF files are written out as BGZF compressed files.
# --recode-INFO-all	Keep all INFO values in the original file.



echo
echo 'Done!'

# remove auxiliary files									    
rm -f $SNIC_TMP/_*
