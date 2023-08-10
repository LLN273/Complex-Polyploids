#!/bin/bash
#
#SBATCH -J VCF2str
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 06:00:00
#SBATCH -A snic2022-22-909
#SBATCH -M snowy
#SBATCH --mail-user luis.leal@ebc.uu.se
#SBATCH --mail-type=FAIL
ulimit -c unlimited



### Scrip used to convert VCF to STRUCTURE file format


##### load modules
module load bioinfo-tools
module load R/4.1.1
module load GDAL/3.1.0


####################################### paths

#remember current path
SRCDIR=$(pwd)  

#input folder
AA=/crex1/proj/snic2020-6-184/private/Luis/P06_birch_exome_IGA-Sweden_2019/30_STRUCTURE/BWA_unmasked

# input/output subfolder 1
RRsub[1]=exome_WGS

# Number of SNPs
NL="10k"

# VCF list
VCF[1]=MASTER_SingleRuns_ALLsamples_GOODonly_NANA_SILVER_WHITE_HUMILIS_PLATYPHYLLA-filtered_tranche_2_PASS-BIALLELIC_FINAL_NoPrivateAlleles_NONCODING

# input extension
OEXT[1]="_ReducedPEND_ReducedSW_LD_SHORT_${NL}SNP.vcf.gz"

## Sample list
SL[1]=/crex1/proj/snic2017-7-149/private/Luis/P06_birch_exome_IGA-Sweden_2019/30_STRUCTURE/00_Sample_list_NANA_SILVER_WHITE_HUMILIS_PLATYPHYLLA_ReducedPEND_ReducedSW_SHORT.args    # 126 samples

### Ploidy for each sample
SPLOIDY_FOLDER=/crex1/proj/snic2017-7-149/private/Luis/P06_birch_exome_IGA-Sweden_2019/30_STRUCTURE
SPLOIDY_file=00_SamplePloidy_nana.pendula.pubescens.platy_Admixture.PLOTS.only_FORMAS-GeneTREE.txt



for k in `seq 1 1 1`; do 			# dataset (root)	

   GO_AA=${AA}/${RRsub[$k]}
   GO_RR=${GO_AA}
   GO_FILE=${VCF[$k]}${OEXT[$k]}
   SL_GO=${SL[$k]}
   mkdir -p $GO_RR

   echo $GO_AA
   echo $GO_RR
   echo $GO_FILE
   echo 

   cd $GO_RR

   # Output file
   OUTstr=${GO_FILE%.vcf.gz}.str


   ####### Use R scrip to get population code and ploidy for each sample (species based, or with more than one population within species)

   ## get list of samples (from VCF file)
   zcat ${GO_AA}/$GO_FILE | grep '^#' | tail -n 1 | sed 's/\t/\n/g' | tail -n +10 > $SNIC_TMP/_samples.txt

   ## Number of sample
   echo 'Number of samples:'
   wc -l $SNIC_TMP/_samples.txt
   echo


   ## Associate each sample in VCF to population id and ploidy
   echo "creating population and ploidy files..."
   echo
   Rscript --no-save $SRCDIR/02_population_ID.R $SPLOIDY_FOLDER/$SPLOIDY_file \
                                                       $SNIC_TMP/_samples.txt \
                                                       $GO_RR \
                                                       ${GO_FILE%.vcf.gz}.pop \
                                                       ${GO_FILE%.vcf.gz}.ploidy

   ## Check number of entries in population and ploidy files (should match number of samples)
   echo "Number of entries in population and ploidy files:"
   wc -l $GO_RR/${GO_FILE%.vcf.gz}.pop
   wc -l $GO_RR/${GO_FILE%.vcf.gz}.ploidy
   echo

   ## reformat sample names (convert '-' to '_')
   zcat ${GO_AA}/$GO_FILE | grep '^#' | sed 's/-/_/g' > $SNIC_TMP/_GO_VCF.vcf
   zcat ${GO_AA}/$GO_FILE | grep -v '^#' >> $SNIC_TMP/_GO_VCF.vcf

   ## Convert from VCF to STRUCTURE file format
   echo "Creating STRUCTURE files..."
   echo
   
   Rscript --no-save $SRCDIR/02_vcf2structure.R $GO_RR/${GO_FILE%.vcf.gz}.pop \
                                                $GO_RR/${GO_FILE%.vcf.gz}.ploidy \
                                                $SNIC_TMP/_GO_VCF.vcf \
                                                $GO_RR \
                                                ${GO_FILE%.vcf.gz}.str


   ## remove header from STRUCTURE input file (otherwise we get an error when running STRUCTURE in the next step)
   cat ${GO_FILE%.vcf.gz}.str | tail -n +2 > $SNIC_TMP/_GO_STR

   ## rename WGS sample names (sample  names in STRUCTURE can only be up to 11 characters long)
   cat $SNIC_TMP/_GO_STR | sed 's/^Nana_Finland_Enontekio.A002_2/NAN99/'  > ${GO_FILE%.vcf.gz}_CLEAN.str

   rm -f $SNIC_TMP/_samples.txt
   rm -f $SNIC_TMP/_GO_VCF.vcf
   rm -f $SNIC_TMP/_GO_STR
   rm -f $SNIC_TMP/_SL

done


echo
echo "Done!"






