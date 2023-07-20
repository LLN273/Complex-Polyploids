#!/bin/bash
#
#SBATCH -J DBImport
#SBATCH -p core
#SBATCH -n 16
#SBATCH -t 08:00:00
#SBATCH -A snic2022-22-589
##SBATCH -A snic2021-22-727
#SBATCH -M snowy
#SBATCH --mail-user luis.leal@ebc.uu.se
#SBATCH --mail-type=FAIL
ulimit -c unlimited



### Script used to combine gVCF files for specimens associated to the same species


module load bioinfo-tools
module load GATK/4.2.0.0




########### input data files and paths

#remember initial path
SRCDIR_INI=$(pwd)   

# gVCF input folder
AA=${1:?msg}
	
#output file
RR=${2:?msg}

# reference genome
refGenome=${3:?msg}			

# reference genome folder
refGenomefolder=${4:?msg}

# sample list
SPL=${5:?msg}

# Output file name
outFile_root=${6:?msg}

# Intervals list (BED format) exome data 
ANNF=${7:?msg}

# Contig list
GO_CONTIGlist=${8?msg}





echo
echo "Sample list (BAM):" $SPL
echo "Intervals list:" $ANNF
echo "Contig list:" $GO_CONTIGlist
echo "Input folder (BAM):" $AA
echo "Output folder:" $RR
echo "output file (root):" db_${outFile_root}

echo "Reference genome:" $refGenome
echo "Reference genome folder:" $refGenomefolder
echo






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
paste -d'-' $SNIC_TMP/__auxA.txt $SNIC_TMP/__aux3.txt > $SNIC_TMP/__targets_aux.interval_list


### select subset of contigs
cat $GO_CONTIGlist | sort > $SNIC_TMP/__GO_CONTIGlist_sort.txt
grep -f $SNIC_TMP/__GO_CONTIGlist_sort.txt $SNIC_TMP/__targets_aux.interval_list > $SNIC_TMP/__targets.intervals









########### read sample list

i=1

while read -r LINE                                                                      
do
   SAMPLE_LIST[$i]=$LINE
   let "i+=1"
done < ${SPL}


### Number of samples
N_SAMPLES=${#SAMPLE_LIST[@]}





########### copy reference genome to scratch disk

cd $SNIC_TMP
rsync -ah $refGenomefolder/${refGenome%.fa*}* .



########### List of g.VCF files 

declare -A Infile_VCF

for i in `seq 1 1 $N_SAMPLES`; do                  
   Infile_VCF[$i]=${SAMPLE_LIST[$i]}-individualHaplotypeCaller.g.vcf
done




########### Create file with path to VCF files (this can be a simple list with path or a 'map file')
##### The sample map is a tab-delimited text file with sample_name--tab--path_to_sample_vcf per line.
##### Using a sample map saves the tool from having to download the GVCF headers in order to determine the sample names (but using vcfs with incompatible headers may result in silent data corruption.)

rm -f $SNIC_TMP/__vcf.sample_map
touch $SNIC_TMP/__vcf.sample_map

rm -f $SNIC_TMP/__vcf.list
touch $SNIC_TMP/__vcf.list

for i in `seq 1 1 $N_SAMPLES`; do
   echo -e ${AA}/${Infile_VCF[$i]} >> $SNIC_TMP/__vcf.list
   #echo -e "${SAMPLE_LIST[$i]}\t${AA}/${Infile_VCF[$i]}" >> $SNIC_TMP/__vcf.sample_map
done

echo
echo "Input data:"
cat $SNIC_TMP/__vcf.list
cat $SNIC_TMP/__vcf.sample_map
echo




########### Combine gVCF files with GenomicsDBImport

mkdir -p $SNIC_TMP/_TEMPgatk

if [[ $N_SAMPLES -lt 50 ]] ; then
   echo
   echo "Running on 4 cores"
   echo 
   # Use 4 cores if number of samples < 50 
   gatk --java-options "-Xmx25G -Xms25G" \
       GenomicsDBImport \
       --genomicsdb-workspace-path $RR/db_${outFile_root} \
       --variant $SNIC_TMP/__vcf.list \
       --intervals $SNIC_TMP/__targets.intervals \
       --reference $SNIC_TMP/$refGenome \
       --tmp-dir $SNIC_TMP/_TEMPgatk \
       --merge-input-intervals true \
       --genomicsdb-shared-posixfs-optimizations true

else
   echo
   echo "Running on 20 cores (16 on snowy)"
   echo 
   # Use whole node on if number of samples > 50
   gatk --java-options "-Xmx100G -Xms100G" \
       GenomicsDBImport \
       --genomicsdb-workspace-path $RR/db_${outFile_root} \
       --variant $SNIC_TMP/__vcf.list \
       --intervals $SNIC_TMP/__targets.intervals \
       --reference $SNIC_TMP/$refGenome \
       --tmp-dir $SNIC_TMP/_TEMPgatk \
       --merge-input-intervals true \
       --genomicsdb-shared-posixfs-optimizations true


fi




########### clean scratch disk
rm -r $SNIC_TMP/*.fa*
rm -r $SNIC_TMP/__vcf.list
rm -rf $SNIC_TMP/_TEMPgatk
rm -f $SNIC_TMP/__GO_CONTIGlist_sort.txt




