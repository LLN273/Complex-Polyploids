#!/bin/bash
#
#SBATCH -J getfasta
#SBATCH -p core 
#SBATCH -n 1
#SBATCH -t 01:00:00
#SBATCH -A snic2022-22-589
#SBATCH -M snowy
#SBATCH --mail-user luis.leal@ebc.uu.se
#SBATCH --mail-type=FAIL



## Script used to fetch exonic sequences for each sample.

echo
echo "getfasta"
date -u

module load bioinfo-tools
module load BEDTools/2.29.2



#read file name (with path)
READ1=${1:?msg}

#read file name (without path)
InFile1=${2:?msg}		

#output folder (root name)
RR=${3:?msg}

#sample name
sampleName=${4:?msg}

# interval list (BED format) [target exons]
ANNF=${5:?msg}

# genomic feature (exons)
Gfeature=${6:?msg}

#remember initial path
SRCDIR_INI=$(pwd)                                           	



# uncompress fasta
zcat $READ1 > $SNIC_TMP/${InFile1%.gz}


echo
echo $sampleName
echo $READ1
echo $Gfeature
echo "Interval list:" $ANNF
echo


#### Get fasta for each locus listed in bed file
bedtools getfasta -name \
	                -fo $SNIC_TMP/${sampleName}.${Gfeature}_AUX.fa \
                  -fi $SNIC_TMP/${InFile1%.gz} \
	                -bed $ANNF


# notes (bedtools getfasta)
# -name			Use the “name” column in the BED file for the FASTA headers in the output FASTA file.
# -fo			Specify an output file name. By default, output goes to stdout.
# -fi 			<input FASTA>
# -bed 			<BED/GFF/VCF>



# clean sequence headers (issue with some headers in intron files)
sed 's/;zero_length_insertion=True//g' $SNIC_TMP/${sampleName}.${Gfeature}_AUX.fa > $RR/${sampleName}.${Gfeature}.fa


##### Remove auxiliary files
cd $SNIC_TMP
rm -f $SNIC_TMP/${InFile1%.gz}
rm -f $SNIC_TMP/${sampleName}.${Gfeature}_AUX.fa



echo
echo "Done!"
date -u
echo








