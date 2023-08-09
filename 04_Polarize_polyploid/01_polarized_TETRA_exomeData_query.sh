#!/bin/bash
#
#SBATCH -J polar1
#SBATCH -p core 
#SBATCH -n 1
#SBATCH -t 06:00:00
#SBATCH -A snic2022-22-589
#SBATCH -M snowy
#SBATCH --mail-user luis.leal@ebc.uu.se
#SBATCH --mail-type=FAIL
ulimit -c unlimited



## Script used to polarize allotetraploid sequence in the MSA.
## Polarization done separately for each locus.


# load modules
module load bioinfo-tools
module load python3/3.8.7


#### paths and folder names

# remember initial path
SRCDIR_INI=$(pwd)      

# input folder
AA=${1:?msg}

# output folder 
RR=${2:?msg}

# Number of species
NSPEC=${3:?msg}

# Polyploid ID
POLYP=${4:?msg}

# Ref genome ID
REFspecies=${5:?msg}

## gene list
GL=${6:?msg}

## genomic feature
GO_FEATURE=${7:?msg}


echo
echo $GO_FEATURE
echo $AA
echo $RR
echo $GL
echo "Number of species:" $NSPEC
echo "Polyploid:" $POLYP
echo "Reference genome:" $REFspecies
echo



######## read gene list

i=1

while read -r LINE                                                                      
do
   GENE_LIST[$i]=$LINE
   let "i+=1"
done < ${GL}

N_GENES=${#GENE_LIST[@]}



cd $RR


for i in `seq 1 1 $N_GENES`; do

   # input MSA
   GO_MSA=${GENE_LIST[$i]}.${GO_FEATURE}_trimN.fa
   rsync -ah $AA/$GO_MSA $RR/$GO_MSA

   echo $i $GO_MSA

   # polarize sequences (first check whether source MSA is empty)
   if [ -s $RR/$GO_MSA ]
   then
      python3 01s_polarizeTETRA.py $GO_MSA $NSPEC $REFspecies $POLYP            
   else
      touch $RR/${GO_MSA%.fa}-ALT.fasta    	#if MSA file is empty
   fi

   rm -f $SNIC_TMP/AUX*
   rm -f $RR/$GO_MSA

done
