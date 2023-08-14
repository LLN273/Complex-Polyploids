#!/bin/bash
#


########### Run IQ-TREE2 with simulated data
        

#### paths and folder names

#remember initial path
SRCDIR_INI=$(pwd)                                           	 

#input folder
AA=${1:?msg}	

#input file suffix
InFileSUF=${2:?msg}	

#output folder (root name)
RR=${3:?msg}

#number of bootstraps
BS=${4:?msg}

# output file extension
OFE=${5:?msg}

# substitution model
MODEL_IQTREE2=${6:?msg}

# model: List of rate heterogeneity among sites (when using Modelfinder and GHOST)
MODEL_IQTREE2_RATE=${7:?msg}

#number of genes/loci
NLOCUS=${8:?msg}	

# Outgroup sequence
OUTG=${9:?msg}

## Replicate
REP=${10:?msg}

## Iteration 
NALT=${11:?msg}

## ABC subfolder name
#ABCsn=${11:?msg}


#echo
#echo "input file suffix:" $InFileSUF
#echo "Input folder:" $AA
#echo "Replicate:" $REP
#echo 'Number of genes/loci:' $NLOCUS
#echo "Number of bootstraps replicates:" $BS
#echo 'Substitution model:' $MODEL_IQTREE2
#echo 'mrate:' $MODEL_IQTREE2_RATE
#echo 'Outgroup:' $OUTG




########### run IQ-TREE2 with ultrafast bootstrap approximation (--ufboot)

mkdir -p $RR/$REP
mkdir -p $SNIC_TMP/ALT${NALT}/$REP

for i in $( eval echo {0001..${NLOCUS}} ); do		# locus/genes; must include leading zeros (LOCUS_0001, LOCUS_02, etc)

   GO_MSA=data_${i}_${InFileSUF}

   AUX_fa=$GO_MSA
   AUX_prefix_out=${GO_MSA%-*}_UFBoot_${OFE}
   rsync -ah $AA/$REP/$GO_MSA $SNIC_TMP/ALT${NALT}/$REP/$AUX_fa

   rm -f $RR/$REP/${AUX_prefix_out}*

   cd $RR/$REP

   #echo
   #echo $REP $i $AUX_fa

   iqtree-omp --ufboot $BS \
              -s $SNIC_TMP/ALT${NALT}/$REP/$AUX_fa \
              --seqtype DNA \
              -m ${MODEL_IQTREE2} \
              -nt 1 \
              --seed 4545 \
              --prefix ${AUX_prefix_out} \
              --redo \
              -o $OUTG \
              -quiet


done 

# clean auxiliary files
rm -f $SNIC_TMP/ALT${NALT}/$REP/$AUX_fa
rm -f $RR/$REP/*.model.gz
rm -f $RR/$REP/*.mldist
rm -f $RR/$REP/*.ckp.gz
rm -f $RR/$REP/*.bionj
rm -f $RR/$REP/*.ufboot
rm -f $RR/$REP/*.splits.nex



















 








