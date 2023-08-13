#!/bin/bash
#

## Script used to decrease ILS by replacing random genes with genes reflecting true gene phylogeny



# load modules
module load bioinfo-tools
module load R/4.0.0




#### paths and folder names

#remember initial path
SRCDIR_INI=$(pwd)                                           	 

#input folder
GO_AA=${1:?msg}	

#output folder (root name)
GO_RR=${2:?msg}

#Replicate
r_aux=${3:?msg}

# number of loci/genes	
NLOCUS=${4:?msg}

## new gene phylogeny (must replicate species tree if we aim at reducing ILS)
GP_GOLD=${5:?msg}

## proportion of genes to be replaced [0-1]
GP_prop=${6:?msg}

echo
echo 'Replicate:' $r_aux
echo 'Number of loci/genes:' $NLOCUS
echo 'Input folder:' $GO_AA
echo 'Output folder:' $GO_RR
echo



######### Randomize gene list

rm -f $SNIC_TMP/_00_GENE_LIST.txt

for i in $( eval echo {0001..${NLOCUS}} ); do				
   echo $i >> $SNIC_TMP/_00_GENE_LIST.txt
done

cat $SNIC_TMP/_00_GENE_LIST.txt | shuf > $SNIC_TMP/_GENE_LIST.txt



## read randomized gene list

i=1

while read -r LINE                                                                      
do
   GENE_LIST[$i]=$LINE
   let "i+=1"
done < $SNIC_TMP/_GENE_LIST.txt


####################### number of events
GP_TOTAL="$(awk '{print $1*$2}' <<<"$NLOCUS $GP_prop")"


cd $GO_RR

rsync -ah $GO_AA/l_trees.trees .
rsync -ah $GO_AA/s_tree.trees .
rsync -ah $GO_AA/g_*.trees .
rsync -ah $GP_GOLD ./00_reference_GENE.trees


### replace fraction of ILS genes by reference species tree (no ILS)

GP_stop_round="$(printf "%.0f\n" $GP_TOTAL )"	# round to closest integer

for i in `seq 1 1 $GP_stop_round`; do					

   l_aux=${GENE_LIST[$i]}
   echo $i $l_aux "1"

   ## gene tree
   InFile="g_trees${l_aux}.trees"

   #replace gene tree
   cp ./00_reference_GENE.trees $InFile

done


# clean auxiliary files
rm -f $GO_RR/00_reference_GENE.trees

