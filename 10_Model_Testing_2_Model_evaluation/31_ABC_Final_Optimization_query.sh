#!/bin/bash
#
#SBATCH -J ABC
#SBATCH -p core
#SBATCH -n 20
#SBATCH -t 48:00:00
#SBATCH -A snic2022-22-909
#SBATCH --mail-user luis.leal@ebc.uu.se
#SBATCH --mail-type=FAIL
ulimit -c unlimited


#### Script used to perform final optimization of model parameters using ABC


echo
echo "Starting Uppmax jobs ..."
date -u
echo


# load modules
module load bioinfo-tools
module load R/4.0.0
module load python3/3.8.7
module load iqtree/2.0-rc2-omp-mpi
echo





###################################### paths and folder names

# Remember initial path
SRCDIR_INI=$(pwd)   

# input folder
AA=${1:?msg}

# output folder
RR_MAIN=${2:?msg}

# number of replicates
REP=${3:?msg}

# number of loci
NLOCUS=${4:?msg}

# polyploidization model
myMODEL=${5:?msg}

# folder containing prior lists
GO_PF=${6:?msg}
					                        
# index for first set of 100 prior sets                              #              
prior_index_start=${7:?msg}


echo
echo "Input folder:" $AA
echo "Output folder:" $RR_MAIN
echo "MODEL:" $myMODEL
echo "Number of replicates:" $REP
echo "Number of loci:" $NLOCUS
echo "Folder containing priors" $GO_PF
echo "index first prior set:" $prior_index_start



##### Read priors

unset H_PENDPLATYGF
x=0

while read -r LINE                                                                      
do
   H_PENDPLATYGF[$x]=$LINE
   let "x+=1"
done < $GO_PF/platy_pend_gene_flow.txt

N_SAMPLES=${#H_PENDPLATYGF[@]}



unset H_0
x=0

while read -r LINE                                                                      
do
   H_0[$x]=$LINE
   let "x+=1"
done < $GO_PF/H0.txt


   
unset H_1
x=0

while read -r LINE                                                                      
do
   H_1[$x]=$LINE
   let "x+=1"
done < $GO_PF/H1.txt



unset H_2
x=0

while read -r LINE                                                                      
do
   H_2[$x]=$LINE
   let "x+=1"
done < $GO_PF/H2.txt



unset H_3
x=0

while read -r LINE                                                                      
do
   H_3[$x]=$LINE
   let "x+=1"
done < $GO_PF/H3.txt



unset H_4
x=0

while read -r LINE                                                                      
do
   H_4[$x]=$LINE
   let "x+=1"
done < $GO_PF/H4.txt



unset H_5
x=0

while read -r LINE                                                                      
do
   H_5[$x]=$LINE
   let "x+=1"
done < $GO_PF/H5.txt



unset H_6
x=0

while read -r LINE                                                                      
do
   H_6[$x]=$LINE
   let "x+=1"
done < $GO_PF/H6.txt



unset H_7
x=0

while read -r LINE                                                                      
do
   H_7[$x]=$LINE
   let "x+=1"
done < $GO_PF/H7.txt



unset H_8
x=0

while read -r LINE                                                                      
do
   H_8[$x]=$LINE
   let "x+=1"
done < $GO_PF/H8.txt




unset H_9
x=0

while read -r LINE                                                                      
do
   H_9[$x]=$LINE
   let "x+=1"
done < $GO_PF/H9.txt



unset ILS_REDUCTION
x=0

while read -r LINE                                                                      
do
   ILS_REDUCTION[$x]=$LINE
   let "x+=1"
done < $GO_PF/ILSreduction.txt




####### Perform simulations

aux_1=0
prior_set_start="$(awk '{print $1+$2}' <<<"$prior_index_start $aux_1")"
aux_99=99
prior_set_end="$(awk '{print $1+$2}' <<<"$prior_set_start $aux_99")"



for n in `seq $prior_set_start 1 $prior_set_end`; do               # prior sets

   GO_ILS_reduction=${ILS_REDUCTION[$n]}
   GO_AA=${AA}
   GO_RR=$RR_MAIN/$n
   mkdir -p $GO_RR

   rm -f $GO_RR/00_PRIORS.txt
   echo ${ILS_REDUCTION[$n]} >> $GO_RR/00_PRIORS.txt
   echo ${H_PENDPLATYGF[$n]} >> $GO_RR/00_PRIORS.txt
   echo ${H_0[$n]} >> $GO_RR/00_PRIORS.txt
   echo ${H_1[$n]} >> $GO_RR/00_PRIORS.txt
   echo ${H_2[$n]} >> $GO_RR/00_PRIORS.txt
   echo ${H_3[$n]} >> $GO_RR/00_PRIORS.txt
   echo ${H_4[$n]} >> $GO_RR/00_PRIORS.txt
   echo ${H_5[$n]} >> $GO_RR/00_PRIORS.txt
   echo ${H_6[$n]} >> $GO_RR/00_PRIORS.txt
   echo ${H_7[$n]} >> $GO_RR/00_PRIORS.txt
   echo ${H_8[$n]} >> $GO_RR/00_PRIORS.txt
   echo ${H_9[$n]} >> $GO_RR/00_PRIORS.txt
   
   echo
   echo $n
   echo $GO_ILS_reduction
   echo $GO_AA
   echo $GO_RR
   echo ${ILS_REDUCTION[$n]} ${H_PENDPLATYGF[$n]} ${H_0[$n]} ${H_1[$n]} ${H_2[$n]} ${H_3[$n]} ${H_4[$n]} ${H_5[$n]} ${H_6[$n]} ${H_7[$n]} ${H_8[$n]} ${H_9[$n]}
   date -u

                                            									
   ./20_POLARIZATION_IQTREE_PAIRINGprofiles_MAIN.sh $GO_AA \
                                                    $GO_RR \
                                                    $myMODEL \
                                                    $REP \
                                                    $NLOCUS \
                                                    $GO_RR/00_PRIORS.txt


done


echo 
echo "Done!"
date -u

