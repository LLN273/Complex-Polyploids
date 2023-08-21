#!/bin/bash
#
#SBATCH -J ABC_SA
#SBATCH -p core
#SBATCH -n 20
#SBATCH -t 240:00:00
#SBATCH -A naiss2023-22-717
#SBATCH --mail-user luis.leal@ebc.uu.se
#SBATCH --mail-type=FAIL
ulimit -c unlimited



### Script used to perform initial model optimization using simulated annealing.
### Script must be run independently for each polyploidization model, population, and chain 
### There is no cross-talk between chains


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

################## Remember initial path
SRCDIR_INI=$(pwd)   

################## Folder containing simulated (unpolarized) MSAs (root)
AA=${1:?msg}

##################  Output folder (root); chain number is the last digit in the path
RR_MAIN=${2:?msg}

################## Observed pairing profiles (each line in input file correspond to a different polarization geometry, in the following order: B. pendula; B. nana; B. humilis; B. platyphylla)
PP_Obs=${3:?msg}

################## Number of replicates
REP=${4:?msg}

################## Number of loci/genes		
NLOCUS=${5:?msg}

################## ILS prior
ILS_PRIOR=${6:?msg}

################## Simulation conditions (polyploidization model)
myMODEL=${7:?msg}

################## Homoeologous Exchange across subgenomes  (1: allowed;  0: not allowed) >> valid for allopolyploid models only
HomoeologousExchange=${8:?msg}

################## pendula-platyphylla gene flow (1: allowed;  0: not allowed)
geneFlowFlag=${9:?msg}

################## Prior flag
Priorflag=${10:?msg}	

################## Initial priors (when provided by the user)
PP=${11:?msg}



echo
echo "MODEL:" $myMODEL
echo "Input folder:" $AA
echo "Output folder:" $RR_MAIN
echo "Observed pairing profiles:" $PP_Obs
echo "Number of replicates:" $REP
echo "Number of loci:" $NLOCUS
echo "ILS prior:" $ILS_PRIOR
if [[ ${Priorflag} = "random" ]] ; then
   echo "Initial priors: random generated"
fi
if [[ ${Priorflag} = "user" ]] ; then
   echo "Initial priors:" $PP
fi
echo "HomoeologousExchange flag:" $HomoeologousExchange
echo "pendula-platyphylla gene flow flag:" $geneFlowFlag
echo


################# Optimize model

python3 05_ABC_SimAnnealing_optimization.py $AA \
                                            $RR_MAIN \
                                            $PP_Obs \
                                            $REP \
                                            $NLOCUS \
                                            $ILS_PRIOR \
                                            $myMODEL \
                                            $HomoeologousExchange \
                                            $geneFlowFlag \
                                            $Priorflag \
                                            $PP
                                            
											
echo 
echo "ABC_SimAnnealing has finished."
date -u

