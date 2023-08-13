#!/bin/bash
#
#SBATCH -J SA
#SBATCH -p core
#SBATCH -n 20
#SBATCH -t 100:00:00
#SBATCH -A snic2022-22-909
#SBATCH --mail-user luis.leal@ebc.uu.se
#SBATCH --mail-type=FAIL
ulimit -c unlimited



### Script used to perform initial model optimization using simulated annealing.
### Script must be run independently for each polyploidization model, population, and chain 



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

################## Folder containing simulated MSAs (root), prior to polarization
AA=/crex1/proj/snic2020-6-184/nobackup/private/Luis/P11_SIMULATED_FASTA_PHYLOGENY_HYBRIDIZATION_exomeData/01_simphy_ILS


################## Simulation conditions (polyploidization models)
myMODEL="PPPP"    # B. pendula, autopolyploid
#myMODEL="AAAA"    # ancestor of B. pendula/B. platyphylla, autopolyploid
#myMODEL="PPNN"    # B. pendula/B. nana, allopolyploid
#myMODEL="PPHH"    # B. pendula/B. humilis, allopolyploid
#myMODEL="PPPyPy"    # B. pendula/B. platyphylla, allopolyploid
#myMODEL="PPNH"    # B. pendula/B. nana/B. humilis, allopolyploid
#myMODEL="AANN"    # ancestor of B. pendula-B. platyphylla /B. nana, allopolyploid
#myMODEL="AAHH"    # ancestor of B. pendula-B. platyphylla /B. humilis, allopolyploid
#myMODEL="AANH"    # ancestor of B. pendula-B. platyphylla /B. nana/B. humilis, allopolyploid


################## Population
pop="SVsouth"
#pop="Arctic"
#pop="Spain"
#pop="Central_Asia"
#pop="LT"
#pop="UA"
#pop="JOK"


################## Homoeologous Exchange across subgenomes  (1: allowed;  0: not allowed) >> valid for allopolyploid models only
HomoeologousExchange=0

################## Chain number (update as required)
chain_number=1

##################  Output folder (root)

if [[ ${HomoeologousExchange} = 0 ]] ; then
	RR_MAIN=/crex1/proj/snic2020-6-184/nobackup/private/Luis/P11_SIMULATED_FASTA_PHYLOGENY_HYBRIDIZATION_exomeData/02_ABC_Simulated_annealing/${myMODEL}/${pop}/${chain_number}
	echo ""
fi

if [[ ${HomoeologousExchange} = 1 ]] ; then
   RR_MAIN=/crex1/proj/snic2020-6-184/nobackup/private/Luis/P11_SIMULATED_FASTA_PHYLOGENY_HYBRIDIZATION_exomeData/02_ABC_Simulated_annealing/${myMODEL}_HE/${pop}/${chain_number}
   echo ""
fi



################## Observed pairing profiles (each line in input file correspond to a different polarization geometry, in the following order: B. pendula; B. nana; B. humilis; B. platyphylla)
PP_Obs=/crex1/proj/snic2017-7-149/private/Luis/P11_SIMULATED_FASTA_PHYLOGENY_HYBRIDIZATION_exomeData/02_ABC_Simulated_annealing/00_PAIRING_PROFILES_${pop}.txt		

################## Prior flag
#Priorflag="user"		# priors provided by user (PP file, above)
Priorflag="random"		# random generated priors

################## Initial priors (when provided by the user)
#PP=/crex1/proj/snic2017-7-149/private/Luis/P11_SIMULATED_FASTA_PHYLOGENY_HYBRIDIZATION_exomeData/02_ABC_Simulated_annealing/PRIORS_simAnnealing_userDefined/00_PRIORS_ILS_gf_H1_H7_SIMANN_${myMODEL}_${pop}.txt    		# PPPP model
PP=""

################## Number of replicates
REP=25		

################## Number of loci/genes		
NLOCUS=50

################## ILS prior
ILS_PRIOR=40

################## pendula-platyphylla gene flow (1: allowed;  0: not allowed)
geneFlowFlag=1




##################################### Start model optimization

rm -rf $RR_MAIN
mkdir -p $RR_MAIN

if [[ ${PP} = "" ]] ; then
   PP=$SRCDIR_INI/PPempty.txt
   touch $PP
fi


echo
echo "Input folder:" $AA
echo "Output folder:" $RR_MAIN
echo "Observed pairing profiles:" $PP_Obs
echo "MODEL:" $myMODEL
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


python3 05_SimAnnealing_optimization.py $AA \
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
                                            
											
