#!/bin/bash
#

### Script used to perform initial model optimization using simulated annealing.
### Script must be run independently for each polyploidization model, population, and chain 
### There is no cross-talk between chains



###################################### paths and folder names

################## Remember initial path
SRCDIR_INI=$(pwd)   

################## Folder containing simulated (unpolarized) MSAs (root)
AA=/crex1/proj/snic2020-6-184/nobackup/private/Luis/P11_SIMULATED_FASTA_PHYLOGENY_HYBRIDIZATION_exomeData/01_simphy_ILS

################## Simulation conditions (polyploidization models)
myMODEL[1]="PPPP"    # B. pendula, autopolyploid
myMODEL[2]="AAAA"    # ancestor of B. pendula/B. platyphylla, autopolyploid
myMODEL[3]="PPNN"    # B. pendula/B. nana, allopolyploid
myMODEL[4]="PPHH"    # B. pendula/B. humilis, allopolyploid
myMODEL[5]="PPPyPy"    # B. pendula/B. platyphylla, allopolyploid
myMODEL[6]="PPNH"    # B. pendula/B. nana/B. humilis, allopolyploid
myMODEL[7]="AANN"    # ancestor of B. pendula-B. platyphylla /B. nana, allopolyploid
myMODEL[8]="AAHH"    # ancestor of B. pendula-B. platyphylla /B. humilis, allopolyploid
myMODEL[9]="AANH"    # ancestor of B. pendula-B. platyphylla /B. nana/B. humilis, allopolyploid

################## Population
pop[1]="SVsouth"
pop[2]="Arctic"
pop[3]="Spain"
pop[4]="Central_Asia"
pop[5]="LT"
pop[6]="UA"
pop[7]="JOK"

################## Homoeologous Exchange across subgenomes  (1: allowed;  0: not allowed) >> valid for allopolyploid models only
HomoeologousExchange=0

################## Number of chains
N_chains=4

################## Prior flag
#Priorflag="user"		# priors provided by user (PP file, below)
Priorflag="random"		# random generated priors

################## Number of replicates
REP=25		

################## Number of loci/genes		
NLOCUS=50

################## ILS prior
ILS_PRIOR=40

################## pendula-platyphylla gene flow (1: allowed;  0: not allowed)
geneFlowFlag=1







for k in `seq 1 1 1`; do               # models

   for n in `seq 1 1 7`; do               # population

      for j in `seq 1 1 $N_chains`; do               # chain
	  
         GO_myMODEL=${myMODEL[$k]}
         GO_pop=${pop[$n]}
         chain_number=${j}
		 
         ################## Observed pairing profiles (each line in input file correspond to a different polarization geometry, in the following order: B. pendula; B. nana; B. humilis; B. platyphylla)
         PP_Obs=/crex1/proj/snic2017-7-149/private/Luis/P11_SIMULATED_FASTA_PHYLOGENY_HYBRIDIZATION_exomeData/02_ABC_Simulated_annealing/00_PAIRING_PROFILES_${GO_pop}.txt	

         ##################  Output folder (root); chain number is the last digit in the path

         if [[ ${HomoeologousExchange} = 0 ]] ; then
            RR_MAIN=/crex1/proj/snic2020-6-184/nobackup/private/Luis/P11_SIMULATED_FASTA_PHYLOGENY_HYBRIDIZATION_exomeData/02_ABC_Simulated_annealing/${GO_myMODEL}/${GO_pop}/${chain_number}
            echo ""
         fi

         if [[ ${HomoeologousExchange} = 1 ]] ; then
            RR_MAIN=/crex1/proj/snic2020-6-184/nobackup/private/Luis/P11_SIMULATED_FASTA_PHYLOGENY_HYBRIDIZATION_exomeData/02_ABC_Simulated_annealing/${GO_myMODEL}_HE/${GO_pop}/${chain_number}
            echo ""
         fi

         rm -rf $RR_MAIN
         mkdir -p $RR_MAIN
		 

         ################## Initial priors (when provided by the user)
         #PP=/crex1/proj/snic2017-7-149/private/Luis/P11_SIMULATED_FASTA_PHYLOGENY_HYBRIDIZATION_exomeData/02_ABC_Simulated_annealing/PRIORS_simAnnealing_userDefined/00_PRIORS_ILS_gf_H1_H7_SIMANN_${GO_myMODEL}_${GO_pop}.txt    		
         PP=""

         if [[ ${PP} = "" ]] ; then
            PP=$SRCDIR_INI/PPempty.txt
            if [[ -e "$PP" ]]; then
               echo $PP "file exists."
            else
               touch $PP
            fi
         fi


         echo "MODEL:" $GO_myMODEL
         echo "Population:" $GO_pop
         echo "Chain:" $chain_number
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

         sbatch $SRCDIR_INI/05_ABC_SimAnnealing_query.sh $AA \
                                                        $RR_MAIN \
                                                        $PP_Obs \
                                                        $REP \
                                                        $NLOCUS \
                                                        $ILS_PRIOR \
                                                        $GO_myMODEL \
                                                        $HomoeologousExchange \
                                                        $geneFlowFlag \
                                                        $Priorflag \
                                                        $PP
                                            
											
      done
   done
done




