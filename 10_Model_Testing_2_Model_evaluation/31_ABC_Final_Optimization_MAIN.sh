#!/bin/bash


#### Script used to perform final optimization of model parameters using ABC



###################################### paths and folder names

################## Remember initial path
SRCDIR_INI=$(pwd)   

################## Folder containing (unpolarized) MSAs (root)
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

# models that include Homoeologous Exchange across subgenomes
#myMODEL="PPNN_HE"    # B. pendula/B. nana, allopolyploid. 
#myMODEL="PPHH_HE"    # B. pendula/B. humilis, allopolyploid
#myMODEL="PPPyPy_HE"    # B. pendula/B. platyphylla, allopolyploid
#myMODEL="PPNH_HE"    # B. pendula/B. nana/B. humilis, allopolyploid
#myMODEL="AANN_HE"    # ancestor of B. pendula-B. platyphylla /B. nana, allopolyploid
#myMODEL="AAHH_HE"    # ancestor of B. pendula-B. platyphylla /B. humilis, allopolyploid
#myMODEL="AANH_HE"    # ancestor of B. pendula-B. platyphylla /B. nana/B. humilis, allopolyploid


################## Population
pop[1]="SVsouth"
pop[2]="Arctic"
pop[3]="Spain"
pop[4]="Central_Asia"
pop[5]="LT"
pop[6]="UA"
pop[7]="JOK"



##################  Output folder (root)
RR=/crex1/proj/snic2020-6-184/nobackup/private/Luis/P11_SIMULATED_FASTA_PHYLOGENY_HYBRIDIZATION_exomeData/03_ABC_Final_optimization

################## priors folder (generated using 30_ABC_generate_priors_based_on_SA_results_MAIN.py script based on results obtained from simulated annealing optimization)
PP=/crex1/proj/snic2017-7-149/private/Luis/P11_SIMULATED_FASTA_PHYLOGENY_HYBRIDIZATION_exomeData/02_ABC_Simulated_annealing/PRIORS_ABC   		

################## Number of replicates
REP=25		

################## Number of loci/genes		
NLOCUS=50



### Run batches of 100 simulations at a time (each simulation run is informed by a different prior set)

for k in `seq 1 1 7`; do               # populations
   for n in `seq 0 100 999`; do               # prior sets

      RR_MAIN=${RR}/${myMODEL}/${pop[$k]}
      #PP_GO=${PP}/ABC_PRIORS_${myMODEL}_${pop[$k]}
      PP_GO=${PP}
      mkdir -p $RR_MAIN
   
      echo
      echo $n $myMODEL $REP $NLOCUS
      echo $AA
      echo $RR_MAIN
      echo $PP_GO
   
      # Copy priors to output folder
      cd $RR_MAIN
      rsync -ahv $PP_GO .
   
      cd $SRCDIR_INI

      sbatch 31_ABC_Final_Optimization_query.sh $AA \
                                                $RR_MAIN \
                                                $REP \
                                                $NLOCUS \
                                                $myMODEL \
                                                $PP_GO \
                                                $n
                                            
   done											
done






