#!/bin/bash
#
#SBATCH -J L2norm
#SBATCH -p devcore
#SBATCH -n 1
#SBATCH -t 1:00:00
#SBATCH -A snic2022-22-909
#SBATCH -M snowy
#SBATCH --mail-user luis.leal@ebc.uu.se
#SBATCH --mail-type=FAIL
ulimit -c unlimited


### Script used to compute L2 distance for each ABC simulation run

echo
echo "Starting Uppmax jobs ..."
date -u
echo


module load bioinfo-tools
module load python3/3.8.7



###################################### paths and folder names

################## Remember initial path
SRCDIR_INI=$(pwd)   

################## Folder containing results from ABC analysis
AA=/crex1/proj/snic2020-6-184/nobackup/private/Luis/P11_SIMULATED_FASTA_PHYLOGENY_HYBRIDIZATION_exomeData/03_ABC_Final_optimization


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



################## priors folder
PP=/crex1/proj/snic2017-7-149/private/Luis/P11_SIMULATED_FASTA_PHYLOGENY_HYBRIDIZATION_exomeData/02_ABC_Simulated_annealing/PRIORS_ABC 

################# Input file (results from ABC simulations)
Infile=summary_ABC_1000_simulations.txt



################################## Compute L2 distance

for k in `seq 1 1 7`; do               # populations

   AA_MAIN=${AA}/${myMODEL}/${pop[$k]}
   PP_GO=${PP}/ABC_PRIORS_${myMODEL}_${pop[$k]}
   RR_MAIN=$AA_MAIN

   ################## Observed pairing profiles
   PP_Obs=$SRCDIR_INI/00_PAIRING_PROFILES_${pop[$k]}.txt		


   echo
   echo $myMODEL $REP $NLOCUS
   echo $AA_MAIN
   echo $RR_MAIN
   echo $PP_GO
   echo

   python3 33_L2_ABC_results.py $AA_MAIN \
                                $RR_MAIN \
                                $PP_Obs \
                                $PP_GO \
                                $Infile 


   # clean output file (remove square brackets; replace multiple spaces by single space)
   sed 's/[[]//g' $RR_MAIN/L2-norm_ABC_1000_simulations.txt > $SNIC_TMP/_aux1.txt
   sed 's/[]]//g' $SNIC_TMP/_aux1.txt | tr -s ' ' > $RR_MAIN/L2-norm_ABC_1000_simulations.txt
   rm -f $SNIC_TMP/_aux1.txt
                                            								
done


echo
echo "Done!" 
date -u
