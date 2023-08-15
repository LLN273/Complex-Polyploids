#!/bin/bash
#
#SBATCH -J gleanABC
#SBATCH -p devcore
#SBATCH -n 1
#SBATCH -t 1:00:00
#SBATCH -A snic2022-22-909
#SBATCH -M snowy
#SBATCH --mail-user luis.leal@ebc.uu.se
#SBATCH --mail-type=FAIL
ulimit -c unlimited


### Script used to collect results from ABC simulations


echo
echo "Starting Uppmax jobs ..."
date -u
echo



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



################# output file (results from ABC simulations)
OUTfile=summary_ABC_1000_simulations.txt



### Collect results from ABC simulations

for k in `seq 1 1 7`; do               # populations

   AA_MAIN=${AA}/${myMODEL}/${pop[$k]}
   RR_MAIN=$AA_MAIN

   rm -f $RR_MAIN/$OUTfile

   echo "ABC_No" "pendula" "platyphylla" "pend_platy" "populifolia" "pop_out" "nana" "nana_out" "occidentalis" "occ_out" "humilis" "albosinensis" "lenta" "albo_len_hum" "albo_len" "basal" "pendula" "platyphylla" "pend_platy" "populifolia" "pop_out" "nana" "nana_out" "occidentalis" "occ_out" "humilis" "albosinensis" "lenta" "albo_len_hum" "albo_len" "basal" "pendula" "platyphylla" "pend_platy" "populifolia" "pop_out" "nana" "nana_out" "occidentalis" "occ_out" "humilis" "albosinensis" "lenta" "albo_len_hum" "albo_len" "basal" "pendula" "platyphylla" "pend_platy" "populifolia" "pop_out" "nana" "nana_out" "occidentalis" "occ_out" "humilis" "albosinensis" "lenta" "albo_len_hum" "albo_len" "basal" >> $RR_MAIN/$OUTfile

   echo
   echo $myMODEL $REP $NLOCUS
   echo $AA_MAIN
   echo $RR_MAIN
   echo

   for n in `seq 0 1 999`; do               # prior sets

      InFile=$AA_MAIN/$n/03_PAIRINGfrequencies/sister_ID_analysis_priors-JOINT_summary.txt

      sim_pen="$(cat $InFile | head -n1)"
      sim_nana="$(cat $InFile | tail -n3| head -n1)"
      sim_hum="$(cat $InFile | tail -n2| head -n1)"
      sim_platy="$(cat $InFile | tail -n1)"

      echo $n $sim_pen $sim_nana $sim_hum $sim_platy >> $RR_MAIN/$OUTfile
                                  
   done											
done


echo
echo "Done!" 
date -u

