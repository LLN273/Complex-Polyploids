#!/bin/bash

## Script used to fix tree tips so that the two (haploid) individuals belonging to the same species are monophyletic 
## (as these will be recoded as one single diploid individual at a later stage)



#### paths and folder names

#remember initial path
SRCDIR_INI=$(pwd)                                           	 

# Input folder (root) [results from simphy]
AA=/crex1/proj/snic2020-6-184/nobackup/private/Luis/P11_SIMULATED_FASTA_PHYLOGENY_HYBRIDIZATION_exomeData/01_simphy_ILS

# Input/output subfolder (different simulation conditions)
INfolder[1]=A_modILS_birch

## Sample list
SPL=/crex1/proj/snic2017-7-149/private/Luis/P11_SIMULATED_FASTA_PHYLOGENY_HYBRIDIZATION_exomeData/01_simphy_ILS/00_sample_list_ABCD_MASTER.txt

# Number of species
NSPEC[1]=14

## number of replicates
REP=25		

# number of loci/genes		
NLOCUS=50     


echo

for k in `seq 1 1 1`; do 			# different simulation conditions

      echo
      GO_NSPEC=${NSPEC[$k]}
      GO_FOLDER=${INfolder[$k]}
      GO_AA=${AA}/${GO_FOLDER}
      GO_RR=${AA}/${GO_FOLDER}_CLEAN
      mkdir -p $GO_RR

      cd $GO_RR
      rsync -ah $GO_AA/*.command .
      rsync -ah $GO_AA/*.db .
      rsync -ah $GO_AA/*.params .


      for r_aux in $( eval echo {01..${REP}} ); do		

         mkdir -p $GO_RR/$r_aux
         cd $GO_RR/$r_aux
         rsync -ah $SRCDIR_INI/02a_fix_tree_tips_birch.R .

         echo $r_aux $GO_NSPEC $NLOCUS
         echo $GO_AA/$r_aux
         echo $GO_RR/$r_aux
         echo $SPL

         cd $SRCDIR_INI
         sbatch $SRCDIR_INI/02_fix_tree_tips_lowILS_query.sh $GO_AA/$r_aux \
                                                      $GO_RR/$r_aux \
                                                      $r_aux \
                                                      $GO_NSPEC \
                                                      $NLOCUS \
                                                      $SPL
                                          
         echo

      done
done



