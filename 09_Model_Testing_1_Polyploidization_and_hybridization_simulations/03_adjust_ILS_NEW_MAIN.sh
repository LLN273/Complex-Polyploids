#!/bin/bash
#
#SBATCH -J ILSfix
#SBATCH -p core 
#SBATCH -n 1
#SBATCH -t 01:00:00
#SBATCH -A snic2022-22-909
#SBATCH -M snowy
#SBATCH --mail-user luis.leal@ebc.uu.se
#SBATCH --mail-type=FAIL
ulimit -c unlimited



## Script used to decrease ILS by replacing random genes with genes reflecting true gene phylogeny



#### paths and folder names

#remember initial path
SRCDIR_INI=$(pwd)                                           	 

# Input folder (root) [results ASTRAL analysis]
AA=/crex1/proj/snic2020-6-184/nobackup/private/Luis/P11_SIMULATED_FASTA_PHYLOGENY_HYBRIDIZATION_exomeData/01_simphy_ILS

# Input subfolder (different simulation conditions)
INfolder[1]=A_modILS_birch_CLEAN
INfolder[2]=${INfolder[1]}
INfolder[3]=${INfolder[1]}
INfolder[4]=${INfolder[1]}
INfolder[5]=${INfolder[1]}

# Output subfolder (different simulation conditions)
OUTfolder[1]=A_modILS_birch_CLEAN_reducedILS_30perct     # 30% genes replaced
OUTfolder[2]=A_modILS_birch_CLEAN_reducedILS_35perct     # 35% genes replaced
OUTfolder[3]=A_modILS_birch_CLEAN_reducedILS_40perct     # 40% genes replaced
OUTfolder[4]=A_modILS_birch_CLEAN_reducedILS_45perct     # 45% genes replaced
OUTfolder[5]=A_modILS_birch_CLEAN_reducedILS_50perct     # 50% genes replaced

## new gene phylogeny (must replicate species tree if we aim at reducing ILS)
GP_GOLD=/crex1/proj/snic2017-7-149/private/Luis/P11_SIMULATED_FASTA_PHYLOGENY_HYBRIDIZATION_exomeData/01_simphy_ILS/00_reference_GENE_2BpubB.trees


## proportion of genes to be replaced
GP_prop[1]="0.3"
GP_prop[2]="0.35"
GP_prop[3]="0.4"
GP_prop[4]="0.45"
GP_prop[5]="0.5"

## number of replicates
REP=25		

# number of loci/genes		
NLOCUS=50


echo

for k in `seq 1 1 5`; do 			# different simulation conditions

      echo
      GO_FOLDER=${INfolder[$k]}
      GO_OUT=${OUTfolder[$k]}
      GO_AA=${AA}/${GO_FOLDER}
      GO_RR=${AA}/${GO_OUT}
      mkdir -p $GO_RR


      cd $GO_RR
      rsync -ah $GO_AA/*.command .
      rsync -ah $GO_AA/*.db .
      rsync -ah $GO_AA/*.params .


      for r_aux in $( eval echo {01..${REP}} ); do

         mkdir -p $GO_RR/$r_aux

         echo
		     echo
         echo $r_aux $NLOCUS
         echo ${GP_prop[$k]} $GP_GOLD
         echo $GO_AA/$r_aux
         echo $GO_RR/$r_aux
         

         cd $SRCDIR_INI
         $SRCDIR_INI/03_adjust_ILS_NEW_query.sh $GO_AA/$r_aux \
                                                      $GO_RR/$r_aux \
                                                      $r_aux \
                                                      $NLOCUS \
                                                      $GP_GOLD \
                                                      ${GP_prop[$k]}
                                          





         echo

      done
done




