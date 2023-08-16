#!/bin/bash
#
#SBATCH -J fsc_glean
#SBATCH -p devcore
#SBATCH -n 1
#SBATCH -t 00:05:00
#SBATCH -A snic2022-22-909
#SBATCH -M snowy
#SBATCH --mail-user luis.leal@ebc.uu.se
#SBATCH --mail-type=FAIL
ulimit -c unlimited


### Script used to glean fsc simulation results..." (output from 10_fastsimcola26_50X_MAIN.sh)

echo
echo "Gleaning fsc simulation results..." 
date -u


##### load modules
module load bioinfo-tools
module load R/4.1.1
module load R_packages/4.1.1


############### input data files and paths

#remember initial path
SRCDIR_INI=$(pwd)  

### Sample set
sample_set="4_8"

#### input/output subfolder 1
RRsub=exome_WGS

#### Depth filtering
DPfilter="8X"

##### input/output extension
FEXT[1]="PubCentralEurope"
FEXT[2]="PubSVsouth"
FEXT[3]="PubCentralAsia"
FEXT[4]="PubSpain"

# Input folder (observed SFS)
AA=/crex/proj/snic2020-6-184/private/Luis/P06_birch_exome_IGA-Sweden_2019/32c_fastsimcoal2_Bpubescens_DivergenceTime_simulations/BWA_unmasked

######## Evolutionary models
EVOLM[1]=BP1	# ((PEN, PUB), PLA) and PEN->PUB migration only; d=8
EVOLM[2]=BP2	# ((PEN, PLA), PUB) and PEN->PUB migration only; d=8
EVOLM[3]=BP3	# (PEN, (PUB, PLA)) and PEN->PUB migration only; d=8


for i in `seq 1 1 3`; do 			# model
   for k in `seq 1 1 4`; do 			# population	
      
      GO_RR=$AA/${RRsub}/$DPfilter/${sample_set}/model${EVOLM[$i]}/${FEXT[$k]}/C1

      # Output file
      OUTfile=fsm_50X_ALL_RESULTS.txt

      # add header
      cat $AA/${RRsub}/$DPfilter/${sample_set}/model${EVOLM[$i]}/${FEXT[$k]}/C1/1/model${EVOLM[$i]}_${FEXT[$k]}/model${EVOLM[$i]}_${FEXT[$k]}.bestlhoods | head -n1 > $GO_RR/$OUTfile


      ################# Collect fsc results
      for n in `seq 1 1 50`; do 	
	  
	       GO_AA=$AA/${RRsub}/$DPfilter/${sample_set}/model${EVOLM[$i]}/${FEXT[$k]}/C1/${n}/model${EVOLM[$i]}_${FEXT[$k]}	

         echo $n
         cat $GO_AA/model${EVOLM[$i]}_${FEXT[$k]}.bestlhoods | tail -n1 >> $GO_RR/$OUTfile

      done
   
   
   
      ################# compute mean values
   
      OUTfile2="fsm_50X_meanValues.txt"

      Rscript --no-save $SRCDIR_INI/11_compute_meanValue.R $GO_RR/$OUTfile \
                                                           $OUTfile2 \
                                                           $GO_RR


  done
done



echo
echo "Done!"
date -u
echo



