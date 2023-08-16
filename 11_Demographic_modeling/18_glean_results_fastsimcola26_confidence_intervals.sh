#!/bin/bash
#
#SBATCH -J fscCI_glean
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 06:00:00
#SBATCH -A snic2022-22-909
#SBATCH -M snowy
#SBATCH --mail-user luis.leal@ebc.uu.se
#SBATCH --mail-type=FAIL
ulimit -c unlimited


### Script used to glean fsc simulation results..." (output from 16_fastsimcola26_confidence_intervals_MAIN.sh)


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
EVOLM=bmt1





##### best gene flow model (for each population)
pen2pub0[1]=1		
pub2pen0[1]=0			
pen2pla0[1]=1
pla2pen0[1]=1	
pub2pla0[1]=1
pla2pub0[1]=1	
anc2pub1[1]=0
pub2anc1[1]=0

pen2pub0[2]=1		
pub2pen0[2]=0			
pen2pla0[2]=1
pla2pen0[2]=1	
pub2pla0[2]=1
pla2pub0[2]=1	
anc2pub1[2]=0
pub2anc1[2]=0

pen2pub0[3]=1		
pub2pen0[3]=0			
pen2pla0[3]=1
pla2pen0[3]=1	
pub2pla0[3]=1
pla2pub0[3]=1	
anc2pub1[3]=0
pub2anc1[3]=0

pen2pub0[4]=1		
pub2pen0[4]=0			
pen2pla0[4]=1
pla2pen0[4]=1	
pub2pla0[4]=1
pla2pub0[4]=1	
anc2pub1[4]=0
pub2anc1[4]=0



for k in `seq 1 1 4`; do 			# population	

   gf_code=${pen2pub0[$k]}_${pub2pen0[$k]}_${pen2pla0[$k]}_${pla2pen0[$k]}_${pub2pla0[$k]}_${pla2pub0[$k]}_${anc2pub1[$k]}_${pub2anc1[$k]}

   for i in `seq 1 1 100`; do 			# replicates
   
      GO_AA=$AA/${RRsub}/$DPfilter/${sample_set}/model${EVOLM}_confidence_intervals/${FEXT[$k]}/C1/${gf_code}/model${EVOLM}_${FEXT[$k]}/model${EVOLM}_${FEXT[$k]}_${i}
                   
      GO_RR=$AA/${RRsub}/$DPfilter/${sample_set}/model${EVOLM}_confidence_intervals/${FEXT[$k]}/C1/${gf_code}/model${EVOLM}_${FEXT[$k]}/model${EVOLM}_${FEXT[$k]}_${i}

      # Output file
      OUTfile=fsm_50X_AllResults.txt
	  
      # add header
      cat $GO_AA/1/model${EVOLM}_${FEXT[$k]}/model${EVOLM}_${FEXT[$k]}.bestlhoods | head -n1 > $GO_RR/$OUTfile


      ################# Collect fsc results
      for n in `seq 1 1 50`; do 	
	  
         GO_AA2=$GO_AA/${n}/model${EVOLM}_${FEXT[$k]}	

         echo ${FEXT[$k]} $i $n
         cat $GO_AA2/model${EVOLM}_${FEXT[$k]}.bestlhoods | tail -n1 >> $GO_RR/$OUTfile

      done
   
   
   
      ################# compute mean values
   
      OUTfile2="fsm_50X_meanValues.txt"

      Rscript --no-save $SRCDIR_INI/18_compute_meanValue.R $GO_RR/$OUTfile \
                                                           $OUTfile2 \
                                                           $GO_RR


  done
done


echo
echo "Done!"
date -u
echo



