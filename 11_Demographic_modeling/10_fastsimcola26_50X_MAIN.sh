#!/bin/bash
#


#### Script used to run fastsimcoal2 and produce results shown in Table S3
#### Likelihood analysis of three different phylogenies
#### 50 runs per model/population



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
FEXT[1]="PubSVsouth"
FEXT[2]="PubCentralAsia"
FEXT[3]="PubSpain"
FEXT[4]="PubMixed"

# Input folder
AA=/crex/proj/snic2020-6-184/private/Luis/P06_birch_exome_IGA-Sweden_2019/32b_fastsimcoal2_ancestral_state/BWA_unmasked

# Output folder
RR=/crex/proj/snic2020-6-184/private/Luis/P06_birch_exome_IGA-Sweden_2019/32c_fastsimcoal2_Bpubescens_DivergenceTime_simulations/BWA_unmasked

######## Evolutionary models
## all possible phylogenies (d: number of parameters)
EVOLM[1]=BP1	# ((PEN, PUB), PLA) and PEN->PUB migration only; d=8
EVOLM[2]=BP2	# ((PEN, PLA), PUB) and PEN->PUB migration only; d=8
EVOLM[3]=BP3	# (PEN, (PUB, PLA)) and PEN->PUB migration only; d=8



for n in `seq 1 1 3`; do 			# model

   # Template file (scenario)
   TPL=/crex1/proj/snic2017-7-149/private/Luis/P06_birch_exome_IGA-Sweden_2019/32_fastsimcoal2_Bpubescens_DivergenceTime/10_simulation_birch_3pop_divTime_model_${EVOLM[$n]}.tpl		

   # Estimation file
   EST=/crex1/proj/snic2017-7-149/private/Luis/P06_birch_exome_IGA-Sweden_2019/32_fastsimcoal2_Bpubescens_DivergenceTime/10_simulation_birch_3pop_divTime_model_${EVOLM[$n]}.est
   
   
   for k in `seq 1 1 4`; do 			# population	
   
   
      ################# Copy DAF, tpl and est files to output folder
      GO_AA=$AA/${RRsub}/$DPfilter/${sample_set}
      GO_RR=$RR/${RRsub}/$DPfilter/${sample_set}/model${EVOLM[$n]}/${FEXT[$k]}/C1
      GO_TPL=${TPL}
      GO_EST=${EST}
      mkdir -p $GO_RR
	  
      echo
      echo ${EVOLM[$n]} ${FEXT[$k]}
      echo $GO_RR
      echo $GO_TPL
      echo $GO_EST

      sbatch 10_fastsimcola26_50X_query.sh $GO_AA \
                                           $GO_RR \
                                           $GO_TPL \
                                           $GO_EST \
                                           ${EVOLM[$n]} \
                                           ${FEXT[$k]}


   done
done




