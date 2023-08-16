#!/bin/bash
#
#SBATCH -J fscGF_glean
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 01:00:00
#SBATCH -A snic2022-22-589
#SBATCH -M snowy
#SBATCH --mail-user luis.leal@ebc.uu.se
#SBATCH --mail-type=FAIL
ulimit -c unlimited


### Script used to glean fsc simulation results..." (output from 12_fastsimcola26_loopy_MAIN.sh)

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

# Input folder (root)
AA=/crex/proj/snic2020-6-184/private/Luis/P06_birch_exome_IGA-Sweden_2019/32c_fastsimcoal2_Bpubescens_DivergenceTime_simulations/BWA_unmasked



######## Evolutionary models
EVOLM=bmt1



for k in `seq 1 1 4`; do 			# dataset (root)

   #### outputfile
   outFile=model${EVOLM}_${FEXT[$k]}_loopy_OUT.txt
   
   # Output folder
   RR=$AA/${RRsub}/$DPfilter/${sample_set}/model${EVOLM}/${FEXT[$k]}/C1

   rm -f $RR/_aux1.txt
   rm -f $RR/_aux2.txt
   rm -f $RR/_aux3.txt
   rm -f $RR/$outFile
   
   
   ########## no gene flow
   
   gf_code="0_0_0_0_0_0_0_0"
   inFile=model${EVOLM}_${FEXT[$k]}.bestlhoods
   
   
   ##### Collect results for each gene flow permutation (50 runs per permutation)
   
   # Output folder
   GO_RR1=$AA/${RRsub}/$DPfilter/${sample_set}/model${EVOLM}/${FEXT[$k]}/C1/${gf_code}
   
   # Output file
   OUTfile1=fsm_geneflow_50X_ALL_RESULTS.txt
   
   # add header
   cat $AA/${RRsub}/$DPfilter/${sample_set}/model${EVOLM}/${FEXT[$k]}/C1/${gf_code}/1/model${EVOLM}_${FEXT[$k]}/$inFile | head -n1 > $GO_RR1/$OUTfile1
   
   for n in `seq 1 1 50`; do   
      GO_AA1=$AA/${RRsub}/$DPfilter/${sample_set}/model${EVOLM}/${FEXT[$k]}/C1/${gf_code}/${n}/model${EVOLM}_${FEXT[$k]}
      cat $GO_AA1/$inFile | tail -n1 >> $GO_RR1/$OUTfile1 
   done
   
   ##### Compute mean values
   
   OUTfile2="fsm_geneflow_50X_meanValues.txt"

   Rscript --no-save $SRCDIR_INI/14_compute_meanValue_loopy.R $GO_RR1/$OUTfile1 \
                                                              $OUTfile2 \
                                                              $GO_RR1
   

   echo $gf_code 
   echo $gf_code >> $RR/_aux1.txt
   cat $GO_RR1/$OUTfile2 | tail -n 1 >> $RR/_aux2.txt
   
   
   ########### in the presence of gene flow	
   for pen2pub0 in `seq 1 1 1`; do 			
      for pub2pen0 in `seq 0 1 1`; do 			
         for pen2pla0 in `seq 0 1 1`; do 
            for pla2pen0 in `seq 0 1 1`; do 	
               for pub2pla0 in `seq 0 1 1`; do 
                  for pla2pub0 in `seq 0 1 1`; do 	
                     for anc2pub1 in `seq 0 1 1`; do 
                        for pub2anc1 in `seq 0 1 1`; do 
                         
                           gf_code=${pen2pub0}_${pub2pen0}_${pen2pla0}_${pla2pen0}_${pub2pla0}_${pla2pub0}_${anc2pub1}_${pub2anc1}

                           ##### Collect results for each gene flow permutation (50 runs per permutation)
						   
                           # Output folder
                           GO_RR1=$AA/${RRsub}/$DPfilter/${sample_set}/model${EVOLM}/${FEXT[$k]}/C1/${gf_code}
						   
                           # Output file
                           OUTfile1=fsm_geneflow_50X_ALL_RESULTS.txt
						   
                           # add header
                           cat $AA/${RRsub}/$DPfilter/${sample_set}/model${EVOLM}/${FEXT[$k]}/C1/${gf_code}/1/model${EVOLM}_${FEXT[$k]}/$inFile | head -n1 > $GO_RR1/$OUTfile1
						   
                           for n in `seq 1 1 50`; do   
                              GO_AA1=$AA/${RRsub}/$DPfilter/${sample_set}/model${EVOLM}/${FEXT[$k]}/C1/${gf_code}/${n}/model${EVOLM}_${FEXT[$k]}
                              cat $GO_AA1/$inFile | tail -n1 >> $GO_RR1/$OUTfile1 
                           done
						   
                           ##### Compute mean values
						   
                           OUTfile2="fsm_geneflow_50X_meanValues.txt"

                           Rscript --no-save $SRCDIR_INI/14_compute_meanValue_loopy.R $GO_RR1/$OUTfile1 \
                                                                                      $OUTfile2 \
                                                                                      $GO_RR1
						   
   
                           echo $gf_code 
                           echo $gf_code >> $RR/_aux1.txt
                           cat $GO_RR1/$OUTfile2 | tail -n 1 >> $RR/_aux2.txt

                        done
                     done
                  done
               done
            done
         done
      done
   done

   paste $RR/_aux1.txt $RR/_aux2.txt > $RR/_aux3.txt
   sed 's/_/\t/g' $RR/_aux3.txt > $RR/$outFile

   rm -f $RR/_aux1.txt
   rm -f $RR/_aux2.txt
   rm -f $RR/_aux3.txt

done




