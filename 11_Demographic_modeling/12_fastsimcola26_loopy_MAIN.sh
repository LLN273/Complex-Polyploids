#!/bin/bash
#


#### Script used to run fastsimcoal2 and produce results shown in Supplementary_Materials_5.xlsx
#### Likelihood analysis of all possible migration patterns for best model (BP2)
#### 50 runs per population/migration-pattern



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

# Input folder
AA=/crex/proj/snic2020-6-184/private/Luis/P06_birch_exome_IGA-Sweden_2019/32b_fastsimcoal2_ancestral_state/BWA_unmasked

# Output folder
RR=/crex/proj/snic2020-6-184/private/Luis/P06_birch_exome_IGA-Sweden_2019/32c_fastsimcoal2_Bpubescens_DivergenceTime_simulations/BWA_unmasked

######## Evolutionary models
## paper: test different gene flow combinations for BP2 model
EVOLM=bmt1

# Template file (scenario)
TPL=/crex1/proj/snic2017-7-149/private/Luis/P06_birch_exome_IGA-Sweden_2019/32_fastsimcoal2_Bpubescens_DivergenceTime/10_simulation_birch_3pop_divTime_model_${EVOLM}.tpl		

# Estimation file
EST=/crex1/proj/snic2017-7-149/private/Luis/P06_birch_exome_IGA-Sweden_2019/32_fastsimcoal2_Bpubescens_DivergenceTime/10_simulation_birch_3pop_divTime_model_${EVOLM}.est


for k in `seq 1 1 4`; do 			# population	
   for pen2pub0 in `seq 1 1 1`; do 			
      for pub2pen0 in `seq 0 1 1`; do 			
         for pen2pla0 in `seq 0 1 1`; do 
            for pla2pen0 in `seq 0 1 1`; do 	
               for pub2pla0 in `seq 0 1 1`; do 
                  for pla2pub0 in `seq 0 1 1`; do 	
                     for anc2pub1 in `seq 0 1 1`; do 
                        for pub2anc1 in `seq 0 1 1`; do 

                           gf_code=${pen2pub0}_${pub2pen0}_${pen2pla0}_${pla2pen0}_${pub2pla0}_${pla2pub0}_${anc2pub1}_${pub2anc1}

                           ################# Copy DAF + based tpl and est files to output folder
                           GO_AA=$AA/${RRsub}/$DPfilter/${sample_set}
                           GO_RR=$RR/${RRsub}/$DPfilter/${sample_set}/model${EVOLM}/${FEXT[$k]}/C1/$gf_code
   
                           GO_TPL=${TPL}
                           GO_EST=${EST}
                           mkdir -p $GO_RR

                           rsync -ahv ${GO_AA}/3Pop_${FEXT[$k]}_jointDAFpop1_0.obs $GO_RR/model${EVOLM}_${FEXT[$k]}_jointDAFpop1_0.obs
                           rsync -ahv ${GO_AA}/3Pop_${FEXT[$k]}_jointDAFpop2_0.obs $GO_RR/model${EVOLM}_${FEXT[$k]}_jointDAFpop2_0.obs
                           rsync -ahv ${GO_AA}/3Pop_${FEXT[$k]}_jointDAFpop2_1.obs $GO_RR/model${EVOLM}_${FEXT[$k]}_jointDAFpop2_1.obs
                           rsync -ahv $GO_TPL $GO_RR/model${EVOLM}_${FEXT[$k]}.tpl
                           rsync -ahv $GO_EST $GO_RR/model${EVOLM}_${FEXT[$k]}.est


                           ################# modify 'est' file according to gene flow conditions

                           if [[ ${pen2pub0} = 0 ]] ; then
                              sed 's/0 MIG0_pubpen logunif 1e-7 1e-4 output/0 MIG0_pubpen logunif 0 0 output/' $GO_RR/model${EVOLM}_${FEXT[$k]}.est > $GO_RR/_AUX.est
                              mv $GO_RR/_AUX.est $GO_RR/model${EVOLM}_${FEXT[$k]}.est
                           fi

                           if [[ ${pub2pen0} = 0 ]] ; then
                              sed 's/0 MIG0_penpub logunif 1e-7 1e-4 output/0 MIG0_penpub logunif 0 0 output/' $GO_RR/model${EVOLM}_${FEXT[$k]}.est > $GO_RR/_AUX.est
                              mv $GO_RR/_AUX.est $GO_RR/model${EVOLM}_${FEXT[$k]}.est
                           fi

                           if [[ ${pen2pla0} = 0 ]] ; then
                              sed 's/0 MIG0_plapen logunif 1e-7 1e-4 output/0 MIG0_plapen logunif 0 0 output/' $GO_RR/model${EVOLM}_${FEXT[$k]}.est > $GO_RR/_AUX.est
                              mv $GO_RR/_AUX.est $GO_RR/model${EVOLM}_${FEXT[$k]}.est
                           fi

                           if [[ ${pla2pen0} = 0 ]] ; then
                              sed 's/0 MIG0_penpla logunif 1e-7 1e-4 output/0 MIG0_penpla logunif 0 0 output/' $GO_RR/model${EVOLM}_${FEXT[$k]}.est > $GO_RR/_AUX.est
                              mv $GO_RR/_AUX.est $GO_RR/model${EVOLM}_${FEXT[$k]}.est
                           fi

                           if [[ ${pub2pla0} = 0 ]] ; then
                              sed 's/0 MIG0_plapub logunif 1e-7 1e-4 output/0 MIG0_plapub logunif 0 0 output/' $GO_RR/model${EVOLM}_${FEXT[$k]}.est > $GO_RR/_AUX.est
                              mv $GO_RR/_AUX.est $GO_RR/model${EVOLM}_${FEXT[$k]}.est
                           fi

                           if [[ ${pla2pub0} = 0 ]] ; then
                              sed 's/0 MIG0_pubpla logunif 1e-7 1e-4 output/0 MIG0_pubpla logunif 0 0 output/' $GO_RR/model${EVOLM}_${FEXT[$k]}.est > $GO_RR/_AUX.est
                              mv $GO_RR/_AUX.est $GO_RR/model${EVOLM}_${FEXT[$k]}.est
                           fi

                           if [[ ${anc2pub1} = 0 ]] ; then
                              sed 's/0 MIG1_pubpen logunif 1e-7 1e-4 output/0 MIG1_pubpen logunif 0 0 output/' $GO_RR/model${EVOLM}_${FEXT[$k]}.est > $GO_RR/_AUX.est
                              mv $GO_RR/_AUX.est $GO_RR/model${EVOLM}_${FEXT[$k]}.est
                           fi

                           if [[ ${pub2anc1} = 0 ]] ; then
                              sed 's/0 MIG1_penpub logunif 1e-7 1e-4 output/0 MIG1_penpub logunif 0 0 output/' $GO_RR/model${EVOLM}_${FEXT[$k]}.est > $GO_RR/_AUX.est
                              mv $GO_RR/_AUX.est $GO_RR/model${EVOLM}_${FEXT[$k]}.est
                           fi

                           ################# Run fastsimcoal2
                           sbatch $SRCDIR_INI/12_fastsimcola26_loopy_query.sh $GO_RR \
                                                                              $EVOLM \
                                                                              ${FEXT[$k]}


                        done
                     done
                  done
               done
            done
         done
      done
   done
done



# Notes: 	
# -t		template file
# -e 		Estimation file (parameter prior definition file). Parameters drawn from specified distributions are substituted into template file.
# -n		Number of simulations (minimum)
# -N		Maximum number of simulations to estimate the expected SFS
# -d 		computes derived site frequency spectrum (for SNP or DNA as SNP (-s) data only).
# -M		perform parameter estimation by max lhood from SFS values between iterations
# -l		minimum number of  loops (ECM cycles) for which the lhood is computed on both monomorphic and polymorphic sites if REFERENCE parameter is defined. Default is 2.
# -L 		maximum number of loops (ECM cycles) to perform during lhood maximization. Default is 20
# -k	    Number of simulated polymorphic sites to keep in memory before writing them to temporary files. (default: 200,000)
# --multiSFS	generate or use multidimensional SFS
# -m		when using folded SFS 
# -E 		number of draws from parameter priors (Listed parameter values are substituted in template file)
# -c 		number of openMP threads for parameter estimation
# -B		max. no. of batches for multi-threaded runs (default=12)
# -0 		do not use information on the number of monomorphic sites
# -C  		minimum observed SFS entry count taken into account in likelihood computation (default = 1, but value can be < 1. e.g  0.5). Entries of the observed SFS with lower observations wil be collapsed into a single entry.
#			This means that entries with less than XXX SNPs are pooled together. This option is useful when there are many entries in the observed SFS with few SNPs and with a limited number of SNPS to avoid overfitting.
# -q		quiet mode


##### From manual, p. 45: http://cmpg.unibe.ch/software/fastsimcoal2/man/fastsimcoal25.pdf
## OBSERVED SFS FILE NAMES
# Note that the name of the observed SFS file was not specified on the command line. This is because it is assumed to have the same name as the prefix of the template file (here 1PopBot20Mb) and a given suffix, 
# which exact definition depends on the number of population samples and on the type of SFS.
# Note also that all the observed SFS files should only contain a single observed SFS.

# ONE OBSERVED SAMPLE
# If there is a single observed sample in the model, the suffix will be:
# -_DAFpop0.obs if it is a file listing the derived allele SFS (unfolded spectrum)
# -_MAFpop0.obs if it is a file listing the minor allele SFS (folded spectrum)

#TWO OBSERVED SAMPLES
# If there are two observed samples in the model (0 and 1), one would need a file with the following suffix
# -_jointDAFpop1_0.obs if it is a file listing the derived allele SFS (unfolded spectrum)
# -_jointMAFpop1_0.obs if it is a file listing the minor allele SFS (folded spectrum)

#MORE THAN TWO OBSERVED SAMPLES
# If there are more than two observed samples in the model (say 0, 1, and 2), one would need three separate files with the following suffix
# -_jointDAFpop1_0.obs, _jointDAFpop2_1.obs, _jointDAFpop2_0.obs
# For the folded spectrum, the name would begin by _jointMAF

#MULTIDIMENSIONAL SFS
#It is also possible to tell fsc27 to use another format for observed SFS using the command line -multiSFS. In that case, fsc27 expects the observed SFS to be in a single file, even when more than one population sample is specified, with the following suffix:
#- _DSFS.obs

