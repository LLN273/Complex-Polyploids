#!/bin/bash
#
#SBATCH -J fsc_BS
#SBATCH -p devel
#SBATCH -n 16
#SBATCH -t 01:00:00
#SBATCH -A snic2022-22-589
#SBATCH -M snowy
#SBATCH --mail-user luis.leal@ebc.uu.se
#SBATCH --mail-type=FAIL
ulimit -c unlimited


#### Script used to run fastsimcoal2 and produce results shown in Table S4
#### Estimate confidence intervals using parametric bootstrap SFS

echo
echo "Starting fsc simulations: estimate confidence intervals using parametric bootstrap SFS"
date -u


##### load modules
module load bioinfo-tools

# local installation of fastSimcoal2 (ver 2.7.0.9)
fastSimcoal2=/crex1/proj/snic2017-7-149/private/Luis/z_APPS/fastsimcoal2/fsc27_linux64/fsc27093

# CITATION
# The following citations should be used for fsc26:
# Excoffier, L. and M. Foll. 2011. fastsimcoal: a continuous-time coalescent simulator of genomic diversity under arbitrarily complex evolutionary scenarios. Bioinformatics 27: 1332-1334.
# Excoffier, L., Dupanloup, I., Huerta-Sánchez, E., and M. Foll (2013) Robust demographic inference from genomic and SNP data. PLOS Genetics 9(10):e1003905



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
AA=/crex/proj/snic2020-6-184/private/Luis/P06_birch_exome_IGA-Sweden_2019/32b_fastsimcoal2_ancestral_state/BWA_unmasked

# Input file (fastsimcoal2 results best model)
BB=/crex/proj/snic2020-6-184/private/Luis/P06_birch_exome_IGA-Sweden_2019/32c_fastsimcoal2_Bpubescens_DivergenceTime_simulations/BWA_unmasked

# Output folder
RR=/crex/proj/snic2020-6-184/private/Luis/P06_birch_exome_IGA-Sweden_2019/32c_fastsimcoal2_Bpubescens_DivergenceTime_simulations/BWA_unmasked

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


#### Best replicate (for each model)
RREP[1]=47
RREP[2]=4
RREP[3]=14
RREP[4]=7


### Estimating the confidence intervals
### See fsc manual, p. 58-60 (http://cmpg.unibe.ch/software/fastsimcoal2/man/fastsimcoal27.pdf)
### See also https://github.com/OB-lab/James_et_al._2021-MBE
### "We calculated the confidence intervals for the parameters of the [best] model (...) using parametric bootstrap. 
### This approach simulates DNA sequences, and their corresponding SFS, given the chosen model and the parameter values of its best run. Then, it recalculates the parameter values from the simulated SFS. This process was done 100 times.
### For simulating the DNA and SFS of the chosen model, the parameter values of its best run should be specified in a parameter (*.par) input file. This file can be generated editing the maximum likelihood parameter file (*_maxL.par) output by fastsimcoal2" 
### (and replacing the last three sections, as done below).


for k in `seq 1 1 4`; do 			# population	
           
   gf_code=${pen2pub0[$k]}_${pub2pen0[$k]}_${pen2pla0[$k]}_${pla2pen0[$k]}_${pub2pla0[$k]}_${pla2pub0[$k]}_${anc2pub1[$k]}_${pub2anc1[$k]}

   ################# Copy DAF + based tpl and est files to output folder
   GO_AA=$AA/${RRsub}/$DPfilter/${sample_set}
   GO_BB=$BB/${RRsub}/$DPfilter/${sample_set}/model${EVOLM}/${FEXT[$k]}/C1/$gf_code/${RREP[$k]}
   GO_RR=$RR/${RRsub}/$DPfilter/${sample_set}/model${EVOLM}_confidence_intervals/${FEXT[$k]}/C1/${gf_code}
   mkdir -p $GO_RR

   ################# Estimation file: original est file for BP2 model
   EST=${GO_BB}/model${EVOLM}_${FEXT[$k]}.est

   ################# Template file: original tpl file for BP2 model
   TPL=${GO_BB}/model${EVOLM}_${FEXT[$k]}.tpl

   ################# Parameter file (par) produced for best model	
   PAR=${GO_BB}/model${EVOLM}_${FEXT[$k]}/model${EVOLM}_${FEXT[$k]}_maxL.par

   ################# Parameter values file (pv) produced for best model
   PV=${GO_BB}/model${EVOLM}_${FEXT[$k]}/model${EVOLM}_${FEXT[$k]}.pv

   GO_EST=${EST}   
   GO_TPL=${TPL}
   GO_PAR=${PAR}
   GO_PV=${PV}

   rsync -ahv ${GO_AA}/3Pop_${FEXT[$k]}_jointDAFpop1_0.obs $GO_RR/model${EVOLM}_${FEXT[$k]}_jointDAFpop1_0.obs
   rsync -ahv ${GO_AA}/3Pop_${FEXT[$k]}_jointDAFpop2_0.obs $GO_RR/model${EVOLM}_${FEXT[$k]}_jointDAFpop2_0.obs
   rsync -ahv ${GO_AA}/3Pop_${FEXT[$k]}_jointDAFpop2_1.obs $GO_RR/model${EVOLM}_${FEXT[$k]}_jointDAFpop2_1.obs
   rsync -ahv $GO_EST $GO_RR/model${EVOLM}_${FEXT[$k]}.est
   rsync -ahv $GO_TPL $GO_RR/model${EVOLM}_${FEXT[$k]}.tpl
   rsync -ahv $GO_PAR $GO_RR/model${EVOLM}_${FEXT[$k]}_AUX.par

   ################# modify par file: specify that a DNA sequence, representing a given number of independent loci (10,000) of a particular length (100 bp), should be simulated. 

   ### Remove last 7 lines from par file
   head -n -7  $GO_RR/model${EVOLM}_${FEXT[$k]}_AUX.par > $GO_RR/_AUX1.txt

   ### New file ending, required to generate DNA sequence data (100 non-recombining segments of 100,000 bp)
   # //Number of independent loci [chromosome] 
   # 100 0
   # //Per chromosome: Number of linkage blocks
   # 1
   # //per Block: data type, num loci, rec. rate and mut rate + optional parameters
   # DNA 100000 0 9.5e-9 OUTEXP

   PAR_TAIL=${SRCDIR_INI}/16_PAR_ending.txt
      
   ### Create new par file
   cat $GO_RR/_AUX1.txt $PAR_TAIL > $GO_RR/model${EVOLM}_${FEXT[$k]}.par

   rm -f $GO_RR/_AUX1.txt


   ################# Make 100 SFS simulations under evolutionary scenario defined in parameter file
   echo
   cd $GO_RR

   $fastSimcoal2 -i $GO_RR/model${EVOLM}_${FEXT[$k]}.par \
                 -j \
                 -d \
                 -s0 \
                 -x \
                 -I \
                 -q \
                 -c $SLURM_NTASKS \
                 -B $SLURM_NTASKS \
                 -n100 \
                 -k 50000000


   ################# Use fastsimcoal to estimate model parameters based on each of the 100 pseudo-observed data sets (newly simulated jointDAF files)

   for n in `seq 1 1 100`; do 			

      echo
      echo "Pseudo data sets:" $n

      ## copy par, est, tpl and pv files to newly created BS folders
      rsync -ahv $GO_EST $GO_RR/model${EVOLM}_${FEXT[$k]}/model${EVOLM}_${FEXT[$k]}_${n}/model${EVOLM}_${FEXT[$k]}.est
      rsync -ahv $GO_TPL $GO_RR/model${EVOLM}_${FEXT[$k]}/model${EVOLM}_${FEXT[$k]}_${n}/model${EVOLM}_${FEXT[$k]}.tpl
      rsync -ahv $GO_PV $GO_RR/model${EVOLM}_${FEXT[$k]}/model${EVOLM}_${FEXT[$k]}_${n}/model${EVOLM}_${FEXT[$k]}.pv

      ## run fsc
      cd $GO_RR/model${EVOLM}_${FEXT[$k]}/model${EVOLM}_${FEXT[$k]}_${n}
	  
      rm -f model${EVOLM}_${FEXT[$k]}_DAF*
      rm -f model${EVOLM}_${FEXT[$k]}_numPolymSites.obs
	  
      cd $SRCDIR_INI
      sbatch ./16_fastsimcola26_confidence_intervals_query.sh $GO_RR/model${EVOLM}_${FEXT[$k]}/model${EVOLM}_${FEXT[$k]}_${n} \
                                                              ${EVOLM} \
                                                              ${FEXT[$k]}
	  

	  
   done
done



# Notes: 	
# -t		template file
# -e 		Estimation file (parameter prior definition file). Parameters drawn from specified distributions are substituted into template file.
# -i --ifile 	Name of parameter file
# -n		Number of simulations (minimum)
# -N		Maximum number of simulations to estimate the expected SFS
# -d 		computes derived site frequency spectrum (for SNP or DNA as SNP (-s) data only).
# -M		perform parameter estimation by max lhood from SFS values between iterations
# -l		minimum number of  loops (ECM cycles) for which the lhood is computed on both monomorphic and polymorphic sites if REFERENCE parameter is defined. Default is 2.
# -L 		maximum number of loops (ECM cycles) to perform during lhood maximization. Default is 20
# -k	    	Number of simulated polymorphic sites to keep in memory before writing them to temporary files. (default: 200,000)
# --multiSFS	generate or use multidimensional SFS
# -m		when using folded SFS 
# -E 		number of draws from parameter priors (Listed parameter values are substituted in template file)
# -c 		number of openMP threads for parameter estimation
# -B		max. no. of batches for multi-threaded runs (default=12)
# -0 		do not use information on the number of monomorphic sites
# -C  		minimum observed SFS entry count taken into account in likelihood computation (default = 1, but value can be < 1. e.g  0.5). Entries of the observed SFS with lower observations wil be collapsed into a single entry.
#			This means that entries with less than XXX SNPs are pooled together. This option is useful when there are many entries in the observed SFS with few SNPs and with a limited number of SNPS to avoid overfitting.
# -q		quiet mode
# -j  		output one simulated or bootstrapped SFS per file in a separate directory for easier analysis (requires -d or -m and -s0 options)
# -s  		Output DNA as SNP data, with a given maximum number to output (use 0 to output all SNPs in the DNA sequence(s)).
# -x  		Does not generate Arlequin output files
# -I 		Generates DNA mutations according to an infinite site (IS) mutation model. Under this model, each mutation is supposed to occur at a different but random site on the DNA sequence. Under
#		the IS model, if different mutations are allocated to the same DNA sequence position, they are generated independently, but marked to have occurred at the same position in the Arlequin output .arp file
# --initValues 	Specifies a file (*.pv) containing initial parameter values for parameter estimation. This is especially useful to reduce the number of runs necessary to estimate parameters when estimating confidence intervals by bootstrap


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




echo
echo "Done!"
date -u
echo



