#!/bin/bash
#
#SBATCH -J fsc_BS
#SBATCH -p core
#SBATCH -n 16
#SBATCH -t 24:00:00
#SBATCH -A snic2022-22-589
#SBATCH -M snowy
#SBATCH --mail-user luis.leal@ebc.uu.se
#SBATCH --mail-type=FAIL
ulimit -c unlimited


#### Script used to run fastsimcoal2 and produce results shown in Table S4
#### Estimate confidence intervals using parametric bootstrap SFS

echo
echo "Starting fsc simulations"
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

# output folder
RR=${1:?msg}

# Evolutionary models
EVOLM=${2:?msg}

# Population extension
FEXT=${3:?msg}


### Estimating the confidence intervals
### See fsc manual, p. 58-60 (http://cmpg.unibe.ch/software/fastsimcoal2/man/fastsimcoal27.pdf)
### See also https://github.com/OB-lab/James_et_al._2021-MBE
### "We calculated the confidence intervals for the parameters of the [best] model (...) using parametric bootstrap. 
### This approach simulates DNA sequences, and their corresponding SFS, given the chosen model and the parameter values of its best run. Then, it recalculates the parameter values from the simulated SFS. This process was done 100 times.
### For simulating the DNA and SFS of the chosen model, the parameter values of its best run should be specified in a parameter (*.par) input file. This file can be generated editing the maximum likelihood parameter file (*_maxL.par) output by fastsimcoal2" 
### (and replacing the last three sections, as done below).


for k in `seq 1 1 50`; do 			# runs	
                      
   echo 
   echo $k

   GO_RR=${RR}/${k}
   mkdir -p ${GO_RR}
   cd ${GO_RR}

   rsync -ahv ${RR}/*.pv .
   rsync -ahv ${RR}/*.tpl .
   rsync -ahv ${RR}/*.est .
   
   # copy parametric bootstrap SFS to folder
   rsync -ahv ${RR}/*.obs .


   ## run fsc
   cd $GO_RR
	  
   rm -f model${EVOLM}_${FEXT}_DAF*
   rm -f model${EVOLM}_${FEXT}_numPolymSites.obs
	  
   $fastSimcoal2 -t model${EVOLM}_${FEXT}.tpl \
                 -e model${EVOLM}_${FEXT}.est \
                 -n100000 \
                 -d \
                 -M \
                 -L30 \
                 -0 \
                 –initvalues model${EVOLM}_${FEXT}.pv \
                 -c16 \
                 -B16 \
                 -q
	  

	
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



