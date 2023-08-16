#!/bin/bash
#
#SBATCH -J fsc
#SBATCH -p core
#SBATCH -n 16
#SBATCH -t 24:00:00
#SBATCH -A snic2022-22-909
#SBATCH -M snowy
#SBATCH --mail-user luis.leal@ebc.uu.se
#SBATCH --mail-type=FAIL
ulimit -c unlimited


#### Script used to run fastsimcoal2 and produce results shown in Table S3
#### Likelihood analysis of three different phylogenies
#### 50 runs per model/population



##### load modules
module load bioinfo-tools


# local installation of fastSimcoal2 (ver 2.7.0.9)
fastSimcoal2=/crex1/proj/snic2017-7-149/private/Luis/z_APPS/fastsimcoal2/fsc27_linux64/fsc27093

# CITATION
# The following citations should be used for fsc26:
# Excoffier, L. and M. Foll. 2011. fastsimcoal: a continuous-time coalescent simulator of genomic diversity under arbitrarily complex evolutionary scenarios. Bioinformatics 27: 1332-1334.
# Excoffier, L., Dupanloup, I., Huerta-Sánchez, E., and M. Foll (2013) Robust demographic inference from genomic and SNP data. PLOS Genetics 9(10):e1003905

## Running fsc:
#  Simulate data under an evolutionary scenario with parameter values randomly drawn from priors:
#
#  $fsc26 -t test.tpl -n 10 -e test.est -E 100
#
#  fsc26 will use the scenario defined in the template file test.tpl and 
#  generate 100 sets of parameter values by randomly drawing these values from the priors defined in the file test.est. 
#  10 simulations will be done for each sets of randomly drawn parameter values.



############### input data files and paths

#remember initial path
SRCDIR_INI=$(pwd)  

# input folder
GO_AA=${1:?msg}

# output folder
RR=${2:?msg}

# template file
GO_TPL=${3:?msg}

# Estimation file
GO_EST=${4:?msg}

# Evolutionary models
EVOLM=${5:?msg}

# Population extension
FEXT=${6:?msg}




for i in `seq 1 1 50`; do 			# runs	

   GO_RR=${RR}/${i}
   mkdir -p ${GO_RR}

   rsync -ahv ${GO_AA}/3Pop_${FEXT}_jointDAFpop1_0.obs $GO_RR/model${EVOLM}_${FEXT}_jointDAFpop1_0.obs
   rsync -ahv ${GO_AA}/3Pop_${FEXT}_jointDAFpop2_0.obs $GO_RR/model${EVOLM}_${FEXT}_jointDAFpop2_0.obs
   rsync -ahv ${GO_AA}/3Pop_${FEXT}_jointDAFpop2_1.obs $GO_RR/model${EVOLM}_${FEXT}_jointDAFpop2_1.obs
   rsync -ahv $GO_TPL $GO_RR/model${EVOLM}_${FEXT[$k]}.tpl
   rsync -ahv $GO_EST $GO_RR/model${EVOLM}_${FEXT[$k]}.est

   ################# Run fastsimcoal2
   cd $GO_RR

   $fastSimcoal2 -t model${EVOLM}_${FEXT[$k]}.tpl \
                 -e model${EVOLM}_${FEXT[$k]}.est \
                 -M \
                 -n100000 \
                 -d \
                 -l10 \
                 -L40 \
                 -0 \
                 -c $SLURM_NTASKS \
                 -B $SLURM_NTASKS 

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

