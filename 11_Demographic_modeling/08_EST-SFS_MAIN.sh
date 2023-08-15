#!/bin/bash
#
#SBATCH -J EST_SFS
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 12:00:00
#SBATCH -A snic2022-22-909
#SBATCH -M snowy
#SBATCH --mail-user luis.leal@ebc.uu.se
#SBATCH --mail-type=FAIL
ulimit -c unlimited


# Script used to run EST-SFS
# Focal species: B. pendula
# OG1: B. populifolia
# OG2: B. nana
# OG3: B. occidentalis



##### load modules
module load bioinfo-tools

# Path to local install of EST_SFS, v. 2.04
# Citation:
# Keightley, P. D and Jackson B. C. (2018). Inferring the probability of the derived versus the ancestral allelic state at a polymorphic site. Genetics 209: 897-906.
EST_SFS=/crex1/proj/snic2017-7-149/private/Luis/z_APPS/est-sfs-release-2.04/est-sfs



####################################### paths

#remember current path
SRCDIR=$(pwd)  

### Sample set
sample_set="4_8"

### Path to input data files
AA=/crex1/proj/snic2020-6-184/private/Luis/P06_birch_exome_IGA-Sweden_2019/32b_fastsimcoal2_ancestral_state/BWA_unmasked

### input/output subfolder 1
RRsub=exome_WGS

#### Depth filtering
DPfilter="8X"

### Input data file (root name)
DataFileRoot=EST_input_data

##### input/output extension
FEXT[1]="PubSpain"
FEXT[2]="PubCentralEurope"
FEXT[3]="PubSVsouth"
FEXT[4]="PubCentralAsia"


#### EST-SFS configuration file
configFile=/crex1/proj/snic2017-7-149/private/Luis/P06_birch_exome_IGA-Sweden_2019/32_fastsimcoal2_Bpubescens_DivergenceTime/08_config-rate6.txt

## Note: 
# n_outgroup [1, 2, or 3]
# model [0, 1 or 2]
# nrandom [0 or positive integer]

# model:
# 0 = Jukes-Cantor model
# 1 = Kimura 2-parameter model
# 2 = Rate-6 model (see Keightley and Jackson 2018 for details)

# nrandom:
# If nrandom = 0, the ML search algorithm starts with pre-set parameter starting values. If nrandom >= 1, the algorithm returns the highest log likelihood found in nrandom runs using starting values for the parameters
# randomly sampled from wide ranges. The MLs of the individual runs are reported. In the case of the Rate-6 model, it is recommended that at least 10 random starting value runs are carried out to ensure convergence.

#### EST-SFS seed file
#### This text file contains a single positive integer. It is overwritten by a new random number seed at the end of the run.
seedfile=/crex1/proj/snic2017-7-149/private/Luis/P06_birch_exome_IGA-Sweden_2019/32_fastsimcoal2_Bpubescens_DivergenceTime/08_seedfile.txt

#### EST-SFS output file 1: uSFS (root name)
output_file_sfs=EST-SFS_output-file-sfs

#### EST-SFS output file 2: p-values (root name)
output_file_pvalues=EST-SFS_output-file-pvalues

# Notes about file with output p-values:
# The file begins with several lines started by “0” containing various self-explanatory outputs from the program. There then follow numbered lines, containing fields as follows:
# 1. Line number. This corresponds to the line number of the data file.
# 2. Configuration index. Not of interest to the user.
# 3. The probability of the major allele being ancestral.
# 4-7. For two outgroups: the probabilities of the first internal node (Fig. 1) having states A, C, G or T.
# 4-19. For three outgroups: the probabilities of the first and second internal nodes (Fig. 1) having states [A,A], [A,C], [A,G], [A,T], [C,A], [C,C], [C,G], [C,T], etc.


for k in `seq 1 1 4`; do 			# dataset (root)	

   GO_AA=${AA}/${RRsub}/$DPfilter/${sample_set}
   GO_RR=${GO_AA}

   GO_DATA_FILE=${DataFileRoot}_${FEXT[$k]}_CLEAN.txt
   GO_output_file_sfs=${output_file_sfs}_${FEXT[$k]}.txt
   GO_output_file_pvalues=${output_file_pvalues}_${FEXT[$k]}.txt

   echo
   echo $GO_AA
   echo $GO_RR
   echo $GO_DATA_FILE
   echo $GO_output_file_sfs
   echo $GO_output_file_pvalues
   echo 

   ##################### Run EST-SFS
   cd $GO_RR
   $EST_SFS $configFile \
            $GO_AA/$GO_DATA_FILE \
            $seedfile \
            $GO_RR/$GO_output_file_sfs \
            $GO_RR/$GO_output_file_pvalues
				   
done



echo
echo 'Done!'




