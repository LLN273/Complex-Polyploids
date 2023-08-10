#!/bin/bash
#
#SBATCH -J str
#SBATCH -p core
#SBATCH -n 16
#SBATCH -t 100:00:00                          
#SBATCH -A snic2022-22-909
#SBATCH -M snowy
#SBATCH --mail-user luis.leal@ebc.uu.se
#SBATCH --mail-type=FAIL
ulimit -c unlimited


### Run STRUCTURE


##### load modules
module load bioinfo-tools
module load python3/3.7.2

# path to structure threader
structure_threader=/home/luisleal/.local/bin/structure_threader

# path to STRUCTURE 2.3.4.
STRUCTURE=/home/luisleal/.local/bin/structure

# Citation (structure_threader):
# Pina-Martins, F., Silva, D. N., Fino, J., & Paulo, O. S. (2017). Structure_threader: An improved method for automation and parallelization of programs structure, fastStructure and MavericK on multicore CPU systems. Molecular Ecology Resources, n/a-n/a. doi:10.1111/1755-0998.12702

# STRUCTURE:
# Pritchard JK, Stephens M, Donnelly P (2000) Inference of population structure using multilocus genotype data. Genetics, 155, 945–959.

# structure threader homepage:
# https://github.com/StuntsPT/Structure_threader
# https://structure-threader.readthedocs.io/en/latest/
# 


####################################### paths

#remember current path
SRCDIR=$(pwd)  

#input folder
AA=/crex1/proj/snic2020-6-184/private/Luis/P06_birch_exome_IGA-Sweden_2019/30_STRUCTURE/BWA_unmasked

# input subfolder
RRsub[1]=exome_WGS

# extra extension
EXT[1]=15_NANA_SILVER_WHITE_HUMILIS_PLATYPHYLLA_ReducedPEND_ReducedSW_10kSNP_SHORT

# Structure input file
STR[1]=MASTER_SingleRuns_ALLsamples_GOODonly_NANA_SILVER_WHITE_HUMILIS_PLATYPHYLLA-filtered_tranche_2_PASS-BIALLELIC_FINAL_NoPrivateAlleles_NONCODING

# input extension
OEXT[1]="ReducedPEND_ReducedSW_LD_SHORT_10kSNP_CLEAN.str"


### Number of replicates 
GO_replicates=16  



for k in `seq 1 1 1`; do 			# dataset (root)	

   GO_AA=${AA}
   GO_RR=${AA}/${EXT[$k]}
   GO_infile=${STR[$k]}_${OEXT[$k]}
   GO_mainparams=00_mainparams_${EXT[$k]}
   GO_extraparams=00_extraparams_${EXT[$k]}
   mkdir -p $GO_RR

   echo $GO_infile
   echo $GO_mainparams
   echo $GO_extraparams
   echo $GO_AA
   echo $GO_RR
   echo 

   cd $GO_RR

   ## Copy "mainparams" and "extraparams" files to output folder
   cd $GO_RR
   rsync -ahv $SRCDIR/$GO_mainparams ./mainparams
   rsync -ahv $SRCDIR/$GO_extraparams ./extraparams
   rsync -ahv $GO_AA/$GO_infile .

   ###### when using numeric labels for each sample 
   ## prepare structure input file (replace spaces by tabs; add empty line at the end of the file
   cat $GO_RR/$GO_infile | cut -f2- | sed 's/ /\t/g'  > $SNIC_TMP/_auxB.str
   
   NROWS="$(cat $SNIC_TMP/_auxB.str | wc -l)"
   aux_4="4"
   NSAMPLES="$(awk '{print $2/$1}' <<<"$aux_4 $NROWS")"
   echo "NSAMPLES:" $NSAMPLES
   
   seq 1 $NSAMPLES > $SNIC_TMP/_auxA.str
   cat $SNIC_TMP/_auxA.str $SNIC_TMP/_auxA.str $SNIC_TMP/_auxA.str $SNIC_TMP/_auxA.str | sort -n > $SNIC_TMP/_auxA2.str

   paste -d'\t' $SNIC_TMP/_auxA2.str $SNIC_TMP/_auxB.str > $GO_RR/strdataset		
   

   ## Run strucure in parallel

   echo
   cd $GO_RR
   $structure_threader run -Klist {2..9} \
                           -R $GO_replicates \
                           -i $GO_RR/strdataset \
                           -o $GO_RR/STRUCTURE_RESULTS \
                           -t $SLURM_NTASKS \
                           -st $STRUCTURE


   ################ usage:
   # ./structure_threader run -K Ks -R replicates -i infile -o outpath -t num_of_threads -st path_to_structure

   # Where:
   # -K 		number of "Ks" to run (from 1 to K) 
   # -Klist {2..8}	  # runs K values selected by user
   # -R 		number of replicate runs for each value of "K"
   # -i 		input file for STRUCTURE
   # -o 		directory where the output results should be stored
   # -t 		number of threads to use
   # -st 		path for the STRUCTURE binary.
   # --seed	define a random seed starting value
   # The program should be run in the same directory where the files "mainparams" and "extraparams" for your STRUCTURE run are placed.


done


