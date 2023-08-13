#!/bin/bash
#
#SBATCH -J INDELIble
#SBATCH -p core 
#SBATCH -n 4
#SBATCH -t 06:00:00
#SBATCH -A snic2022-22-909
#SBATCH -M snowy
#SBATCH --mail-user luis.leal@ebc.uu.se
#SBATCH --mail-type=FAIL
ulimit -c unlimited
set -eo pipefail


## Script used to create fasta sequences (simulated MSA) with INDELIble, based on simulated gene phylogenies.
## For wrapper's usage, see simphy manual: https://github.com/adamallo/SimPhy/wiki/Manual#521-input-files-newick-tree-format


# load modules
module load bioinfo-tools
module load perl/5.26.2


# path to INDELIble
INDELIBLE="/crex1/proj/snic2017-7-149/private/Luis/z_APPS/INDELible/INDELibleV1.03/src/"
export PATH="$INDELIBLE:$PATH"

# simphy's INDELIble wrapper
IDW=/crex1/proj/snic2017-7-149/private/Luis/z_APPS/SimPhy/SimPhy_1.0.2/scripts/INDELIble_wrapper.pl


#### paths and folder names

# remember initial path
SRCDIR_INI=$(pwd)      

# Input folder
#AA=/crex1/proj/snic2020-6-184/nobackup/private/Luis/P11_SIMULATED_FASTA_PHYLOGENY_HYBRIDIZATION_exomeData/01_simphy_ILS/A_modILS_birch_CLEAN_reducedILS_30perct	# reduced ILS (30%); 1500bp read length
#AA=/crex1/proj/snic2020-6-184/nobackup/private/Luis/P11_SIMULATED_FASTA_PHYLOGENY_HYBRIDIZATION_exomeData/01_simphy_ILS/A_modILS_birch_CLEAN_reducedILS_35perct	# reduced ILS (35%); 1500bp read length
AA=/crex1/proj/snic2020-6-184/nobackup/private/Luis/P11_SIMULATED_FASTA_PHYLOGENY_HYBRIDIZATION_exomeData/01_simphy_ILS/A_modILS_birch_CLEAN_reducedILS_40perct	# reduced ILS (40%); 1500bp read length
#AA=/crex1/proj/snic2020-6-184/nobackup/private/Luis/P11_SIMULATED_FASTA_PHYLOGENY_HYBRIDIZATION_exomeData/01_simphy_ILS/A_modILS_birch_CLEAN_reducedILS_45perct	# reduced ILS (45%); 1500bp read length
#AA=/crex1/proj/snic2020-6-184/nobackup/private/Luis/P11_SIMULATED_FASTA_PHYLOGENY_HYBRIDIZATION_exomeData/01_simphy_ILS/A_modILS_birch_CLEAN_reducedILS_50perct	# reduced ILS (50%); 1500bp read length

# Configuration file
CONFIG=${SRCDIR_INI}/INDELible_complex_1500.txt											# read length reduced to (mean=1500,st=100); used in P11 sims
		

# Run perl wrapper
# Run as: INDELIble_wrapper.pl directory input_config seed number-of-cores

perl $IDW \
     $AA \
     $CONFIG \
     22 \
     $SLURM_NTASKS
     

