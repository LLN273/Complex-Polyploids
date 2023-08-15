#!/bin/bash
#
#SBATCH -J preEST
#SBATCH -p core
#SBATCH -n 4
#SBATCH -t 12:00:00
#SBATCH -A snic2022-22-909
#SBATCH -M snowy
#SBATCH --mail-user luis.leal@ebc.uu.se
#SBATCH --mail-type=FAIL
ulimit -c unlimited


# Script used to prepare data required to run EST-SFS (step 2)
# Collects VCF for focal species (B. pendula) and the three outgroups (B. nana, B. populifolia, and B. occidentalis)
# Prepares EST-SFS input file


##### load modules
module load bioinfo-tools
module load python3



####################################### paths

#remember current path
SRCDIR=$(pwd)  

### Sample set
sample_set="4_8"

### Path to VCF files associated to focal species
AA=/crex1/proj/snic2020-6-184/private/Luis/P06_birch_exome_IGA-Sweden_2019/32a_fastsimcoal2_Bpubescens_DivergenceTime_VCFfiles/BWA_unmasked

### input/output subfolder 1
RRsub=exome_WGS

#### Depth filtering
DPfilter="8X"

### Input VCf (focal species; root name)
VCF=MASTER_SingleRuns_ALLsamples_GOODonly_SILVER_WHITE_PLATYPHYLLA

##### input/output extension
FEXT[1]="PubSpain"
FEXT[2]="PubCentralEurope"
FEXT[3]="PubSVsouth"
FEXT[4]="PubCentralAsia"

### Folder containing VCF files associated to outgroup species
BB=/crex1/proj/snic2020-6-184/private/Luis/P06_birch_exome_IGA-Sweden_2019/32b_fastsimcoal2_ancestral_state/BWA_unmasked

## VCF files associated to outgroup species (species 1 is closer to focal species)
VCF_OUTG[1]=Populifolia-Canada.A005-11_snps-CLEAN_READY.vcf					# B. populifolia
VCF_OUTG[2]=Nana-Finland-Enontekio.A002-2_snps-CLEAN_READY.vcf				# B. nana
VCF_OUTG[3]=Occidentalis-Canada-Alberta.A009-17_snps-CLEAN_READY.vcf		# B. occidentalis

### output folder
RR=/crex1/proj/snic2020-6-184/private/Luis/P06_birch_exome_IGA-Sweden_2019/32b_fastsimcoal2_ancestral_state/BWA_unmasked

### Focal species: sample list
if [[ ${sample_set} = "4_8" ]] ; then
   SPL=/crex1/proj/snic2017-7-149/private/Luis/P06_birch_exome_IGA-Sweden_2019/32_fastsimcoal2_Bpubescens_DivergenceTime/00_samples_ancestral_state_PENDULA_8.txt
fi



for k in `seq 1 1 4`; do 			# dataset (root)	

   GO_AA=${AA}/${RRsub}/$DPfilter/${sample_set}
   GO_BB=$BB/${RRsub}
   GO_RR=${RR}/${RRsub}/$DPfilter/${sample_set}
   GO_SL=$SPL
   GO_FILE_focal=${VCF}_${FEXT[$k]}_intronic_and_4FD_DPfilter.vcf
   GO_out_file=EST_input_data_${FEXT[$k]}.txt
   mkdir -p $GO_RR

   echo
   echo $GO_AA
   echo $GO_RR
   echo $GO_FILE_focal
   echo $GO_SL
   echo 


   ##################### unzip vcf files
   #zcat $GO_AA/${GO_FILE_focal}.gz > $GO_RR/${GO_FILE_focal}
   #zcat $GO_BB/${VCF_OUTG[1]}.gz > $GO_BB/${VCF_OUTG[1]}
   #zcat $GO_BB/${VCF_OUTG[2]}.gz > $GO_BB/${VCF_OUTG[2]}
   #zcat $GO_BB/${VCF_OUTG[3]}.gz > $GO_BB/${VCF_OUTG[3]}


   ##################### Run python script
   cd $GO_RR
   python3 $SRCDIR/07s_EST_SFS_Create_InputFile.py $GO_RR/${GO_FILE_focal} \
                                                   $SPL \
                                                   $GO_BB/${VCF_OUTG[1]} \
                                                   $GO_BB/${VCF_OUTG[2]} \
                                                   $GO_BB/${VCF_OUTG[3]} \
                                                   $GO_out_file


   ### Clean output file
   sed 's/[][]//g' $GO_RR/$GO_out_file > $SNIC_TMP/_AUX_1.txt
   sed 's/ //g' $SNIC_TMP/_AUX_1.txt > $GO_RR/${GO_out_file%.txt}_CLEAN.txt
   
   
   rm -f $SNIC_TMP/_*
								   
done



## clean auxilary files
rm -f $SNIC_TMP/*


echo
echo 'Done!'





