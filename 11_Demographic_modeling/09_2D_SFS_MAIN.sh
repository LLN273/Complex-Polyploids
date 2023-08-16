#!/bin/bash
#
#SBATCH -J 2D_SFS
#SBATCH -p devcore
#SBATCH -n 1
#SBATCH -t 00:10:00
#SBATCH -A snic2022-22-909
#SBATCH -M snowy
#SBATCH --mail-user luis.leal@ebc.uu.se
#SBATCH --mail-type=FAIL
ulimit -c unlimited



# Script used to create 2D-SFS
# Populations:
# Pop0: B. pubescens
# Pop1: B. pendula
# Pop2: B. platyphylla




##### load modules
module load bioinfo-tools
module load python3


####################################### paths


#remember current path
SRCDIR=$(pwd)  

### Sample set
sample_set="4_8"


#### Path to input data files (VCF files)
AA=/crex1/proj/snic2020-6-184/private/Luis/P06_birch_exome_IGA-Sweden_2019/32a_fastsimcoal2_Bpubescens_DivergenceTime_VCFfiles/BWA_unmasked

#### input/output subfolder 1
RRsub=exome_WGS

#### Depth filtering
DPfilter="8X"

##### input/output extension
FEXT[1]="PubSpain"
FEXT[2]="PubCentralEurope"
FEXT[3]="PubSVsouth"
FEXT[4]="PubCentralAsia"

#### Input VCF file, root name (joint VCF)
VCF=MASTER_SingleRuns_ALLsamples_GOODonly_SILVER_WHITE_PLATYPHYLLA

#### Ancestral allele info (EST_SFS pvalues output)
BB=/crex1/proj/snic2020-6-184/private/Luis/P06_birch_exome_IGA-Sweden_2019/32b_fastsimcoal2_ancestral_state/BWA_unmasked
	
#### EST-SFS file name (root name)
#### Assumed to be output for data consisting of 1 focal species and 3 outgroups
EST_SFS_pvalues=EST-SFS_output-file-pvalues

### Population ploidy
ploidy_pop0=4  	# B. pubescens
ploidy_pop1=2	# B. pendula
ploidy_pop2=2	# B. platyphylla

### Population code (SFS output file_)
code_pop0=0  	# B. pubescens
code_pop1=1		# B. pendula
code_pop2=2		# B. platyphylla


if [[ ${sample_set} = "4_8" ]] ; then
	samples_pub[1]=/crex1/proj/snic2017-7-149/private/Luis/P06_birch_exome_IGA-Sweden_2019/32_fastsimcoal2_Bpubescens_DivergenceTime/00_samples_ancestral_state_PUBESCENS_PubSpain_4.txt
	samples_pub[2]=/crex1/proj/snic2017-7-149/private/Luis/P06_birch_exome_IGA-Sweden_2019/32_fastsimcoal2_Bpubescens_DivergenceTime/00_samples_ancestral_state_PUBESCENS_PubCentralEurope_4.txt
	samples_pub[3]=/crex1/proj/snic2017-7-149/private/Luis/P06_birch_exome_IGA-Sweden_2019/32_fastsimcoal2_Bpubescens_DivergenceTime/00_samples_ancestral_state_PUBESCENS_PubSVsouth_4.txt
	samples_pub[4]=/crex1/proj/snic2017-7-149/private/Luis/P06_birch_exome_IGA-Sweden_2019/32_fastsimcoal2_Bpubescens_DivergenceTime/00_samples_ancestral_state_PUBESCENS_PubCentralAsia_4.txt
	samples_pend=/crex1/proj/snic2017-7-149/private/Luis/P06_birch_exome_IGA-Sweden_2019/32_fastsimcoal2_Bpubescens_DivergenceTime/00_samples_ancestral_state_PENDULA_8.txt
	samples_platy=/crex1/proj/snic2017-7-149/private/Luis/P06_birch_exome_IGA-Sweden_2019/32_fastsimcoal2_Bpubescens_DivergenceTime/00_samples_ancestral_state_PLATYPHYLA_8.txt
fi



for k in `seq 1 1 4`; do 			# dataset (root)	

   GO_AA=${AA}/${RRsub}/$DPfilter/${sample_set}
   GO_BB=${BB}/${RRsub}/$DPfilter/${sample_set}
   GO_RR=${GO_BB}

   GO_VCF_FILE=${VCF}_${FEXT[$k]}_intronic_and_4FD_DPfilter.vcf
   GO_EST_SFS_pvalues=${EST_SFS_pvalues}_${FEXT[$k]}.txt
   GO_SAMPLES_0=${samples_pub[$k]}
   GO_SAMPLES_1=$samples_pend
   GO_SAMPLES_2=$samples_platy
   
   echo
   echo $GO_AA/$GO_VCF_FILE
   echo $GO_BB/$GO_EST_SFS_pvalues
   echo $GO_RR
   echo $GO_SAMPLES_0
   echo $GO_SAMPLES_1
   echo $GO_SAMPLES_2
   echo "Ploidies:" $ploidy_pop0 $ploidy_pop1 $ploidy_pop2
   echo "Population codes:" $code_pop0 $code_pop1 $code_pop2
   echo 

   #unzip input VCF file
   zcat $GO_AA/${GO_VCF_FILE}.gz > $SNIC_TMP/${GO_VCF_FILE}



   ##################### Create SFS
   cd $GO_RR
   python3 $SRCDIR/09s_2D_SFS_Creator.py $SNIC_TMP/${GO_VCF_FILE} \
                                         $GO_BB/$GO_EST_SFS_pvalues \
                                         $GO_SAMPLES_0 \
                                         $GO_SAMPLES_1 \
                                         $ploidy_pop0 \
                                         $ploidy_pop1 \
                                         $code_pop0 \
                                         $code_pop1 \
                                         $GO_RR/3Pop_${FEXT[$k]}_jointDAFpop1_0.obs
   
   python3 $SRCDIR/09s_2D_SFS_Creator.py $SNIC_TMP/${GO_VCF_FILE} \
                                         $GO_BB/$GO_EST_SFS_pvalues \
                                         $GO_SAMPLES_0 \
                                         $GO_SAMPLES_2 \
                                         $ploidy_pop0 \
                                         $ploidy_pop2 \
                                         $code_pop0 \
                                         $code_pop2 \
                                         $GO_RR/3Pop_${FEXT[$k]}_jointDAFpop2_0.obs
   
   python3 $SRCDIR/09s_2D_SFS_Creator.py $SNIC_TMP/${GO_VCF_FILE} \
                                         $GO_BB/$GO_EST_SFS_pvalues \
                                         $GO_SAMPLES_1 \
                                         $GO_SAMPLES_2 \
                                         $ploidy_pop1 \
                                         $ploidy_pop2 \
                                         $code_pop1 \
                                         $code_pop2 \
                                         $GO_RR/3Pop_${FEXT[$k]}_jointDAFpop2_1.obs
				 

   ## Add header to DAF files
   echo "1 observations" > $SNIC_TMP/_auxFinal.obs
   cat $SNIC_TMP/_auxFinal.obs $GO_RR/3Pop_${FEXT[$k]}_jointDAFpop1_0.obs > $SNIC_TMP/_auxFinal2.obs
   rsync -ahv $SNIC_TMP/_auxFinal2.obs $GO_RR/3Pop_${FEXT[$k]}_jointDAFpop1_0.obs
   
   cat $SNIC_TMP/_auxFinal.obs $GO_RR/3Pop_${FEXT[$k]}_jointDAFpop2_0.obs > $SNIC_TMP/_auxFinal2.obs
   rsync -ahv $SNIC_TMP/_auxFinal2.obs $GO_RR/3Pop_${FEXT[$k]}_jointDAFpop2_0.obs
   
   cat $SNIC_TMP/_auxFinal.obs $GO_RR/3Pop_${FEXT[$k]}_jointDAFpop2_1.obs > $SNIC_TMP/_auxFinal2.obs
   rsync -ahv $SNIC_TMP/_auxFinal2.obs $GO_RR/3Pop_${FEXT[$k]}_jointDAFpop2_1.obs
						   
done


## clean auxiliary files
rm -f $SNIC_TMP/*

echo
echo 'Done!'






