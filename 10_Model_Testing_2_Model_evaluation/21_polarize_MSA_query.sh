#!/bin/bash
#
ulimit -c unlimited
set -eo pipefail


## Polarize polyploid sequence

# load modules
#module load bioinfo-tools
#module load python3/3.8.7
#module load R/4.0.0


#### paths and folder names

# remember initial path
SRCDIR_INI=$(pwd)      

# input folder
AA=${1:?msg}

# output folder 
RR=${2:?msg}

# Number of species
NSPEC=${3:?msg}

# subfolder name (it is also the root name for each library)
GO_FOLDER=${4:?msg}

## replicate number
r_aux=${5:?msg}

## number of loci
NLOCUS=${6:?msg}

## PARENT 1
anc_1=${7:?msg}

## PARENT 2
anc_2=${8:?msg}

## PARENT 3
anc_3=${9:?msg}

## PARENT 4
anc_4=${10:?msg}

## Reference sequence ID1
REFspecies_1=${11:?msg}

## Iteration 
NALT=${12:?msg}

## Gene flow to pendula, species
GF_2pend_species=${13:?msg}

## Gene flow to pendula, fraction of pendula genes replaced by platyphylla
HYB_FRACTION_2pen=${14:?msg}

## HYBR0, species
HYBR_0[1]=${15:?msg}
HYBR_0[2]=${16:?msg}
HYBR_0[3]=${17:?msg}
HYBR_0[4]=${18:?msg}

## HYBR0, fraction of genes replaced
H_0=${19:?msg}

## Hybridization Blocks 1 to 7, species
HYBR_1[1]=${20:?msg}
HYBR_1[2]=${21:?msg}
HYBR_1[3]=${22:?msg}
HYBR_1[4]=${23:?msg}

HYBR_2[1]=${24:?msg}
HYBR_2[2]=${25:?msg}
HYBR_2[3]=${26:?msg}
HYBR_2[4]=${27:?msg}

HYBR_3[1]=${28:?msg}
HYBR_3[2]=${29:?msg}
HYBR_3[3]=${30:?msg}
HYBR_3[4]=${31:?msg}

HYBR_4[1]=${32:?msg}
HYBR_4[2]=${33:?msg}
HYBR_4[3]=${34:?msg}
HYBR_4[4]=${35:?msg}

HYBR_5[1]=${36:?msg}
HYBR_5[2]=${37:?msg}
HYBR_5[3]=${38:?msg}
HYBR_5[4]=${39:?msg}

HYBR_6[1]=${40:?msg}
HYBR_6[2]=${41:?msg}
HYBR_6[3]=${42:?msg}
HYBR_6[4]=${43:?msg}

HYBR_7[1]=${44:?msg}
HYBR_7[2]=${45:?msg}
HYBR_7[3]=${46:?msg}
HYBR_7[4]=${47:?msg}

HYBR_8[1]=${48:?msg}
HYBR_8[2]=${49:?msg}
HYBR_8[3]=${50:?msg}
HYBR_8[4]=${51:?msg}

# Hybridization fraction, Blocks 1 to 7
HYBRY_FRACTION_1=${52:?msg}
HYBRY_FRACTION_2=${53:?msg}
HYBRY_FRACTION_3=${54:?msg}
HYBRY_FRACTION_4=${55:?msg}
HYBRY_FRACTION_5=${56:?msg}
HYBRY_FRACTION_6=${57:?msg}
HYBRY_FRACTION_7=${58:?msg}
HYBRY_FRACTION_8=${59:?msg}

### Extra introgression species
Extra_introgr_species=${60:?msg}

# Hybridization levels for B. pendula
GO_HYBRY_EXTRA=${61:?msg}

###### Simulation conditions (polyploidization model)
myMODEL=${62:?msg}

## ABC subfolder name
#ABCsn=${55:?msg}

#echo
#echo $GO_FOLDER
#echo $AA
#echo $RR
#echo "Number of species:" $NSPEC
#echo "Number of loci:" $NLOCUS
#echo "Iteration:" $NALT
#echo 'Reference sequence ID 1:' $REFspecies_1
#echo 'Parental species IDs:' $anc_1 $anc_2 $anc_3 $anc_4
#echo 'Gene flow to pendula, species:' $GF_2pend_species
#echo 'Gene flow to pendula, fraction of pendula genes replaced [0-1]:' $HYB_FRACTION_2pen
#echo 'Block 0:' "${HYBR_0[@]}" $H_0
#echo 'Block 1:' "${HYBR_1[@]}" ${HYBRY_FRACTION_1}
#echo 'Block 2:' "${HYBR_2[@]}" ${HYBRY_FRACTION_2}
#echo 'Block 3:' "${HYBR_3[@]}" ${HYBRY_FRACTION_3}
#echo 'Block 4:' "${HYBR_4[@]}" ${HYBRY_FRACTION_4}
#echo 'Block 5:' "${HYBR_5[@]}" ${HYBRY_FRACTION_5}
#echo 'Block 6:' "${HYBR_6[@]}" ${HYBRY_FRACTION_6}
#echo 'Block 7:' "${HYBR_7[@]}" ${HYBRY_FRACTION_7}
#echo 'Extra introgression species:' $Extra_introgr_species
#echo 'Hybridization levels for extra species:' $GO_HYBRY_EXTRA
#echo


####################### hybridization, number of events
GF_P2P="$(awk '{print $1*$2}' <<<"$NLOCUS $HYB_FRACTION_2pen")"
HGT_TOTAL_0="$(awk '{print $1*$2}' <<<"$NLOCUS $H_0")"
HGT_TOTAL_1="$(awk '{print $1*$2}' <<<"$NLOCUS $HYBRY_FRACTION_1")"
HGT_TOTAL_2="$(awk '{print $1*$2}' <<<"$NLOCUS $HYBRY_FRACTION_2")"
HGT_TOTAL_3="$(awk '{print $1*$2}' <<<"$NLOCUS $HYBRY_FRACTION_3")"
HGT_TOTAL_4="$(awk '{print $1*$2}' <<<"$NLOCUS $HYBRY_FRACTION_4")"
HGT_TOTAL_5="$(awk '{print $1*$2}' <<<"$NLOCUS $HYBRY_FRACTION_5")"
HGT_TOTAL_6="$(awk '{print $1*$2}' <<<"$NLOCUS $HYBRY_FRACTION_6")"
HGT_TOTAL_7="$(awk '{print $1*$2}' <<<"$NLOCUS $HYBRY_FRACTION_7")"
HGT_TOTAL_8="$(awk '{print $1*$2}' <<<"$NLOCUS $HYBRY_FRACTION_8")"

HGT0_stop_round="$(printf "%.0f\n" $HGT_TOTAL_0 )"	# round to closest integer
HGT1_stop_round="$(printf "%.0f\n" $HGT_TOTAL_1 )"	# round to closest integer
HGT2_stop_round="$(printf "%.0f\n" $HGT_TOTAL_2 )"	# round to closest integer
HGT3_stop_round="$(printf "%.0f\n" $HGT_TOTAL_3 )"	# round to closest integer
HGT4_stop_round="$(printf "%.0f\n" $HGT_TOTAL_4 )"	# round to closest integer
HGT5_stop_round="$(printf "%.0f\n" $HGT_TOTAL_5 )"	# round to closest integer
HGT6_stop_round="$(printf "%.0f\n" $HGT_TOTAL_6 )"	# round to closest integer
HGT7_stop_round="$(printf "%.0f\n" $HGT_TOTAL_7 )"	# round to closest integer
HGT8_stop_round="$(printf "%.0f\n" $HGT_TOTAL_8 )"	# round to closest integer


##### Extra introgression (B. pendula)
HGT_TOTAL_extra_aux="$(awk '{print $1*$2}' <<<"$NLOCUS $GO_HYBRY_EXTRA")"
HGT_TOTAL_extra_aux_stop_round="$(printf "%.0f\n" $HGT_TOTAL_extra_aux )"	# round to closest integer


if [ $HGT_TOTAL_extra_aux_stop_round -eq "0" ] ; then
   HGT_TOTAL_extra=0 
   
elif [ ${myMODEL} = "PPPP" ] || [ ${myMODEL} = "AAAA" ] ; then
   aux_HE="$(awk '{print $1+$2+$3+$4}' <<<"$HGT0_stop_round $HGT1_stop_round $HGT2_stop_round $HGT3_stop_round")"
   if [ $aux_HE -gt $HGT_TOTAL_extra_aux_stop_round ] ; then
	  HGT_TOTAL_extra=$HGT_TOTAL_extra_aux_stop_round
   else
      HGT_TOTAL_extra="$(awk '{print $1+$2+$3}' <<<"$HGT_TOTAL_extra_aux_stop_round $HGT5_stop_round $HGT6_stop_round")"  
   fi
	   
elif [ ${myMODEL} = "PPNN" ] || [ ${myMODEL} = "AANN" ] || [ ${myMODEL} = "PPHH" ] || [ ${myMODEL} = "AAHH" ] || [ ${myMODEL} = "PPNN_HE" ] || [ ${myMODEL} = "AANN_HE" ] || [ ${myMODEL} = "PPHH_HE" ] || [ ${myMODEL} = "AAHH_HE" ] ; then
   aux_HE="$(awk '{print $1+$2+$3}' <<<"$HGT1_stop_round $HGT2_stop_round $HGT3_stop_round")"
   if [ $aux_HE -gt $HGT_TOTAL_extra_aux_stop_round ] ; then
	  HGT_TOTAL_extra="$(awk '{print $1+$2}' <<<"$HGT_TOTAL_extra_aux_stop_round $HGT0_stop_round")"
   else
      HGT_TOTAL_extra="$(awk '{print $1+$2+$3+$4+$5+$6}' <<<"$HGT_TOTAL_extra_aux_stop_round $HGT0_stop_round $HGT4_stop_round $HGT5_stop_round $HGT6_stop_round $HGT7_stop_round")"
   fi
   
elif [ ${myMODEL} = "PPPyPy" ] || [ ${myMODEL} = "PPPyPy_HE" ] ; then
   aux_HE="$(awk '{print $1+$2+$3+$4}' <<<"$HGT0_stop_round $HGT1_stop_round $HGT2_stop_round $HGT3_stop_round")"
   aux_HE2="$(awk '{print $1+$2+$3+$4+$5}' <<<"$HGT0_stop_round $HGT1_stop_round $HGT2_stop_round $HGT3_stop_round $HGT7_stop_round")"
   if [ $aux_HE -gt $HGT_TOTAL_extra_aux_stop_round ] ; then
      HGT_TOTAL_extra=$HGT_TOTAL_extra_aux_stop_round
   elif [ $aux_HE2 -gt $HGT_TOTAL_extra_aux_stop_round ] ; then
      HGT_TOTAL_extra="$(awk '{print $1+$2+$3+$4}' <<<"$HGT_TOTAL_extra_aux_stop_round $HGT4_stop_round $HGT5_stop_round $HGT6_stop_round")"
   else
      HGT_TOTAL_extra="$(awk '{print $1+$2+$3+$4+$5}' <<<"$HGT_TOTAL_extra_aux_stop_round $HGT4_stop_round $HGT5_stop_round $HGT6_stop_round $HGT8_stop_round")"
   fi
   
elif [ ${myMODEL} = "PPNH" ] || [ ${myMODEL} = "AANH" ] || [ ${myMODEL} = "PPNH_HE" ] || [ ${myMODEL} = "AANH_HE" ] ; then
   aux_HE="$(awk '{print $1+$2+$3}' <<<"$HGT1_stop_round $HGT2_stop_round $HGT3_stop_round")"
   if [ $aux_HE -gt $HGT_TOTAL_extra_aux_stop_round ] ; then
      HGT_TOTAL_extra="$(awk '{print $1+$2}' <<<"$HGT_TOTAL_extra_aux_stop_round $HGT0_stop_round")"
   else   
      HGT_TOTAL_extra="$(awk '{print $1+$2+$3+$4+$5+$6+$7}' <<<"$HGT_TOTAL_extra_aux_stop_round $HGT0_stop_round $HGT4_stop_round $HGT5_stop_round $HGT6_stop_round $HGT7_stop_round $HGT8_stop_round")"
   fi
fi


if [ $HGT_TOTAL_extra -gt $NLOCUS ] ; then 
	HGT_TOTAL_extra=$NLOCUS
fi





mkdir -p $SNIC_TMP/ALT${NALT}/$r_aux
mkdir -p $RR/$r_aux
cd $RR/$r_aux
rsync -ah $SRCDIR_INI/21s_PLATY2PEND_SIMPLE_INTERSPECIFIC_HYBRIDIZATION.py .
rsync -ah $SRCDIR_INI/21s_polarizeMSA_fullCODED_ALT_TETRA_HYBRIDIZATION_BOTH_SinglePolarRef_2HYBRIDS_maxHETEROZY.py .
rsync -ah $SRCDIR_INI/21s_polarizeMSA_fullCODED_ALT_TETRA_HYBRIDIZATION_BOTH_SinglePolarRef_2HYBRIDS_maxHETEROZY_step_2.py .



## Randomize gene list
cat $SRCDIR_INI/00_GENE_LIST.txt | head -${NLOCUS} | shuf > $SNIC_TMP/ALT${NALT}/_GENE_LIST_${r_aux}.txt


## read randomized gene list

unset GENE_LIST
i=1

while read -r LINE                                                                      
do
   GENE_LIST[$i]=$LINE
   let "i+=1"
done < $SNIC_TMP/ALT${NALT}/_GENE_LIST_${r_aux}.txt

N_GENES=${#GENE_LIST[@]}


######################################################### Platyphylla to pendula gene flow

#echo
#echo "Gene flow from platyphylla to pendula"
#echo

#cd $SNIC_TMP/ALT${NALT}
#mkdir -p $SNIC_TMP/ALT${NALT}/$r_aux
cd $SNIC_TMP/ALT${NALT}/$r_aux
rsync -ah $SRCDIR_INI/21s_PLATY2PEND_SIMPLE_INTERSPECIFIC_HYBRIDIZATION.py .

GF_stop_round="$(printf "%.0f\n" $GF_P2P )"	# round to closest integer

if [[ $GF_stop_round -gt "0" ]] ; then
   for n in `seq 1 1 $GF_stop_round`; do

      l_aux=${GENE_LIST[$n]}
      GF_FLAG="1"   #equivalent to HYBRID_CODE="10"
      #echo $n $l_aux "PP_GF" $GF_FLAG

      # get simulated sequences (before polarization); remove spaces and blank lines
      cat $AA/$r_aux/data_${l_aux}_TRUE.fasta | tr -d "[:blank:]" | grep -v '^$' > $SNIC_TMP/ALT${NALT}/$r_aux/data_${l_aux}_TRUE_CLEAN.fasta

      # implement hybridization
      python3 21s_PLATY2PEND_SIMPLE_INTERSPECIFIC_HYBRIDIZATION.py data_${l_aux}_TRUE_CLEAN.fasta $NSPEC $anc_1 $anc_2 $anc_3 $anc_4 $GF_2pend_species $GF_FLAG
      
      mv data_${l_aux}_TRUE_CLEAN-ALT.fasta $RR/$r_aux/data_${l_aux}_TRUE_CLEAN.fasta

   done
else
   GF_stop_round=0
fi



## no gene flow from platyphylla to pendula
aux_1=1
NO_GF_start="$(awk '{print $1+$2}' <<<"$GF_stop_round $aux_1")"	
   
if [[ $GF_stop_round -lt $NLOCUS ]] ; then
   for n in `seq $NO_GF_start 1 $NLOCUS`; do
   
      l_aux=${GENE_LIST[$n]}
      GF_FLAG="0"   #equivalent to HYBRID_CODE="00"
      #echo $n $l_aux "0" $GF_FLAG

      # get simulated sequences (before polarization); remove spaces and blank lines
      cat $AA/$r_aux/data_${l_aux}_TRUE.fasta | tr -d "[:blank:]" | grep -v '^$' > $SNIC_TMP/ALT${NALT}/$r_aux/data_${l_aux}_TRUE_CLEAN.fasta

      # polarize sequences
      python3 21s_PLATY2PEND_SIMPLE_INTERSPECIFIC_HYBRIDIZATION.py data_${l_aux}_TRUE_CLEAN.fasta $NSPEC $anc_1 $anc_2 $anc_3 $anc_4 $GF_2pend_species $GF_FLAG

      mv data_${l_aux}_TRUE_CLEAN-ALT.fasta $RR/$r_aux/data_${l_aux}_TRUE_CLEAN.fasta

   done
fi



####### Re-randomize gene list

cat $SRCDIR_INI/00_GENE_LIST.txt | head -${NLOCUS} | shuf > $SNIC_TMP/ALT${NALT}/_GENE_LIST_AUX_${r_aux}.txt

# Extend gene list (three exact copies, in sequence) >> required to deal with cases such as Hhum=0.7 & Hnana=0.4 if sites do not overlap (in fact, in this case they will overlap in 10% of cases)
cat $SNIC_TMP/ALT${NALT}/_GENE_LIST_AUX_${r_aux}.txt $SNIC_TMP/ALT${NALT}/_GENE_LIST_AUX_${r_aux}.txt $SNIC_TMP/ALT${NALT}/_GENE_LIST_AUX_${r_aux}.txt $SNIC_TMP/ALT${NALT}/_GENE_LIST_AUX_${r_aux}.txt > $SNIC_TMP/ALT${NALT}/_GENE_LIST_${r_aux}.txt

## read randomized gene list

unset GENE_LIST
i=1

while read -r LINE                                                                      
do
   GENE_LIST[$i]=$LINE
   let "i+=1"
done < $SNIC_TMP/ALT${NALT}/_GENE_LIST_${r_aux}.txt

N_GENES=${#GENE_LIST[@]}



##################################################################################
################################################################################## Hybridization + Polarization
##################################################################################

let "NSPEC+=2"		# 2x2 sequences added (equivalent to twice _0_0 _0_1)
cd $RR/$r_aux
#echo

   
######################################################### Special hybridization event affecting only the first haplotype in polyploid (used to simulate recent hybridization with B. pendula)

HGTextra_stop_round="$(printf "%.0f\n" $HGT_TOTAL_extra )"	# round to closest integer

Hextra_start=1
Hextra_stop=$HGTextra_stop_round

if [[ $Hextra_stop -gt 0 ]] ; then
   for n in `seq $Hextra_start 1 $Hextra_stop`; do

      l_aux=${GENE_LIST[$n]}

      # hybrid code
      HYBRID_CODE="1000"		# hybridize only first haplotype
      #echo $n $l_aux "Hextra" $HYBRID_CODE

      # polarize sequences
      python3 21s_polarizeMSA_fullCODED_ALT_TETRA_HYBRIDIZATION_BOTH_SinglePolarRef_2HYBRIDS_maxHETEROZY.py data_${l_aux}_TRUE_CLEAN.fasta $NSPEC $Extra_introgr_species $Extra_introgr_species $Extra_introgr_species $Extra_introgr_species ${HYBRID_CODE}
      
      mv data_${l_aux}_TRUE_CLEAN-ALT.fasta $RR/$r_aux/data_${l_aux}_TRUE_CLEAN.fasta
	  
   done
fi



######################################################### Hybridization 0

#HGT0_stop_round="$(printf "%.0f\n" $HGT_TOTAL_0 )"	# round to closest integer

#aux_1=1
#HO_start="$(awk '{print $1+$2}' <<<"$GF_stop_round $aux_1")"
#HO_stop="$(awk '{print $1+$2}' <<<"$GF_stop_round $HGT0_stop_round")"
HO_start=1
HO_stop=$HGT0_stop_round

#echo "HO" $HGT0_stop_round $HO_stop
#echo

if [[ $HGT0_stop_round -gt 0 ]] ; then
   for n in `seq $HO_start 1 $HO_stop`; do

      l_aux=${GENE_LIST[$n]}

      # hybrid code
      HYBRID_CODE="0111" 	# note: first haplotype has been reserved for introgression from B. pendula
 
      #echo $n $l_aux "H0" $HYBRID_CODE

      # polarize sequences
      python3 21s_polarizeMSA_fullCODED_ALT_TETRA_HYBRIDIZATION_BOTH_SinglePolarRef_2HYBRIDS_maxHETEROZY_step_2.py data_${l_aux}_TRUE_CLEAN.fasta $NSPEC $REFspecies_1 "${HYBR_0[@]}" ${HYBRID_CODE} "PUBESCENS_N" "PUBESCENS_PEN" "PUBESCENS_H" "PUBESCENS_PLA" "PUBESCENS_ANC"
      
   done

else

   #aux_1=1
   #HO_stop="$(awk '{print $1-$2}' <<<"$HO_start $aux_1")"
   HO_stop=0

fi



######################################################### Hybridization 1

#HGT1_stop_round="$(printf "%.0f\n" $HGT_TOTAL_1 )"	# round to closest integer
aux_1=1
H1_start="$(awk '{print $1+$2}' <<<"$HO_stop $aux_1")"
H1_stop="$(awk '{print $1+$2}' <<<"$HO_stop $HGT1_stop_round")"


if [[ $HGT1_stop_round -gt 0 ]] ; then
   for n in `seq $H1_start 1 $H1_stop`; do

      l_aux=${GENE_LIST[$n]}

      # hybrid code
      HYBRID_CODE="0111" 

      #echo $n $l_aux "H1" $HYBRID_CODE

      # polarize sequences
      python3 21s_polarizeMSA_fullCODED_ALT_TETRA_HYBRIDIZATION_BOTH_SinglePolarRef_2HYBRIDS_maxHETEROZY_step_2.py data_${l_aux}_TRUE_CLEAN.fasta $NSPEC $REFspecies_1 "${HYBR_1[@]}" ${HYBRID_CODE} "PUBESCENS_N" "PUBESCENS_PEN" "PUBESCENS_H" "PUBESCENS_PLA" "PUBESCENS_ANC"
      
   done

else

   aux_1=1
   H1_stop="$(awk '{print $1-$2}' <<<"$H1_start $aux_1")"

fi



######################################################### Hybridization 2

#HGT2_stop_round="$(printf "%.0f\n" $HGT_TOTAL_2 )"	# round to closest integer
aux_1=1
H2_start="$(awk '{print $1+$2}' <<<"$H1_stop $aux_1")"
H2_stop="$(awk '{print $1+$2}' <<<"$H1_stop $HGT2_stop_round")"

if [[ $HGT2_stop_round -gt 0 ]] ; then
   for n in `seq $H2_start 1 $H2_stop`; do

   l_aux=${GENE_LIST[$n]}

   # hybrid code
   HYBRID_CODE="0111" 

   #echo $n $l_aux "H2" $HYBRID_CODE

   # polarize sequences
   python3 21s_polarizeMSA_fullCODED_ALT_TETRA_HYBRIDIZATION_BOTH_SinglePolarRef_2HYBRIDS_maxHETEROZY_step_2.py data_${l_aux}_TRUE_CLEAN.fasta $NSPEC $REFspecies_1 "${HYBR_2[@]}" ${HYBRID_CODE} "PUBESCENS_N" "PUBESCENS_PEN" "PUBESCENS_H" "PUBESCENS_PLA" "PUBESCENS_ANC"
      
   done

else

   aux_1=1
   H2_stop="$(awk '{print $1-$2}' <<<"$H2_start $aux_1")"

fi



######################################################### Hybridization 3

#HGT3_stop_round="$(printf "%.0f\n" $HGT_TOTAL_3 )"	# round to closest integer
aux_1=1
H3_start="$(awk '{print $1+$2}' <<<"$H2_stop $aux_1")"
H3_stop="$(awk '{print $1+$2}' <<<"$H2_stop $HGT3_stop_round")"

if [[ $HGT3_stop_round -gt 0 ]] ; then
   for n in `seq $H3_start 1 $H3_stop`; do

      l_aux=${GENE_LIST[$n]}

      # hybrid code
      HYBRID_CODE="0111" 

      #echo $n $l_aux "H3" $HYBRID_CODE

      # polarize sequences
      python3 21s_polarizeMSA_fullCODED_ALT_TETRA_HYBRIDIZATION_BOTH_SinglePolarRef_2HYBRIDS_maxHETEROZY_step_2.py data_${l_aux}_TRUE_CLEAN.fasta $NSPEC $REFspecies_1 "${HYBR_3[@]}" ${HYBRID_CODE} "PUBESCENS_N" "PUBESCENS_PEN" "PUBESCENS_H" "PUBESCENS_PLA" "PUBESCENS_ANC"
      
   done

else

   aux_1=1
   H3_stop="$(awk '{print $1-$2}' <<<"$H3_start $aux_1")"

fi



######################################################### Hybridization 4

#HGT4_stop_round="$(printf "%.0f\n" $HGT_TOTAL_4 )"	# round to closest integer
aux_1=1
H4_start="$(awk '{print $1+$2}' <<<"$H3_stop $aux_1")"
H4_stop="$(awk '{print $1+$2}' <<<"$H3_stop $HGT4_stop_round")"

if [[ $HGT4_stop_round -gt 0 ]] ; then
   for n in `seq $H4_start 1 $H4_stop`; do

      l_aux=${GENE_LIST[$n]}

      # hybrid code
      HYBRID_CODE="1111" 

      #echo $n $l_aux "H4" $HYBRID_CODE

      # polarize sequences
      python3 21s_polarizeMSA_fullCODED_ALT_TETRA_HYBRIDIZATION_BOTH_SinglePolarRef_2HYBRIDS_maxHETEROZY_step_2.py data_${l_aux}_TRUE_CLEAN.fasta $NSPEC $REFspecies_1 "${HYBR_4[@]}" ${HYBRID_CODE} "PUBESCENS_N" "PUBESCENS_PEN" "PUBESCENS_H" "PUBESCENS_PLA" "PUBESCENS_ANC"
      
   done

else

   aux_1=1
   H4_stop="$(awk '{print $1-$2}' <<<"$H4_start $aux_1")"

fi




######################################################### Hybridization 5

#HGT5_stop_round="$(printf "%.0f\n" $HGT_TOTAL_5 )"	# round to closest integer
aux_1=1
H5_start="$(awk '{print $1+$2}' <<<"$H4_stop $aux_1")"
H5_stop="$(awk '{print $1+$2}' <<<"$H4_stop $HGT5_stop_round")"

if [[ $HGT5_stop_round -gt 0 ]] ; then
   for n in `seq $H5_start 1 $H5_stop`; do

      l_aux=${GENE_LIST[$n]}

      # hybrid code
      HYBRID_CODE="1111" 

      #echo $n $l_aux "H5" $HYBRID_CODE

      # polarize sequences
      python3 21s_polarizeMSA_fullCODED_ALT_TETRA_HYBRIDIZATION_BOTH_SinglePolarRef_2HYBRIDS_maxHETEROZY_step_2.py data_${l_aux}_TRUE_CLEAN.fasta $NSPEC $REFspecies_1 "${HYBR_5[@]}" ${HYBRID_CODE} "PUBESCENS_N" "PUBESCENS_PEN" "PUBESCENS_H" "PUBESCENS_PLA" "PUBESCENS_ANC"
      
   done

else

   aux_1=1
   H5_stop="$(awk '{print $1-$2}' <<<"$H5_start $aux_1")"

fi




######################################################### Hybridization 6

#HGT6_stop_round="$(printf "%.0f\n" $HGT_TOTAL_6 )"	# round to closest integer
aux_1=1
H6_start="$(awk '{print $1+$2}' <<<"$H5_stop $aux_1")"
H6_stop="$(awk '{print $1+$2}' <<<"$H5_stop $HGT6_stop_round")"

if [[ $HGT6_stop_round -gt 0 ]] ; then
   for n in `seq $H6_start 1 $H6_stop`; do

      l_aux=${GENE_LIST[$n]}

      # hybrid code
      HYBRID_CODE="1111" 

      #echo $n $l_aux "H6" $HYBRID_CODE


      # polarize sequences
      python3 21s_polarizeMSA_fullCODED_ALT_TETRA_HYBRIDIZATION_BOTH_SinglePolarRef_2HYBRIDS_maxHETEROZY_step_2.py data_${l_aux}_TRUE_CLEAN.fasta $NSPEC $REFspecies_1 "${HYBR_6[@]}" ${HYBRID_CODE} "PUBESCENS_N" "PUBESCENS_PEN" "PUBESCENS_H" "PUBESCENS_PLA" "PUBESCENS_ANC"
      
   done

else

   aux_1=1
   H6_stop="$(awk '{print $1-$2}' <<<"$H6_start $aux_1")"

fi




######################################################### Hybridization 7

#HGT7_stop_round="$(printf "%.0f\n" $HGT_TOTAL_7 )"	# round to closest integer
aux_1=1
H7_start="$(awk '{print $1+$2}' <<<"$H6_stop $aux_1")"
H7_stop="$(awk '{print $1+$2}' <<<"$H6_stop $HGT7_stop_round")"

if [[ $HGT7_stop_round -gt 0 ]] ; then
   for n in `seq $H7_start 1 $H7_stop`; do

      l_aux=${GENE_LIST[$n]}

      # hybrid code
      HYBRID_CODE="1111" 

      #echo $n $l_aux "H7" $HYBRID_CODE

      # polarize sequences
      python3 21s_polarizeMSA_fullCODED_ALT_TETRA_HYBRIDIZATION_BOTH_SinglePolarRef_2HYBRIDS_maxHETEROZY_step_2.py data_${l_aux}_TRUE_CLEAN.fasta $NSPEC $REFspecies_1 "${HYBR_7[@]}" ${HYBRID_CODE} "PUBESCENS_N" "PUBESCENS_PEN" "PUBESCENS_H" "PUBESCENS_PLA" "PUBESCENS_ANC"
      
   done

else

   aux_1=1
   H7_stop="$(awk '{print $1-$2}' <<<"$H7_start $aux_1")"

fi



######################################################### Hybridization 8
#HGT8_stop_round="$(printf "%.0f\n" $HGT_TOTAL_8 )"	# round to closest integer

aux_1=1
H8_start="$(awk '{print $1+$2}' <<<"$H7_stop $aux_1")"
H8_stop="$(awk '{print $1+$2}' <<<"$H7_stop $HGT8_stop_round")"

if [[ $HGT8_stop_round -gt 0 ]] ; then
   for n in `seq $H8_start 1 $H8_stop`; do

      l_aux=${GENE_LIST[$n]}

      # hybrid code
      HYBRID_CODE="1111" 

      #echo $n $l_aux "H8" $HYBRID_CODE

      # polarize sequences
      python3 21s_polarizeMSA_fullCODED_ALT_TETRA_HYBRIDIZATION_BOTH_SinglePolarRef_2HYBRIDS_maxHETEROZY_step_2.py data_${l_aux}_TRUE_CLEAN.fasta $NSPEC $REFspecies_1 "${HYBR_8[@]}" ${HYBRID_CODE} "PUBESCENS_N" "PUBESCENS_PEN" "PUBESCENS_H" "PUBESCENS_PLA" "PUBESCENS_ANC"
      
   done

else

   aux_1=1
   H8_stop="$(awk '{print $1-$2}' <<<"$H8_start $aux_1")"

fi






######################################################### loci with no HGT (no hybridization)
aux_1=1
NO_HYBRI_start="$(awk '{print $1+$2}' <<<"$H8_stop $aux_1")"	
   
if [[ $H8_stop -lt $NLOCUS ]] ; then
   for n in `seq $NO_HYBRI_start 1 $NLOCUS`; do

      l_aux=${GENE_LIST[$n]}
      HYBRID_CODE="0000"
      #echo $n $l_aux "0" $HYBRID_CODE

      # polarize sequences
      python3 21s_polarizeMSA_fullCODED_ALT_TETRA_HYBRIDIZATION_BOTH_SinglePolarRef_2HYBRIDS_maxHETEROZY_step_2.py data_${l_aux}_TRUE_CLEAN.fasta $NSPEC $REFspecies_1 $anc_1 $anc_2 $anc_3 $anc_4 ${HYBRID_CODE} "PUBESCENS_N" "PUBESCENS_PEN" "PUBESCENS_H" "PUBESCENS_PLA" "PUBESCENS_ANC"

   done
fi





###

let "NSPEC-=2"

# delete auxiliary
#rm -f $SNIC_TMP/ALT${NALT}/_*
#rm -f $SNIC_TMP/ALT${NALT}/data*
#rm -f $SNIC_TMP/ALT${NALT}/*.py
rm -rf $SNIC_TMP/ALT${NALT}/$r_aux
rm -f $RR/$r_aux/data_*_TRUE_CLEAN.fasta
rm -f $RR/$r_aux/data_*_TRUE_CLEAN.fasta


#echo
#echo 'Done!'

