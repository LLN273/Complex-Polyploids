#!/bin/bash
#


## Glean gene tree counts from outputs produced by 24_IQTREE_gene_tree_DISTANCE_ALT.R


#### paths and folder names

#remember initial path
SRCDIR_INI=$(pwd)                                           	 

# input folder 
GO_AA=${1:?msg}

# input file
In_file=${2:?msg}




#echo
GO_RR=${GO_AA}

## output file
outFile=${In_file%.txt}_summary.txt
rm -f $GO_RR/$outFile
echo "pendula" "platyphylla" "pend_platy" "populifolia" "pop_out" "nana" "nana_out" "occidentalis" "occ_out" "humilis" "albosinensis" "lenta" "albo_len_hum" "albo_len" "basal" > $GO_RR/$outFile


# get counts (%)
cd $GO_AA
TC_PEND="$(grep -w 'PENDULA_SIM' ${In_file} | cut -f3)"
TC_PLATY="$(grep -w 'PLATYPHYLLA_SIM' ${In_file} | cut -f3)"
TC_PEND_PLATY="$(grep -w 'PENDULA_SIM_PLATYPHYLLA_SIM' ${In_file} | cut -f3)"
TC_POP="$(grep -w 'POPULIFOLIA_SIM' ${In_file} | cut -f3)"
TC_POP_OUT="$(grep -w 'PENDULA_SIM_PLATYPHYLLA_SIM_POPULIFOLIA_SIM' ${In_file} | cut -f3)"
TC_NANA="$(grep -w 'NANA_SIM' ${In_file} | cut -f3)"
TC_NANA_OUT="$(grep -w 'NANA_SIM_PENDULA_SIM_PLATYPHYLLA_SIM_POPULIFOLIA_SIM' ${In_file} | cut -f3)"
TC_OCC="$(grep -w 'OCCIDENTALIS_SIM' ${In_file} | cut -f3)"
TC_OCC_OUT="$(grep -w 'NANA_SIM_OCCIDENTALIS_SIM_PENDULA_SIM_PLATYPHYLLA_SIM_POPULIFOLIA_SIM' ${In_file} | cut -f3)"
TC_HUM="$(grep -w 'HUMILIS_SIM' ${In_file} | cut -f3)"
TC_ALBOS="$(grep -w 'ALBOSINENSIS_SIM' ${In_file} | cut -f3)"
TC_LEN="$(grep -w 'LENTA_SIM' ${In_file} | cut -f3)"
TC_HUMALBLEN="$(grep -w 'ALBOSINENSIS_SIM_HUMILIS_SIM_LENTA_SIM' ${In_file} | cut -f3)"
TC_ALBLEN="$(grep -w 'ALBOSINENSIS_SIM_LENTA_SIM' ${In_file} | cut -f3)"
TC_BASAL="$(grep -w 'ALBOSINENSIS_SIM_ALNUS_SIM_HUMILIS_SIM_LENTA_SIM_NANA_SIM_OCCIDENTALIS_SIM_PENDULA_SIM_PLATYPHYLLA_SIM_POPULIFOLIA_SIM' ${In_file} | cut -f3)"

if [[ $TC_PEND = "" ]]; then TC_PEND=0; fi
if [[ $TC_PLATY = "" ]]; then TC_PLATY=0; fi
if [[ $TC_PEND_PLATY = "" ]]; then TC_PEND_PLATY=0; fi
if [[ $TC_POP = "" ]]; then TC_POP=0; fi
if [[ $TC_POP_OUT = "" ]]; then TC_POP_OUT=0; fi
if [[ $TC_NANA = "" ]]; then TC_NANA=0; fi
if [[ $TC_NANA_OUT = "" ]]; then TC_NANA_OUT=0; fi
if [[ $TC_OCC = "" ]]; then TC_OCC=0; fi
if [[ $TC_OCC_OUT = "" ]]; then TC_OCC_OUT=0; fi
if [[ $TC_HUM = "" ]]; then TC_HUM=0; fi
if [[ $TC_ALBOS = "" ]]; then TC_ALBOS=0; fi
if [[ $TC_LEN = "" ]]; then TC_LEN=0; fi
if [[ $TC_HUMALBLEN = "" ]]; then TC_HUMALBLEN=0; fi
if [[ $TC_ALBLEN = "" ]]; then TC_ALBLEN=0; fi
if [[ $TC_BASAL = "" ]]; then TC_BASAL=0; fi

echo $TC_PEND $TC_PLATY $TC_PEND_PLATY $TC_POP $TC_POP_OUT $TC_NANA $TC_NANA_OUT $TC_OCC $TC_OCC_OUT $TC_HUM $TC_ALBOS $TC_LEN $TC_HUMALBLEN $TC_ALBLEN $TC_BASAL >> $GO_RR/$outFile

