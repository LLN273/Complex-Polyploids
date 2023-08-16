#!/bin/bash
#
#SBATCH -J infog
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 24:00:00
#SBATCH -A snic2022-22-589
#SBATCH -M snowy
#SBATCH --mail-user luis.leal@ebc.uu.se
#SBATCH --mail-type=FAIL
ulimit -c unlimited

## Get list of parsimony-informative sites for each gene
## Run separately for each sample
 
echo
echo "Starting Uppmax jobs ..."
echo


# load modules
module load bioinfo-tools


#### paths and folder names

#remember initial path
SRCDIR_INI=$(pwd)                                           	 


# input folder (root)
AA[1]=/crex1/proj/snic2020-6-184/private/Luis/P06_birch_exome_IGA-Sweden_2019/43b_polarizationProtocol_IQ-TREE2_IUPAC_ALT1_SLIM/BWA_unmasked
AA[2]=/crex1/proj/snic2020-6-184/private/Luis/P06_birch_exome_IGA-Sweden_2019/44b_polarizationProtocol_IQ-TREE2_IUPAC_ALT2_SLIM/BWA_unmasked
AA[3]=/crex1/proj/snic2020-6-184/private/Luis/P06_birch_exome_IGA-Sweden_2019/46b_polarizationProtocol_IQ-TREE2_IUPAC_ALT3_SLIM/BWA_unmasked
AA[4]=/crex1/proj/snic2020-6-184/private/Luis/P06_birch_exome_IGA-Sweden_2019/45b_polarizationProtocol_IQ-TREE2_IUPAC_ALT4_SLIM/BWA_unmasked

# input subfolders
SUBFOLDER[1]=SINGLE_TRANSCRIPTS

### Output folder
RR[1]=${SRCDIR_INI}/INFORMATIVE_SITES_exomeData_white_TETRA_100_SLIM_ALT1
RR[2]=${SRCDIR_INI}/INFORMATIVE_SITES_exomeData_white_TETRA_100_SLIM_ALT2
RR[3]=${SRCDIR_INI}/INFORMATIVE_SITES_exomeData_white_TETRA_100_SLIM_ALT3
RR[4]=${SRCDIR_INI}/INFORMATIVE_SITES_exomeData_white_TETRA_100_SLIM_ALT4

# Gene list
GL=/crex1/proj/snic2017-7-149/private/birch_phylogeny_scripts/P03_MS_alignment/00_gene_list.txt		
	
# MAIN sample list 
SPL=/crex1/proj/snic2017-7-149/private/birch_phylogeny_scripts/P03_MS_alignment/00_white_ploidyAnalysis_REDUX.txt

# reference sequence ID [can be any sequence in the MSA]
REFGEN[1]="PENDULA_FINLAND_LOIMAA.LOIMAA"
REFGEN[2]="NANA_FINLAND_ENONTEKIO.A002_2"
REFGEN[3]="PLATYPHYLLA_RUSSIA.A001_1"
REFGEN[4]="HUM02"

# Number of species
NSPEC=11

# FEATURE: GENE EXON INTRON
FT[1]=exons

# model name
MODEL_outfilename=MFP_ModelFinder



######## read sample list (MAIN)

declare -A SAMPLE_LIST

i=1

while read -r LINE                                                                      
do
   SAMPLE_LIST[$i]=$LINE
   SAMPLE_LIST_1[$i]=$LINE
   let "i+=1"
done < ${SPL}


### Number of samples
N_SAMPLES=${#SAMPLE_LIST_1[@]}




######## read gene list

i=1

while read -r LINE                                                                      
do
   GENE_LIST[$i]=$LINE
   let "i+=1"
done < ${GL}

N_GENES=${#GENE_LIST[@]}




for n in `seq 1 1 1`; do		# genomic feature (exons)

   GO_FEATURE=${FT[$n]}

   for j in `seq 4 1 4`; do  			# polarization

      for m in `seq 1 1 ${N_SAMPLES}`; do 			# polyploid ID

         GO_SAMPLE=${SAMPLE_LIST[$m]}
         old='-'
         new='_'
         GO_POLYP="${GO_SAMPLE//${old}/${new}}"	# Polyploid ID ('-' replaced by '_')
         GO_POLYP=${GO_POLYP^^}			# uppercase
         GO_REF=${REFGEN[$j]}

         GO_AA=${AA[$j]}/REFGEN_${GO_REF}/${SUBFOLDER[$n]}/${GO_SAMPLE}
         GO_RR=${RR[$j]}/REFGEN_${GO_REF}
         mkdir -p $GO_RR

         AUX_prefix_out="_SLIM-NoMAN02NoNAN01"

         #output file
         OUTfile=INFORMATIVE_SITES_03.${GO_FEATURE}_UFBoot_${MODEL_outfilename}_${GO_POLYP}${AUX_prefix_out}.txt
         rm -f ${GO_RR}/$OUTfile

         echo
         echo $OUTfile
         echo $GO_AA
         echo $GO_RR

         for i in `seq 1 1 $N_GENES`; do                 # cycle through genes

            GO_GENE=${GENE_LIST[$i]}

            #input file
            InFileSUF=${GO_GENE}.${GO_FEATURE}_trimN-ALT.fasta_UFBoot_${MODEL_outfilename}${AUX_prefix_out}.log

            cd $GO_AA

            # get number of informative sites
            NPS_AUX="$(cat ${InFileSUF} | grep 'parsimony-informative' | cut -f1 -d' ')"
            NPS[i]=$NPS_AUX

            if [[ "${NPS[$i]}" -ge 100 ]]
            then                  
               NPSF[i]=OK
            else
               NPSF[i]=FAIL
            fi

            #echo ${NPS[$i]} ${NPSF[$i]}

            # get number of sequences in alignment
            NSEQ[i]="$(cat ${InFileSUF} | grep 'Alignment has' | cut -f3 -d' ')"

            if [[ "${NSEQ[$i]}" = ${NSPEC} ]]
            then                  
               NSEQF[i]=OK
            else
               NSEQF[i]=FAIL
            fi

            # check whether there are sequences that contain more than 50% gaps/ambiguity
            GA_AUX="$(grep 'sequences contain more than 50% gaps/ambiguity' ${InFileSUF})"
         
            if [[ "${GA_AUX}" = "" ]]
            then                  
               GA[i]=OK
            else
               GA[i]=FAIL
            fi
 
            # check whether there are sequences with identical nucleotide composition
            INC_AUX="$(grep 'is identical to' ${InFileSUF})"

            if [[ "${INC_AUX}" = "" ]]
            then                  
               INC[i]=OK
            else
               INC[i]=FAIL
            fi

            # check whether bs analysis converged
            BS_AUX="$(grep 'WARNING: bootstrap analysis did not converge. You should rerun with higher number of iterations (-nm option)' ${InFileSUF})"

            if [[ "${BS_AUX}" = "" ]]
            then                  
               BS[i]=OK
            else
               BS[i]=FAIL
            fi


         done

         # save to file
         rm -f $GO_RR/_AUX_${OUTfile}
         for i in `seq 1 1 $N_GENES`; do
           echo ${GENE_LIST[$i]} ${NPS[$i]} ${NPSF[$i]} ${GA[$i]} ${INC[$i]} ${BS[$i]} ${NSEQ[$i]} ${NSEQF[$i]} >> $GO_RR/_AUX_${OUTfile}
         done

         # sort output file
         cat $GO_RR/_AUX_${OUTfile} | sort -k2 -nr > $GO_RR/${OUTfile}
         rm -f $GO_RR/_AUX_${OUTfile}

      done
   done
done


echo "Done!"




