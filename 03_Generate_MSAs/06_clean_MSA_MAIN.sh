#!/bin/bash
#
#SBATCH -J filterN
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 48:00:00
#SBATCH -A snic2022-22-589
#SBATCH -M snowy
#SBATCH --mail-user luis.leal@ebc.uu.se
#SBATCH --mail-type=FAIL
ulimit -c unlimited


## Scrip used for removing regions with too many masked positions (N) in the MSAs (site removed if # masked sites > 2)


# load modules
module load bioinfo-tools
module load trimAl/1.4.1


##### Paths and folders

#remember initial path
SRCDIR_INI=$(pwd) 

# Path to folder containing input data
AA=/crex1/proj/snic2020-6-184/private/Luis/P06_birch_exome_IGA-Sweden_2019/42_polarizationProtocol_MS_alignment_IUPAC_JG/BWA_unmasked

# input/output subfolder 1
RRsub[1]=exome_WGS
RRsub[2]=exome_WGS

# MAIN sample list 
SPL[1]=/crex1/proj/snic2017-7-149/private/birch_phylogeny_scripts/P03_MS_alignment/00_silver_ploidyAnalysis_REDUX.txt			
SPL[2]=/crex1/proj/snic2017-7-149/private/birch_phylogeny_scripts/P03_MS_alignment/00_white_ploidyAnalysis_REDUX.txt			

#file containing gene list (or transcript list)
GL=/crex1/proj/snic2017-7-149/private/birch_phylogeny_scripts/P03_MS_alignment/00_gene_list.txt	

# gap threshold
gapTH=0.857 		# maximum number of masked sites allowed is 2 [2/14=0.143]  >>>>>>> 1-0.143=0.857   

# Output file containing stats (root)
OUTF=Glean_09c_filter-masked



######## read sample list (MAIN)

declare -A SAMPLE_LIST

i=1

while read -r LINE                                                                      
do
   SAMPLE_LIST[1,$i]=$LINE
   SAMPLE_LIST_1[$i]=$LINE
   let "i+=1"
done < ${SPL[1]}


### Number of samples
N_SAMPLES[1]=${#SAMPLE_LIST_1[@]}



i=1

while read -r LINE                                                                      
do
   SAMPLE_LIST[2,$i]=$LINE
   SAMPLE_LIST_2[$i]=$LINE
   let "i+=1"
done < ${SPL[2]}


### Number of samples
N_SAMPLES[2]=${#SAMPLE_LIST_2[@]}




######## read gene list

i=1

while read -r LINE                                                                      
do
   GENE_LIST[$i]=$LINE
   let "i+=1"
done < ${GL}

N_GENES=${#GENE_LIST[@]}



# subfolders containing input data (gene or exon) 
SUBFOLDER[1]=SINGLE_GENES
SUBFOLDER[2]=SINGLE_TRANSCRIPTS
SUBFOLDER[3]=SINGLE_INTRONS

# genomeic feature
GF[1]=gene
GF[2]=exons
GF[3]=introns



cd $SRCDIR_INI

for k in `seq 1 1 2`; do 				# sample set (MAIN)

   for p in `seq 1 1 ${N_SAMPLES[$k]}`; do  				# loop through samples (MAIN)

      GO_SAMPLE=${SAMPLE_LIST[$k,$p]}
      GO_AA=${AA}/${RRsub[$k]}/${GO_SAMPLE}
      GO_RR=${GO_AA}
      GO_gapTH=${gapTH}

      for n in `seq 2 1 2`; do		# genomic feature (gene, exon, or intron)

               GO_FEATURE=${GF[$n]}

               echo $GO_AA
               echo $GO_RR
               echo $GO_gapTH
               echo

               OUTF_GO=${OUTF}_${GO_SAMPLE}_${GO_FEATURE}.txt
               rm -f $GO_RR/$OUTF_GO
               echo "gene" "Masked_before" "Masked_after" "Percentage_masked_before" "Percentage_masked_after"  > $GO_RR/$OUTF_GO
			
               for i in `seq 1 1 $N_GENES`; do                 # cycle through genes

                  GO_GENE=${GENE_LIST[$i]}
          	      REF_FASTA=${GO_GENE}.${GO_FEATURE}.fa
               
 	                # convert sequence to big caps
                  tr a-z A-Z < $GO_AA/$REF_FASTA > $SNIC_TMP/__${GO_GENE}_AUX_1.fa

                  # in sequence headers, convert '-' to '_'
                  sed '/^>/ s/-/_/g' $SNIC_TMP/__${GO_GENE}_AUX_1.fa > $SNIC_TMP/__${GO_GENE}_AUX_2.fa

                  # convert gaps (-) to X
                  sed 's/-/X/g' $SNIC_TMP/__${GO_GENE}_AUX_2.fa > $SNIC_TMP/__${GO_GENE}_AUX_3.fa
	 
                  # convert masked nucleotides (N) to gaps (-)
                  sed 's/N/-/g' $SNIC_TMP/__${GO_GENE}_AUX_3.fa > $SNIC_TMP/__${GO_GENE}_AUX_4.fa

                  # remove sites with too many missing sites (use trimAl, with missing sites coded as gaps)
                  trimal -in $SNIC_TMP/__${GO_GENE}_AUX_4.fa \
                         -out $SNIC_TMP/__${GO_GENE}_AUX_5.fa \
                         -gapthreshold $GO_gapTH \
                         -cons 0

    	            # -gapthreshold <n>    1 - (fraction of sequences with a gap allowed). Range: [0 - 1]
                  # -gapthreshold 0.904 >>>>> maximum number of masked sites allowed is 2 [2/21=0.096]  >>>>>>> 1-0.096=0.904
        		      # -gapthreshold 0.93 >>>>> maximum number of masked sites allowed is 2 [2/30=0.0667]  >>>>>>> 1-0.07=0.93
                  # -cons <n>                Minimum percentage of the positions in the original alignment to conserve. Range: [0 - 100]

                  # convert gaps (-) to masked nucleotides (N)
                  sed 's/-/N/g' $SNIC_TMP/__${GO_GENE}_AUX_5.fa > $SNIC_TMP/__${GO_GENE}_AUX_6.fa

                  # convert X sites to gaps (-)
                  sed 's/X/-/g' $SNIC_TMP/__${GO_GENE}_AUX_6.fa > $SNIC_TMP/__${GO_GENE}_AUX_7.fa
	   	   
  	              # convert sequence to big caps
                  tr a-z A-Z < $SNIC_TMP/__${GO_GENE}_AUX_7.fa > $SNIC_TMP/__${GO_GENE}_AUX_8.fa

                  # convert to single line fasta
                  awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' < $SNIC_TMP/__${GO_GENE}_AUX_8.fa > $SNIC_TMP/__${GO_GENE}_AUX_9.fa
                  cat $SNIC_TMP/__${GO_GENE}_AUX_9.fa | tail -n +2 > $GO_RR/${REF_FASTA%.fa}_trimN.fa


                  ### STATS

                  # Number of masked sites before trimming masked
                  N_before="$(cat $GO_AA/$REF_FASTA | grep -v '^>' | grep -io [N] | wc -l )"

                  # Number of masked sites after trimming masked
                  N_after="$(cat $GO_RR/${REF_FASTA%.fa}_trimN.fa | grep -v '^>' | grep -io [N] | wc -l )"

          	      # Total number of sites before trimming masked
                  T_before="$(cat $GO_AA/$REF_FASTA | grep -v '^>' | grep -io [A-T] | wc -l )"

                  # Total number of sites after trimming masked
                  T_after="$(cat $GO_RR/${REF_FASTA%.fa}_trimN.fa | grep -v '^>' | grep -io [A-T] | wc -l )"
		   
                  # Percentage masked, before trimming
                  aux_100="100"   
                  if [ -z "$N_before" ]; then
                     PERCT_before="0"
                  else
                     PERCT_before="$(awk '{print $1*$2/$3}' <<<"$aux_100 $N_before $T_before")"		   
                  fi

                  # Percentage masked, after trimming   
                  if [ -z "$N_after" ]; then
                     PERCT_after="0"
                  else
                     PERCT_after="$(awk '{print $1*$2/$3}' <<<"$aux_100 $N_after $T_after")"
                  fi	

                  echo $GO_GENE $N_before $N_after $PERCT_before $PERCT_after >> $GO_RR/$OUTF_GO
		   
	                # remove auxiliary files
                  rm -f $SNIC_TMP/__*.fa	   

               done
      done
   done
done


echo
echo done!
echo

