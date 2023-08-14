#!/bin/bash
#

## Read unpolarized MSAs and priors; generate polarized MSAs and apply polyploidization and hybridization model; run IQTREE2 and generate pairing profiles.


###################################### paths and folder names

###### remember initial path
SRCDIR_INI=$(pwd)                                           	 

###### input folder (MSAs, before polarization)
AA=${1:?msg}

###### output folder (root)
RR_MAIN=${2:?msg}

###### Simulation conditions (polyploidization models)
myMODEL=${3:?msg}

######  Number of replicates	
REP=${4:?msg}

######  Number of loci/genes		
NLOCUS=${5:?msg}

###### PRIORS
PP=${6:?msg}




###### Gene list
GL_large=/crex1/proj/snic2017-7-149/private/Luis/P11_SIMULATED_FASTA_PHYLOGENY_HYBRIDIZATION_exomeData/02_ABC_Simulated_annealing/00_GENE_LIST.txt




######################################  IQTree parameters

# Number of bootstraps
BS=1000

#### Evolutionary model
MODEL_IQTREE2=MFP		# MFP >>> ModelFinder used to determine best substitution model (but no correction for heterotachy)

####### MODEL RATE: List of rate heterogeneity among sites (when using Modelfinder and GHOST)
MODEL_IQTREE2_RATE="F*H"				

####### model name to be included in output file(s)
MODEL_outfilename=MFP_ModelFinder



######################################  Polyploidization and hybridization model

if [[ ${myMODEL} = "PPPP" ]] ; then
	
   ############## PPPP (pendula autopolyploid)		
   #
   # Parental species:
   GO_PAR[1]="PUBESCENS_PEN"	# pendula
   GO_PAR[2]="PUBESCENS_PEN"	# pendula
   GO_PAR[3]="PUBESCENS_PEN"	# pendula
   GO_PAR[4]="PUBESCENS_PEN"	# pendula
   ## Simulate pendula-platyphylla ancestral gene flow
   GF_2pend_species="PLATYPHYLLA_SIM"
   ### Hybridization blocks
   read -r GO_HYBR0[1] GO_HYBR0[2] GO_HYBR0[3] GO_HYBR0[4] <<< "PUBESCENS_PEN PUBESCENS_PEN PUBESCENS_PLA PUBESCENS_PLA"	# block 0
   read -r GO_HYBR1[1] GO_HYBR1[2] GO_HYBR1[3] GO_HYBR1[4] <<< "PUBESCENS_PEN PUBESCENS_PEN PUBESCENS_H PUBESCENS_N"	# block 1	
   read -r GO_HYBR2[1] GO_HYBR2[2] GO_HYBR2[3] GO_HYBR2[4] <<< "PUBESCENS_PEN PUBESCENS_PEN PUBESCENS_H PUBESCENS_PEN"	# block 2	
   read -r GO_HYBR3[1] GO_HYBR3[2] GO_HYBR3[3] GO_HYBR3[4] <<< "PUBESCENS_PEN PUBESCENS_PEN PUBESCENS_PEN PUBESCENS_N"	# block 3	
   read -r GO_HYBR4[1] GO_HYBR4[2] GO_HYBR4[3] GO_HYBR4[4] <<< "PUBESCENS_PEN PUBESCENS_PEN PUBESCENS_PEN PUBESCENS_PEN"	# block 4	
   read -r GO_HYBR5[1] GO_HYBR5[2] GO_HYBR5[3] GO_HYBR5[4] <<< "PUBESCENS_H PUBESCENS_N PUBESCENS_H PUBESCENS_N"	# block 5	
   read -r GO_HYBR6[1] GO_HYBR6[2] GO_HYBR6[3] GO_HYBR6[4] <<< "PUBESCENS_PEN PUBESCENS_PEN PUBESCENS_PEN PUBESCENS_PEN"	# block 6
   read -r GO_HYBR7[1] GO_HYBR7[2] GO_HYBR7[3] GO_HYBR7[4] <<< "PUBESCENS_PEN PUBESCENS_PEN PUBESCENS_PEN PUBESCENS_PEN"	# block 7	
   read -r GO_HYBR8[1] GO_HYBR8[2] GO_HYBR8[3] GO_HYBR8[4] <<< "PUBESCENS_PEN PUBESCENS_PEN PUBESCENS_PEN PUBESCENS_PEN"	# block 8
   
   ### Extra introgression 
   Extra_introgr_species="PUBESCENS_PEN"

fi


if [[ ${myMODEL} = "AAAA" ]] ; then
	
   ############## AAAA (ancestor of B. pendula/B. platyphylla, autopolyploid)		
   #
   # Parental species:
   GO_PAR[1]="PUBESCENS_ANC"	# ancestor of B. pendula/B. platyphylla
   GO_PAR[2]="PUBESCENS_ANC"	# ancestor of B. pendula/B. platyphylla
   GO_PAR[3]="PUBESCENS_ANC"	# ancestor of B. pendula/B. platyphylla
   GO_PAR[4]="PUBESCENS_ANC"	# ancestor of B. pendula/B. platyphylla
   ## Simulate pendula-platyphylla ancestral gene flow
   GF_2pend_species="PLATYPHYLLA_SIM"
   ### Hybridization blocks
   read -r GO_HYBR0[1] GO_HYBR0[2] GO_HYBR0[3] GO_HYBR0[4] <<< "PUBESCENS_ANC PUBESCENS_ANC PUBESCENS_PLA PUBESCENS_PLA"	# block 0
   read -r GO_HYBR1[1] GO_HYBR1[2] GO_HYBR1[3] GO_HYBR1[4] <<< "PUBESCENS_ANC PUBESCENS_ANC PUBESCENS_H PUBESCENS_N"	# block 1	
   read -r GO_HYBR2[1] GO_HYBR2[2] GO_HYBR2[3] GO_HYBR2[4] <<< "PUBESCENS_ANC PUBESCENS_ANC PUBESCENS_H PUBESCENS_ANC"	# block 2	
   read -r GO_HYBR3[1] GO_HYBR3[2] GO_HYBR3[3] GO_HYBR3[4] <<< "PUBESCENS_ANC PUBESCENS_ANC PUBESCENS_ANC PUBESCENS_N"	# block 3	
   read -r GO_HYBR4[1] GO_HYBR4[2] GO_HYBR4[3] GO_HYBR4[4] <<< "PUBESCENS_ANC PUBESCENS_ANC PUBESCENS_ANC PUBESCENS_ANC"	# block 4	
   read -r GO_HYBR5[1] GO_HYBR5[2] GO_HYBR5[3] GO_HYBR5[4] <<< "PUBESCENS_H PUBESCENS_N PUBESCENS_H PUBESCENS_N"	# block 5	
   read -r GO_HYBR6[1] GO_HYBR6[2] GO_HYBR6[3] GO_HYBR6[4] <<< "PUBESCENS_PEN PUBESCENS_PEN PUBESCENS_PEN PUBESCENS_PEN"	# block 6
   read -r GO_HYBR7[1] GO_HYBR7[2] GO_HYBR7[3] GO_HYBR7[4] <<< "PUBESCENS_ANC PUBESCENS_ANC PUBESCENS_ANC PUBESCENS_ANC"	# block 7	
   read -r GO_HYBR8[1] GO_HYBR8[2] GO_HYBR8[3] GO_HYBR8[4] <<< "PUBESCENS_ANC PUBESCENS_ANC PUBESCENS_ANC PUBESCENS_ANC"	# block 8
   
   ### Extra introgression
   Extra_introgr_species="PUBESCENS_PEN"

fi


if [ ${myMODEL} = "PPNN" ] || [ ${myMODEL} = "PPNN_HE" ] ; then
	
   ############## PPNN (pendula/nana allopolyploid)		
   #
   # Parental species:
   GO_PAR[1]="PUBESCENS_PEN"	# pendula
   GO_PAR[2]="PUBESCENS_PEN"	# pendula
   GO_PAR[3]="PUBESCENS_N"	# nana
   GO_PAR[4]="PUBESCENS_N"	# nana
   ## Simulate pendula-platyphylla ancestral gene flow
   GF_2pend_species="PLATYPHYLLA_SIM"
   ### Hybridization blocks
   read -r GO_HYBR0[1] GO_HYBR0[2] GO_HYBR0[3] GO_HYBR0[4] <<< "PUBESCENS_PLA PUBESCENS_PLA PUBESCENS_N PUBESCENS_N"	# block 0
   read -r GO_HYBR1[1] GO_HYBR1[2] GO_HYBR1[3] GO_HYBR1[4] <<< "PUBESCENS_PEN PUBESCENS_PEN PUBESCENS_H PUBESCENS_H"	# block 1	
   read -r GO_HYBR2[1] GO_HYBR2[2] GO_HYBR2[3] GO_HYBR2[4] <<< "PUBESCENS_PEN PUBESCENS_PEN PUBESCENS_H PUBESCENS_N"	# block 2	
   read -r GO_HYBR3[1] GO_HYBR3[2] GO_HYBR3[3] GO_HYBR3[4] <<< "PUBESCENS_PEN PUBESCENS_PEN PUBESCENS_N PUBESCENS_H"	# block 3	
   read -r GO_HYBR4[1] GO_HYBR4[2] GO_HYBR4[3] GO_HYBR4[4] <<< "PUBESCENS_H PUBESCENS_H PUBESCENS_N PUBESCENS_N"	# block 4	
   read -r GO_HYBR5[1] GO_HYBR5[2] GO_HYBR5[3] GO_HYBR5[4] <<< "PUBESCENS_PEN PUBESCENS_PEN PUBESCENS_PEN PUBESCENS_PEN"	# block 5	
   read -r GO_HYBR6[1] GO_HYBR6[2] GO_HYBR6[3] GO_HYBR6[4] <<< "PUBESCENS_N PUBESCENS_N PUBESCENS_N PUBESCENS_N"	# block 6
   read -r GO_HYBR7[1] GO_HYBR7[2] GO_HYBR7[3] GO_HYBR7[4] <<< "PUBESCENS_PEN PUBESCENS_PEN PUBESCENS_N PUBESCENS_N"	# block 7	
   read -r GO_HYBR8[1] GO_HYBR8[2] GO_HYBR8[3] GO_HYBR8[4] <<< "PUBESCENS_PEN PUBESCENS_PEN PUBESCENS_N PUBESCENS_N"	# block 8
   
   ### Extra introgression
   Extra_introgr_species="PUBESCENS_PEN"

fi


if [ ${myMODEL} = "PPHH" ] || [ ${myMODEL} = "PPHH_HE" ] ; then
	
   ############## PPHH (pendula/humilis allopolyploid)		
   #
   # Parental species:
   GO_PAR[1]="PUBESCENS_PEN"	# pendula
   GO_PAR[2]="PUBESCENS_PEN"	# pendula
   GO_PAR[3]="PUBESCENS_H"	# humilis
   GO_PAR[4]="PUBESCENS_H"	# humilis
   ## Simulate pendula-platyphylla ancestral gene flow
   GF_2pend_species="PLATYPHYLLA_SIM"
   ### Hybridization blocks
   read -r GO_HYBR0[1] GO_HYBR0[2] GO_HYBR0[3] GO_HYBR0[4] <<< "PUBESCENS_PLA PUBESCENS_PLA PUBESCENS_H PUBESCENS_H"	# block 0
   read -r GO_HYBR1[1] GO_HYBR1[2] GO_HYBR1[3] GO_HYBR1[4] <<< "PUBESCENS_PEN PUBESCENS_PEN PUBESCENS_N PUBESCENS_N"	# block 1	
   read -r GO_HYBR2[1] GO_HYBR2[2] GO_HYBR2[3] GO_HYBR2[4] <<< "PUBESCENS_PEN PUBESCENS_PEN PUBESCENS_N PUBESCENS_H"	# block 2	
   read -r GO_HYBR3[1] GO_HYBR3[2] GO_HYBR3[3] GO_HYBR3[4] <<< "PUBESCENS_PEN PUBESCENS_PEN PUBESCENS_H PUBESCENS_N"	# block 3	
   read -r GO_HYBR4[1] GO_HYBR4[2] GO_HYBR4[3] GO_HYBR4[4] <<< "PUBESCENS_N PUBESCENS_N PUBESCENS_H PUBESCENS_H"	# block 4	
   read -r GO_HYBR5[1] GO_HYBR5[2] GO_HYBR5[3] GO_HYBR5[4] <<< "PUBESCENS_PEN PUBESCENS_PEN PUBESCENS_PEN PUBESCENS_PEN"	# block 5	
   read -r GO_HYBR6[1] GO_HYBR6[2] GO_HYBR6[3] GO_HYBR6[4] <<< "PUBESCENS_H PUBESCENS_H PUBESCENS_H PUBESCENS_H"	# block 6
   read -r GO_HYBR7[1] GO_HYBR7[2] GO_HYBR7[3] GO_HYBR7[4] <<< "PUBESCENS_PEN PUBESCENS_PEN PUBESCENS_H PUBESCENS_H"	# block 7	
   read -r GO_HYBR8[1] GO_HYBR8[2] GO_HYBR8[3] GO_HYBR8[4] <<< "PUBESCENS_PEN PUBESCENS_PEN PUBESCENS_H PUBESCENS_H"	# block 8
   
   ### Extra introgression 
   Extra_introgr_species="PUBESCENS_PEN"

fi


if [ ${myMODEL} = "PPPyPy" ] || [ ${myMODEL} = "PPPyPy_HE" ] ; then
	
   ############## PPPyPy (pendula/platyphylla allopolyploid)		
   #
   # Parental species:
   GO_PAR[1]="PUBESCENS_PEN"	# pendula
   GO_PAR[2]="PUBESCENS_PEN"	# pendula
   GO_PAR[3]="PUBESCENS_PLA"	# platyphylla
   GO_PAR[4]="PUBESCENS_PLA"	# platyphylla
   ## Simulate pendula-platyphylla ancestral gene flow
   GF_2pend_species="PLATYPHYLLA_SIM"
   ### Hybridization blocks
   read -r GO_HYBR0[1] GO_HYBR0[2] GO_HYBR0[3] GO_HYBR0[4] <<< "PUBESCENS_PEN PUBESCENS_PEN PUBESCENS_H PUBESCENS_H"	# block 0
   read -r GO_HYBR1[1] GO_HYBR1[2] GO_HYBR1[3] GO_HYBR1[4] <<< "PUBESCENS_PEN PUBESCENS_PEN PUBESCENS_H PUBESCENS_N"	# block 1	
   read -r GO_HYBR2[1] GO_HYBR2[2] GO_HYBR2[3] GO_HYBR2[4] <<< "PUBESCENS_PEN PUBESCENS_PEN PUBESCENS_H PUBESCENS_PLA"	# block 2	
   read -r GO_HYBR3[1] GO_HYBR3[2] GO_HYBR3[3] GO_HYBR3[4] <<< "PUBESCENS_PEN PUBESCENS_PEN PUBESCENS_PLA PUBESCENS_N"	# block 3	
   read -r GO_HYBR4[1] GO_HYBR4[2] GO_HYBR4[3] GO_HYBR4[4] <<< "PUBESCENS_H PUBESCENS_N PUBESCENS_PLA PUBESCENS_PLA"	# block 4	
   read -r GO_HYBR5[1] GO_HYBR5[2] GO_HYBR5[3] GO_HYBR5[4] <<< "PUBESCENS_PEN PUBESCENS_PEN PUBESCENS_PEN PUBESCENS_PEN"	# block 5	
   read -r GO_HYBR6[1] GO_HYBR6[2] GO_HYBR6[3] GO_HYBR6[4] <<< "PUBESCENS_H PUBESCENS_N PUBESCENS_H PUBESCENS_N"	# block 6
   read -r GO_HYBR7[1] GO_HYBR7[2] GO_HYBR7[3] GO_HYBR7[4] <<< "PUBESCENS_PEN PUBESCENS_PEN PUBESCENS_N PUBESCENS_N"	# block 7	
   read -r GO_HYBR8[1] GO_HYBR8[2] GO_HYBR8[3] GO_HYBR8[4] <<< "PUBESCENS_PLA PUBESCENS_PLA PUBESCENS_PLA PUBESCENS_PLA"	# block 8
   
   ### Extra introgression 
   Extra_introgr_species="PUBESCENS_PEN"

fi


if [ ${myMODEL} = "PPNH" ] || [ ${myMODEL} = "PPNH_HE" ] ; then
	
   ############## PPNH (pendula/nana-humilis allopolyploid)		
   #
   # Parental species:
   GO_PAR[1]="PUBESCENS_PEN"	# pendula
   GO_PAR[2]="PUBESCENS_PEN"	# pendula
   GO_PAR[3]="PUBESCENS_N"	# nana
   GO_PAR[4]="PUBESCENS_H"	# humilis
   ## Simulate pendula-platyphylla ancestral gene flow
   GF_2pend_species="PLATYPHYLLA_SIM"
   ### Hybridization blocks
   read -r GO_HYBR0[1] GO_HYBR0[2] GO_HYBR0[3] GO_HYBR0[4] <<< "PUBESCENS_PLA PUBESCENS_PLA PUBESCENS_N PUBESCENS_H"	# block 0
   read -r GO_HYBR1[1] GO_HYBR1[2] GO_HYBR1[3] GO_HYBR1[4] <<< "PUBESCENS_PEN PUBESCENS_PEN PUBESCENS_PEN PUBESCENS_PEN"	# block 1	
   read -r GO_HYBR2[1] GO_HYBR2[2] GO_HYBR2[3] GO_HYBR2[4] <<< "PUBESCENS_PEN PUBESCENS_PEN PUBESCENS_PEN PUBESCENS_H"	# block 2	
   read -r GO_HYBR3[1] GO_HYBR3[2] GO_HYBR3[3] GO_HYBR3[4] <<< "PUBESCENS_PEN PUBESCENS_PEN PUBESCENS_N PUBESCENS_PEN"	# block 3	
   read -r GO_HYBR4[1] GO_HYBR4[2] GO_HYBR4[3] GO_HYBR4[4] <<< "PUBESCENS_N PUBESCENS_H PUBESCENS_N PUBESCENS_H"	# block 4	
   read -r GO_HYBR5[1] GO_HYBR5[2] GO_HYBR5[3] GO_HYBR5[4] <<< "PUBESCENS_N PUBESCENS_N PUBESCENS_N PUBESCENS_N"	# block 5	
   read -r GO_HYBR6[1] GO_HYBR6[2] GO_HYBR6[3] GO_HYBR6[4] <<< "PUBESCENS_H PUBESCENS_H PUBESCENS_H PUBESCENS_H"	# block 6
   read -r GO_HYBR7[1] GO_HYBR7[2] GO_HYBR7[3] GO_HYBR7[4] <<< "PUBESCENS_PEN PUBESCENS_PEN PUBESCENS_N PUBESCENS_H"	# block 7	
   read -r GO_HYBR8[1] GO_HYBR8[2] GO_HYBR8[3] GO_HYBR8[4] <<< "PUBESCENS_PEN PUBESCENS_PEN PUBESCENS_N PUBESCENS_H"	# block 8
   
   ### Extra introgression 
   Extra_introgr_species="PUBESCENS_PEN"

fi


if [ ${myMODEL} = "AANN" ] || [ ${myMODEL} = "AANN_HE" ] ; then
	
   ############## AANN (pendula-platyphylla ancestor /nana allopolyploid)		
   #
   # Parental species:
   GO_PAR[1]="PUBESCENS_ANC"	# ancestor of B. pendula/B. platyphylla
   GO_PAR[2]="PUBESCENS_ANC"	# ancestor of B. pendula/B. platyphylla
   GO_PAR[3]="PUBESCENS_N"	# nana
   GO_PAR[4]="PUBESCENS_N"	# nana
   ## Simulate pendula-platyphylla ancestral gene flow
   GF_2pend_species="PLATYPHYLLA_SIM"
   ### Hybridization blocks
   read -r GO_HYBR0[1] GO_HYBR0[2] GO_HYBR0[3] GO_HYBR0[4] <<< "PUBESCENS_PLA PUBESCENS_PLA PUBESCENS_N PUBESCENS_N"	# block 0
   read -r GO_HYBR1[1] GO_HYBR1[2] GO_HYBR1[3] GO_HYBR1[4] <<< "PUBESCENS_ANC PUBESCENS_ANC PUBESCENS_H PUBESCENS_H"	# block 1	
   read -r GO_HYBR2[1] GO_HYBR2[2] GO_HYBR2[3] GO_HYBR2[4] <<< "PUBESCENS_ANC PUBESCENS_ANC PUBESCENS_H PUBESCENS_N"	# block 2	
   read -r GO_HYBR3[1] GO_HYBR3[2] GO_HYBR3[3] GO_HYBR3[4] <<< "PUBESCENS_ANC PUBESCENS_ANC PUBESCENS_N PUBESCENS_H"	# block 3	
   read -r GO_HYBR4[1] GO_HYBR4[2] GO_HYBR4[3] GO_HYBR4[4] <<< "PUBESCENS_H PUBESCENS_H PUBESCENS_N PUBESCENS_N"	# block 4	
   read -r GO_HYBR5[1] GO_HYBR5[2] GO_HYBR5[3] GO_HYBR5[4] <<< "PUBESCENS_ANC PUBESCENS_ANC PUBESCENS_ANC PUBESCENS_ANC"	# block 5	
   read -r GO_HYBR6[1] GO_HYBR6[2] GO_HYBR6[3] GO_HYBR6[4] <<< "PUBESCENS_N PUBESCENS_N PUBESCENS_N PUBESCENS_N"	# block 6
   read -r GO_HYBR7[1] GO_HYBR7[2] GO_HYBR7[3] GO_HYBR7[4] <<< "PUBESCENS_PEN PUBESCENS_PEN PUBESCENS_PEN PUBESCENS_PEN"	# block 7	
   read -r GO_HYBR8[1] GO_HYBR8[2] GO_HYBR8[3] GO_HYBR8[4] <<< "PUBESCENS_ANC PUBESCENS_ANC PUBESCENS_N PUBESCENS_N"	# block 8
   
   ### Extra introgression 
   Extra_introgr_species="PUBESCENS_PEN"

fi


if [ ${myMODEL} = "AAHH" ] || [ ${myMODEL} = "AAHH_HE" ] ; then
	
   ############## AAHH (pendula-platyphylla ancestor /humilis allopolyploid)		
   #
   # Parental species:
   GO_PAR[1]="PUBESCENS_ANC"	# ancestor of B. pendula/B. platyphylla
   GO_PAR[2]="PUBESCENS_ANC"	# ancestor of B. pendula/B. platyphylla
   GO_PAR[3]="PUBESCENS_H"	# humilis
   GO_PAR[4]="PUBESCENS_H"	# humilis
   ## Simulate pendula-platyphylla ancestral gene flow
   GF_2pend_species="PLATYPHYLLA_SIM"
   ### Hybridization blocks
   read -r GO_HYBR0[1] GO_HYBR0[2] GO_HYBR0[3] GO_HYBR0[4] <<< "PUBESCENS_PLA PUBESCENS_PLA PUBESCENS_H PUBESCENS_H"	# block 0
   read -r GO_HYBR1[1] GO_HYBR1[2] GO_HYBR1[3] GO_HYBR1[4] <<< "PUBESCENS_ANC PUBESCENS_ANC PUBESCENS_N PUBESCENS_N"	# block 1	
   read -r GO_HYBR2[1] GO_HYBR2[2] GO_HYBR2[3] GO_HYBR2[4] <<< "PUBESCENS_ANC PUBESCENS_ANC PUBESCENS_N PUBESCENS_H"	# block 2	
   read -r GO_HYBR3[1] GO_HYBR3[2] GO_HYBR3[3] GO_HYBR3[4] <<< "PUBESCENS_ANC PUBESCENS_ANC PUBESCENS_H PUBESCENS_N"	# block 3	
   read -r GO_HYBR4[1] GO_HYBR4[2] GO_HYBR4[3] GO_HYBR4[4] <<< "PUBESCENS_N PUBESCENS_N PUBESCENS_H PUBESCENS_H"	# block 4	
   read -r GO_HYBR5[1] GO_HYBR5[2] GO_HYBR5[3] GO_HYBR5[4] <<< "PUBESCENS_ANC PUBESCENS_ANC PUBESCENS_ANC PUBESCENS_ANC"	# block 5	
   read -r GO_HYBR6[1] GO_HYBR6[2] GO_HYBR6[3] GO_HYBR6[4] <<< "PUBESCENS_H PUBESCENS_H PUBESCENS_H PUBESCENS_H"	# block 6
   read -r GO_HYBR7[1] GO_HYBR7[2] GO_HYBR7[3] GO_HYBR7[4] <<< "PUBESCENS_PEN PUBESCENS_PEN PUBESCENS_PEN PUBESCENS_PEN"	# block 7	
   read -r GO_HYBR8[1] GO_HYBR8[2] GO_HYBR8[3] GO_HYBR8[4] <<< "PUBESCENS_ANC PUBESCENS_ANC PUBESCENS_H PUBESCENS_H"	# block 8
   
   ### Extra introgression 
   Extra_introgr_species="PUBESCENS_PEN"

fi


if [ ${myMODEL} = "AANH" ] || [ ${myMODEL} = "AANH_HE" ] ; then
	
   ############## PPNH (pendula-platyphylla ancestor / nana-humilis allopolyploid)		
   #
   # Parental species:
   GO_PAR[1]="PUBESCENS_ANC"	# ancestor of B. pendula/B. platyphylla
   GO_PAR[2]="PUBESCENS_ANC"	# ancestor of B. pendula/B. platyphylla
   GO_PAR[3]="PUBESCENS_N"	# nana
   GO_PAR[4]="PUBESCENS_H"	# humilis
   ## Simulate pendula-platyphylla ancestral gene flow
   GF_2pend_species="PLATYPHYLLA_SIM"
   ### Hybridization blocks
   read -r GO_HYBR0[1] GO_HYBR0[2] GO_HYBR0[3] GO_HYBR0[4] <<< "PUBESCENS_PLA PUBESCENS_PLA PUBESCENS_N PUBESCENS_H"	# block 0
   read -r GO_HYBR1[1] GO_HYBR1[2] GO_HYBR1[3] GO_HYBR1[4] <<< "PUBESCENS_ANC PUBESCENS_ANC PUBESCENS_N PUBESCENS_H"	# block 1	
   read -r GO_HYBR2[1] GO_HYBR2[2] GO_HYBR2[3] GO_HYBR2[4] <<< "PUBESCENS_ANC PUBESCENS_ANC PUBESCENS_ANC PUBESCENS_H"	# block 2	
   read -r GO_HYBR3[1] GO_HYBR3[2] GO_HYBR3[3] GO_HYBR3[4] <<< "PUBESCENS_ANC PUBESCENS_ANC PUBESCENS_N PUBESCENS_ANC"	# block 3	
   read -r GO_HYBR4[1] GO_HYBR4[2] GO_HYBR4[3] GO_HYBR4[4] <<< "PUBESCENS_N PUBESCENS_H PUBESCENS_N PUBESCENS_H"	# block 4	
   read -r GO_HYBR5[1] GO_HYBR5[2] GO_HYBR5[3] GO_HYBR5[4] <<< "PUBESCENS_N PUBESCENS_N PUBESCENS_N PUBESCENS_N"	# block 5	
   read -r GO_HYBR6[1] GO_HYBR6[2] GO_HYBR6[3] GO_HYBR6[4] <<< "PUBESCENS_H PUBESCENS_H PUBESCENS_H PUBESCENS_H"	# block 6
   read -r GO_HYBR7[1] GO_HYBR7[2] GO_HYBR7[3] GO_HYBR7[4] <<< "PUBESCENS_PEN PUBESCENS_PEN PUBESCENS_PEN PUBESCENS_PEN"	# block 7	
   read -r GO_HYBR8[1] GO_HYBR8[2] GO_HYBR8[3] GO_HYBR8[4] <<< "PUBESCENS_ANC PUBESCENS_ANC PUBESCENS_ANC PUBESCENS_ANC"	# block 8
   
   ### Extra introgression 
   Extra_introgr_species="PUBESCENS_PEN"

fi



############################################################################  START
   
################### Read priors

unset PRIOR_LIST
x=1

while read -r LINE                                                                      
do
   PRIOR_LIST[$x]=$LINE
   let "x+=1"
done < ${PP}

N_PRIORS=${#PRIOR_LIST[@]}



################### Loci of interest
cat $GL_large | head -${NLOCUS} > $SNIC_TMP/00_gene_list_TRIM_MAX-00000.txt
GO_GL=$SNIC_TMP/00_gene_list_TRIM_MAX-00000.txt



###################
################### POLARIZATION
###################
echo
echo "1. POLARIZATION..." 
date -u

# ILS 
GO_ILS_reduction=${PRIOR_LIST[1]}
GO_FOLDER=A_modILS_birch_CLEAN_reducedILS_${GO_ILS_reduction}perct
	  
# input folder
GO_AA=${AA}/${GO_FOLDER}
	  
# Number of species (including root species)
GO_NSPEC=14	
	  
# platyphylla to pendula gene flow [0-1]
GO_HYB_FRACTION_2pen=${PRIOR_LIST[2]}
	 
# Hybridization levels for each block
GO_HYBRY_FRACTION0=${PRIOR_LIST[3]}
GO_HYBRY_FRACTION1=${PRIOR_LIST[4]}
GO_HYBRY_FRACTION2=${PRIOR_LIST[5]}
GO_HYBRY_FRACTION3=${PRIOR_LIST[6]}
GO_HYBRY_FRACTION4=${PRIOR_LIST[7]}
GO_HYBRY_FRACTION5=${PRIOR_LIST[8]}
GO_HYBRY_FRACTION6=${PRIOR_LIST[9]}
GO_HYBRY_FRACTION7=${PRIOR_LIST[10]}
GO_HYBRY_FRACTION8=${PRIOR_LIST[11]}
	  
# Hybridization levels for extra species
GO_HYBRY_EXTRA=${PRIOR_LIST[12]}


######## POLAR 1: Reference sequence: B. pendula
#echo "Polarization 1..."
	  
# Polarizing reference sequence ID
GO_REFspecies_1="PENDULA_SIM"
		 
# polarization label
NALT=1

# output folder
GO_RR=$SNIC_TMP/01_POLAR/ALT${NALT}
rm -rf $GO_RR
mkdir -p $GO_RR

#run in parallel for different replicates 
for r_aux in $( eval echo {01..${REP}} ); do		# replicates; must include leading zeros (REPLICATE_01, REPLICATE_02, etc)
            $SRCDIR_INI/21_polarize_MSA_query.sh $GO_AA $GO_RR $GO_NSPEC $GO_FOLDER $r_aux $NLOCUS "${GO_PAR[@]}" $GO_REFspecies_1 $NALT $GF_2pend_species $GO_HYB_FRACTION_2pen "${GO_HYBR0[@]}" ${GO_HYBRY_FRACTION0} \
                                  "${GO_HYBR1[@]}" "${GO_HYBR2[@]}" "${GO_HYBR3[@]}" "${GO_HYBR4[@]}" "${GO_HYBR5[@]}" "${GO_HYBR6[@]}" "${GO_HYBR7[@]}" "${GO_HYBR8[@]}" \
								${GO_HYBRY_FRACTION1} ${GO_HYBRY_FRACTION2} ${GO_HYBRY_FRACTION3} ${GO_HYBRY_FRACTION4} ${GO_HYBRY_FRACTION5} ${GO_HYBRY_FRACTION6} ${GO_HYBRY_FRACTION7} ${GO_HYBRY_FRACTION8} \
								${Extra_introgr_species} ${GO_HYBRY_EXTRA} ${myMODEL} &
done
		 
#wait


######## POLAR 2: Reference sequence: B. nana
#echo "Polarization 2..."
	  
# Polarizing reference sequence ID
GO_REFspecies_1="NANA_SIM"
	 
# polarization label
NALT=2

# output folder
GO_RR=$SNIC_TMP/01_POLAR/ALT${NALT}
rm -rf $GO_RR
mkdir -p $GO_RR

#run in parallel for different replicates 
for r_aux in $( eval echo {01..${REP}} ); do		# replicates; must include leading zeros (REPLICATE_01, REPLICATE_02, etc)
            $SRCDIR_INI/21_polarize_MSA_query.sh $GO_AA $GO_RR $GO_NSPEC $GO_FOLDER $r_aux $NLOCUS "${GO_PAR[@]}" $GO_REFspecies_1 $NALT $GF_2pend_species $GO_HYB_FRACTION_2pen "${GO_HYBR0[@]}" ${GO_HYBRY_FRACTION0} \
                               "${GO_HYBR1[@]}" "${GO_HYBR2[@]}" "${GO_HYBR3[@]}" "${GO_HYBR4[@]}" "${GO_HYBR5[@]}" "${GO_HYBR6[@]}" "${GO_HYBR7[@]}" "${GO_HYBR8[@]}" \
							${GO_HYBRY_FRACTION1} ${GO_HYBRY_FRACTION2} ${GO_HYBRY_FRACTION3} ${GO_HYBRY_FRACTION4} ${GO_HYBRY_FRACTION5} ${GO_HYBRY_FRACTION6} ${GO_HYBRY_FRACTION7} ${GO_HYBRY_FRACTION8} \
							${Extra_introgr_species} ${GO_HYBRY_EXTRA} ${myMODEL} &
done
		 
#wait
	  
	  
######## POLAR 3: Reference sequence: B. humilis
#echo "Polarization 3..."
	  
# Polarizing reference sequence ID
GO_REFspecies_1="HUMILIS_SIM"
 
# polarization label
NALT=3

# output folder
GO_RR=$SNIC_TMP/01_POLAR/ALT${NALT}
rm -rf $GO_RR
mkdir -p $GO_RR

#run in parallel for different replicates 
		for r_aux in $( eval echo {01..${REP}} ); do		# replicates; must include leading zeros (REPLICATE_01, REPLICATE_02, etc)
            $SRCDIR_INI/21_polarize_MSA_query.sh $GO_AA $GO_RR $GO_NSPEC $GO_FOLDER $r_aux $NLOCUS "${GO_PAR[@]}" $GO_REFspecies_1 $NALT $GF_2pend_species $GO_HYB_FRACTION_2pen "${GO_HYBR0[@]}" ${GO_HYBRY_FRACTION0} \
                               "${GO_HYBR1[@]}" "${GO_HYBR2[@]}" "${GO_HYBR3[@]}" "${GO_HYBR4[@]}" "${GO_HYBR5[@]}" "${GO_HYBR6[@]}" "${GO_HYBR7[@]}" "${GO_HYBR8[@]}" \
							${GO_HYBRY_FRACTION1} ${GO_HYBRY_FRACTION2} ${GO_HYBRY_FRACTION3} ${GO_HYBRY_FRACTION4} ${GO_HYBRY_FRACTION5} ${GO_HYBRY_FRACTION6} ${GO_HYBRY_FRACTION7} ${GO_HYBRY_FRACTION8} \
							${Extra_introgr_species} ${GO_HYBRY_EXTRA} ${myMODEL} &
done
		 
#wait
  

######## POLAR 4: Reference sequence: B. platyphylla
#echo "Polarization 4..."

# Polarizing reference sequence ID
GO_REFspecies_1="PLATYPHYLLA_SIM"

# polarization label
NALT=4

# output folder
GO_RR=$SNIC_TMP/01_POLAR/ALT${NALT}
rm -rf $GO_RR
mkdir -p $GO_RR

#run in parallel for different replicates 
for r_aux in $( eval echo {01..${REP}} ); do		# replicates; must include leading zeros (REPLICATE_01, REPLICATE_02, etc)
            $SRCDIR_INI/21_polarize_MSA_query.sh $GO_AA $GO_RR $GO_NSPEC $GO_FOLDER $r_aux $NLOCUS "${GO_PAR[@]}" $GO_REFspecies_1 $NALT $GF_2pend_species $GO_HYB_FRACTION_2pen "${GO_HYBR0[@]}" ${GO_HYBRY_FRACTION0} \
                               "${GO_HYBR1[@]}" "${GO_HYBR2[@]}" "${GO_HYBR3[@]}" "${GO_HYBR4[@]}" "${GO_HYBR5[@]}" "${GO_HYBR6[@]}" "${GO_HYBR7[@]}" "${GO_HYBR8[@]}" \
							${GO_HYBRY_FRACTION1} ${GO_HYBRY_FRACTION2} ${GO_HYBRY_FRACTION3} ${GO_HYBRY_FRACTION4} ${GO_HYBRY_FRACTION5} ${GO_HYBRY_FRACTION6} ${GO_HYBRY_FRACTION7} ${GO_HYBRY_FRACTION8} \
							${Extra_introgr_species} ${GO_HYBRY_EXTRA} ${myMODEL} &
done
         
wait




###################
################### IQ-TREE
###################
echo
echo "2. IQ-TREE..." 
date -u
		  
AA2=$SNIC_TMP/01_POLAR
RR2=$SNIC_TMP/02_IQTREE
	  
#GO_FOLDER=${INfolder[$k]}
#GO_POLAR=ALT
#GO_NLOCUS=${NLOCUS}
	  
#input file suffix
InFileSUF=TRUE_CLEAN-ALT.fasta
		  
# Outgroup sequence
OUTG="ALNUS_SIM"
	  
######## POLAR 1: Reference sequence: B. pendula
#echo "Polarization 1..."
		  
# polarization label
NALT=1

GO_AA2=${AA2}/ALT${NALT}
GO_RR2=${RR2}/ALT${NALT}
rm -rf $GO_RR2
mkdir -p $GO_RR2
			   
#run in parallel for different replicates 
for r_aux in $( eval echo {01..${REP}} ); do		# replicates; must include leading zeros (REPLICATE_01, REPLICATE_02, etc)
	  		      $SRCDIR_INI/22_IQ-TREE2_SINGLEGenes_ultrafastBS_SimData_query.sh $GO_AA2 $InFileSUF $GO_RR2 $BS ${MODEL_outfilename} $MODEL_IQTREE2 $MODEL_IQTREE2_RATE $NLOCUS $OUTG $r_aux $NALT &
done
			   
#wait
			
			
			   
######## POLAR 2: Reference sequence: B. nana
#echo "Polarization 2..."
	  
# polarization label
NALT=2

GO_AA2=${AA2}/ALT${NALT}
GO_RR2=${RR2}/ALT${NALT}
rm -rf $GO_RR2
mkdir -p $GO_RR2
   
#run in parallel for different replicates 
for r_aux in $( eval echo {01..${REP}} ); do		# replicates; must include leading zeros (REPLICATE_01, REPLICATE_02, etc)
		    $SRCDIR_INI/22_IQ-TREE2_SINGLEGenes_ultrafastBS_SimData_query.sh $GO_AA2 $InFileSUF $GO_RR2 $BS ${MODEL_outfilename} $MODEL_IQTREE2 $MODEL_IQTREE2_RATE $NLOCUS $OUTG $r_aux $NALT &
done
		   
#wait  
			   
			   
			   
######## POLAR 3: Reference sequence: B. humilis
#echo "Polarization 3..."
	  
# polarization label
NALT=3

GO_AA2=${AA2}/ALT${NALT}
GO_RR2=${RR2}/ALT${NALT}
rm -rf $GO_RR2
mkdir -p $GO_RR2

#run in parallel for different replicates 
for r_aux in $( eval echo {01..${REP}} ); do		# replicates; must include leading zeros (REPLICATE_01, REPLICATE_02, etc)
		    $SRCDIR_INI/22_IQ-TREE2_SINGLEGenes_ultrafastBS_SimData_query.sh $GO_AA2 $InFileSUF $GO_RR2 $BS ${MODEL_outfilename} $MODEL_IQTREE2 $MODEL_IQTREE2_RATE $NLOCUS $OUTG $r_aux $NALT &
done
		 
#wait    
			   
			   
######## POLAR 4: Reference sequence: B. platyphylla
#echo "Polarization 4..."
	  
# polarization label
NALT=4

GO_AA2=${AA2}/ALT${NALT}
GO_RR2=${RR2}/ALT${NALT}
rm -rf $GO_RR2
mkdir -p $GO_RR2

#run in parallel for different replicates 
for r_aux in $( eval echo {01..${REP}} ); do		# replicates; must include leading zeros (REPLICATE_01, REPLICATE_02, etc)
		    $SRCDIR_INI/22_IQ-TREE2_SINGLEGenes_ultrafastBS_SimData_query.sh $GO_AA2 $InFileSUF $GO_RR2 $BS ${MODEL_outfilename} $MODEL_IQTREE2 $MODEL_IQTREE2_RATE $NLOCUS $OUTG $r_aux $NALT &
done
		 
wait

   
               
			   
###################
################### ## Create file with IQ-TEST2 consensus tree obtained for each gene | Use either maximum number of genes available (MAX) or random gene lists previously created
###################			   
echo
echo "3. Glean IQ-TREE gene trees..." 
date -u
	  
AA3=$SNIC_TMP/02_IQTREE
RR3=${RR_MAIN}/03_PAIRINGfrequencies
		 
# input file suffix
#file_aux=""
InFileSUF=TRUE_CLEAN_UFBoot_${MODEL_outfilename}.iqtree
			   
### read gene list

unset GENE_LIST
x=1

while read -r LINE                                                                      
do
   GENE_LIST[$x]=$LINE
   let "x+=1"
 done < ${GO_GL}

N_GENES=${#GENE_LIST[@]}
		 
		 
######## POLAR 1: Reference sequence: B. pendula
#echo "Polarization 1..."
	  
# polarization label
NALT=1
		 
GO_AA3=${AA3}/ALT${NALT}
GO_RR3=${RR3}/ALT${NALT}
rm -rf $GO_RR3
		 
for r_aux in $( eval echo {01..${REP}} ); do		# replicates; must include leading zeros (REPLICATE_01, REPLICATE_02, etc)

   mkdir -p $GO_RR3/$r_aux
  
   ########output file
   OUTfile=GENE_TREES_${r_aux}_UFBoot_MAX-${N_GENES}.newick
   rm -f ${GO_RR3}/$r_aux/$OUTfile
			   
   for i in `seq 1 1 $N_GENES`; do                 # cycle through genes
			   
                  #input file
                  GO_GENE=${GENE_LIST[$i]}
                  InFile=data_${GO_GENE}_${InFileSUF}
                  #InFile=data_${GO_GENE}_TRUE_CLEAN_UFBoot_${MODEL_outfilename}.iqtree
			
                  # get IQ-TEST2 consensus tree for each gene           
                  cat $GO_AA3/$r_aux/${InFile} | grep -A 2 'Consensus tree in newick format:' | tail -1 >> ${GO_RR3}/$r_aux/$OUTfile
   done
done
			
# get genes from all replicates
cat ${GO_RR3}/*/GENE_TREES_*_UFBoot_MAX-*.newick > ${GO_RR3}/GENE_TREES_ALL_REPLICATES.newick
		 


######## POLAR 2: Reference sequence: B. nana
#echo "Polarization 2..."
	  
# polarization label
NALT=2
	 
GO_AA3=${AA3}/ALT${NALT}
GO_RR3=${RR3}/ALT${NALT}
rm -rf $GO_RR3
	 
for r_aux in $( eval echo {01..${REP}} ); do		# replicates; must include leading zeros (REPLICATE_01, REPLICATE_02, etc)

   mkdir -p $GO_RR3/$r_aux

   ########output file
   OUTfile=GENE_TREES_${r_aux}_UFBoot_MAX-${N_GENES}.newick
   rm -f ${GO_RR3}/$r_aux/$OUTfile
		    
   for i in `seq 1 1 $N_GENES`; do                 # cycle through genes
		   
               #input file
               GO_GENE=${GENE_LIST[$i]}
               InFile=data_${GO_GENE}_${InFileSUF}
               #InFile=data_${GO_GENE}_TRUE_CLEAN_UFBoot_${MODEL_outfilename}.iqtree
		
               # get IQ-TEST2 consensus tree for each gene           
               cat $GO_AA3/$r_aux/${InFile} | grep -A 2 'Consensus tree in newick format:' | tail -1 >> ${GO_RR3}/$r_aux/$OUTfile
   done
done
		
# get genes from all replicates
cat ${GO_RR3}/*/GENE_TREES_*_UFBoot_MAX-*.newick > ${GO_RR3}/GENE_TREES_ALL_REPLICATES.newick
	 
		   
		   
######## POLAR 3: Reference sequence: B. humilis
#echo "Polarization 3..."
	  
# polarization label
NALT=3
	 
GO_AA3=${AA3}/ALT${NALT}
GO_RR3=${RR3}/ALT${NALT}
rm -rf $GO_RR3
	 
for r_aux in $( eval echo {01..${REP}} ); do		# replicates; must include leading zeros (REPLICATE_01, REPLICATE_02, etc)

   mkdir -p $GO_RR3/$r_aux

   ########output file
   OUTfile=GENE_TREES_${r_aux}_UFBoot_MAX-${N_GENES}.newick
   rm -f ${GO_RR3}/$r_aux/$OUTfile
		      
   for i in `seq 1 1 $N_GENES`; do                 # cycle through genes
		   
               #input file
               GO_GENE=${GENE_LIST[$i]}
               InFile=data_${GO_GENE}_${InFileSUF}
               #InFile=data_${GO_GENE}_TRUE_CLEAN_UFBoot_${MODEL_outfilename}.iqtree
		
               # get IQ-TEST2 consensus tree for each gene           
               cat $GO_AA3/$r_aux/${InFile} | grep -A 2 'Consensus tree in newick format:' | tail -1 >> ${GO_RR3}/$r_aux/$OUTfile
   done
done
		
# get genes from all replicates
cat ${GO_RR3}/*/GENE_TREES_*_UFBoot_MAX-*.newick > ${GO_RR3}/GENE_TREES_ALL_REPLICATES.newick
	 

    
######## POLAR 4: Reference sequence: B. platyphylla
#echo "Polarization 4..."
	  
# polarization label
NALT=4
	 
GO_AA3=${AA3}/ALT${NALT}
GO_RR3=${RR3}/ALT${NALT}
rm -rf $GO_RR3
	 
for r_aux in $( eval echo {01..${REP}} ); do		# replicates; must include leading zeros (REPLICATE_01, REPLICATE_02, etc)

   mkdir -p $GO_RR3/$r_aux

   ########output file
   OUTfile=GENE_TREES_${r_aux}_UFBoot_MAX-${N_GENES}.newick
   rm -f ${GO_RR3}/$r_aux/$OUTfile
		   
   for i in `seq 1 1 $N_GENES`; do                 # cycle through genes
		   
               #input file
               GO_GENE=${GENE_LIST[$i]}
               InFile=data_${GO_GENE}_${InFileSUF}
               #InFile=data_${GO_GENE}_TRUE_CLEAN_UFBoot_${MODEL_outfilename}.iqtree
		
               # get IQ-TEST2 consensus tree for each gene           
               cat $GO_AA3/$r_aux/${InFile} | grep -A 2 'Consensus tree in newick format:' | tail -1 >> ${GO_RR3}/$r_aux/$OUTfile
   done
done
		
# get genes from all replicates
cat ${GO_RR3}/*/GENE_TREES_*_UFBoot_MAX-*.newick > ${GO_RR3}/GENE_TREES_ALL_REPLICATES.newick
	 


  
###################
################### For each gene family, ID species closer to polyploid
###################	
echo
echo "4. Pairing frequencies..." 
date -u

AA4=${RR_MAIN}/03_PAIRINGfrequencies
	  
## tetraploid IDs
TETRAPLOID="TETRAPLOID_SIM"
	  
## outgroup
outgroup_ID="ALNUS_SIM"
		 

######## POLAR 1: Reference sequence: B. pendula
#echo "Polarization 1..."
	  
# polarization label
NALT=1
		  
GO_AA4=${AA4}/ALT${NALT}
GO_RR4=$GO_AA4
		  
## output file
OUTFile="sister_ID_analysis_priors-ALT${NALT}.txt"
		 
Rscript --no-save /crex1/proj/snic2017-7-149/private/Luis/P11_SIMULATED_FASTA_PHYLOGENY_HYBRIDIZATION_exomeData/02_ABC_Simulated_annealing/24_IQTREE_gene_tree_DISTANCE_ALT.R \
				                                                     ${GO_AA4}/GENE_TREES_ALL_REPLICATES.newick \
                                                                     $TETRAPLOID \
                                                                     $outgroup_ID \
                                                                     $GO_RR4 \
                                                                     $OUTFile
		  
$SRCDIR_INI/23_pairing_profile_summary.sh $GO_RR4 $OUTFile
		  

######## POLAR 2: Reference sequence: B. nana
#echo "Polarization 2..."
		  
# polarization label
NALT=2
		  
GO_AA4=${AA4}/ALT${NALT}
GO_RR4=$GO_AA4
		  
## output file
OUTFile="sister_ID_analysis_priors-ALT${NALT}.txt"
		 
Rscript --no-save /crex1/proj/snic2017-7-149/private/Luis/P11_SIMULATED_FASTA_PHYLOGENY_HYBRIDIZATION_exomeData/02_ABC_Simulated_annealing/24_IQTREE_gene_tree_DISTANCE_ALT.R \
				                                                     ${GO_AA4}/GENE_TREES_ALL_REPLICATES.newick \
                                                                     $TETRAPLOID \
                                                                     $outgroup_ID \
                                                                     $GO_RR4 \
                                                                     $OUTFile
			 
$SRCDIR_INI/23_pairing_profile_summary.sh $GO_RR4 $OUTFile
		  
		  
######## POLAR 3: Reference sequence: B. humilis
#echo "Polarization 3..."
		  
# polarization label
NALT=3
		  
GO_AA4=${AA4}/ALT${NALT}
GO_RR4=$GO_AA4
		  
## output file
OUTFile="sister_ID_analysis_priors-ALT${NALT}.txt"
		 
Rscript --no-save /crex1/proj/snic2017-7-149/private/Luis/P11_SIMULATED_FASTA_PHYLOGENY_HYBRIDIZATION_exomeData/02_ABC_Simulated_annealing/24_IQTREE_gene_tree_DISTANCE_ALT.R \
				                                                     ${GO_AA4}/GENE_TREES_ALL_REPLICATES.newick \
                                                                     $TETRAPLOID \
                                                                     $outgroup_ID \
                                                                     $GO_RR4 \
                                                                     $OUTFile
			 
$SRCDIR_INI/23_pairing_profile_summary.sh $GO_RR4 $OUTFile
		  
		  
######## POLAR 4: Reference sequence: B. platyphylla
#echo "Polarization 4..."
		  
# polarization label
NALT=4
		  
GO_AA4=${AA4}/ALT${NALT}
GO_RR4=$GO_AA4
		  
## output file
OUTFile="sister_ID_analysis_priors-ALT${NALT}.txt"
		 
Rscript --no-save /crex1/proj/snic2017-7-149/private/Luis/P11_SIMULATED_FASTA_PHYLOGENY_HYBRIDIZATION_exomeData/02_ABC_Simulated_annealing/24_IQTREE_gene_tree_DISTANCE_ALT.R \
				                                                     ${GO_AA4}/GENE_TREES_ALL_REPLICATES.newick \
                                                                     $TETRAPLOID \
                                                                     $outgroup_ID \
                                                                     $GO_RR4 \
                                                                     $OUTFile
			 
$SRCDIR_INI/23_pairing_profile_summary.sh $GO_RR4 $OUTFile
		  
		  
####### Joint pairing profile (all polarization geometries)

cat $AA4/ALT1/sister_ID_analysis_priors-ALT1_summary.txt $AA4/ALT2/sister_ID_analysis_priors-ALT2_summary.txt $AA4/ALT3/sister_ID_analysis_priors-ALT3_summary.txt $AA4/ALT4/sister_ID_analysis_priors-ALT4_summary.txt | grep -v '^p'> $AA4/sister_ID_analysis_priors-JOINT_summary.txt		  



####### Clean scratch disk
rm -f $SNIC_TMP/_*
rm -f $SNIC_TMP/data*
rm -f $SNIC_TMP/*.py

rm -f $SNIC_TMP/*_UFBoot_*
rm -f $SNIC_TMP/*fa
rm -f $SNIC_TMP/*fasta
rm -f $SNIC_TMP/AUX*
rm -rf $SNIC_TMP/ALT*
rm -rf $SNIC_TMP/0*

######################################################### END

#echo
#echo "Phylo analysis done!"
#date -u






