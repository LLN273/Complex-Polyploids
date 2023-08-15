#  
# Luis Leal (2021)
#
# Written in Python 3
# UPPMAX: module load python3/3.6.0



### Script used to polarize tetraploid sequence in simulated MSAs 
### Can also be used to simulate introgression events into tetraploid. 



######################################################### LOAD STANDARD MODULES

import sys
import re                                                 
import ast                                                
import os                                                 
import time
from collections import Counter
from pdb import set_trace as bp                           
import random						  






######################################################## OPEN INPUT FILES


try:
    fhand = open(sys.argv[1], 'r')                  			# input fasta file
except:
    print('\n Error: input file missing.')
    exit() 


try:
    Nspecies = sys.argv[2]                  				# Number of species
except:
    print('\n Error: provide the following arguments:')
    print('21s_polarizeMSA_fullCODED_ALT_TETRA_HYBRIDIZATION_BOTH_SinglePolarRef_2HYBRIDS_maxHETEROZY_step_2.py [input.fasta] [No. species] [refSequence ID1] [HYB species 1] [HYB species 2] [HYB species 3] [HYB species 4] [hybridization code] [ancestral ID1] [ancestral ID2] [ancestral ID3] [ancestral ID4] [ancestral ID5]  \n')
    exit() 


try:
    RefGenID = sys.argv[3]                  				# Reference sequence ID
except:
    print('\n Error: provide the following arguments:')
    print('21s_polarizeMSA_fullCODED_ALT_TETRA_HYBRIDIZATION_BOTH_SinglePolarRef_2HYBRIDS_maxHETEROZY_step_2.py [input.fasta] [No. species] [refSequence ID1] [HYB species 1] [HYB species 2] [HYB species 3] [HYB species 4] [hybridization code] [ancestral ID1] [ancestral ID2] [ancestral ID3] [ancestral ID4] [ancestral ID5] \n')
    exit() 


try:
    HYBRID_1 = sys.argv[4]                  				# hybridizing sequence 1 ID
except:
    print('\n Error: provide the following arguments:')
    print('21s_polarizeMSA_fullCODED_ALT_TETRA_HYBRIDIZATION_BOTH_SinglePolarRef_2HYBRIDS_maxHETEROZY_step_2.py [input.fasta] [No. species] [refSequence ID1] [HYB species 1] [HYB species 2] [HYB species 3] [HYB species 4] [hybridization code] [ancestral ID1] [ancestral ID2] [ancestral ID3] [ancestral ID4] [ancestral ID5] \n')
    exit() 


try:
    HYBRID_2 = sys.argv[5]                  				# hybridizing sequence 2 ID
except:
    print('\n Error: provide the following arguments:')
    print('21s_polarizeMSA_fullCODED_ALT_TETRA_HYBRIDIZATION_BOTH_SinglePolarRef_2HYBRIDS_maxHETEROZY_step_2.py [input.fasta] [No. species] [refSequence ID1] [HYB species 1] [HYB species 2] [HYB species 3] [HYB species 4] [hybridization code] [ancestral ID1] [ancestral ID2] [ancestral ID3] [ancestral ID4] [ancestral ID5] \n')
    exit() 


try:
    HYBRID_3 = sys.argv[6]                  				# hybridizing sequence 3 ID
except:
    print('\n Error: provide the following arguments:')
    print('21s_polarizeMSA_fullCODED_ALT_TETRA_HYBRIDIZATION_BOTH_SinglePolarRef_2HYBRIDS_maxHETEROZY_step_2.py [input.fasta] [No. species] [refSequence ID1] [HYB species 1] [HYB species 2] [HYB species 3] [HYB species 4] [hybridization code] [ancestral ID1] [ancestral ID2] [ancestral ID3] [ancestral ID4] [ancestral ID5] \n')
    exit() 


try:
    HYBRID_4 = sys.argv[7]                  				# hybridizing sequence 4 ID
except:
    print('\n Error: provide the following arguments:')
    print('21s_polarizeMSA_fullCODED_ALT_TETRA_HYBRIDIZATION_BOTH_SinglePolarRef_2HYBRIDS_maxHETEROZY_step_2.py [input.fasta] [No. species] [refSequence ID1] [HYB species 1] [HYB species 2] [HYB species 3] [HYB species 4] [hybridization code] [ancestral ID1] [ancestral ID2] [ancestral ID3] [ancestral ID4] [ancestral ID5] \n')
    exit() 


try:
    HYBRID_CODE = sys.argv[8]                  				# hybridization coding (eg '0110' specifies that for this locus there are 2 gene copies from 'hybridizing species 1')
except:
    print('\n Error: provide the following arguments:')
    print('21s_polarizeMSA_fullCODED_ALT_TETRA_HYBRIDIZATION_BOTH_SinglePolarRef_2HYBRIDS_maxHETEROZY_step_2.py [input.fasta] [No. species] [refSequence ID1] [HYB species 1] [HYB species 2] [HYB species 3] [HYB species 4] [hybridization code] [ancestral ID1] [ancestral ID2] [ancestral ID3] [ancestral ID4] [ancestral ID5] \n')
    exit() 


try:
    DELSEQ = sys.argv[9]                  				# Sequence to be removed from output phylogeny  
except:
    print('\n Error: provide the following arguments:')
    print('21s_polarizeMSA_fullCODED_ALT_TETRA_HYBRIDIZATION_BOTH_SinglePolarRef_2HYBRIDS_maxHETEROZY_step_2.py [input.fasta] [No. species] [refSequence ID1] [HYB species 1] [HYB species 2] [HYB species 3] [HYB species 4] [hybridization code] [ancestral ID1] [ancestral ID2] [ancestral ID3] [ancestral ID4] [ancestral ID5] \n')
    exit() 


try:
    DELSEQ2 = sys.argv[10]                  				# Sequence to be removed from output phylogeny  
except:
    print('\n Error: provide the following arguments:')
    print('21s_polarizeMSA_fullCODED_ALT_TETRA_HYBRIDIZATION_BOTH_SinglePolarRef_2HYBRIDS_maxHETEROZY_step_2.py [input.fasta] [No. species] [refSequence ID1] [HYB species 1] [HYB species 2] [HYB species 3] [HYB species 4] [hybridization code] [ancestral ID1] [ancestral ID2] [ancestral ID3] [ancestral ID4] [ancestral ID5] \n')
    exit() 


try:
    DELSEQ3 = sys.argv[11]                  				# Sequence to be removed from output phylogeny  
except:
    print('\n Error: provide the following arguments:')
    print('21s_polarizeMSA_fullCODED_ALT_TETRA_HYBRIDIZATION_BOTH_SinglePolarRef_2HYBRIDS_maxHETEROZY_step_2.py [input.fasta] [No. species] [refSequence ID1] [HYB species 1] [HYB species 2] [HYB species 3] [HYB species 4] [hybridization code] [ancestral ID1] [ancestral ID2] [ancestral ID3] [ancestral ID4] [ancestral ID5] \n')
    exit() 


try:
    DELSEQ4 = sys.argv[12]                  				# Sequence to be removed from output phylogeny  
except:
    print('\n Error: provide the following arguments:')
    print('21s_polarizeMSA_fullCODED_ALT_TETRA_HYBRIDIZATION_BOTH_SinglePolarRef_2HYBRIDS_maxHETEROZY_step_2.py [input.fasta] [No. species] [refSequence ID1] [HYB species 1] [HYB species 2] [HYB species 3] [HYB species 4] [hybridization code] [ancestral ID1] [ancestral ID2] [ancestral ID3] [ancestral ID4] [ancestral ID5] \n')
    exit() 


try:
    DELSEQ5 = sys.argv[13]                  				# Sequence to be removed from output phylogeny  
except:
    print('\n Error: provide the following arguments:')
    print('21s_polarizeMSA_fullCODED_ALT_TETRA_HYBRIDIZATION_BOTH_SinglePolarRef_2HYBRIDS_maxHETEROZY_step_2.py [input.fasta] [No. species] [refSequence ID1] [HYB species 1] [HYB species 2] [HYB species 3] [HYB species 4] [hybridization code] [ancestral ID1] [ancestral ID2] [ancestral ID3] [ancestral ID4] [ancestral ID5] \n')
    exit() 




### AUXILIARY FUNCTION: code nucleotide; extended coding for heterozygotic sites according to IUPAC nomenclature:
### See also: http://www.iqtree.org/doc/Frequently-Asked-Questions#how-does-iq-tree-treat-gapmissingambiguous-characters

def extendedCODE(input_list):
   species_nuc_uniq = list(set(input_list))  # get unique nucleotides
   species_nuc_uniq.sort()		     # sort list alphabetically

   outvalue = ''
   if species_nuc_uniq == ['A'] : outvalue = 'A'
   if species_nuc_uniq == ['C'] : outvalue = 'C'
   if species_nuc_uniq == ['G'] : outvalue = 'G'
   if species_nuc_uniq == ['T'] : outvalue = 'T'
   if species_nuc_uniq == ['A','C'] : outvalue = 'M'
   if species_nuc_uniq == ['A','G'] : outvalue = 'R'
   if species_nuc_uniq == ['A','T'] : outvalue = 'W'
   if species_nuc_uniq == ['C','G'] : outvalue = 'S'
   if species_nuc_uniq == ['C','T'] : outvalue = 'Y'
   if species_nuc_uniq == ['G','T'] : outvalue = 'K'
   if species_nuc_uniq == ['A','C','G'] : outvalue = 'V'
   if species_nuc_uniq == ['A','C','T'] : outvalue = 'H'
   if species_nuc_uniq == ['A','G','T'] : outvalue = 'D'
   if species_nuc_uniq == ['C','G','T'] : outvalue = 'B'
   if species_nuc_uniq == ['A','C','G','T'] : outvalue = 'N'

   return outvalue

### end function




### AUXILIARY FUNCTION: expand IUPAC code

def unfoldCODE(species_nuc):

   outvalue = ''
   if species_nuc == 'A' : outvalue = ['A']
   if species_nuc == 'C' : outvalue = ['C']
   if species_nuc == 'G' : outvalue = ['G']
   if species_nuc == 'T' : outvalue = ['T']
   if species_nuc == 'M' : outvalue = ['A','C']
   if species_nuc == 'R' : outvalue = ['A','G']
   if species_nuc == 'W' : outvalue = ['A','T']
   if species_nuc == 'S' : outvalue = ['C','G']
   if species_nuc == 'Y' : outvalue = ['C','T']
   if species_nuc == 'K' : outvalue = ['G','T']
   if species_nuc == 'V' : outvalue = ['A','C','G']
   if species_nuc == 'H' : outvalue = ['A','C','T']
   if species_nuc == 'D' : outvalue = ['A','G','T']
   if species_nuc == 'B' : outvalue = ['C','G','T']
   if species_nuc == 'N' : outvalue = ['N']
   if outvalue == '' :
      print('\n Error: sequences contain non-IUPAC nucleotide character\n')
      exit()

   return outvalue

### end function










######################################################## Read fasta file

INfasta_species = list()
INfasta_seq = list()
FLAG_aux = 0
for line in fhand:
    aux_NC = line.replace("\n", "")
    if FLAG_aux == 0 :
      INfasta_species.append(aux_NC)
      FLAG_aux = 1
    else :
      INfasta_seq.append(aux_NC)
      FLAG_aux = 0

# remove empty entries ("")
INfasta_species_CLEAN = INfasta_species
INfasta_seq_CLEAN = INfasta_seq


# save species/haplotype and associated sequence to dictionary
INfasta_dict = dict()
for i in range(len(INfasta_species_CLEAN)) :
   INfasta_dict[INfasta_species_CLEAN[i]] = INfasta_seq_CLEAN[i]

#for i in range(len(INfasta_species_CLEAN)) :
#   print(INfasta_species_CLEAN[i], INfasta_dict[INfasta_species_CLEAN[i]])

#print(len(INfasta_dict))



# get species list from fasta file
species_list = list()
for tag in INfasta_species_CLEAN :
   tag = tag.replace(">", "")                   #remove > from sequence name
   species_list.append(tag[:(len(tag) - 4)])	#remove "_0_0" from "14_0_0" ["14" is the species name]

species_list = list(set(species_list))		# unique elements only
#print(species_list)


# Check number of sequences
DIC_LEN=len(INfasta_dict)
if DIC_LEN != (2 * int(Nspecies)) :
   print('DIC_LEN:', DIC_LEN, '; Nspecies:', Nspecies)
   print('\n Error: check number of sequences in fasta file\n')
   exit()

if len(species_list) != (int(Nspecies) + 2) :
   print('len(species_list):', len(species_list), ';Nspecies:', Nspecies)
   print('\n Error: check number of sequences in fasta file\n')
   exit()



# Check all sequences have the same length
FLAG_aux = 0
for i in range(len(INfasta_species_CLEAN)) :
   len_seq = len(INfasta_dict[INfasta_species_CLEAN[i]])
   #print(len_seq)
   if FLAG_aux == 0 :
      len_ref = len_seq
      FLAG_aux = 1
   else :
      if len_seq != len_ref :
         print('\n Error: sequences have different lenghts\n')
         exit()

          





###################################################### Create diploid/tetraploid sequences

####### Create tetraploid sequences, IUPAC coded

### no hybridization
ANC1_ID = '>TETRAPLOID_SIM_1_0_0'
ANC2_ID = '>TETRAPLOID_SIM_2_0_0'
ANC3_ID = '>TETRAPLOID_SIM_3_0_0'
ANC4_ID = '>TETRAPLOID_SIM_4_0_0'


### hybridization: if there is hybridization followed by introgression between HYBRID species (source) and TETRAploid species (target)

if HYBRID_CODE[0]  == "1" :
   ANC1_ID = '>' + HYBRID_1 + '_0_0'


if HYBRID_CODE[1]  == "1" :
   ANC2_ID = '>' + HYBRID_2 + '_0_1'


if HYBRID_CODE[2]  == "1" :
   ANC3_ID = '>' + HYBRID_3 + '_0_0'


if HYBRID_CODE[3]  == "1" :
   ANC4_ID = '>' + HYBRID_4 + '_0_1'


#print(HYBRID_CODE)
#print(ANC1_ID)
#print(ANC2_ID)
#print(ANC3_ID)
#print(ANC4_ID)

ANC1seq = list()
ANC2seq = list()
ANC3seq = list()
ANC4seq = list()

for i in INfasta_dict[ANC1_ID] :
   ANC1seq.append(i)

for i in INfasta_dict[ANC2_ID] :
   ANC2seq.append(i)

for i in INfasta_dict[ANC3_ID] :
   ANC3seq.append(i)

for i in INfasta_dict[ANC4_ID] :
   ANC4seq.append(i)

TETRAseq = list()
TETRAseq_aux = list(zip(ANC1seq,ANC2seq,ANC3seq,ANC4seq))

# reduce each locus to its variants, IUPAC coded (mimicking MSA obtained from WGS data)
for i in TETRAseq_aux :
   species_nuc_uniq = extendedCODE(i)
   TETRAseq.append(species_nuc_uniq)
   #print(i,species_nuc_uniq)

TETRAseq = "".join(TETRAseq)



## Create diploid sequences (IUPAC coded)

DIPdic = INfasta_dict	# dictionary only with diploid sequences (before linking haplotypes for same individual)


# remove sequences not present in final phylogeny
DELSEQ1_aux = '>' + DELSEQ + '_0_0'
DELSEQ2_aux = '>' + DELSEQ + '_0_1'
DELSEQ3_aux = '>' + DELSEQ2 + '_0_0'
DELSEQ4_aux = '>' + DELSEQ2 + '_0_1'
DELSEQ5_aux = '>' + DELSEQ3 + '_0_0'
DELSEQ6_aux = '>' + DELSEQ3 + '_0_1'
DELSEQ7_aux = '>' + DELSEQ4 + '_0_0'
DELSEQ8_aux = '>' + DELSEQ4 + '_0_1'

DELSEQ9_aux = '>TETRAPLOID_SIM_1_0_0'
DELSEQ10_aux = '>TETRAPLOID_SIM_2_0_0'
DELSEQ11_aux = '>TETRAPLOID_SIM_3_0_0'
DELSEQ12_aux = '>TETRAPLOID_SIM_4_0_0'

DELSEQ13_aux = '>' + DELSEQ5 + '_0_0'
DELSEQ14_aux = '>' + DELSEQ5 + '_0_1'


del DIPdic[DELSEQ1_aux]
del DIPdic[DELSEQ2_aux]
del DIPdic[DELSEQ3_aux]
del DIPdic[DELSEQ4_aux]
del DIPdic[DELSEQ5_aux]
del DIPdic[DELSEQ6_aux]
del DIPdic[DELSEQ7_aux]
del DIPdic[DELSEQ8_aux]
del DIPdic[DELSEQ9_aux]
del DIPdic[DELSEQ10_aux]
del DIPdic[DELSEQ11_aux]
del DIPdic[DELSEQ12_aux]

try:                            
   del DIPdic[DELSEQ13_aux]
   del DIPdic[DELSEQ14_aux]
except:
   print("Allopolyploid")



DEP_IDs = list()
for key, value in DIPdic.items():
    tag = key.replace(">", "")                   #remove > from sequence name
    DEP_IDs.append(tag[:(len(tag) - 4)])	#remove "_0_0" from "14_0_0" ["14" is the species name]
    #print(key)

DEP_IDs = list(set(DEP_IDs))		# unique elements only
#DEP_IDs.remove('0')			# remove outgroup


DICT_seq = dict()	# dictionary containing linked haplotypes, all species

for item in DEP_IDs :
   hap_A = '>' + item + '_0_0'
   hap_B = '>' + item + '_0_1'
   seq_A = list()
   seq_B = list()

   for i in INfasta_dict[hap_A] :
      seq_A.append(i)

   for i in INfasta_dict[hap_B] :
      seq_B.append(i)

   hap_AB = list()
   hap_AB_aux = list(zip(seq_A,seq_B))

   for i in hap_AB_aux :
      species_nuc_uniq = extendedCODE(i)
      hap_AB.append(species_nuc_uniq)
      #print(item, i,species_nuc_uniq)

   hap_AB = "".join(hap_AB)

   DICT_seq[item] = hap_AB



## add outgroup to dictionary
#hap_A = '>' + '0' + '_0_0'
#DICT_seq['0'] = INfasta_dict[hap_A]

## add tetraploid to dictionary
POLYPID = 'TETRAPLOID_SIM'
DICT_seq[POLYPID] = TETRAseq

## tag reference sequence
#REFG_ID =  RefGenID + 'R'
#DICT_seq[REFG_ID] = DICT_seq.pop(RefGenID)


#for key, value in DICT_seq.items():
#   print(key,value)

#for key, value in DICT_seq.items():
#   print(key,len(DICT_seq[key]))







###################################################### Unfold tetraploid and ref sequences

## Reference sequence 1
REFseq = list()
for i in DICT_seq[RefGenID] :
   nuc_ste = unfoldCODE(i)
   REFseq.append(nuc_ste)
   #print(i,nuc_ste)



## Polyploid sequence
TETRAseq_u = list()
for i in DICT_seq[POLYPID] :
   nuc_ste = unfoldCODE(i)
   TETRAseq_u.append(nuc_ste)
   #print(i,nuc_ste)










###################################################### Create polarized sequence for ALLOTETRAPLOID

DICT_POLAR_ALT = dict()	        # dictionary containing polarized polyploid sequence (also all other sequences, but which stay as in the input MSA)

###  Add non-polyploid sequences to dictionary
for key, value in DICT_seq.items():
   aux_ID = '>' + key
   DICT_POLAR_ALT[aux_ID] = value
   #print(aux_ID)
   #print(key,DICT_POLAR_ALT[aux_ID])


# remove unpolarized polyploid from dictionary (polarized polyploid will be added in the next step)
POLY_ID = '>' + POLYPID
in_POLY = DICT_POLAR_ALT[POLY_ID]
del DICT_POLAR_ALT[POLY_ID]

#for key, value in DICT_POLAR_ALT.items():
   #print(key,len(value))
   #print(key,value)



##### tetraploid

ALT_list = list()
species_seq = TETRAseq_u[:]
for i in range(len(REFseq)) :		# for each nucleotide in reference sequence
   ref_nuc_ls = REFseq[i]
   species_nuc = list(species_seq[i])
   #print(ref_nuc_ls, species_nuc)
   #
   if species_nuc[0] == 'N' :		# if polyploid is masked
      aux_nuc = 'N' 
      ALT_list.append(aux_nuc)
      #print(ref_nuc_ls, species_nuc, aux_nuc) 
   elif ref_nuc_ls[0] == 'N' :		# if reference is masked
      aux_nuc = 'N' 
      ALT_list.append(aux_nuc)
      #print(ref_nuc_ls, species_nuc, aux_nuc)
   elif (species_nuc == ['A'] or species_nuc == ['C'] or species_nuc == ['G'] or species_nuc == ['T']) :		# polyploid if fixed
      ALT_list.append(species_nuc[0])
      #print(ref_nuc_ls, species_nuc, species_nuc[0])
   elif species_nuc == ref_nuc_ls :			# if reference alleles were found in the polyploid and there are no other alleles present
      aux_nuc = extendedCODE(species_nuc)
      ALT_list.append(aux_nuc)
      #print(ref_nuc_ls, species_nuc, aux_nuc)
   elif all(x0 in species_nuc for x0 in ref_nuc_ls) :		# if all reference alleles were found also in polyploid, but there are other alleles present as well
      aux_nuc0 = [pp for pp in species_nuc if pp not in ref_nuc_ls]		# remove reference alleles from polyploid allele set
      if aux_nuc0 == [] : aux_nuc0 = species_nuc
      aux_nuc = extendedCODE(aux_nuc0)
      ALT_list.append(aux_nuc)
      #print(ref_nuc_ls, species_nuc, aux_nuc)
   elif any(x0 in ref_nuc_ls for x0 in species_nuc) :		# if at least one of the reference alleles is found also in the polyploid
      aux_nuc0 = [pp for pp in species_nuc if pp not in ref_nuc_ls]		# remove reference alleles from polyploid allele set, if present
      if aux_nuc0 == [] : aux_nuc0 = species_nuc
      aux_nuc = extendedCODE(aux_nuc0)
      ALT_list.append(aux_nuc)
      #print(ref_nuc_ls, species_nuc, aux_nuc)
   else :							# reference allele not found in polyploid
      aux_nuc = extendedCODE(species_nuc)
      ALT_list.append(aux_nuc)
      #print(ref_nuc_ls, species_nuc, aux_nuc)






# update dictionary
DICT_POLAR_ALT[POLY_ID]=''.join(ALT_list)

#for key, value in DICT_POLAR_ALT.items():
#   print(key,len(DICT_POLAR_ALT[key]))






## Count differences between polarized and non-polarized polyploid sequence

tetra_RAW = in_POLY
tetra_ALT = DICT_POLAR_ALT[POLY_ID]

counter_ALT_RAW = 0

for i in range(len(tetra_RAW)) :		
   hap_RAW = tetra_RAW[i]
   hap_ALT = tetra_ALT[i]
   if hap_RAW != hap_ALT :
      counter_ALT_RAW += 1
      #print(hap_RAW,hap_ALT)



##### Save polarized MSA to file

#open output file 
straux1=str(sys.argv[1])
straux1 = straux1[:(len(straux1)-6)]
outfile_ALT = straux1 + '-ALT.fasta'
outfile_stats1 = straux1 + '.stats'

outfile2 = open(outfile_ALT, 'w')   
for key, value in DICT_POLAR_ALT.items() :
    outfile2.write(key)
    outfile2.write('\n')
    outfile2.write(value)
    outfile2.write('\n')

outfile2.flush()
outfile2.close()




 
      


