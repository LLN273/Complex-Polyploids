#  
# Luis Leal (2021)
#
# Written in Python 3
# UPPMAX: module load python3/3.6.0



### Script used to simulate gene flow between B. pendula and B. platyphylla


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
    print('21s_PLATY2PEND_SIMPLE_INTERSPECIFIC_HYBRIDIZATION.py [input.fasta] [No. species] [ancestral ID1] [ancestral ID2] [ancestral ID3] [ancestral ID4] [hybridization ID1] [hybridization flag] \n')
    exit() 


try:
    AncestralID_1 = sys.argv[3]                  				# ID first ancestral sequence  (tetraploid)
except:
    print('\n Error: provide the following arguments:')
    print('21s_PLATY2PEND_SIMPLE_INTERSPECIFIC_HYBRIDIZATION.py [input.fasta] [No. species] [ancestral ID1] [ancestral ID2] [ancestral ID3] [ancestral ID4] [hybridization ID1] [hybridization flag] \n')
    exit() 


try:
    AncestralID_2 = sys.argv[4]                  				# ID second ancestral sequence  (tetraploid)
except:
    print('\n Error: provide the following arguments:')
    print('21s_PLATY2PEND_SIMPLE_INTERSPECIFIC_HYBRIDIZATION.py [input.fasta] [No. species] [ancestral ID1] [ancestral ID2] [ancestral ID3] [ancestral ID4] [hybridization ID1] [hybridization flag] \n')
    exit() 


try:
    AncestralID_3 = sys.argv[5]                  				# ID third ancestral sequence  (tetraploid)
except:
    print('\n Error: provide the following arguments:')
    print('21s_PLATY2PEND_SIMPLE_INTERSPECIFIC_HYBRIDIZATION.py [input.fasta] [No. species] [ancestral ID1] [ancestral ID2] [ancestral ID3] [ancestral ID4] [hybridization ID1] [hybridization flag] \n')
    exit() 


try:
    AncestralID_4 = sys.argv[6]                  				# ID fourth ancestral sequence  (tetraploid)
except:
    print('\n Error: provide the following arguments:')
    print('21s_PLATY2PEND_SIMPLE_INTERSPECIFIC_HYBRIDIZATION.py [input.fasta] [No. species] [ancestral ID1] [ancestral ID2] [ancestral ID3] [ancestral ID4] [hybridization ID1] [hybridization flag] \n')
    exit() 


try:
    GF_2pend_species = sys.argv[7]                  				# hybridizing species ID
except:
    print('\n Error: provide the following arguments:')
    print('21s_PLATY2PEND_SIMPLE_INTERSPECIFIC_HYBRIDIZATION.py [input.fasta] [No. species] [ancestral ID1] [ancestral ID2] [ancestral ID3] [ancestral ID4] [hybridization ID1] [hybridization flag] \n')
    exit() 


try:
    GF_FLAG = sys.argv[8]                  				# gene flow flag (0: no gene-flow; 1: gene-flow)
except:
    print('\n Error: provide the following arguments:')
    print('21s_PLATY2PEND_SIMPLE_INTERSPECIFIC_HYBRIDIZATION.py [input.fasta] [No. species] [ancestral ID1] [ancestral ID2] [ancestral ID3] [ancestral ID4] [hybridization ID1] [hybridization flag] \n')
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
   print('\n Error: check number of sequences in fasta file\n')
   exit()

if len(species_list) != int(Nspecies) :
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

          




###################################################### Tetraploid sequences

####### Create tetraploid sequences
TETRAPLOID_SIM_1_0_0 = '>' + AncestralID_1 + '_0_0'
TETRAPLOID_SIM_1_0_1 = '>' + AncestralID_2 + '_0_1'
TETRAPLOID_SIM_2_0_0 = '>' + AncestralID_3 + '_0_0'
TETRAPLOID_SIM_2_0_1 = '>' + AncestralID_4 + '_0_1'

###### gene flow from B. platyphylla to B. pendula 
if GF_FLAG == "1" :
   if GF_2pend_species == "PLATYPHYLLA_SIM" :
      PLATYPHYLLA_SIM_0_0 = '>' + "PENDULA_SIM" + '_0_0'
      PLATYPHYLLA_SIM_0_1 = '>' + "PENDULA_SIM" + '_0_1'

  
#print(INfasta_dict[TETRAPLOID_SIM_1_0_0])
#



### tetraploid nucleotide sequences

INfasta_dict[">TETRAPLOID_SIM_1_0_0"] = INfasta_dict[TETRAPLOID_SIM_1_0_0]
INfasta_dict[">TETRAPLOID_SIM_2_0_0"] = INfasta_dict[TETRAPLOID_SIM_1_0_1]
INfasta_dict[">TETRAPLOID_SIM_3_0_0"] = INfasta_dict[TETRAPLOID_SIM_2_0_0]
INfasta_dict[">TETRAPLOID_SIM_4_0_0"] = INfasta_dict[TETRAPLOID_SIM_2_0_1]


### Update B. platyphylla nucleotide sequence
if GF_FLAG == "1" :
   INfasta_dict[">PLATYPHYLLA_SIM_0_0"] = INfasta_dict[PLATYPHYLLA_SIM_0_0]
   INfasta_dict[">PLATYPHYLLA_SIM_0_1"] = INfasta_dict[PLATYPHYLLA_SIM_0_1]


#for key, value in INfasta_dict.items():
#   print(key,value)

#print(len(INfasta_dict))



##### Save  MSA to file

#open output file 
straux1=str(sys.argv[1])
straux1 = straux1[:(len(straux1)-6)]
outfile_ALT = straux1 + '-ALT.fasta'

outfile2 = open(outfile_ALT, 'w')   
for key, value in INfasta_dict.items() :
    outfile2.write(key)
    outfile2.write('\n')
    outfile2.write(value)
    outfile2.write('\n')


outfile2.close()




 
      


