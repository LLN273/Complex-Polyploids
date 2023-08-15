#  
# Luis Leal (2023)
#
# Written in Python 3



### Script used to create EST-SFS input file
### All species assumed to be diploid
### Input:
### - Joint VCF file that includes samples associated to focal species
### - sample list for focal species 
### - VCF outgroup 1
### - VCF outgroup 2
### - VCF outgroup 3
### - outputfile


print('\n Parsing started ...')


######################################################### LOAD STANDARD MODULES

import sys
import re                                                 
import ast                                                
import os                                                 
import time
import numpy as np
import pandas as pd
from pdb import set_trace as bp                           
import csv 



######################################################## OPEN INPUT FILES


try:
    fhand = open(sys.argv[1], 'r')                  			# Joint VCF file that includes samples associated to focal species
except:
    print('\n Error: input file missing.')
    print('03s_EST_SFS_Create_InputFile.py [input.vcf focal] [sample.list focal] [input.vcf outgroup1] [input.vcf outgroup2] [input.vcf outgroup3] [output.file] \n')
    exit() 
    
    
try:
    shand1 = open(sys.argv[2], 'r')                  			# sample list for focal species 
except:
    print('\n Error: input file missing.')
    print('03s_EST_SFS_Create_InputFile.py [input.vcf focal] [sample.list focal] [input.vcf outgroup1] [input.vcf outgroup2] [input.vcf outgroup3] [output.file] \n')
    exit()     
    
    
try:
    fhand_og1 = open(sys.argv[3], 'r')                  			# VCF outgroup 1
except:
    print('\n Error: input file missing.')
    print('03s_EST_SFS_Create_InputFile.py [input.vcf focal] [sample.list focal] [input.vcf outgroup1] [input.vcf outgroup2] [input.vcf outgroup3] [output.file] \n')
    exit()      
    
    
try:
    fhand_og2 = open(sys.argv[4], 'r')                  			# VCF outgroup 2
except:
    print('\n Error: input file missing.')
    print('03s_EST_SFS_Create_InputFile.py [input.vcf focal] [sample.list focal] [input.vcf outgroup1] [input.vcf outgroup2] [input.vcf outgroup3] [output.file] \n')
    exit()   
    
    
try:
    fhand_og3 = open(sys.argv[5], 'r')                  			# VCF outgroup 3
except:
    print('\n Error: input file missing.')
    print('03s_EST_SFS_Create_InputFile.py [input.vcf focal] [sample.list focal] [input.vcf outgroup1] [input.vcf outgroup2] [input.vcf outgroup3] [output.file] \n')
    exit()     
    
       
try:
    outFile_name = sys.argv[6]                  			# output filename
except:
    print('\n Error: input file missing.')
    print('03s_EST_SFS_Create_InputFile.py [input.vcf focal] [sample.list focal] [input.vcf outgroup1] [input.vcf outgroup2] [input.vcf outgroup3] [output.file] \n')
    exit()   
    


######################################################## Read sample list

## focal species
SamplesList = list()    # store sample list 

for line in shand1:
    line = line.replace('\n','')
    SamplesList.append(line)
    

#print(SamplesList)
#bp()



### AUXILIARY FUNCTION: searches 'Fixed fields' string (VCF) and gets numeric value of a specific parameter

def getValue(rd_aux00, queryString):
    rm_mark0 = rd_aux00.find(queryString)
    if rm_mark0 != -1 :                                             #if parameter present in x[7]
        rd_aux0 = rd_aux00[(rm_mark0+len(queryString)):]
        rm_mark = rd_aux0.find(';')
        if rm_mark != -1 :                                          #when parameter is in the middle of x[7]
            try :
                outvalue = float(rd_aux0[:rm_mark])       		    # get parameter value
            except: outvalue = ''      								# parameter value is a string (e.g., 'NaN')
        else :
            P_aux1 = rd_aux00[(rm_mark0+len(queryString)):]         #when parameter is at the end of x[7]
            try :
                outvalue = float(P_aux1)       		                # get parameter value
            except: outvalue = ''      								# parameter value is a string (e.g., 'NaN')
    else :
        outvalue = ''
    return outvalue


### end function




### AUXILIARY FUNCTION: gets EST codes for outgroups

def ESTrecode(alleles, genotype):
    outcode = [0,0,0,0]
    if (genotype == "0/0" or genotype == "0|0"):
        if alleles[0] == 'A' : outcode[0] = 1
        if alleles[0] == 'C' : outcode[1] = 1
        if alleles[0] == 'G' : outcode[2] = 1
        if alleles[0] == 'T' : outcode[3] = 1
    elif (genotype == "0/1" or genotype == "0|1"):
        if alleles[0] == 'A' : outcode[0] = 1
        if alleles[0] == 'C' : outcode[1] = 1
        if alleles[0] == 'G' : outcode[2] = 1
        if alleles[0] == 'T' : outcode[3] = 1
        if alleles[1] == 'A' : outcode[0] = 1
        if alleles[1] == 'C' : outcode[1] = 1
        if alleles[1] == 'G' : outcode[2] = 1
        if alleles[1] == 'T' : outcode[3] = 1
    elif (genotype == "1/1" or genotype == "1|1"):
        if alleles[1] == 'A' : outcode[0] = 1
        if alleles[1] == 'C' : outcode[1] = 1
        if alleles[1] == 'G' : outcode[2] = 1
        if alleles[1] == 'T' : outcode[3] = 1
    elif (genotype == "1/2" or genotype == "1|2"):
        if alleles[1] == 'A' : outcode[0] = 1
        if alleles[1] == 'C' : outcode[1] = 1
        if alleles[1] == 'G' : outcode[2] = 1
        if alleles[1] == 'T' : outcode[3] = 1
        if alleles[2] == 'A' : outcode[0] = 1
        if alleles[2] == 'C' : outcode[1] = 1
        if alleles[2] == 'G' : outcode[2] = 1
        if alleles[2] == 'T' : outcode[3] = 1
    else :
        print('\n Error: The following genotype was found in one of the outgroup VCF files:')
        print(genotype)
        print('Please update function ESTrecode in order to process this genotype. \n')
        exit()

    return outcode

### end function




### AUXILIARY FUNCTION: gets EST codes for focal species

def ESTrecode_focal(alleles, genotype_list):
   outcode = [0,0,0,0]
   for GT in genotype_list :
       if GT == "0/0" or GT == "0|0": 
           if alleles[0] == 'A' : outcode[0] = outcode[0] + 1
           if alleles[0] == 'C' : outcode[1] = outcode[1] + 1
           if alleles[0] == 'G' : outcode[2] = outcode[2] + 1
           if alleles[0] == 'T' : outcode[3] = outcode[3] + 1
       elif GT == "0/1" or GT == "0|1":
           if alleles[0] == 'A' : outcode[0] = outcode[0] + 1
           if alleles[0] == 'C' : outcode[1] = outcode[1] + 1
           if alleles[0] == 'G' : outcode[2] = outcode[2] + 1
           if alleles[0] == 'T' : outcode[3] = outcode[3] + 1
           if alleles[1] == 'A' : outcode[0] = outcode[0] + 1
           if alleles[1] == 'C' : outcode[1] = outcode[1] + 1
           if alleles[1] == 'G' : outcode[2] = outcode[2] + 1
           if alleles[1] == 'T' : outcode[3] = outcode[3] + 1
       elif GT == "1/1" or GT == "1|1":
           if alleles[1] == 'A' : outcode[0] = outcode[0] + 1
           if alleles[1] == 'C' : outcode[1] = outcode[1] + 1
           if alleles[1] == 'G' : outcode[2] = outcode[2] + 1
           if alleles[1] == 'T' : outcode[3] = outcode[3] + 1
       elif GT == "1/2" or GT == "1|2":
           if alleles[1] == 'A' : outcode[0] = outcode[0] + 1
           if alleles[1] == 'C' : outcode[1] = outcode[1] + 1
           if alleles[1] == 'G' : outcode[2] = outcode[2] + 1
           if alleles[1] == 'T' : outcode[3] = outcode[3] + 1
           if alleles[2] == 'A' : outcode[0] = outcode[0] + 1
           if alleles[2] == 'C' : outcode[1] = outcode[1] + 1
           if alleles[2] == 'G' : outcode[2] = outcode[2] + 1
           if alleles[2] == 'T' : outcode[3] = outcode[3] + 1
       else :
           print('\n Error: The following genotype was found in the focal species VCF file:')
           print(GT)
           print('Please update function ESTrecode_focal in order to process this genotype. \n')
           exit()
    
   return outcode
### end function










######################################################## Get data for outgroups
######################################################## Read genotype for each locus and recode
######################################################## code as A allele as 1,0,0,0
######################################################## code as C allele as 0,1,0,0
######################################################## code as G allele as 0,0,1,0
######################################################## code as T allele as 0,0,0,1
######################################################## code missing as 0,0,0,0




############### Outgroup 1
print('\n Processing OG1 ...')   

OG1_dict = dict()
counterOG1 = 0
 
for line in fhand_og1:
    #    
    # detect comment lines
    xx = re.findall('^#', line)
    #
    if len(xx) == 0 :   # not a comment line
        # variants
        # counter, all sites 
        counterOG1 = counterOG1 + 1
        #
        x = line.split()
        #
        CONTIG = x[0]
        POS = x[1]
        # gather and recode alleles   
        REF = x[3]
        ALT = x[4]
        ALT = ALT.split(',')     
        ALLELES = [REF,ALT]
        ALLELES = [item for sublist in ALLELES for item in sublist]
        if '.' in ALLELES: ALLELES.remove(".")
        # genotype
        GTfield = x[9]
        GTx = GTfield.split(':')[0]     
        # EST code
        ESTcode = ESTrecode(ALLELES, GTx)
        # save to dictionary
        genome_position = '_'.join([CONTIG, POS]) 
        OG1_dict[genome_position] = ESTcode
        #print(CONTIG,POS,REF,ALT,ALLELES,GTx,ESTcode, gen_position)
        #if counterOG1 > 40 :
        #   print(OG1_dict)
        #   bp()
 
  
        
############### Outgroup 2
print('\n Processing OG2 ...')   

OG2_dict = dict()
counterOG2 = 0

for line in fhand_og2:
    #    
    # detect comment lines
    xx = re.findall('^#', line)
    #
    if len(xx) == 0 :   # not a comment line
        # variants
        # counter, all sites 
        counterOG2 = counterOG2 + 1
        #
        x = line.split()
        #
        CONTIG = x[0]
        POS = x[1]
        # gather and recode alleles   
        REF = x[3]
        ALT = x[4]
        ALT = ALT.split(',')     
        ALLELES = [REF,ALT]
        ALLELES = [item for sublist in ALLELES for item in sublist]
        if '.' in ALLELES: ALLELES.remove(".")
        # genotype
        GTfield = x[9]
        GTx = GTfield.split(':')[0]     
        # EST code
        ESTcode = ESTrecode(ALLELES, GTx)
        # save to dictionary
        genome_position = '_'.join([CONTIG, POS]) 
        OG2_dict[genome_position] = ESTcode
        
        
        
############### Outgroup 3
print('\n Processing OG3 ...')   

OG3_dict = dict()
counterOG3 = 0

for line in fhand_og3:
    #    
    # detect comment lines
    xx = re.findall('^#', line)
    #
    if len(xx) == 0 :   # not a comment line
        # variants
        # counter, all sites 
        counterOG3 = counterOG3 + 1
        #
        x = line.split()
        #
        CONTIG = x[0]
        POS = x[1]
        # gather and recode alleles   
        REF = x[3]
        ALT = x[4]
        ALT = ALT.split(',')     
        ALLELES = [REF,ALT]
        ALLELES = [item for sublist in ALLELES for item in sublist]
        if '.' in ALLELES: ALLELES.remove(".")
        # genotype
        GTfield = x[9]
        GTx = GTfield.split(':')[0]     
        # EST code
        ESTcode = ESTrecode(ALLELES, GTx)
        # save to dictionary
        genome_position = '_'.join([CONTIG, POS]) 
        OG3_dict[genome_position] = ESTcode

#print(SamplesList)
#bp()






######################################################## Get data for focal species
print('\n Processing focal species ...')   
focalSpecies_dict = dict()
counterALL = 0          # counts all sites (variant and non-variant)
indexSample1 = list()   # index of samples, as shown in VCF file
EST_list = list()

for line in fhand:
    #    
    # detect comment lines
    xx = re.findall('^#', line)
    #
    if len(xx) > 0 : 
        # comment lines
        flagS = line.startswith("#CHROM")   # look for line in VCF file containing sample names
        if flagS == True :
            x = line.split()
            #
            #indices = 0, 1, 2, 3, 4, 5, 6, 7, 8
            #x = [i for j, i in enumerate(x) if j not in indices]     # remove '[#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT' ] entries from list
            #
            # Find index of samples of interest
            for sample in SamplesList:
                myindex = x.index(sample)
                indexSample1.append(myindex)
            #print(indexSample1)
        #
    else :
        # variants
        # counter, all sites 
        counterALL = counterALL + 1
        #
        x = line.split()
        #
        CONTIG = x[0]
        POS = x[1]
        # gather and recode alleles   
        REF = x[3]
        ALT = x[4]
        ALT = ALT.split(',')     
        ALLELES = [REF,ALT]
        ALLELES = [item for sublist in ALLELES for item in sublist]
        if '.' in ALLELES: ALLELES.remove(".")
        # get genotypes
        GT1 = [x[y] for y in indexSample1]     
        GT1 = [(r.split(':')[0]) for r in GT1]
        # EST code
        ESTcode = ESTrecode_focal(ALLELES, GT1)
        #print(CONTIG,POS,REF,ALT,ALLELES,GT1,ESTcode)
        # Get outgroup EST codes for same genomic position
        genome_position = '_'.join([CONTIG, POS])
        try:
            og1_code = OG1_dict[genome_position]                  			
        except:
            og1_code = [0,0,0,0]
        try:
            og2_code = OG2_dict[genome_position]                  			
        except:
            og2_code = [0,0,0,0]
        try:
            og3_code = OG3_dict[genome_position]                  			
        except:
            og3_code = [0,0,0,0]
        aux_ESTlist = [ESTcode, og1_code, og2_code, og3_code]
        EST_list.append(aux_ESTlist)
        #if counterALL > 20 :
        #   for n in EST_list:
        #      print(n)
        #   bp()






OUTfile = open(outFile_name, 'w', newline ='') 
with OUTfile:     
    write = csv.writer(OUTfile, delimiter='\t') 
    write.writerows(EST_list) 


OUTfile.close()

print('\n Done!')


