#  
# Luis Leal (2023)
#
# Written in Python 3


### Script used to build the unfolded 2D site frequency spectrum (SFS) for two populations
### Input:
### - Joint VCF file
### - Ancestral allele
### - sample list for species 1
### - sample list for species 2
### - ploidy pop 1
### - ploidy pop 2
### - name pop1 (integer)
### - name pop2 (integer)
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



######################################################## OPEN INPUT FILES

try:
    fhand = open(sys.argv[1], 'r')                  			# vcf file name (Joint VCF)
except:
    print('\n Error: input file missing.')
    print('07s_2D_SFS_Creator.py [input.vcf] [ancestral.info] [pop1.list] [pop2.list] [ploidy pop1] [ploidy pop2] [pop1 label] [pop2 label] [output file] \n')
    exit() 
    
    
try:
    fhand_a = open(sys.argv[2], 'r')                  			# ancestral allele info (EST-FSF output)
except:
    print('\n Error: input file missing.')
    print('07s_2D_SFS_Creator.py [input.vcf] [ancestral.info] [pop1.list] [pop2.list] [ploidy pop1] [ploidy pop2] [pop1 label] [pop2 label] [output file] \n')
    exit() 
    
    
try:
    shand1 = open(sys.argv[3], 'r')                  			# list of samples in pop1
except:
    print('\n Error: input file missing.')
    print('07s_2D_SFS_Creator.py [input.vcf] [ancestral.info] [pop1.list] [pop2.list] [ploidy pop1] [ploidy pop2] [pop1 label] [pop2 label] [output file] \n')
    exit()     
    

try:
    shand2 = open(sys.argv[4], 'r')                  			# list of samples in pop2
except:
    print('\n Error: input file missing.')
    print('07s_2D_SFS_Creator.py [input.vcf] [ancestral.info] [pop1.list] [pop2.list] [ploidy pop1] [ploidy pop2] [pop1 label] [pop2 label] [output file] \n')
    exit()        
      
    
try:
    ploidy1 = sys.argv[5]                  			# ploidy individuals in pop1
except:
    print('\n Error: input file missing.')
    print('07s_2D_SFS_Creator.py [input.vcf] [ancestral.info] [pop1.list] [pop2.list] [ploidy pop1] [ploidy pop2] [pop1 label] [pop2 label] [output file] \n')
    exit()  
    
    
try:
    ploidy2 = sys.argv[6]                  			# ploidy individuals in pop2
except:
    print('\n Error: input file missing.')
    print('07s_2D_SFS_Creator.py [input.vcf] [ancestral.info] [pop1.list] [pop2.list] [ploidy pop1] [ploidy pop2] [pop1 label] [pop2 label] [output file] \n')
    exit()        


try:
    popName1 = sys.argv[7]                  			# label for pop1
except:
    print('\n Error: input file missing.')
    print('07s_2D_SFS_Creator.py [input.vcf] [ancestral.info] [pop1.list] [pop2.list] [ploidy pop1] [ploidy pop2] [pop1 label] [pop2 label] [output file] \n')
    exit()    
    

try:
    popName2 = sys.argv[8]                  			# label for pop2
except:
    print('\n Error: input file missing.')
    print('07s_2D_SFS_Creator.py [input.vcf] [ancestral.info] [pop1.list] [pop2.list] [ploidy pop1] [ploidy pop2] [pop1 label] [pop2 label] [output file] \n')
    exit()      
    

try:
    outFile = sys.argv[9]                  			# output filename
except:
    print('\n Error: input file missing.')
    print('07s_2D_SFS_Creator.py [input.vcf] [ancestral.info] [pop1.list] [pop2.list] [ploidy pop1] [ploidy pop2] [pop1 label] [pop2 label] [output file] \n')
    exit()   







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




######################################################## Read sample lists


## pop1

SamplesPop1List = list()    # store sample list for pop1

for line in shand1:
    line = line.replace('\n','')
    SamplesPop1List.append(line)
    
#print(SamplesPop1List)
#bp()
    
## pop2

SamplesPop2List = list()    # store sample list for pop2

for line in shand2:
    line = line.replace('\n','')
    SamplesPop2List.append(line)    

#print(SamplesPop2List)
#bp()



######################################################## Ancestral allele


ancestralAllele = list()    # store ancestral alleles

for line in fhand_a:
    #
    # detect comment lines
    xx = re.findall('^0', line)
    #
    if len(xx) == 0 :   # not a comment line
       x = line.split(" ")
       aux_allele = x[3:19]                     # probabilities of the first and second internal nodes having states [A,A], [A,C], [A,G], [A,T], [C,A], [C,C], [C,G], [C,T], etc.
       #print(aux_allele)
       max_value = max(aux_allele)              # find highest probability
       max_index = aux_allele.index(max_value)  # position for which max probability was observed
       if (max_index >= 0 and max_index <= 4): ancestral = "A" 
       if (max_index >= 5 and max_index <= 8): ancestral = "C"
       if (max_index >= 9 and max_index <= 12): ancestral = "G"
       if (max_index >= 13 and max_index <= 16): ancestral = "T"
       ancestralAllele.append(ancestral)
       #print(ancestral)

#print(len(ancestralAllele))
#bp()


######################################################## Open VCF file and check number of alternate alleles for each sample

counterALL = 0          # main counter
indexSample1 = list()   # index of pop1 samples, as shown in VCF file
indexSample2 = list()   # index of pop2 samples, as shown in VCF file
SFS = np.zeros([(len(SamplesPop1List)*int(ploidy1)+1), (len(SamplesPop2List)*int(ploidy2)+1)], dtype=int)
#aux = SFS.shape
#print(aux)

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
            # pop1 samples
            for sample in SamplesPop1List:
                myindex = x.index(sample)
                indexSample1.append(myindex)
            #print(indexSample1)
            # pop2 samples
            for sample in SamplesPop2List:
                myindex = x.index(sample)
                indexSample2.append(myindex)
            #print(indexSample2)
        #
    else :
        # variants
        # counter, all sites 
        counterALL = counterALL + 1
        #
        x = line.split()
        #
        # get minor allele frequency for this site
        #INFOf = x[7]                	        # INFO field
        #AF = float(getValue(INFOf, 'AF='))
        #
        # get genotypes
        GT1 = [x[y] for y in indexSample1]     # pop1
        GT1 = [(r.split(':')[0]) for r in GT1]
        GT2 = [x[y] for y in indexSample2]     # pop2
        GT2 = [(r.split(':')[0]) for r in GT2]
        #print(GT1)
        #print(GT2)
        #
        # Get reference and alternate alleles
        REF = x[3]
        ALT = x[4]
        # Get ancestral allele
        AA = ancestralAllele[counterALL-1]
        #print(REF, ALT, AA)
        #
        # Recode REF/ALL alleles if reference allele does not match ancestral allele 
        # (note: cases where ancestral allele is not one of REF or ALT alleles are ignored)
        if (AA == REF or AA == ALT) : 
            if AA == ALT :
                GT1 = [sub.replace('0', 'A').replace('1', '0').replace('A','1') for sub in GT1]
                GT2 = [sub.replace('0', 'A').replace('1', '0').replace('A','1') for sub in GT2]
            # for each population, count number of times alternative allele is observed
            if ('./.' not in GT1) and ('./.' not in GT2) and ('./././.' not in GT1) and ('./././.' not in GT2) and ('.|.' not in GT1) and ('.|.' not in GT2) and ('.|.|.|.' not in GT1) and ('.|.|.|.' not in GT2):        # no missing values
                GT1_hold = GT1
                GT2_hold = GT2
                GT1 = list(''.join(GT1))
                GT1 = list(filter(lambda a: a != '/', GT1))
                GT1 = list(filter(lambda a: a != '|', GT1))
                GT2 = list(''.join(GT2))
                GT2 = list(filter(lambda a: a != '/', GT2))
                GT2 = list(filter(lambda a: a != '|', GT2))
                GT1 = [int(i) for i in GT1]
                GT2 = [int(i) for i in GT2]
                ALT_count_1 = sum(GT1)
                ALT_count_2 = sum(GT2)
                #print(GT1_hold)
                #print(GT2_hold)
                #print(ALT_count_1, ALT_count_2)
                #if counterALL > 20 :
                #    bp()
                #
                # Update 2D-SFS
                SFS[ALT_count_1][ALT_count_2] = SFS[ALT_count_1][ALT_count_2] + 1
      



#print(SFS)
#bp()

# Check VCF length matches total number of ancestral alleles
if counterALL != len(ancestralAllele) :
    print('\n Error: the total number of ancestral alleles does not match VCF length.')
    print('Number of ancestral alleles:', len(ancestralAllele))
    print('VCF length:', counterALL ,'\n')
    exit()   

# convert unfolded SFS to dataframe
SFS_DF = pd.DataFrame(SFS)

# dataframe row names (pop1)
rnList = list()
for i in range(len(SamplesPop1List)*int(ploidy1)+1) :
    aux_rn = 'd' + str(popName1) + '_' + str(i)
    rnList.append(aux_rn)

SFS_DF.index = rnList


# dataframe column names (pop2)
cnList = list()
for i in range(len(SamplesPop2List)*int(ploidy2)+1) :
    aux_cn = 'd' + str(popName2) + '_' + str(i)
    cnList.append(aux_cn)

SFS_DF.columns = cnList


# transpose df
SFS_DF = SFS_DF.transpose()



#open output file 
outfile1_name = outFile
outfile1 = open(outfile1_name, 'w')                                

# save results to file
aux_header = '1 observation'
SFS_DF.to_csv(outfile1, sep ='\t')

outfile1.close()




### Allele count per population
### format: numpy.delete(arr, obj, axis=None)
AL_pop1 = np.delete(SFS, 0, 0)      # delete first (obj=0) row (axis=0) 
AL_pop2 = np.delete(SFS, 0, 1)      # delete first (obj=0) column (axis=1)  
AL_SFS = np.sum(SFS)
AL_pop1 = np.sum(AL_pop1)
AL_pop2 = np.sum(AL_pop2)
print("Allele count Pop", popName1, ":", AL_pop1)
print("Allele count Pop", popName2, ":", AL_pop2)
print("Allele count SFS:", AL_SFS)

print('\n Done!')


