#  
# Luis Leal (2022)
#
# Written in Python 3


# Script generates ABC prior values (for 1,000 simulation runs).
# Priors are sampled from a truncated distribution centered on the best model parameters determined using simulated annealing.


### Usage:
### python3 30_ABC_generate_priors_MAIN.py [priors parameter file] [N.samples]
###
### [priors parameter file]: priors min and max value for each parameter 
### [N.samples]: sample size (1000)


### Parameter file example (in actual file, only the first line is commented with the '#' sign): 


# parameter	mean	sd	distribution	sameAs	
#ILSreduction	40	40	discrete_uniform	0
#platy_pend_gene_flow	0.0	0.0	continuous_uniform	0
#H0	0.0	0.0	continuous_uniform	0
#H1	0.10	0.03	normal_dist	0
#H2	0.20	0.03	normal_dist	0
#H3	0.29	0.03	normal_dist	0
#H4	0.0	0.0	continuous_uniform	0
#H5	0.04	0.025	normal_dist	0
#H6	0.13	0.03	normal_dist	0
#H7	0.0	0.0	continuous_uniform	0
#H8	0.0	0.0	continuous_uniform	0
#H9	0.66	0.05	normal_dist	0




### Notes:
### 1. Continuous uniform distributions will be forced to [0,1] range. 
### 2. The "ILSreduction" distribution can only accept the following values: {0, 5, 10, 15, 20, ..., 90, 95, 100}; required files must be prepared in advanced
### 3. 'sameAs' forces priors to be identical to those computed for a different variable
### 4. Columns must be tab spaces

### Output:
### ABC priors will be saved in "ABC_PRIORS" folder.







######################################################### LOAD STANDARD MODULES

import sys
import re                                                 
import ast                                                
import os                                                 
import time
from collections import Counter
from pdb import set_trace as bp                           
import random						  
import itertools
import numpy as np 
import scipy.stats as stats				  




######################################################## OPEN INPUT FILES

try:
    fhand = open(sys.argv[1], 'r')                  			
except:
    print('\n Error: Cannot find ABC_PARAMETERS.txt input file.\n')
    exit() 



######################################################### Read sample size

try:
    Nsamples = sys.argv[2]                  				# Number of species
except:
    Nsamples = 1000



######################################################## READ INPUT FILE; CHECK FILE INTEGRITY

counterALL = 0          # counts number of parameters in input file
inPAR = list()

for line in fhand:
    #
    # detect comment lines
    xx = re.findall('^#', line)
    #
    if len(xx) == 0 : 
        #
        x = line.split()
        inPAR.append(x)
        counterALL = counterALL + 1

### Check number of parameters 
aux_len = len(inPAR)
#print(aux_len)
if aux_len != 12:
    print('\n Error: Input parameters file not formatted properly (1).')
    print('        Provide the following parameters: ILSreduction, platy_pend_gene_flow, H0, H1, H2, H3, H4, H5, H6, H7, H8, H9.')
    print('        Format: parameter_name	min	max	distribution_type	sameAs.\n')
    exit() 

flat_inPAR = list(itertools.chain(*inPAR))
#print(flat_inPAR)
aux_len2 = len(flat_inPAR)
if aux_len2 != 60:
    print('\n Error: Input parameters file not formatted properly (2).')
    print('        Provide the following parameters: ILSreduction, platy_pend_gene_flow, H0, H1, H2, H3, H4, H5, H6, H7, H8, H9.')
    print('        Format: parameter_name	min	max	distribution_type	sameAs.\n')
    exit() 


### Check whether any parameter is missing 
CP_parname = list()
CP_min = list()
CP_max = list()
CP_dtype = list()
CP_sameAs = list()

try:
    for i in inPAR:
        CP_parname.append(i[0])
        CP_min.append(i[1])
        CP_max.append(i[2])
        CP_dtype.append(i[3])
        CP_sameAs.append(i[4])

except:
    print('\n Error: Input parameters file not formatted properly (3).')
    print('        Format: parameter_name	min	max	distribution_type   sameAs.\n')
    exit()
        

#print(CP_parname)
#print(CP_min)
#print(CP_max)
#print(CP_dtype)


parname_aux = sorted(set(['ILSreduction', 'platy_pend_gene_flow', 'H0', 'H1', 'H2', 'H3', 'H4', 'H5', 'H6', 'H7', 'H8', 'H9']))
CP_parname_aux = sorted(set(CP_parname))
#print(CP_parname_aux)

if parname_aux!=CP_parname_aux:
    print('\n Error: Input parameters file not formatted properly (4).')
    print('        Provide the following parameters: ILSreduction, platy_pend_gene_flow, H0, H1, H2, H3, H4, H5, H6, H7, H8, H9.')
    print('        Format: parameter_name	min	max	distribution_type   sameAs.\n')
    exit() 


### Check numerical parameters
try:
    CP_min_f = [float(i) for i in CP_min]
    CP_max_f = [float(i) for i in CP_max]
    #print(CP_min_f)
    #print(CP_max_f)
except:
    print('\n Error: Input parameters file not formatted properly (5).')
    print('        Minimum and maximum values must be numeric.\n')
    exit()



######################################################## GENERATE PRIORS

CP_min = [float(i) for i in CP_min]
CP_max = [float(i) for i in CP_max]


OUTfolder = "ABC_PRIORS"
try :
    os.mkdir(OUTfolder)
except:
    print()

os.chdir(OUTfolder)

PRIORS_dict = dict()

for i in range(len(CP_parname)):
    GO_parname = CP_parname[i][:]
    if GO_parname != "ILSreduction":
        #
        # Check min,max range (must be within [0,1])
        if CP_min[i] < 0 :
            GO_CP_min = 0.0
        else :
            GO_CP_min = CP_min[i]
        #
        if CP_max[i] > 1 :
            GO_CP_max = 1.0
        else :
            GO_CP_max = CP_max[i]
        #
        #print(GO_parname, GO_CP_min, GO_CP_max)
        #
        ### Generate priors 
        # continuous uniform distribution
        if CP_dtype[i] == "continuous_uniform" :
           if GO_CP_min > GO_CP_max : GO_CP_min = GO_CP_max
           GO_priors = np.random.uniform(GO_CP_min, GO_CP_max, size = int(Nsamples))
           #print(GO_priors)
        #
        # truncated normal distr.
        if CP_dtype[i] == "normal_dist" :
           mu, sigma = GO_CP_min, GO_CP_max    	# mean and standard deviation
           lower_l = mu - 2*sigma	       	# lower limit
           upper_l = mu + 2*sigma          	# upper limit
           if lower_l < 0 : lower_l = 0
           if upper_l > 1 : upper_l = 1
           print(i,GO_CP_min, GO_CP_max, mu, sigma, lower_l, upper_l)
           GO_priors = stats.truncnorm.rvs((lower_l - mu)/sigma,(upper_l - mu)/sigma, loc=mu, scale=sigma, size=int(Nsamples))
           #GO_priors = stats.truncnorm.rvs(lower_l, upper_l, loc=mu, scale=sigma, size=int(Nsamples))
           #GO_priors = np.random.normal(mu, sigma, int(Nsamples))       
        #
        # Check if priors must be identical to priors computed for a different variable 
        if CP_sameAs[i] != '0':
           GO_priors = PRIORS_dict[CP_sameAs[i]]
        #
	# save to dictionary
        PRIORS_dict[GO_parname] = GO_priors
        #
        # Save priors to file
        my_outfile = GO_parname + '.txt'
        outfile1 = open(my_outfile, 'w') 
        for n in GO_priors:
            outfile1.write(str(n))
            outfile1.write('\n')
        #
        outfile1.close()
    else :
        # Generate ILSreduction priors:
        # Check min,max range (must be within [0,100])
        if CP_min[i] < 0 :
            GO_CP_min = 0
        else :
            #GO_CP_min = int(round(CP_min[i]/10))
            GO_CP_min = CP_min[i]/10
        #
        if CP_max[i] > 100 :
            GO_CP_max = 10
        else :
            #GO_CP_max = int(round(CP_max[i]/10))
            GO_CP_max = CP_max[i]/10
        #
        if GO_CP_min > GO_CP_max : GO_CP_min = GO_CP_max
        #print(GO_CP_min,GO_CP_max)
        #
        # Generate priors (discrete uniform distribution)
        GO_CP_min_r = round(2*GO_CP_min)/2
        GO_CP_max_r = round(2*GO_CP_max)/2
        aux_set = list()
        aux_set.append(GO_CP_min_r)
        aux_set.append(GO_CP_max_r)
        FLAGstop = 0
        aux_bas = GO_CP_min_r
        while FLAGstop == 0:
            aux_bas += 0.5
            if aux_bas < GO_CP_max_r :
                aux_set.append(aux_bas)
            else:
                FLAGstop = 1
        #
        GO_priors2 = list()
        for j in range(int(Nsamples)):
            rd = int(10*random.choice(aux_set))
            GO_priors2.append(rd)
        #
        #print(GO_priors2)
        #
        # Save priors to file
        my_outfile = GO_parname + '.txt'
        outfile1 = open(my_outfile, 'w') 
        for n in GO_priors2:
            outfile1.write(str(n))
            outfile1.write('\n')
        #
        outfile1.close()
        


print('Done!\n')


