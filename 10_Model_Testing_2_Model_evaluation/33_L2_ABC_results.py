

# Luis Leal (2023)
### Script used to compute L2 distance for each ABC simulation run


######################################################### libraries (general)

import sys
import re                                                 
import ast                                                
import os                                                 
import time
import datetime
from collections import Counter
from pdb import set_trace as bp                           
import random						  
import itertools
import numpy as np 
import subprocess



######################################################### Auxiliary functions: L2 norm distance    
     
def L2norm_function(obs, sim, w) :
   ### obs: observed pairing profiles, all polarization geometries
   ### sim: simulated pairing profiles, all polarization geometries
   ### w: L2 weights
    
   L2 = np.sqrt(np.sum(np.square(np.absolute(np.multiply(np.subtract(obs,sim),w)))))
    
   return(L2)    
    
    


######################################################### Initialize

################# print options (single line)
np.set_printoptions(linewidth=np.inf)

################# Remember initial path
SRCDIR_INI = os.path.abspath(os.getcwd())

################# Read input and out-folder paths
try:
    AA = sys.argv[1]                  				    # input folder
except:
    print('\n Error: provide the following arguments:')
    print('33_L2_ABC_results.py [input folder] [output folder] [Observed_pairing_patterns.txt] [ABC priors folder] [ABC results file] \n')
    exit() 
    
try:
    RR_MAIN = sys.argv[2]                  				# output folder (MAIN)
except:
    print('\n Error: provide the following arguments:')
    print('33_L2_ABC_results.py [input folder] [output folder] [Observed_pairing_patterns.txt] [ABC priors folder] [ABC results file] \n')
    exit() 




################# Read observed pairing profiles 

PP_obs = list()

try:
    fhand = open(sys.argv[3], 'r')                  			
except:
    print('\n Error: input file missing.')
    print('33_L2_ABC_results.py [input folder] [output folder] [Observed_pairing_patterns.txt] [ABC priors folder] [ABC results file] \n')
    exit() 

for line in fhand:
    x = line.split()
    x = [float(i) for i in x]
    PP_obs.append(x)

PP_obs = np.array(PP_obs)
#print(PP_obs)




################# Read priors (1000 prior sets)

try:
   PP = sys.argv[4]  
   #
   ILSreduction_file = PP + "/ILSreduction.txt"    
   fhand_ILSreduction = open(ILSreduction_file, 'r') 
   #
   platy_pend_gene_flow_file = PP + "/platy_pend_gene_flow.txt"            				    
   fhand_platy_pend_gene_flow = open(platy_pend_gene_flow_file, 'r') 
   #
   H0_file = PP + "/H0.txt"            				    
   fhand_H0 = open(H0_file, 'r') 	
   #
   H1_file = PP + "/H1.txt"            				    
   fhand_H1 = open(H1_file, 'r')
   #
   H2_file = PP + "/H2.txt"            				    
   fhand_H2 = open(H2_file, 'r')
   #
   H3_file = PP + "/H3.txt"            				    
   fhand_H3 = open(H3_file, 'r')
   #
   H4_file = PP + "/H4.txt"            				    
   fhand_H4 = open(H4_file, 'r')
   #
   H5_file = PP + "/H5.txt"            				    
   fhand_H5 = open(H5_file, 'r')
   #
   H6_file = PP + "/H6.txt"            				    
   fhand_H6 = open(H6_file, 'r')
   #
   H7_file = PP + "/H7.txt"            				    
   fhand_H7 = open(H7_file, 'r')
   #
   H8_file = PP + "/H8.txt"            				    
   fhand_H8 = open(H8_file, 'r')
   #
   H9_file = PP + "/H9.txt"            				    
   fhand_H9 = open(H9_file, 'r')	
except:
   print('\n Error: Cannot find prior files. \n')
   exit() 


ILSreduction = list()
platy_pend_gene_flow = list()
H0 = list()
H1 = list()
H2 = list()
H3 = list()
H4 = list()
H5 = list()
H6 = list()
H7 = list()
H8 = list()
H9 = list()

for line in fhand_ILSreduction :
   x = line.split()
   ILSreduction.append(int(x[0]))

for line in fhand_platy_pend_gene_flow :
   x = line.split()
   platy_pend_gene_flow.append(float(x[0]))

for line in fhand_H0 :
   x = line.split()
   H0.append(float(x[0]))

for line in fhand_H1 :
   x = line.split()
   H1.append(float(x[0]))

for line in fhand_H2 :
   x = line.split()
   H2.append(float(x[0]))

for line in fhand_H3 :
   x = line.split()
   H3.append(float(x[0]))

for line in fhand_H4 :
   x = line.split()
   H4.append(float(x[0]))

for line in fhand_H5 :
   x = line.split()
   H5.append(float(x[0]))

for line in fhand_H6 :
   x = line.split()
   H6.append(float(x[0]))

for line in fhand_H7 :
   x = line.split()
   H7.append(float(x[0]))

for line in fhand_H8 :
   x = line.split()
   H8.append(float(x[0]))

for line in fhand_H9 :
   x = line.split()
   H9.append(float(x[0]))


  

################# Read simulated pairing patterns (ABC results, 1000 runs)

try:
    ABC_results_file = sys.argv[5]
    ABC_results_path = AA + '/' + ABC_results_file
    fhand_ABC = open(ABC_results_path, 'r')                  			
except:
    print('\n Error: input file missing.')
    print('33_L2_ABC_results.py [input folder] [output folder] [Observed_pairing_patterns.txt] [ABC priors folder] [ABC results file] \n')
    exit() 


ABC_results = list()
for line in fhand_ABC:
    xx = re.findall('^ABC', line)
    if len(xx) > 0 : 
        # header
        ABC_header = line
        #
    else :
       x = line.split()
       x = [float(i) for i in x]
       ABC_results.append(x)

ABC_results = np.array(ABC_results)




################# Compute L2 distance for each prior set (compare observed and estimated pairing patterns)

## weights used to compute L2 (w_max value applies to pendula, platyphylla, nana and humilis peaks; w_min value applied to all other peaks)
w_max = 5
w_min = 1
w_all = np.array(4*[w_max,w_max,w_max,w_min,w_min,w_max,w_min,w_min,w_min,w_max,w_min,w_min,w_min,w_min,w_min])

PPobs = np.concatenate(PP_obs)
L2dist = list()

for i in range(1000):
   PPsim = ABC_results[i][1:]
   L2_ref = -1
   L2_ref = L2norm_function(PPobs, PPsim, w_all)
   L2dist.append(L2_ref)




################# Save results to file
os.chdir(RR_MAIN)
fileOUT = "L2-norm_ABC_1000_simulations.txt"
outfile2 = open(fileOUT, 'w')   

# header
myheader = "ABC_No" + " " + "ILSreduction" + " " + "platy_pend_gene_flow" + " " + "H0"  + " " + "H1"  + " " + "H2"  + " " + "H3"  + " " + "H4"  + " " + "H5"  + " " + "H6"  + " " + "H7"  + " " + "H8"  + " " + "H9" + " " + "pendula" + " " + "platyphylla" + " " + "pend_platy" + " " + "populifolia" + " " + "pop_out" + " " + "nana" + " " + "nana_out" + " " + "occidentalis" + " " + "occ_out" + " " + "humilis" + " " + "albosinensis" + " " + "lenta" + " " + "albo_len_hum" + " " + "albo_len" + " " + "basal" + " " + "pendula" + " " + "platyphylla" + " " + "pend_platy" + " " + "populifolia" + " " + "pop_out" + " " + "nana" + " " + "nana_out" + " " + "occidentalis" + " " + "occ_out" + " " + "humilis" + " " + "albosinensis" + " " + "lenta" + " " + "albo_len_hum" + " " + "albo_len" + " " + "basal" + " " + "pendula" + " " + "platyphylla" + " " + "pend_platy" + " " + "populifolia" + " " + "pop_out" + " " + "nana" + " " + "nana_out" + " " + "occidentalis" + " " + "occ_out" + " " + "humilis" + " " + "albosinensis" + " " + "lenta" + " " + "albo_len_hum" + " " + "albo_len" + " " + "basal" + " " + "pendula" + " " + "platyphylla" + " " + "pend_platy" + " " + "populifolia" + " " + "pop_out" + " " + "nana" + " " + "nana_out" + " " + "occidentalis" + " " + "occ_out" + " " + "humilis" + " " + "albosinensis" + " " + "lenta" + " " + "albo_len_hum" + " " + "albo_len" + " " + "basal" + " " + "L2"

outfile2.write(myheader)
outfile2.write('\n')


# observed
observed = "Observed" + " " + "-" + " " + "-" + " " + "-"  + " " + "-"  + " " + "-"  + " " + "-"  + " " + "-"  + " " + "-"  + " " + "-"  + " " + "-"  + " " + "-"  + " " + "-" + " " + str(PPobs) + " " + "0"
outfile2.write(observed)
outfile2.write('\n')

# ABC simulations
for i in range(1000):
   counter = i + 1
   aux1 = str(counter)  + " " + str(ILSreduction[i]) + " " + str(platy_pend_gene_flow[i]) + " " + str(H0[i])  + " " + str(H1[i])  + " " + str(H2[i])  + " " + str(H3[i])  + " " + str(H4[i])  + " " + str(H5[i])  + " " + str(H6[i])  + " " + str(H7[i])  + " " + str(H8[i]) + " " + str(H9[i]) + " " + str(ABC_results[i][1:]) + " " + str(L2dist[i])
   outfile2.write(aux1)
   outfile2.write('\n')


outfile2.close()
os.chdir(SRCDIR_INI)


#### END






