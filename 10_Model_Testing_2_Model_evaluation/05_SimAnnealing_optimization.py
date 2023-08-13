
# Luis Leal (2023)
# Optimize polyploidization and hybridization parameters using simulated annealing


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






######################################################### Simulated annealing function


def simann_control(myPRIORS, L2_ref, HomoeologousExchange, w_all, geneFlowFlag) :

    ### Performs simulated annealing

    # initialize variables
    temperature = 5E3
    REFERENCE_TEMP = 1E-3
    COOLING_FACTOR = 0.999  
    PRIORS_best = np.copy(myPRIORS)
    L2_best = L2_ref
    step_size = 0.1                     # maximum step size during optimization
    step_counter = 1
    FlagNotBest = 0
    
    while temperature > REFERENCE_TEMP :

        step_counter = step_counter + 1
        
        ####### adjust step size used during optimization (decreases as temperature decreases)
        if (L2_ref < 0.55 and FlagNotBest < 50) :
           step_size = 0.05
        elif L2_ref > 2 :
           step_size = 0.20 
        elif L2_ref > 1 :
           step_size = 0.15
        else :     
           step_size = 0.1
        
        ####### generate new priors
        PRIORS_new = UpdatePriors_function(polymodel, myPRIORS, step_size, HomoeologousExchange, geneFlowFlag)            
                    
        #print("np.sum(PRIORS_new[2:]):", np.sum(PRIORS_new[2:]))
        #print("PRIORS_new_2:", PRIORS_new)
        #aux = np.subtract(PRIORS_new, myPRIORS)
        #print("difference:", aux.astype(float))
        #bp()

        ####### Save priors to file
        RRpp = RR_MAIN
        os.chdir(RRpp)
        PPnew = "00_PRIORS_ABC_SimAnnealing.txt"
        outfile1 = open(PPnew, 'w')   
        for k in PRIORS_new :
           outfile1.write(str(k))
           outfile1.write('\n')
        #
        outfile1.close()
        os.chdir(SRCDIR_INI)
    

        
        ################# Get simulated pairing patterns based on new prior set
        PPnew_path = RR_MAIN + "/" + PPnew
        cmd = ['./20_POLARIZATION_IQTREE_PAIRINGprofiles_MAIN.sh' + ' ' +  AA + ' ' + RR_MAIN + ' ' + polymodel + ' ' + str(REP) + ' ' + str(NLOCUS) + ' ' + PPnew_path]
        process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        res = process.communicate()
        
        process.wait()


        ################# Read simulated pairing patterns

        RRpp = RR_MAIN + "/03_PAIRINGfrequencies"
        os.chdir(RRpp)
        try:
           fhand_pp = open("sister_ID_analysis_priors-JOINT_summary.txt", 'r')                  			
        except:
           print('\n Error: simulated pairing profiles missing.')
           exit()

        os.chdir(SRCDIR_INI)


        PP_sim = list()
        for line in fhand_pp:
            x = line.split()
            x = [float(i) for i in x]
            PP_sim.append(x)

        PP_sim = np.array(PP_sim)
        PPsim = np.concatenate(PP_sim)
        
        
        
        ################# Compute new L2 distance (compare observed and estimated pairing patterns)
        L2_new = -1
        L2_new = L2norm_function(PPobs, PPsim, w_all)
          

        #compare distances
        dist_change = L2_new - L2_ref
        
        #print("L2_new:", L2_new)
        #print("dist_change:", dist_change)
        #bp()
        
        print("\n")
        print("Iteration:", step_counter, "; Temperature:", temperature)
        print("PRIORS:", PRIORS_new)
        print("Sub-routine process returncode (0 = no errors):", process.returncode)   # 0 return code indicates there were no errors
        print("L2 =", L2_new)

        #cost assessment: accept solution if better than previous, or probabilistically using the metropolis acceptance criterion
        # As temperature goes down, poor solutions have a lower probability of being accepted.
        if dist_change < 0 or (np.exp(-25000*dist_change/temperature) > np.random.uniform(0,1) and FlagNotBest > 50) :
            myPRIORS = np.copy(PRIORS_new)
            L2_ref = L2_new
            FlagNotBest = FlagNotBest + 1
            if L2_new < L2_best :
                PRIORS_best = np.copy(PRIORS_new)
                L2_best = L2_new
                FlagNotBest = 0
                print(">> NEW OPTIMAL SOLUTION.")
                #
                # Store results folder
                cmd = ['rm' + ' ' + '-rf' + ' ' + RR_MAIN + '/03_PAIRINGfrequencies_BEST']
                process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                res = process.communicate()
                #
                cmd = ['mv' + ' ' + RR_MAIN + '/03_PAIRINGfrequencies' + ' ' + RR_MAIN + '/03_PAIRINGfrequencies_BEST']
                process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                res = process.communicate()
                #
                cmd = ['mv' + ' ' + RR_MAIN + '/' + PPnew + ' ' + RR_MAIN + '/03_PAIRINGfrequencies_BEST' + '/' + PPnew]
                process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                res = process.communicate()

        # Update log file
        os.chdir(RR_MAIN)
        fileOUT = "ABC_SimAnnealing_log.txt"
        outfile2 = open(fileOUT, 'a')   
        now = datetime.datetime.now()
        out1 = str(now.strftime("%Y-%m-%d %H:%M:%S")) + '\t' + str(step_counter) + '\t' + str(PRIORS_new) + '\t' + str(L2_new)
        outfile2.write(out1)
        outfile2.write('\n')
        out1 = str(now.strftime("%Y-%m-%d %H:%M:%S")) + '\t' + str(step_counter) + '\t' + str(PRIORS_best) + '\t' + str(L2_best) + '\t' + 'Best.' 
        outfile2.write(out1)
        outfile2.write('\n')
        outfile2.close()
        os.chdir(SRCDIR_INI)        
        
        # decrease temperature
        temperature = temperature*COOLING_FACTOR
  
        now = datetime.datetime.now()
        print(now.strftime("%Y-%m-%d %H:%M:%S"))

    return PRIORS_best, L2_best, step_counter, RRpp
    
    
    
    
    
    
######################################################### Auxiliary functions: L2 norm distance    
    
    
def L2norm_function(obs, sim, w) :
   ### obs: observed pairing profiles, all polarization geometries
   ### sim: simulated pairing profiles, all polarization geometries
   ### w: L2 weights
    
   L2 = np.sqrt(np.sum(np.square(np.absolute(np.multiply(np.subtract(obs,sim),w)))))
    
   return(L2)    
    
    

def GenerateRandomPriors_function(polymodel, ILS_PRIOR_init, HomoeologousExchange, geneFlowFlag) :
   # generate random set of priors; number of model parameters varies with polyploidization/hybridization model
   #    
   myPRIORS = np.zeros(12)
   myPRIORS[0] = ILS_PRIOR_init                            # ILS prior
   if geneFlowFlag == "1" :
      myPRIORS[1] = np.random.rand(1,1)                       # platy-pendula gene flow
   else :
      myPRIORS[1] = 0
   #
   if polymodel == "PPPP" :			## PPPP model: H4 = H6 = H7 = H8 = 0; sum(H0 to H8) <= 1
      flagP = 0
      while flagP == 0:
         auxP = np.random.rand(1,4)
         if np.sum(auxP) <= 0.5:        # for the initial priors (only), keep total levels of hybridization between 0 and 50% (this value is allowed to change freely during subsequent iterations, from 0 to 100%).
            flagP = 1
      # 
      #myPRIORS[2] = auxP[0][0]   
      myPRIORS[3] = auxP[0][0]  
      myPRIORS[4] = auxP[0][1]
      myPRIORS[5] = auxP[0][2]
      myPRIORS[7] = auxP[0][3]
      #
   elif polymodel == "AAAA" :			## AAAA model: H4 = H6 = H7 = H8 = 0; sum(H0 to H8) <= 1
      myPRIORS[11] = np.random.rand(1,1)                       # introgression from B. pendula
      flagP = 0
      while flagP == 0:
         auxP = np.random.rand(1,5)
         if np.sum(auxP) <= 0.5:        # for the intial priors (only), keep total levels of hybridization between 0 and 50% (this value is allowed to change freely during subsequent iterations, from 0 to 100%).
            flagP = 1
      # 
      #myPRIORS[2] = auxP[0][0]   
      myPRIORS[3] = auxP[0][0]  
      myPRIORS[4] = auxP[0][1]
      myPRIORS[5] = auxP[0][2]
      myPRIORS[7] = auxP[0][3]  
      myPRIORS[8] = auxP[0][4] 
      #
   elif polymodel == "PPNN" or polymodel == "PPHH" :			## PPNN or PPHH models
      flagP = 0
      while flagP == 0:
         if HomoeologousExchange == "0" :      # no  homoeologous exchange allowed
            auxP = np.random.rand(1,3)
         elif HomoeologousExchange == "1" :   # homoeologous exchange allowed
            auxP = np.random.rand(1,5)
         #
         if np.sum(auxP) <= 0.5:        # for the initial priors (only), keep total levels of hybridization between 0 and 50% (this value is allowed to change freely during subsequent iterations, from 0 to 100%).
            flagP = 1
         # 
         #myPRIORS[2] = auxP[0][0]   
         myPRIORS[3] = auxP[0][0]  
         myPRIORS[4] = auxP[0][1]
         myPRIORS[5] = auxP[0][1]   # H2=H3
         myPRIORS[6] = auxP[0][2]
         if HomoeologousExchange == "1" :
             myPRIORS[7] = auxP[0][3]
             myPRIORS[8] = auxP[0][4]
   elif polymodel == "PPPyPy" :			## PPPyPy model
      flagP = 0
      while flagP == 0:
         if HomoeologousExchange == "0" :      # no  homoeologous exchange allowed
            auxP = np.random.rand(1,7)
         elif HomoeologousExchange ==  "1" :   # homoeologous exchange allowed
            auxP = np.random.rand(1,9)
         #
         if np.sum(auxP) <= 1.0:        # keep at 1.0 as condition is difficult to satisfy for 0.5 limit when there are 9 independent variables
            flagP = 1
         # 
         #myPRIORS[2] = auxP[0][0]   
         myPRIORS[2] = auxP[0][0]
         myPRIORS[3] = auxP[0][1]  
         myPRIORS[4] = auxP[0][2]
         myPRIORS[5] = auxP[0][3]
         myPRIORS[6] = auxP[0][4]
         myPRIORS[8] = auxP[0][5]
         myPRIORS[9] = auxP[0][6]
         if HomoeologousExchange == "1" :
             myPRIORS[7] = auxP[0][7]
             myPRIORS[10] = auxP[0][8]
   elif polymodel == "PPNH" :			## PPNH model
      flagP = 0
      while flagP == 0:
         if HomoeologousExchange == "0" :      # no  homoeologous exchange allowed
            auxP = np.random.rand(1,3)
         elif HomoeologousExchange ==  "1" :   # homoeologous exchange allowed
            auxP = np.random.rand(1,6)
         #
         if np.sum(auxP) <= 0.5:        # for the initial priors (only), keep total levels of hybridization between 0 and 50% (this value is allowed to change freely during subsequent iterations, from 0 to 100%).
            flagP = 1
         # 
         #myPRIORS[2] = auxP[0][0]   
         #myPRIORS[3] = auxP[0][0]  
         myPRIORS[4] = auxP[0][0]
         myPRIORS[5] = auxP[0][1]
         myPRIORS[6] = auxP[0][2]
         if HomoeologousExchange == "1" :
             myPRIORS[3] = auxP[0][3]
             myPRIORS[7] = auxP[0][4]
             myPRIORS[8] = auxP[0][5]
   elif polymodel == "AANN" or polymodel == "AAHH" :			## AANN or AAHH models
      myPRIORS[11] = np.random.rand(1,1)                       # introgression from B. pendula
      flagP = 0
      while flagP == 0:
         if HomoeologousExchange == "0" :      # no  homoeologous exchange allowed
            auxP = np.random.rand(1,3)
         elif HomoeologousExchange == "1" :   # homoeologous exchange allowed
            auxP = np.random.rand(1,6)
         #
         if np.sum(auxP) <= 0.5:        # for the initial priors (only), keep total levels of hybridization between 0 and 50% (this value is allowed to change freely during subsequent iterations, from 0 to 100%).
            flagP = 1
         # 
         #myPRIORS[2] = auxP[0][0]   
         myPRIORS[3] = auxP[0][0]  
         myPRIORS[4] = auxP[0][1]
         myPRIORS[5] = auxP[0][1]   # H2=H3
         myPRIORS[6] = auxP[0][2]
         if HomoeologousExchange == "1" :
             myPRIORS[7] = auxP[0][3]
             myPRIORS[8] = auxP[0][4]
             myPRIORS[9] = auxP[0][5]    
   elif polymodel == "AANH" :			## AANH model
      myPRIORS[11] = np.random.rand(1,1)                       # introgression from B. pendula
      flagP = 0
      while flagP == 0:
         if HomoeologousExchange == "0" :      # no  homoeologous exchange allowed
            auxP = np.random.rand(1,3)
         elif HomoeologousExchange ==  "1" :   # homoeologous exchange allowed
            auxP = np.random.rand(1,7)
         #
         if np.sum(auxP) <= 0.5:        # for the initial priors (only), keep total levels of hybridization between 0 and 50% (this value is allowed to change freely during subsequent iterations, from 0 to 100%).
            flagP = 1
         # 
         #myPRIORS[2] = auxP[0][0]   
         #myPRIORS[3] = auxP[0][0]  
         myPRIORS[4] = auxP[0][0]
         myPRIORS[5] = auxP[0][1]
         myPRIORS[6] = auxP[0][2]
         if HomoeologousExchange == "1" :
             #myPRIORS[3] = auxP[0][3] 
             myPRIORS[7] = auxP[0][3]
             myPRIORS[8] = auxP[0][4]   
             myPRIORS[9] = auxP[0][5]
             myPRIORS[10] = auxP[0][6]
      #
   else :
       print('\n Error: Cannot generate random priors -- polyploidization model unknown.')
       print('Check model name or update funcion GenerateRandomPriors_function() \n')
       exit()
   #
   #       
   myPRIORS = myPRIORS.tolist()                   
   myPRIORS[0] = int(myPRIORS[0])
   
   return(myPRIORS)





def UpdatePriors_function(polymodel, myPRIORS, step_size, HomoeologousExchange, geneFlowFlag) :
    
   PRIORS_new = np.copy(myPRIORS)
   flagP = 0
   #
   if polymodel == "PPPP" :            ## PPPP model: H4 = H6 = H7 = H8 = 0; sum(H0 to H8) <= 1
      while flagP == 0:
         aux_np = np.random.uniform(-1, 1, 5)*step_size
         PRIORS_new[1] = myPRIORS[1] + aux_np[0]
         PRIORS_new[3] = myPRIORS[3] + aux_np[1]
         PRIORS_new[4] = myPRIORS[4] + aux_np[2]
         PRIORS_new[5] = myPRIORS[5] + aux_np[3]
         PRIORS_new[7] = myPRIORS[7] + aux_np[4]
         if (np.sum(PRIORS_new[2:10]) <= 1) and (( PRIORS_new[1:] >= 0) & (PRIORS_new[1:] <= 1)).all() :     # check priors are within expected ranges
            flagP = 1
   #
   elif polymodel == "AAAA" :       ## AAAA model: H4 = H6 = H7 = H8 = 0; sum(H0 to H8) <= 1
      while flagP == 0:
         aux_np = np.random.uniform(-1, 1, 7)*step_size
         PRIORS_new[1] = myPRIORS[1] + aux_np[0]
         PRIORS_new[3] = myPRIORS[3] + aux_np[1]
         PRIORS_new[4] = myPRIORS[4] + aux_np[2]
         PRIORS_new[5] = myPRIORS[5] + aux_np[3]
         PRIORS_new[7] = myPRIORS[7] + aux_np[4]
         PRIORS_new[8] = myPRIORS[8] + aux_np[5]
         PRIORS_new[11] = myPRIORS[11] + aux_np[6]
         if (np.sum(PRIORS_new[2:10]) <= 1) and (( PRIORS_new[1:] >= 0) & (PRIORS_new[1:] <= 1)).all() :     # check priors are within expected ranges
            flagP = 1  
   elif polymodel == "PPNN" or polymodel == "PPHH" :			## PPNN or PPHH models  
      while flagP == 0:
         if HomoeologousExchange == "0" :      # no  homoeologous exchange allowed
            aux_np = np.random.uniform(-1, 1, 4)*step_size
         elif HomoeologousExchange ==  "1" :   # homoeologous exchange allowed
            aux_np = np.random.uniform(-1, 1, 6)*step_size
         #    
         PRIORS_new[1] = myPRIORS[1] + aux_np[0]
         PRIORS_new[3] = myPRIORS[3] + aux_np[1]
         PRIORS_new[4] = myPRIORS[4] + aux_np[2]
         PRIORS_new[5] = myPRIORS[5] + aux_np[2]    # H2=H3
         PRIORS_new[6] = myPRIORS[6] + aux_np[3]
         if HomoeologousExchange == "1" :
             PRIORS_new[7] = myPRIORS[7] + aux_np[4]
             PRIORS_new[8] = myPRIORS[8] + aux_np[5]
         if (np.sum(PRIORS_new[2:10]) <= 1) and (( PRIORS_new[1:] >= 0) & (PRIORS_new[1:] <= 1)).all() :     # check priors are within expected ranges
            flagP = 1
   elif polymodel == "PPPyPy" :			## PPPyPy model
      while flagP == 0:
         if HomoeologousExchange == "0" :      # no  homoeologous exchange allowed
            aux_np = np.random.uniform(-1, 1, 8)*step_size
         elif HomoeologousExchange ==  "1" :   # homoeologous exchange allowed
            aux_np = np.random.uniform(-1, 1, 10)*step_size
         #    
         PRIORS_new[1] = myPRIORS[1] + aux_np[0]
         PRIORS_new[2] = myPRIORS[2] + aux_np[1]
         PRIORS_new[3] = myPRIORS[3] + aux_np[2]
         PRIORS_new[4] = myPRIORS[4] + aux_np[3]
         PRIORS_new[5] = myPRIORS[5] + aux_np[4]
         PRIORS_new[6] = myPRIORS[6] + aux_np[5]
         PRIORS_new[8] = myPRIORS[8] + aux_np[6]
         PRIORS_new[9] = myPRIORS[9] + aux_np[7]
         if HomoeologousExchange == "1" :
             PRIORS_new[7] = myPRIORS[7] + aux_np[8]
             PRIORS_new[10] = myPRIORS[10] + aux_np[9]
         if (np.sum(PRIORS_new[2:10]) <= 1) and (( PRIORS_new[1:] >= 0) & (PRIORS_new[1:] <= 1)).all() :     # check priors are within expected ranges
            flagP = 1
   elif polymodel == "PPNH" :			## PPNH model
      while flagP == 0:
         if HomoeologousExchange == "0" :      # no  homoeologous exchange allowed
            aux_np = np.random.uniform(-1, 1, 4)*step_size
         elif HomoeologousExchange ==  "1" :   # homoeologous exchange allowed
            aux_np = np.random.uniform(-1, 1, 7)*step_size
         #    
         PRIORS_new[1] = myPRIORS[1] + aux_np[0]
         #PRIORS_new[3] = myPRIORS[3] + aux_np[1]
         PRIORS_new[4] = myPRIORS[4] + aux_np[1]
         PRIORS_new[5] = myPRIORS[5] + aux_np[2]
         PRIORS_new[6] = myPRIORS[6] + aux_np[3]
         if HomoeologousExchange == "1" :
             PRIORS_new[3] = myPRIORS[3] + aux_np[4]
             PRIORS_new[7] = myPRIORS[7] + aux_np[5]
             PRIORS_new[8] = myPRIORS[8] + aux_np[6]
         if (np.sum(PRIORS_new[2:10]) <= 1) and (( PRIORS_new[1:] >= 0) & (PRIORS_new[1:] <= 1)).all() :     # check priors are within expected ranges
            flagP = 1
   elif polymodel == "AANN" or polymodel == "AAHH" :			## AANN or AAHH models  
      while flagP == 0:
         if HomoeologousExchange == "0" :      # no  homoeologous exchange allowed
            aux_np = np.random.uniform(-1, 1, 5)*step_size
         elif HomoeologousExchange ==  "1" :   # homoeologous exchange allowed
            aux_np = np.random.uniform(-1, 1, 8)*step_size
         #    
         PRIORS_new[1] = myPRIORS[1] + aux_np[0]
         PRIORS_new[3] = myPRIORS[3] + aux_np[1]
         PRIORS_new[4] = myPRIORS[4] + aux_np[2]
         PRIORS_new[5] = myPRIORS[5] + aux_np[2]    # H2=H3
         PRIORS_new[6] = myPRIORS[6] + aux_np[3]
         PRIORS_new[11] = myPRIORS[11] + aux_np[4]
         if HomoeologousExchange == "1" :
             PRIORS_new[7] = myPRIORS[7] + aux_np[5]
             PRIORS_new[8] = myPRIORS[8] + aux_np[6]
             PRIORS_new[9] = myPRIORS[9] + aux_np[7]
         if (np.sum(PRIORS_new[2:10]) <= 1) and (( PRIORS_new[1:] >= 0) & (PRIORS_new[1:] <= 1)).all() :     # check priors are within expected ranges
            flagP = 1
   elif polymodel == "AANH" :			## AANH model
      while flagP == 0:
         if HomoeologousExchange == "0" :      # no  homoeologous exchange allowed
            aux_np = np.random.uniform(-1, 1, 5)*step_size
         elif HomoeologousExchange ==  "1" :   # homoeologous exchange allowed
            aux_np = np.random.uniform(-1, 1, 9)*step_size
         #    
         PRIORS_new[1] = myPRIORS[1] + aux_np[0]
         #PRIORS_new[3] = myPRIORS[3] + aux_np[1]
         PRIORS_new[4] = myPRIORS[4] + aux_np[1]
         PRIORS_new[5] = myPRIORS[5] + aux_np[2]
         PRIORS_new[6] = myPRIORS[6] + aux_np[3]
         PRIORS_new[11] = myPRIORS[11] + aux_np[4]
         if HomoeologousExchange == "1" :
             #PRIORS_new[3] = myPRIORS[3] + aux_np[5]
             PRIORS_new[7] = myPRIORS[7] + aux_np[5]
             PRIORS_new[8] = myPRIORS[8] + aux_np[6]
             PRIORS_new[9] = myPRIORS[9] + aux_np[7]
             PRIORS_new[10] = myPRIORS[10] + aux_np[8]
         if (np.sum(PRIORS_new[2:10]) <= 1) and (( PRIORS_new[1:] >= 0) & (PRIORS_new[1:] <= 1)).all() :     # check priors are within expected ranges
            flagP = 1
      #
   else :
       print('\n Error: Cannot update random priors -- polyploidization model unknown.')
       print('Check model name or update funcion UpdatePriors_function() \n')
       exit()
   
   if geneFlowFlag == "0" :
      PRIORS_new[1] = 0
         
   PRIORS_new = PRIORS_new.tolist()                   
   PRIORS_new[0] = int(PRIORS_new[0])
            
   return(PRIORS_new)

    
    
    
######################################################### Initialize
#########################################################


################# print options (simgle line)
np.set_printoptions(linewidth=np.inf)


################# Remember initial path
SRCDIR_INI = os.path.abspath(os.getcwd())


################# Read input and out folder paths
try:
    AA = sys.argv[1]                  				    # input folder
except:
    print('\n Error: provide the following arguments:')
    print('05_ABC_SimAnnealing_optimization.py [input folder] [output folder] [Observed_pairing_patterns.txt] [No. replicates] [No genes per replicate] [ILS reduction prior] [polyploidization model] [homoeologous exchange flag] [gene flow flag] [prior flag] [initial priors] \n')
    exit() 
    
try:
    RR_MAIN = sys.argv[2]                  				# output folder (MAIN)
except:
    print('\n Error: provide the following arguments:')
    print('05_ABC_SimAnnealing_optimization.py [input folder] [output folder] [Observed_pairing_patterns.txt] [No. replicates] [No genes per replicate] [ILS reduction prior] [polyploidization model] [homoeologous exchange flag] [gene flow flag] [prior flag] [initial priors] \n')
    exit() 


################# Read observed pairing pattern 
PP_obs = list()

try:
    fhand = open(sys.argv[3], 'r')                  			
except:
    print('\n Error: input file missing.')
    print('05_ABC_SimAnnealing_optimization.py [input folder] [output folder] [Observed_pairing_patterns.txt] [No. replicates] [No genes per replicate] [ILS reduction prior] [polyploidization model] [homoeologous exchange flag] [gene flow flag] [prior flag] [initial priors] \n')
    exit() 

for line in fhand:
    x = line.split()
    x = [float(i) for i in x]
    PP_obs.append(x)

PP_obs = np.array(PP_obs)



################# Read initial simulation conditions

try:
    REP = int(sys.argv[4])                  				# Number of replicates
except:
    print('\n Error: provide the following arguments:')
    print('05_ABC_SimAnnealing_optimization.py [Observed_pairing_patterns.txt] [No. replicates] [No. genes per replicate] [ILS reduction prior] [polyploidization model] [prior flag] [initial priors] [input folder] [output folder] \n')
    exit() 
    
try:
    NLOCUS = int(sys.argv[5])                  				# Number of genes per replicate
except:
    print('\n Error: provide the following arguments:')
    print('05_ABC_SimAnnealing_optimization.py [input folder] [output folder] [Observed_pairing_patterns.txt] [No. replicates] [No genes per replicate] [ILS reduction prior] [polyploidization model] [homoeologous exchange flag] [gene flow flag] [prior flag] [initial priors] \n')
    exit()     

try:
    ILS_PRIOR_init = int(sys.argv[6])                  				# ILS attenuation (prior)
except:
    print('\n Error: provide the following arguments:')
    print('05_ABC_SimAnnealing_optimization.py [input folder] [output folder] [Observed_pairing_patterns.txt] [No. replicates] [No genes per replicate] [ILS reduction prior] [polyploidization model] [homoeologous exchange flag] [gene flow flag] [prior flag] [initial priors] \n')
    exit() 



################# Read polyploidization model name
try:
    polymodel = sys.argv[7]                  				# polyploidization model name
except:
    print('\n Error: provide the following arguments:')
    print('05_ABC_SimAnnealing_optimization.py [input folder] [output folder] [Observed_pairing_patterns.txt] [No. replicates] [No genes per replicate] [ILS reduction prior] [polyploidization model] [homoeologous exchange flag] [gene flow flag] [prior flag] [initial priors] \n')
    exit() 



################# Read Homoeologous Exchange flag
try:
    HomoeologousExchange = sys.argv[8]                  				# Homoeologous Exchange flag
except:
    print('\n Error: provide the following arguments:')
    print('05_ABC_SimAnnealing_optimization.py [input folder] [output folder] [Observed_pairing_patterns.txt] [No. replicates] [No genes per replicate] [ILS reduction prior] [polyploidization model] [homoeologous exchange flag] [gene flow flag] [prior flag] [initial priors] \n')
    exit() 


################# pendula-platyphylla gene flow flag
try:
    geneFlowFlag = sys.argv[9]                  				# gene flow flag
except:
    print('\n Error: provide the following arguments:')
    print('05_ABC_SimAnnealing_optimization.py [input folder] [output folder] [Observed_pairing_patterns.txt] [No. replicates] [No genes per replicate] [ILS reduction prior] [polyploidization model] [homoeologous exchange flag] [gene flow flag] [prior flag] [initial priors] \n')
    exit() 


################# Initial prior set

## Read priors flag
try:
    Priorflag = sys.argv[10]                  				# priors flag
except:
    print('\n Error: provide the following arguments:')
    print('05_ABC_SimAnnealing_optimization.py [input folder] [output folder] [Observed_pairing_patterns.txt] [No. replicates] [No genes per replicate] [ILS reduction prior] [polyploidization model] [homoeologous exchange flag] [gene flow flag] [prior flag] [initial priors] \n')
    exit() 
 


if Priorflag == "user" :			# starting prior set provided by user
   myPRIORS = list()
   try:
      PP = sys.argv[11]                  				    
      fhand_pm = open(PP, 'r')                  			
   except:
      print('\n Error: provide the following arguments:')
      print('05_ABC_SimAnnealing_optimization.py [input folder] [output folder] [Observed_pairing_patterns.txt] [No. replicates] [No genes per replicate] [ILS reduction prior] [polyploidization model] [homoeologous exchange flag] [gene flow flag] [prior flag] [initial priors] \n')
      exit() 
   #
   for line in fhand_pm :
      x = line.split()
      myPRIORS.append(float(x[0]))
   #
   myPRIORS[0] = int(myPRIORS[0])
   #
elif Priorflag == "random" :			# random generated starting prior set, apart from ILS value, which is provided by the user
   myPRIORS = GenerateRandomPriors_function(polymodel, ILS_PRIOR_init, HomoeologousExchange, geneFlowFlag)                   
   #
   # Save priors to file
   os.chdir(RR_MAIN)
   PPrand = "00_PRIORS_ABC_SimAnnealing.txt"
   outfile1 = open(PPrand, 'w')   
   for k in myPRIORS :
      outfile1.write(str(k))
      outfile1.write('\n')
   #
   outfile1.close()
   os.chdir(SRCDIR_INI)
   PP = RR_MAIN + "/" + PPrand
   #
else: 
   print('\n Error: prior flag not found. Please provide the following arguments:')
   print('05_ABC_SimAnnealing_optimization.py [input folder] [output folder] [Observed_pairing_patterns.txt] [No. replicates] [No genes per replicate] [ILS reduction prior] [polyploidization model] [homoeologous exchange flag] [gene flow flag] [prior flag] [initial priors] \n')

#print("\n myPRIORS:", myPRIORS)
#print("np.sum(myPRIORS[2:])", np.sum(myPRIORS[2:]))
#print("PP=", PP)
#bp()   





    
################# Get simulated pairing patterns based on initial prior set
cmd = ['./20_POLARIZATION_IQTREE_PAIRINGprofiles_MAIN.sh' + ' ' +  AA + ' ' + RR_MAIN + ' ' + polymodel + ' ' + str(REP) + ' ' + str(NLOCUS) + ' ' + PP]
process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
res = process.communicate()

print("res =", res)
print("stderr =", res[1])
for line in res[0].decode(encoding='utf-8').split('\n'):
   print(line)

process.wait()

print("\n")
print("Iteration: 1")
print("PRIORS:", myPRIORS)
print("Sub-routine process returncode (0 = no errors):", process.returncode)   # 0 return code indicates there were no errors 

#now = datetime.datetime.now()
#print(now.strftime("%Y-%m-%d %H:%M:%S"))



################# Read simulated pairing patterns

RRpp = RR_MAIN + "/03_PAIRINGfrequencies"
os.chdir(RRpp)
try:
   fhand_pp = open("sister_ID_analysis_priors-JOINT_summary.txt", 'r')                  			
except:
   print('\n Error: simulated pairing profiles missing.')
   exit()

os.chdir(SRCDIR_INI)


PP_sim = list()
for line in fhand_pp:
    x = line.split()
    x = [float(i) for i in x]
    PP_sim.append(x)

PP_sim = np.array(PP_sim)

#print('\n PP_sim:', PP_sim)


################# Compute initial L2 distance (compare observed and estimated pairing patterns)
PPobs = np.concatenate(PP_obs)
PPsim = np.concatenate(PP_sim)

#print('\n PPobs:', PPobs)
#print('PPsim:', PPsim)

## weights used to compute L2 (w_max value applies to pendula, platyphylla, nana and humilis peaks; w_min value applied to all other peaks)
w_max = 5
w_min = 1
w_all = np.array(4*[w_max,w_max,w_max,w_min,w_min,w_max,w_min,w_min,w_min,w_max,w_min,w_min,w_min,w_min,w_min])

## L2 distance
L2_ref = -1
L2_ref = L2norm_function(PPobs, PPsim, w_all)


# Store results folder
cmd = ['rm' + ' ' + '-rf' + ' ' + RR_MAIN + '/03_PAIRINGfrequencies_BEST']
process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
res = process.communicate()

cmd = ['mv' + ' ' + RR_MAIN + '/03_PAIRINGfrequencies' + ' ' + RR_MAIN + '/03_PAIRINGfrequencies_BEST']
process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
res = process.communicate()

if Priorflag == "user" :
   cmd = ['cp' + ' ' + PP + ' ' + RR_MAIN + '/03_PAIRINGfrequencies_BEST' + '/00_PRIORS_ABC_SimAnnealing.txt']
   process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
   res = process.communicate()
elif Priorflag == "random" :
   cmd = ['mv' + ' ' + PP + ' ' + RR_MAIN + '/03_PAIRINGfrequencies_BEST' + '/00_PRIORS_ABC_SimAnnealing.txt']
   process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
   res = process.communicate()

# Print initial L2
print("L2 =", L2_ref)
now = datetime.datetime.now()
print(now.strftime("%Y-%m-%d %H:%M:%S"))



################# Open log file
os.chdir(RR_MAIN)
fileOUT = "ABC_SimAnnealing_log.txt"
outfile2 = open(fileOUT, 'w')   

out1 = "ABC_SimAnnealing log file"
outfile2.write(out1)
outfile2.write('\n')
out1 = 'Date' + '\t' + 'Iteration' + '\t' + 'Priors' + '\t' + 'L2'
outfile2.write(out1)
outfile2.write('\n')
now = datetime.datetime.now()
out1 = str(now.strftime("%Y-%m-%d %H:%M:%S")) + '\t' + '1' + '\t' + str(myPRIORS) + '\t' + str(L2_ref)
outfile2.write(out1)
outfile2.write('\n')

outfile2.close()
os.chdir(SRCDIR_INI)





######################################################### Start priors optimization using simulated annealing
#########################################################

opt_OUT = simann_control(np.array(myPRIORS), L2_ref, HomoeologousExchange, w_all, geneFlowFlag)




######################################################### Optimized solution:
#########################################################
print('/n')
print('/n')
print("ABC-SimAnnealing has finished.")
print("Results were saved to:", opt_OUT[3])
print("Number of iterations:", opt_OUT[2])
print("Optimized solution:")
print("PRIORS:", opt_OUT[0])
print("L2 =", opt_OUT[1])


####### Save to file
os.chdir(opt_OUT[3])
fileOUT = "ABC_SimAnnealing_optimal_conditions.txt"
outfile3 = open(fileOUT, 'w')   
out1 = "ABC_SimAnnealing: Optimized results"
outfile3.write(out1)
outfile3.write('\n')
outfile3.write('\n')
out1 = "Number of iterations:" + ' ' +str(opt_OUT[2])
outfile3.write(out1)
outfile3.write('\n')
out1 = "Optimized solution:"
outfile3.write(out1)
outfile3.write('\n')
out1 = "PRIORS:" + ' ' + str(opt_OUT[0])
outfile3.write(out1)
outfile3.write('\n')
out1 = "L2 =" + ' ' + str(opt_OUT[1])
outfile3.write(out1)
outfile3.write('\n')
outfile3.write('\n')
now = datetime.datetime.now()
outfile3.write(str(now.strftime("%Y-%m-%d %H:%M:%S")))
outfile3.write('\n')

outfile3.close()
os.chdir(SRCDIR_INI)



now = datetime.datetime.now()
print('/n')
print(now.strftime("%Y-%m-%d %H:%M:%S"))


#### END





