# Polyploidization models evaluation


## 1. Initial model parameter optimization using simulated annealing

a) Download content of this folder to local sever.
b) Run 05_SimAnnealing_MAIN.sh script 
- remember to update paths; must also update bash file header and paths to iqtree, R, and python3.
- you must first generate simulated MSAs using scripts located in the 09_Model_Testing_1_Polyploidization_and_hybridization_simulations folder [here]([https://github.com/LLN273/PolyAncestor/blob/main/01s_polarizeTETRA.py](https://github.com/LLN273/Complex-Polyploids/tree/main/09_Model_Testing_1_Polyploidization_and_hybridization_simulations))
- this step is computationally heavy; it is recommended you run this script in an cluster environment and reserve a whole node with twenty or more cores.
c) Monitor ABC_SimAnnealing_log.txt file in results folder to check for convergence. For best results, you should let the algorithm run for several hundred or even thousands of iterations. Using multiple chains (using separate bash jobs) also helps to avoid local minima.

## 2. Generate ABC priors based on results from simulated annealing results

a) Run script 30_ABC_generate_priors_based_on_SA_results_MAIN.py
b) You must provide a parameter file. See 30_ABC_PARAMETERS_PPPP_SVsouth.txt for an example of a suitable parameter file.

## 3. Final parameter optimization using ABC
a) Run script 31_ABC_Final_Optimization_MAIN.sh
- Remember to update paths in script 31_ABC_Final_Optimization_query.sh; must also update bash file header and paths to iqtree, R, and python3.
- this step is computationally heavy; it is recommended you run this script in an cluster environment and reserve a whole node with twenty or more cores.
b) Once job is done, you can run the following scripts to get the L2 distance for each ABC simulation run:
- 32_glean_results_ABC_Final_Optimization.sh
- 33_compute_L2-norm_ABC_Final_Optimization.sh


