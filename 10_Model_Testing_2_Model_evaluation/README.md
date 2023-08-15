# Polyploidization models evaluation


## 1. Initial model parameter optimization using simulated annealing

1. Download contents of this folder to local sever.
2. Run 05_SimAnnealing_MAIN.sh script 
- remember to update paths; you must also update the bash file header and paths to iqtree, R, and python3.
- you must first generate simulated MSAs using the scripts located in the 09_Model_Testing_1_Polyploidization_and_hybridization_simulations folder [here](https://github.com/LLN273/Complex-Polyploids/tree/main/09_Model_Testing_1_Polyploidization_and_hybridization_simulations). You can also use the synthetic MSAs generated for this paper, available at dryad [here](https://doi.org/10.5061/dryad.5tb2rbp9f). (See folder 07_SIMULATION_DATA_AND_POLYPLOIDIZATION_MODEL_TESTING/01_simphy_ILS/)
- this step is computationally heavy; it is recommended that you run this script in a cluster environment and reserve a whole node with twenty or more cores.
3. Monitor ABC_SimAnnealing_log.txt file in results folder to check for convergence. For best results, you should let the algorithm run for several hundred or even thousands of iterations. Using multiple chains (using separate bash jobs) also helps to avoid local minima.

## 2. Generate ABC priors based on results from simulated annealing results

1. Run script 30_ABC_generate_priors_based_on_SA_results_MAIN.py
2. You must provide a parameter file. See 30_ABC_PARAMETERS_PPPP_SVsouth.txt for an example of a suitable parameter file.

## 3. Final parameter optimization using ABC
1. Run script 31_ABC_Final_Optimization_MAIN.sh
- Remember to update paths in script 31_ABC_Final_Optimization_query.sh; you must also update the bash file header and paths to iqtree, R, and python3.
- this step is computationally heavy; it is recommended that you run this script in a cluster environment and reserve a whole node with twenty or more cores.
2. Once job is done, you can run the following scripts to get the L2 distance for each ABC simulation run:
- 32_glean_results_ABC_Final_Optimization.sh
- 33_compute_L2-norm_ABC_Final_Optimization.sh


