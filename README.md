# NUS_Singapore_Sci
# RESCUE

Please refer to the zipped folder NUS_Singapore_Sci-master_zipped.zip 

Molecular dynamics:
Bash script with the steps used to run GROMACS in the Linux platform. Note that the script does not include PDB and MDP files required to run the simulation.


gRNA creation:
Python script that allows users to generate gBlock DNA sequences for gRNA required by the Cas13b-editing system. The user inputd his/her target sequence containing the spacer region, and the program will generate the optimal DNA sequence containing the Kozak sequence for U6 promoter, U6 terminator and flanking BbsI restriction sites.


Least Square Regression:
Python script that allows for least square regression on a series of data. The user has to key in the datapoints by hand, and by running the script with different starting positions (coefficients_guess), one can obtain different approximated solutions, where their accuracy is reflected in their associated Square Distance (sum00).


Ensemble Kalman Filter:
Filter for a simple chemistry model for RESCUE
For documentation, see RESCUE_EnKF_wiki.pdf

Running script: run_idealized_enkf_test.py

Results from 20000 ensemble runs before EnKF are available in prior_ens.txt

Function libraries:
1) funclib_enkf.py	----- Library of functions for the EnKF
2) funclib_model.py	----- Library of functions for running a simple RESCUE model
3) funclib_stats.py	----- Library of functions for making histogram plots
