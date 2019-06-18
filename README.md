# neutrino_floor
Compute the neutrino floor for noble liquid based detectors
Andrea Gaspert 06/2019

###How to compute neutrino floor###

Simply need to modify and run go_oak_multi.py

Minimally, the detector material ('Ar' or 'Xe'), the mode ('cs' or 'exp') and the DM masses at which to compute the floor (list of floats) need to be defined in this file (see the go_oak_multi.py section for details and other parameters)

###file content####

###constants.py

Contains all constants relevant to different detectors, to neutrino fluxes and dark matter distributions
from constants import * to use all constants in code

###speed_dist.py

To pre-compute DM speed distribution in lab frame. Only runs if using SHM++ parameters in Monte-Carlo and if speed_beta_eta.npy
does not already exist in directory (numpy file with lab frames distributions for 64 differents combinations of beta and eta)

Running this code takes multiple hours.

###dm_events.py

Contains all functions to evaluate DM recoil dist. in detector
int_new computes differential event dist (dR/dE) for SHM with a given v0 and vesc, for a given DM mass
int_irreg computes differential event dist (dR/dE) for SHM ++ with given v0, vesc and two integers between 0-7 corresponding
to values of 0.1<=eta<=0.8 and 0.82 <= beta <=0.97. 

###neutrino_dist.py

Computes differential event dist (dR/dE) for neutrino-induced NR and ER in a given detector. Only runs if .npy files containing 
this information not found in directory (dist_neutrino_argon_long.npy and dist_neutrino_xenon_long.npy)

Only need to run once, takes multiple hours. Requires directory named data containing all text files with neutrino flux distibutions.

###neutrino_events.py

Contains all the necessary functions to compute neutrino-induced NR and ER in detectors. Necessary to run neutrino_dist.py.


###optim.py

Contains all the function to find neutrino floor cross-section for a given DM mass and exposure OR the necessary exposure to reach
the neutrino floor for a given DM mass and cross-section

###run_floor_multi.py

Called by the go file to distribute tasks amongst the available cpus

###go_oak_multi.py

File to compute the neutrino floor for a given list of DM masses. Task modified by changing this file. This file is designed to run the
code on the oak server (TRIUMF). If want to run on your own device, add line batch_mode = False.

You should also change the name of the output file to fit your own task (job_name / log_name). 

ARGS is a dictionary where one can define parameters of the detectors /  of how to compute the neutrino floor. If nothing specified, different parameters have certain 
standard values

###Parameters:

- typ : 'Ar' or 'Xe' 
To specify wether the code is computed for an Argon or a Xenon detector. This needs to be defined or else error will be raised
- mod : 'cs' or 'exp' 
'cs' if you want to find the cross section of the floor for a given DM mass and exposure. 'exp' if you want to find the exposure necessary for a given DM mass/ cross-section combination to be below floor. This needs to be defined or else error will be raised
- m: list
List of DM masses at which to compute the neutrino floor (in GeV). Need to be defined of else error will be raised
- E_min : float
Threshold energy of the detector in keVnr. If not defined, ideal threshold will be used by default (1 keV for Xenon and 15 keV for Argon)
- E_max: float
Maximal energy of ROI zone of detector in keVnr. If not defined, ideal maximal energy used by default (90 keV for Xenon and  120 keV for Argon)
- exp : float 
Exposure at which to compute floor in ton-years. Only defined if working in 'cs' mode. If not defined, exposure at which systematic errors dominate used by default (1e4 t-y for Xenon and 1e6 t-y for argon)
- cs : float 
Cross-section at which to compute floor exposure (cm^2). Only defined if working in 'exp' mode. If not defined, error will be raised.
- ER : bool
If True, neutrino-induced ER are considered in the computation of neutrino floor.
- er_rej : float 
ER rejection efficiency. Only defined if ER=True. If not defined, ideal values used by default (1e3 for Xenon and 1e8 for Argon)
- n_bin : int 
Number of energy bins considered. Energy bins are logarithmically-spaced between the threshold and maximal energy. Default is 50 if not defined. This allows the bins to be larger than the energy resolution (~8%)
- n_mc : int 
Number of pseudo-experiments to run to estimate the significance of discovery for every mass/cross-section pair. Default is 100 if not defined.
- n_dicho : int 
number of iteration before the dichotomic search for the cross-section / exposure a the floor stops. Default is 10 if not defined. This allows a precision of 2%.
- speed_dist: 'reg' or 'pp' or 'other' 
The DM speed distribution to use. If 'reg', regular SHM is used. If 'pp', SHM++ is used and algorithm also profiles over values of the sausage fraction and anisotropy. If 'other' is used, a function defining another DM speed distribution IN THE LAB FRAME is used.
- dist: function_name 
If 'other' in speed_dist, need to specify name of the function where the speed distribution of DM IN THE LAB FRAME is defined. The function also needs to be imported or defined in the optim.py file.

###extract_file.py

Code to find log files in a directory containing neutrino floor information and extracting this information in the form of a .npy 
Need to change filename to analize desired files

###figure_plot.py

Code to plot neutrino floors. If want to show exclusion curves from experiments, exp_show = True
Need to modify to load .npy files containing the results you want to show.

