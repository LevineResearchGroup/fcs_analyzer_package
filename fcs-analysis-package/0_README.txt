
########################Symphotime instructions########################
See powerpoint for more detailed instructions
In Symphotime use the following settings to create FCS curves:

FRET-FCS
A = Ch2, Al488, 0-50 ns
B = Ch1, Al594, 0-50 ns

PIE-FCCS
A = Ch2, Al488, 0-50 ns
B = Ch1, Al594, 50-100 ns

Export immediately after calculating correlation curves (rather than after fitting - i.e., do not include fits in export).

Export as .dat ascii file.
Name .dat in a standard manner, so that they can be easily sorted by the python helper script.

################Python instructions##################################
Requirements: lmfit can be installed from conda using
>conda install lmfit
Requires uncertainties package
>conda install uncertainties

Put .dat ascii files in Data folder (or change python script to point towards were you prefer to keep them).
Use FCS_calibration.py to determine calibration constants.

Use Batch_fitting_FCS.py to fit FCS experiments in batch mode.
	Batch_fitting_FCS will populate figures in Figures folder and summary results in Results folder.
	
	
###############Using the scripts###############################

Most of the user input is in the following lines of code


name = 'monomericAB_FRET_FCS_grouped' #name of .dat file in the folder ./Data


# Possible key values are DD (autocorrelation Donor Donor), AA (auto, accepptor acceptor), DxA (donor acceptor cross correlation)
key = 'DD'

set_kappa = 6.18 #from calibration
td_ref = ufloat(0.0289, 0.00033) #from calibration (ms)
D_ref = ufloat(400, 10) #from literature, for calibration (um^2/s)
temperature_ref = ufloat(22, 0.5) # temperature at which reference D was taken (celsius)
temperature_lab = ufloat(22,0.5) #our labs temeprature (celsius)

tau_diff2_fix_value = 0.036 #from 1-comp fit to monomer (ms)
	
#####################Results################

Summary of each curves fit is output to terminal
First the fitting summary for one component, then two component. See powerpoint for explanation of some variables, and how to choose each model.
See lmfit help, and statistics/model fitting literature for further explanation if necessary. 

e.g.,

Results of measurement at time t = 135.00 min
Guessing tau = 0.031698 ms

List of fitted parameters for Model(diffusion_3d): 

Name           Value     Stderr        Min        Max
A0            0.7374   0.002515       0.01        inf
Ginf      -1.688e-05  3.294e-05       -inf        inf
kappa           6.18          0       -inf        inf
tau_diff     0.03607  0.0002085      1e-06        inf
[[Model]]
    Model(diffusion_3d)
[[Fit Statistics]]
    # fitting method   = least_squares
    # function evals   = 5
    # data points      = 213
    # variables        = 3
    chi-square         = 237.465394
    reduced chi-square = 1.13078759
    Akaike info crit   = 29.1594343
    Bayesian info crit = 39.2433108
[[Variables]]
    tau_diff:  0.03607116 +/- 2.0854e-04 (0.58%) (init = 0.03169848)
    A0:        0.73736230 +/- 0.00251525 (0.34%) (init = 0.7442337)
    Ginf:     -1.6881e-05 +/- 3.2936e-05 (195.11%) (init = 0)
    kappa:     6.18 (fixed)
[[Correlations]] (unreported correlations are < 0.100)
    C(tau_diff, A0)   = -0.872
    C(tau_diff, Ginf) = -0.144
Guessing tau = 0.031698 ms

List of fitted parameters for Model(twocomp_diff3d): 

Name            Value     Stderr        Min        Max
A0             0.7373   0.002603       0.01        inf
Ginf       -1.961e-05  3.435e-05       -inf        inf
kappa            6.18          0       -inf        inf
p1           0.001434    0.01068          0          1
tau_diff1      0.1403     0.9024      1e-06        inf
tau_diff2       0.036          0      1e-06        inf
[[Model]]
    Model(twocomp_diff3d)
[[Fit Statistics]]
    # fitting method   = least_squares
    # function evals   = 49
    # data points      = 213
    # variables        = 4
    chi-square         = 237.160575
    reduced chi-square = 1.13473959
    Akaike info crit   = 30.8858444
    Bayesian info crit = 44.3310131
[[Variables]]
    p1:         0.00143370 +/- 0.01067509 (744.58%) (init = 0.5)
    tau_diff1:  0.14029246 +/- 0.90242754 (643.25%) (init = 0.03169848)
    tau_diff2:  0.036 (fixed)
    A0:         0.73726608 +/- 0.00260284 (0.35%) (init = 0.7442337)
    Ginf:      -1.9606e-05 +/- 3.4347e-05 (175.19%) (init = 0)
    kappa:      6.18 (fixed)
[[Correlations]] (unreported correlations are < 0.100)
    C(p1, tau_diff1)   = -0.976
    C(p1, A0)          = -0.804
    C(tau_diff1, A0)   =  0.707
    C(tau_diff1, Ginf) = -0.234
    C(p1, Ginf)        =  0.184
    C(A0, Ginf)        = -0.110



This is followed by a quick summary of a measurement group:

Quick summary
Analyzing DD curves
Mean 1-comp diffusion time = 0.0368 +/- 0.00066 ms
Mean 1-comp diffusion coeff slow = 313.8159 +/- 5.56288 um^2 s^-1
Mean 1-comp Rh slow = 0.7218 +/- 0.01285 nm
Mean 2-comp slow diffusion time = 0.1033 +/- 0.05002 ms
Mean 2-comp diffusion coeff slow = 147.6770 +/- 88.99894 um^2 s^-1
Mean 2-comp Rh slow = 2.0234 +/- 0.97979 nm



####################Further analysis######################

MetaPlots.py has some examples for reading the .dat files output by Batch_fitting_FCS.py in the result folder.

You could add different outputs using the template below, replacing Rh1c with a different variable of interest, e.g., p1

with open('./Results/' + name + '_Rh1c_' + key + '.dat', "w" ) as f:
    f.write('t \t R \t err \n')
    for t,R,err in zip(t_measurement, Rh1c,err_Rh1c):
        f.write('%.3f \t %.3f \t %.3f \n' %(t,R,err))    

