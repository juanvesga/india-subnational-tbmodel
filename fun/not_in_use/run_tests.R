rm(list = ls())
##-----------------------------------------------------------------
##- Run calibrations
##-----------------------------------------------------------------
# This step calls a calibration algorithm, and returns a best fit parameter set
# results saved as ..res/bestset and printed in console
# The function takes data input as targets foir calibration, and district name
#(it can be called without the string, but as default it used Bihar_urb

source("fun/Calibration.R")

# Some data to fit the model to

data<-c(300, # prevalence per 100k in 2017 
        270, # incidence per 100k in 2017	
        450, # prevalence per 100k in 2017 among high risk
        405, # incidence per 100k in 2017	in high risk
        100, # Notified cases per 100k in 2017	
        0.05,	# Fraction of notificatioons from private sector	
        0.04)# Fraction of notif. that are MDR	


# Call calibration function 

Calibration(data,'Bihar_urb')





##-----------------------------------------------------------------
##- Run interventions
##-----------------------------------------------------------------

#This step requires an existing model calibration. The function "Simulate" 
#takes as input a set of parameters for the best fitted simulation, 
#the initial values of State variables (sfin) , the simulation
#(i.e. intervention) conditions and the district name 
#Returns an object with simulation output, print sim output to console and
# results also saved in ..res/sim..
source("fun/Simulate.R")

# This parameter set is the result of calibration (same as in res/bestset...)
bestset<-c(
  12.21262429
  ,7.84303776
  ,4.17630773
  ,0.04995892
  ,1.10423528
  ,2.35512514
  ,2.35984532
  ,1.13814334
  ,0.06434389
  ,1.71717568
  ,0.54973477
  ,0.10621163
  ,0.10823185
  ,0.45737501
  ,0.04317104
  ,0.71967630
  ,0.66289460
  ,0.97275453
  ,0.31176922)


# Load previously saved values to initialize State variables
f <- paste("res/","output","_","Bihar_urb","_","mle", sep = "")
runs <- readRDS(f)
sfin<-runs$sfin

sims<-
Simulate (
  bestset,   # set of best fit parameters
  sfin,      # Values to initialise State Variables 
  0.95,      # FL completion rates
  0.75,      # SL completion rates
  0.5,       # probability xpert upfront
  0.95,      # probability xpert follow-up
  0.95,      # fraction of private providers engaged
  0.8,       # proibability of siagnosis in engaged providers
  0.8,       # prob of xpert in engaged providers
  1,         # 0.1 for childen, 1 for al population
  4,         # number of household contacts to e screened
  0.6,       # Expected completion rates
  0.2,       # Prportion of slum population reached
  0,         # ACF screening by Xray =1 ; Verbally =0;
  1,         # Confirmation by Xpert =1 ; smear =0;
  0.5,       # Operational losses on expected yield
  'Bihar_urb'# string with name of the district , set Urban Bihar as defaut
)


#--- Others

# Plot simulation results
source("fun/plot_sim_intv.R")
plot_sim_intv(sims)



# Plot model fits
source("fun/plot_targets.R")
plot_targets(runs, data, 'Bihar_urb')

