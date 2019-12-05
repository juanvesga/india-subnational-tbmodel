
rm(list = ls())
#  LOAD LIBRARIES AND SOURCE FUNCTIONS
#__________________________________________________________________________
library(Matrix)
library(readxl)
library(deSolve)
library(pracma)
library(neldermead)
library(profvis)
library(Rcpp)
library(abind)
library(openxlsx)

source("fun/get_addresses.R")
source("fun/get_objective.R")
source("fun/goveqs_basis.R")
source("fun/scale_up.R")
source("fun/allocate_parameters.R")
source("fun/make_model.R")
source("fun/return_output.R")
source("fun/get_input.R")
source("fun/sim_pack_intervention.R")
source("fun/Calibration.R")
source("fun/Simulate.R")
sourceCpp("fun/compute_dx.cpp")
sourceCpp("fun/scale_matrix.cpp")

#__________________________________________________________________________
#  Load data 
#__________________________________________________________________________

# Load Data : This is required for calibration and simulation steps
fulldata      <-as.data.frame(read_excel("data/Single_row.xlsx"))



##-----------------------------------------------------------------
## 1) Run calibrations
##-----------------------------------------------------------------

for (ii in 1:nrow(fulldata)){
  
  # Get state and district
  state<-fulldata[ii,1]
  district<-fulldata[ii,2]
  
  # Get dto inputs
  input_data<-fulldata[ii,-c(1,2)]
  
  # Call calibration function 
  x<-Calibration(input_data, state, district)
  
  
  # Save parameters 
  xin<-data.frame(t(x$x))
  wb<-loadWorkbook("data/Single_row.xlsx")
  writeData(wb, "Y", xin, startCol = 3, startRow = 2, colNames = FALSE)
  
  # Save workbook
  saveWorkbook(wb, "data/Single_row.xlsx", overwrite = TRUE)
  
}


##-----------------------------------------------------------------
##- 2) Run interventions
##-----------------------------------------------------------------
# Load previously saved values (Y) to initialize State variables
x0<-as.data.frame(read_excel("data/Single_row.xlsx", sheet = "Y"))

for (ii in 1:nrow(fulldata)){
  state<-fulldata[ii,1]
  district<-fulldata[ii,2]
  bestset<-as.numeric(x0[ii,-c(1,2)])
  
  # Get dto inputs
  input_data<-fulldata[ii,-c(1,2)]
  
  # Call intervention function "Simulate"
  sims<-
    Simulate (
      input_data, #  District base characteristics
      bestset,   # set of best fit parameters
      c(),       # Values to initialise State Variables (leave empty-optional) 
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
      0.5       # Operational losses on expected yield
    )
  
  # Save results (Z)
  tmp<-c(sims$inc_itv,sims$mort_itv,sims$ic_fl,sims$ic_sl,sims$ic_sm,sims$ic_xp,sims$ic_xr,sims$icr_all)
  results<-data.frame(t(tmp))
  wb<-loadWorkbook("data/Single_row.xlsx")
  writeData(wb, "Z", results, startCol = 3, startRow = 2, colNames = FALSE)
  saveWorkbook(wb, "data/Single_row.xlsx", overwrite = TRUE)
  
}



