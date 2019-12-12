
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
library(ggplot2)
library(gridExtra)
library(abind)
library(openxlsx)
library(foreach)
library(doParallel)
library(myTBpack)
library(progress)
# library(doSNOW)
# library(bigstatsr)

source("fun/get_addresses.R")
source("fun/get_distribution_fns.R")
source("fun/get_objective.R")
source("fun/goveqs_basis.R")
source("fun/scale_up.R")
source("fun/allocate_parameters.R")
source("fun/make_model.R")
source("fun/return_output.R")
source("fun/plot_targets.R")
source("fun/get_input.R")
source("fun/sim_pack_intervention.R")
source("fun/plot_sim_intv.R")
source("fun/Calibration.R")
source("fun/Simulate.R")
source("fun/plot_sim_intv.R")

# sourceCpp("fun/compute_dx.cpp")
# sourceCpp("fun/scale_matrix.cpp")

#__________________________________________________________________________
#  Load data 
#__________________________________________________________________________
file<-"data/district_data.xlsx"
# Load Data
fulldata      <-as.data.frame(read_excel(file))
niter<-nrow(fulldata)


##-----------------------------------------------------------------
##- Run calibrations
##-----------------------------------------------------------------
m <- matrix(-99, ncol = 19, nrow = nrow(fulldata))
q <- matrix(-99, ncol = 106, nrow = nrow(fulldata))
z <- matrix(0,nrow(fulldata),22 ) #matrix(-99, ncol = 22, nrow = nrow(fulldata))

## Loading required package: iterators
cl <- makeCluster(11)
registerDoParallel(cl)


collect1 <- function(resultx=NULL, resultsfin=NULL)
{
  me <- list(
    resultx = resultx,
    resultsfin = resultsfin
    
  )
  
  # Set the name for the class
  class(me) <- append(class(me),"collect1")
  return(me)
}



 oper1<-foreach(ii = 1:niter,  .packages = c("deSolve","neldermead","pracma","Matrix","readxl","Rcpp", "RcppArmadillo","abind","myTBpack") )%dopar% {
  # system.time(
  #  for (ii in 1:1){ 
  
  # Get state and district
  state<-fulldata[ii,1]
  district<-fulldata[ii,2]
  
  # Get dto inputs
  input_data<-fulldata[ii,-c(1,2)]
  
  # Call calibration function 
  
  x<-Calibration(input_data, state, district)
  results<-collect1()
  results$resultx   <-as.numeric(x$x)
  results$resultsfin<-as.numeric(x$sfin)
  
  return(results)
  
  
}
# )

stopCluster(cl)

foreach(i=1:niter) %do% {
  m[i,]<-oper1[[i]]$resultx
  q[i,]<-oper1[[i]]$resultsfin
  
}

# Save parameters 
wb<-loadWorkbook(file)
writeData(wb, "Y", m, startCol = 3, startRow = 2, colNames = FALSE)

# SAve final values of state variables
writeData(wb, "sfin", q, startCol = 3, startRow = 2, colNames = FALSE)

# Save workbook
saveWorkbook(wb, file, overwrite = TRUE)









##-----------------------------------------------------------------
##- Run interventions
##-----------------------------------------------------------------
# Load previously saved values to initialize State variables
x0<-as.data.frame(read_excel(file, sheet = "Y"))
sf0<-as.data.frame(read_excel(file, sheet = "sfin"))


## Loading required package: iterators
cl <- makeCluster(10)
registerDoParallel(cl)


collect2 <- function(resultout=NULL)
{
  me <- list(
    resultout = resultout
  )
  
  # Set the name for the class
  class(me) <- append(class(me),"collect2")
  return(me)
}






# oper2<-foreach(ii = 1:niter,  .packages = c("deSolve","pracma","Matrix","readxl","Rcpp", "RcppArmadillo","abind","myTBpack") )%dopar% {


for(ii in 1:10){
  state<-fulldata[ii,1]
  print(state)
  district<-fulldata[ii,2]
  sfin<-list(sf0[ii,-c(1,2)])
  # Get dto inputs
  input_data<-fulldata[ii,-c(1,2)]
  bestset<-as.numeric(x0[ii,-c(1,2)])
  
  # Call intervention function "Simulate"
  
  sims<-
    Simulate (
      input_data, #  District base characteristics
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
      0.5       # Operational losses on expected yield
    )
  results<-collect2()
  results$resultout <- c(sims$inc_itv,sims$mort_itv,sims$ic_fl,sims$ic_sl,sims$ic_sm,sims$ic_xp,sims$ic_xr,sims$icr_all)
  return(results)
  
}

foreach(i=1:niter) %do% {
  z[i,]<-oper2[[i]]$resultout
}

stopCluster(cl)


wb<-loadWorkbook(file)
writeData(wb, "Z", z, startCol = 3, startRow = 2, colNames = FALSE)
saveWorkbook(wb, file, overwrite = TRUE)


