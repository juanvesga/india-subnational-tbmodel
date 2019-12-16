
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
m <- matrix(0, ncol = 19, nrow = niter)
z <- matrix(0,niter,24 ) #matrix(-99, ncol = 22, nrow = nrow(fulldata))

## Loading required package: iterators
cl <- makeCluster(12)
registerDoParallel(cl)


collect1 <- function(resultx=NULL)
{
  me <- list(
    resultx = resultx
    
  )
  
  # Set the name for the class
  class(me) <- append(class(me),"collect1")
  return(me)
}


oper1<-foreach(ii = 1:niter,  .packages = c("deSolve","pracma","Matrix","readxl","Rcpp", "RcppArmadillo","abind","myTBpack") )%dopar% {
 
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
  
  return(results)
  
  
}


stopCluster(cl)

foreach(i=1:niter) %do% {
  m[i,]<-oper1[[i]]$resultx
  
}

# Save parameters 
wb<-loadWorkbook(file)
writeData(wb, "Y", m, startCol = 3, startRow = 2, colNames = FALSE)

# Save workbook
saveWorkbook(wb, file, overwrite = TRUE)









##-----------------------------------------------------------------
##- Run interventions
##-----------------------------------------------------------------
# Load previously saved values to initialize State variables
Y<-as.data.frame(read_excel(file, sheet = "Y"))

## Loading required package: iterators
cl <- makeCluster(12)
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



  oper2<-foreach(ii = 1:niter,  .packages = c("deSolve","pracma","Matrix","readxl","Rcpp", "RcppArmadillo","abind","myTBpack") )%dopar% {


 # for(ii in 1:niter){

   # Get dto inputs
   X<-fulldata[ii,-c(1,2)] 
  
   #calibration params  
   Y1<-as.numeric(Y[ii,-c(1,2)])
  

  
 
  # Get dto inputs
  X<-fulldata[ii,-c(1,2)]
  Y1<-as.numeric(Y[ii,-c(1,2)])
  
  # Call intervention function "Simulate"
  
  sims<-
    Simulate (
      X, #  District base characteristics
      Y1,   # set of best fit parameters
      
      # ALL Values belwo correspont to I
      0.95,      # Target for first-line treatment success
      0.75,      # Target for second-line treatment success
      0.5,       # Proportion of diagnoses being done by CB-NAAT vs smear (see change in order vs.) 
      0.95,      # Proportion of diagnoses having DST result within two weeks of diagnosis
      
      0.95,      # Target proportion of private providers to be engaged
      0.8,       # Proportion of diagnoses using Xpert in private engaged
      0.85,      # Target treatment outcome in privatre engaged
      
      1,          #Population  (< 5 only or all household contacts):  0.1 for childen, 1 for al population
      4,         # Number of contacts envisaged to be screened per index case
      0.6,       # Anticipated regimen completion rates
      
      0.2,       # Proportion of risk group to be screened per year
      0,         # Screening algorithm (verbally, or also with X-ray):  Xray =1 ; Verbally =0;
      1,         # Confirmatory test (smear or Xpert): Xpert =1 ; Smear =0;
      0.8,       # Proportion of presumptives successfully linked to microbiological testing’
      0.8        # Proportion of diagnosed cases successfully initiating treatment 
    )
  
  results<-collect2()
  results$resultout<-c(sims$inc_itv,
                       sims$mort_itv,
                       sims$ic_fl,
                       sims$ic_sl,
                       sims$ic_sm,
                       sims$ic_xp,
                       sims$ic_xr,
                       sims$icr_all,
                       sims$inc_av_n,
                       sims$mort_av_n)
  
  return(results)
  
}

stopCluster(cl) 
 
foreach(i=1:niter) %do% {
  z[i,]<-oper2[[i]]$resultout
}




wb<-loadWorkbook(file)
writeData(wb, "Z", z, startCol = 3, startRow = 2, colNames = FALSE)
saveWorkbook(wb, file, overwrite = TRUE)


