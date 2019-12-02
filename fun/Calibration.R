#cxalibration functioin: takes data as input and fits a TB model.
#Returns calibrated parameter set and model output

Calibration <- function(input_data , # Vector with data points, in this order:
                        #prev_all	: prevalence in 2017
                        #inc_all  : incidence estimation in 2017
                        #prev_hi  : prevalence in the high risk group
                        #inc_hi	  : incidence estimation i  high risk
                        #notif_all: cases notified 2017
                        #prpe	    : prportion of notification coming from private sector
                        #pr_mdr   : fraction of notification being Drug resistant TB
                        
                        district = 'Bihar_urb'  # string with name of the district , set Urban Bihar as defaut
) {
  #__________________________________________________________________________
  #
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
  
  source("fun/get_addresses.R")
  source("fun/get_objective.R")
  source("fun/goveqs_basis.R")
  source("fun/scale_up.R")
  source("fun/allocate_parameters.R")
  source("fun/make_model.R")
  source("fun/return_output.R")
  source("fun/get_input.R")
  
  
  sourceCpp("fun/compute_dx.cpp")
  sourceCpp("fun/scale_matrix.cpp")
  
  
  # Call all model fixed parameters required
  
  input <- get_input()
  prm <- input$prm
  ref <- input$ref
  sel <- input$sel
  agg <- input$agg
  gps <- input$gps
  
  #__________________________________________________________________________
  #  Create instances of model functions
  #__________________________________________________________________________
  #Function handle for returning output
  obj       <-  function(x)
    get_objective(x, prm, ref, sel, agg, gps, input_data, FALSE)
  
  #Function handle for calibration with simplex
  obj_spx   <-  function(x)
     get_objective(x, prm, ref, sel, agg, gps, input_data)
  
  
  #Random draw of parameters: pick a parameter set starting point randomly
  xu = rep(0, length(input$prm$xnames))
  for (u in 1:length(input$prm$xnames)) {
    a <- input$prm$bds[1, u]
    b <- input$prm$bds[2, u]
    xu[u] <- (b - a) * runif(1) + a
  }
  
  # If a previous best parameter set exists, use it as starting point
  #otherwise, use a random start 
  f <- paste("res/", "bestset", "_", district, "_", "mle", sep = "")
  if (file.exists(f)) {
    x0 = readRDS(f)
  } else{
    x0 = xu
  }
  
  
  #__________________________________________________________________________
  #  Execute calibration
  #__________________________________________________________________________
  fminsearch <- neldermead::fminsearch
  opt <- optimset(TolX = 1.e-2, MaxIter=400, Display = "iter")
  x1 <- fminsearch(fun = obj_spx,
                   x0,
                   options = opt)
  x  <- x1$optbase$xopt
  
  #save best parameter set in folder../res
  file <-
    paste("res/", "bestset", "_", district, "_", "mle", sep = "")
  saveRDS(x, file)
  
  # print best parameter set
  print(x)
  
  # Run an instance of the model with the current best
  #set of parameters and save results
  
  runs <- return_output(obj, x, ref, district, c())
  
  
}

#
