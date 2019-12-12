#cxalibration functioin: takes data as input and fits a TB model.
#Returns calibrated parameter set and model output

Calibration <- function(input_data , # Vector with data targets points, in this order:
                        #prev_all	: prevalence in 2017
                        #inc_all  : incidence estimation in 2017
                        #notif_all: cases notified 2017
                        #prpe	    : prportion of notification coming from private sector
                        #pr_mdr   : fraction of notification being Drug resistant TB
                        state,
                        district) {
  #__________________________________________________________________________
  #
  #  Get data into targets or model input
  #__________________________________________________________________________
  target_data<-as.numeric(c(input_data$prev_all,
                            input_data$inc_all,
                            input_data$notif_all,
                            input_data$prpe,
                            input_data$pr_mdr))
  
  
  
  # ------------------
  # Call all model fixed parameters required
  
  input <- get_input(input_data)
  prm <- input$prm
  ref <- input$ref
  sel <- input$sel
  agg <- input$agg
  gps <- input$gps
  
  #__________________________________________________________________________
  #  Create instances of model functions
  #__________________________________________________________________________
  
  avgpar<-colSums(prm$bds)/2
  ratio<-avgpar[1]/avgpar
  ratio<-ratio*200
  
  #Function handle for returning output
  obj       <-  function(x)
    get_objective(x, prm, ref, sel, agg, gps, target_data, FALSE)
  
  
  
  #Random draw of parameters: pick a parameter set starting point randomly
  xu = rep(0, length(input$prm$xnames))
  for (u in 1:length(input$prm$xnames)) {
    a <- input$prm$bds[1, u]
    b <- input$prm$bds[2, u]
    xu[u] <- (b - a) * runif(1) + a
  }
  
  # If a previous best parameter set exists, use it as starting point
  #otherwise, use a random start 
  dat<-as.data.frame(read_excel("data/district_data.xlsx",sheet = "Y"))
  tmp<- dat[dat$state== state & dat$district == district,  ]
  if (tmp$beta != -99) {
    x0 <- as.numeric(tmp[c(3:length(tmp))])
  } else{
    x0 <- xu
  }
  rm(dat)
  
  #__________________________________________________________________________
  #  Execute calibration
  #__________________________________________________________________________
  
  # Function handle for calibration with simplex
  obj_spx   <-  function(x)
    get_objective(x/ratio, prm, ref, sel, agg, gps, target_data)

  fminsearch <- neldermead::fminsearch
  fminbnd<-neldermead::fminbnd
  # opt <- optimset(TolX = 2.e-2, TolFun=0.1 ,MaxIter=100, Display = "iter")
  opt <- optimset(TolFun=1e-2,TolX=1e-2,MaxIter=200, Display = "iter")

  x1 <- fminbnd(fun = obj_spx,
                   x0*ratio,
                   xmin=prm$bds[1,]*ratio,
                   xmax=prm$bds[2,]*ratio,
                   options = opt)


  x  <- x1$optbase$xopt/ratio
  
  #__________________________________________________________________________
  #  Execute calibration
  #   #__________________________________________________________________________

  # # Function handle for calibration with simplex
  # obj_spx   <-  function(x)
  #   get_objective(x, prm, ref, sel, agg, gps, target_data)
  # 
  #  x1 <- optim(x0,obj_spx, 
  #              method='L-BFGS-B',
  #              lower=prm$bds[1,],
  #              upper=prm$bds[2,],
  #             control = list(maxit = 100,
  #                            trace=TRUE,
  #                            parscale = ratio))
  # x  <- x1$par
  #__________________________________________________________________________
  #  Execute calibration
  #   #__________________________________________________________________________
  # obj_spx   <-  function(x)
  # get_objective(x/ratio, prm, ref, sel, agg, gps, target_data)
  # 
  # system.time(
  # x1 <- nlm (obj_spx,x0*ratio,
  #            steptol=1e-2,
  #            gradtol=1e-2,
  #            iterlim=100,
  #            print.level=1,
  #            fscale=0)
  # )
  # 
  # x  <- x1$estimate/ratio
  # 
  #   browser()
  
  #save best parameter set in folder../res
  # file <-
  #   paste("res/", "bestset", "_", district, "_", "mle", sep = "")
  # saveRDS(x, file)
  
  # print best parameter set
  # print(x)
  
  # Run an instance of the model with the current best
  #set of parameters and save results
  
  runs <- return_output(obj, x, ref)
  
  return(list(x=x, sfin=runs$sfin))
}

#
