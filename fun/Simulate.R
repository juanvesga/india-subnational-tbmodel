
Simulate <- function(input_data,# characteristics of the district 
                     params,    # set of best fit parameters
                     sfin,      # Values to initialise State Variables 
                     pdef,      # FL completion rates
                     pdef2,     # SL completion rates
                     xpert_uf,  # probability xpert upfront
                     xpert_fu,  # probability xpert follow-up
                     pse,       # fraction of private providers engaged
                     pDx_pe,    # proibability of siagnosis in engaged providers
                     xpert_fu_pe, # prob of xpert in engaged providers
                     ptreach,   # 0.1 for childen, 1 for al population
                     hcscreen,  # number of household contacts to e screened
                     ptdef,     # Expected completion rates
                     acfhi,     # Prportion of slum population reached
                     xray_acf,  # ACF screening by Xray =1 ; Verbally =0;
                     xpert_acf, # Confirmation by Xpert =1 ; smear =0;
                     acf_loss  # Operational losses on expected yield
) {
  

  #__________________________________________________________________________
  #  Select Intervention parameters
  #__________________________________________________________________________
  
  itv <- list()
  #---------- 1)  Improve first-line treatment outcomes
  itv$pdef <- 1-pdef
  
  #---------- 2) Improve second-line treatment outcomes
  itv$pdef2 <- 1-pdef2
  
  #---------- 3) Universal DST in routine, public TB services
  # Increase probability of being diagnosed with xpert up front
  # and getting a DST (xpert) test done
  
  itv$xpert_uf <- xpert_uf
  itv$xpert_fu <- xpert_fu
  
  
  
  #---------- 4) Private sector engagement
  # Engae a fraction of the private providers
  # and improve probability of Dx
  itv$pse <- pse
  itv$pDx_pe <- pDx_pe
  itv$xpert_fu_pe <- xpert_fu_pe
  
  # ----------5)Preventive therapy amongst household contacts
  itv$ptreach  <- ptreach    # 0.1 for childen, 1 for al population
  itv$hcscreen <- hcscreen   #number of household contacts to e screened
  itv$ptdef    <- 1-ptdef #Expected completion rates
  
  
  # ----------5)Active case finding (ACF)
  itv$acfhi     <- acfhi # Prportion of slum population reached
  itv$xray_acf  <- xray_acf   # ACF screening by Xray =1 ; Verbally =0;
  itv$xpert_acf <- xpert_acf   # Confirmation by Xpert =1 ; smear =0;
  itv$acf_loss  <- acf_loss # Operational losses on expected yield
  
  
  # Call all model fixed parameters required
  input <- get_input(input_data)
  prm <- input$prm
  ref <- input$ref
  prm$itv <- itv
  
  sel <- input$sel
  agg <- input$agg
  gps <- input$gps
  
  
  #__________________________________________________________________________
  #  Execute simulation of intervention
  #__________________________________________________________________________
  itvs <-sim_pack_intervention(params, sfin, prm, ref, sel, agg, gps, c())
  
  # save results
  # file <- paste("res/", "sim", "_", district,sep ="")
  # saveRDS(itvs, file) 
  
  
  #Print out results
  # yrs<-c(2018:2025)
  # incidence=itvs$inc
  # mortality=itvs$mort
  # out1<-as.data.frame(cbind(yrs,incidence,mortality))
  # print(out1)
  # 
  # icr_fl<-itvs$ic_fl
  # icr_sl<-itvs$ic_sl
  # icr_sm<-itvs$ic_sm
  # icr_xpert<-itvs$ic_xp
  # icr_xray<-itvs$ic_xr
  # icr_total<-itvs$icr_all
  # 
  # out2<-as.data.frame(cbind(icr_fl,icr_sl,icr_sm,icr_xpert,icr_xray,icr_total))
  # print(out2)
  # 
  # return object with simulation output
  return(itvs) 
  
}
