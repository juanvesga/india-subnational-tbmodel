
Simulate <- function(input_data,  # characteristics of the district 
                     params,      # set of best fit parameters
                     tx_out_fl,   # Target for first-line treatment success
                     tx_out_sl,   # Target for second-line treatment success
                     xpert_uf,    # 1-Proportion of diagnoses being done by CB-NAAT vs smear (see change in order vs.) 
                     xpert_fu,    # Proportion of diagnoses having DST result within two weeks of diagnosis
                     
                     pse,         # Target proportion of private providers to be engaged
                     xpert_pe,    # Proportion of diagnoses using Xpert in private engaged
                     tx_out_pe,   # Target treatment outcome in privatre engaged
                     
                     ptreach,     # Population  (< 5 only or all household contacts):  0.1 for childen, 1 for al population
                     hcscreen,    # Number of contacts envisaged to be screened per index case
                     ptdef,       # Anticipated regimen completion rates
                     
                     acfhi,       # Proportion of risk group to be screened per year
                     xray_acf,    # Screening algorithm (verbally, or also with X-ray):  Xray =1 ; Verbally =0;
                     xpert_acf,   # Confirmatory test (smear or Xpert): Xpert =1 ; Smear =0;
                     acf_linked,  # Proportion of presumptives successfully linked to microbiological testing’
                     acf_tx       # Proportion of diagnosed cases successfully initiating treatment 
) {
  

  #__________________________________________________________________________
  #  Select Intervention parameters
  #__________________________________________________________________________
  
  itv <- list()
  #---------- 1)  Improve first-line treatment outcomes
  itv$pdef <- 1-tx_out_fl
  
  #---------- 2) Improve second-line treatment outcomes
  itv$pdef2 <- 1-tx_out_sl
  
  #---------- 3) Universal DST in routine, public TB services
  # Increase probability of being diagnosed with xpert up front
  # and getting a DST (xpert) test done
  
  itv$xpert_uf <- xpert_uf
  itv$xpert_fu <- xpert_fu
  
  
  
  #---------- 4) Private sector engagement
  # Engae a fraction of the private providers
  # and improve probability of Dx
  itv$pse     <- pse
  itv$xpert_pe<- xpert_pe
  itv$pdef_pe  <- 1-tx_out_pe
  
  # ----------5)Preventive therapy amongst household contacts
  itv$ptreach  <- ptreach    # 0.1 for childen, 1 for al population
  itv$hcscreen <- hcscreen   #number of household contacts to e screened
  itv$ptdef    <- 1-ptdef #Expected completion rates
  
  
  # ----------5)Active case finding (ACF)
  itv$acfhi       <- acfhi # Prportion of slum population reached
  itv$xray_acf    <- xray_acf   # ACF screening by Xray =1 ; Verbally =0;
  itv$xpert_acf   <- xpert_acf   # Confirmation by Xpert =1 ; smear =0;
  itv$acf_linked  <- acf_linked # Operational losses on expected yield
  itv$acf_tx      <- acf_tx # Operational losses on expected yield
  
  
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
  itvs <-sim_pack_intervention(params, prm, ref, sel, agg, gps)
  
  
  return(itvs) 
  
}
