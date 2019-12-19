
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
file <-"data/district_data.xlsx"

# Load Data : This is required for calibration and simulation steps
fulldata      <-as.data.frame(read_excel(file))

##-----------------------------------------------------------------
## 1) Run calibrations
##-----------------------------------------------------------------

# Select state and district
state<-"bihar"
district<-"darbhanga"

# Get dto inputs
tmp<- fulldata[fulldata$state== state & fulldata$district == district,  ]
X_Y0<-tmp[,-c(1,2)]

# Call calibration function (C)
# Note: If a prevoius calibration (Y0) exists function C will retrieve it and 
# use it as starting parameters
C<-Calibration(X_Y0, state, district)


# Save parameters 
Y1<-data.frame(t(C$x))
wb<-loadWorkbook(file)
writeData(wb, "Y", Y1, startCol = 3, startRow = 2, colNames = FALSE)

# Save workbook
saveWorkbook(wb, file, overwrite = TRUE)




##-----------------------------------------------------------------
##- 2) Run interventions
##-----------------------------------------------------------------
# Load previously saved values (Y) to initialize State variables
Y<-as.data.frame(read_excel(file, sheet = "Y"))
tmp<- Y[Y$state== state & Y$district == district,  ]

Y1 <- as.numeric(tmp[c(3:length(tmp))])

# Get dto inputs
tmp<- fulldata[fulldata$state== state & fulldata$district == district,  ]
X<-tmp[,-c(1,2)]

# Call intervention function "Simulate"
sims<-
  Simulate (
    X, #  District base characteristics
    Y1,   # set of best fit parameters
    
    # ALL Values belwo correspont to I
    0.95,      # Target for first-line treatment success
    
    0.75,      # Target for second-line treatment success
    
    0.5,       # 1-Proportion of diagnoses being done by CB-NAAT vs smear (see change in order vs.) 
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

# Save results (Z)
tmp<-c(sims$inc_itv,
       sims$mort_itv,
       sims$ic_fl,
       sims$ic_sl,
       sims$ic_sm,
       sims$ic_xp,
       sims$ic_xr,
       sims$icr_all,
       sims$inc_av_n,
       sims$mort_av_n)

Z<-data.frame(t(tmp))
wb<-loadWorkbook(file)
writeData(wb, "Z", Z, startCol = 3, startRow = 2, colNames = FALSE)
saveWorkbook(wb, file, overwrite = TRUE)






