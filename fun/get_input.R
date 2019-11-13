
get_input<-function(location){
  
  
  
  #__________________________________________________________________________
  #  Load data & Call distributions for lhd
  #__________________________________________________________________________
  
  # Load Data
  tmp      <-as.data.frame(read_excel("data/data.xlsx"))
  rownames(tmp)<-tmp[,1]  
  
  # All data including targets and input data 
  data_raw<-tmp[,-1]
  
  # Model fitting data only
  drops <- c('stdard_pr',	'size_hirisk',	'tsr_fl',	'tsr_sl', 'popN')
  data<-data_raw[location , !(names(data_raw) %in% drops)]
  datapoints<-data[location,]
  
  # Create likelihood distributions for datapoints
  lhd<-Make_distr_fns(datapoints,1)
  
  #__________________________________________________________________________
  #  Build model structure
  #__________________________________________________________________________
  
  # Model states
  gps<-list( sectors=c("pu","pr","pe"),
             strain =c("ds","mdr"),
             risk   =c("lo","hi"))
  
  states1 = c('U')
  states2 = c('Lf','Ls','Pt','Ia','Is','E','Rlo','Rhi','R')
  states3 = c('Dx','Tx','Tx2')
  states4 = c('S')
  states5 = c('SDx','STx')
  
  
  
  i<-list()
  s<-list()
  groups<-list(states1,gps$risk)
  ref_all <- get_addresses(groups, i, s, 0)
  groups<-list(states2,gps$risk,gps$strain)
  ref_all<-  get_addresses(groups, ref_all$i, ref_all$s, ref_all$i$nstates)
  groups<-list(states3,gps$risk,gps$strain, gps$sectors)
  ref_all<-  get_addresses(groups, ref_all$i, ref_all$s, ref_all$i$nstates)
  groups<-list(states4,gps$risk)
  ref_all<-  get_addresses(groups, ref_all$i, ref_all$s, ref_all$i$nstates)
  groups<-list(states5,gps$risk,gps$sectors)
  ref_all<-  get_addresses(groups, ref_all$i, ref_all$s, ref_all$i$nstates)
  
  #For indexing without NTS
  i<-list()
  s<-list()
  groups<-list(states1,gps$risk)
  ref <- get_addresses(groups, i, s, 0)
  groups<-list(states2,gps$risk,gps$strain)
  ref<-  get_addresses(groups, ref$i, ref$s, ref$i$nstates)
  groups<-list(states3,gps$risk,gps$strain, gps$sectors)
  ref<-  get_addresses(groups, ref$i, ref$s, ref$i$nstates)
  
  
  # Include the auxiliaries
  
  auxnames <- c('inc', 'notif' , 'mort' ,'remo' ,'pt', 'txfp','dx','pmo','acf' )
  auxinds  <- c(  3      , 3    ,  1     ,3     ,1   ,  1    , 3  , 2   ,1)
  ref$i$aux<-list()
  lim_all<-ref_all$i$nstates
  for (ii in 1:length(auxnames)){
    inds <- lim_all + (1:auxinds[ii])
    ref$i$aux[[auxnames[ii]]] <- inds
    lim_all <- inds[length(inds)]
  }
  ref$i$nx <- lim_all
  
  
  #States
  ref$s$infectious <- c(ref$s$Ia, ref$s$Is, ref$s$E, ref$s$Dx)
  ref$s$prevalent  <- c(ref$s$infectious, ref$s$Tx , ref$s$Tx2)
  ref$s$TBmort     <- ref$s$infectious
  ref$s$sought     <- c(ref$s$Dx, ref$s$E )
  ref$s$treated_pr <- intersect(c(ref$s$Tx, ref$s$Tx2) ,ref$s$pr)
  ref$s$treated_pu <- intersect(c(ref$s$Tx, ref$s$Tx2) ,ref$s$pu)
  ref$s$treated_pe <- intersect(c(ref$s$Tx, ref$s$Tx2) ,ref$s$pe)
  ref$s$treated_pay <- intersect(c(ref$s$Tx, ref$s$Tx2) ,c(ref$s$pe,ref$s$pu))
  ref$s$postdx     <-c(ref$s$Dx, ref$s$Tx, ref$s$Tx2, ref$s$Rlo, ref$s$Rhi, ref$s$R)
  ref$s$respSymp   <-c(ref_all$s$S, ref_all$s$SDx, ref_all$s$STx)
  ref$s$notbtreated_pay <-intersect(ref_all$s$STx ,c(ref_all$s$pe,ref_all$s$pu))
  
  ref$s$nstates<-(1:ref$i$nstates)
  
  #__________________________________________________________________________
  #  Prepare selectors and aggregartors of model output
  #__________________________________________________________________________
  
  tmp <- matrix(1,ref_all$i$nstates,ref_all$i$nstates)
  tmp[ref_all$s$lo,ref_all$s$hi]<-0
  tmp[ref_all$s$hi,ref_all$s$lo]<-0
  check<-tmp - diag(diag(tmp))
  
  #  Incidence: 1.Total 2.Hi 3. MDR
  tmp <- matrix(0,3,ref_all$i$nstates)
  tmp[1,ref$s$Ia] <- 1
  tmp[2,intersect(ref$s$Ia,ref$s$hi) ] <- 1
  tmp[3,intersect(ref$s$Ia,ref$s$mdr)] <- 1
  agg <-list(inc= tmp)
  
  tmp <-matrix(0,ref_all$i$nstates,ref_all$i$nstates)
  tmp[ref$s$Ia,] <- 1
  tmp<-tmp*check
  sel<-list(inc= tmp - diag(diag(tmp)))
  
  # DR acquisition selectors
  tmp <-matrix(0,ref_all$i$nstates,ref_all$i$nstates)
  tmp[intersect(ref$s$Tx,ref$s$mdr),intersect(ref$s$Tx,ref$s$ds)] <- 1
  tmp<-tmp*check
  sel$acqu <- tmp - diag(diag(tmp))
  
  tmp <-matrix(0,ref_all$i$nstates,ref_all$i$nstates)
  tmp [ref$s$Ia,ref$s$Ls] <- 1
  tmp <-tmp*check
  sel$remo <- tmp - diag(diag(tmp))
  
  # TB notification in RNTCP
  tmp <- matrix(0,3,ref_all$i$nstates)
  tmp[1,ref$s$treated_pay]   <- 1
  tmp[2,intersect(ref$s$mdr,ref$s$treated_pay)]  <- 1
  tmp[3,ref$s$treated_pe]  <- 1
  agg$notif <-tmp
  
  tmp <-matrix(0,ref_all$i$nstates,ref_all$i$nstates)
  tmp[ref$s$treated_pay,c(ref$s$Ia, ref$s$Is, ref$s$E, ref$s$Dx)] <- 1
  tmp<-tmp*check
  sel$notif <- tmp - diag(diag(tmp))
  
  
  # Selectors mortality
  tmp <- matrix(0,1,ref_all$i$nstates)
  tmp[1,ref$s$infectious] <- 1
  agg$mort <-tmp
  
  # Preventive therapy
  tmp <- matrix(0,1,ref_all$i$nstates)
  tmp[1,intersect(ref$s$Pt,c(ref$s$pu, ref$s$pe))] <- 1
  agg$pt <-tmp
  
  tmp <-matrix(0,ref_all$i$nstates,ref_all$i$nstates)
  tmp[ref$s$Pt,] <- 1
  tmp<-tmp*check
  sel$pt<-tmp - diag(diag(tmp))
  
  
  
  # Selectors for false positive treatment (Tx)
  tmp <- matrix(0,1,ref_all$i$nstates)
  tmp[1,ref$s$notbtreated_pay] <- 1
  agg$txfp <-tmp
  
  tmp <-matrix(0,ref_all$i$nstates,ref_all$i$nstates)
  tmp[ref_all$s$STx,] <- 1
  tmp<-tmp*check
  sel$txfp<- tmp - diag(diag(tmp))
  
  #Aggregator for counting Dx : Will apply to matrices for smear, xpert and Xray
  tmp <- matrix(0,1,ref_all$i$nstates)
  tmp[1,c(ref$s$notbtreated_pay,ref$s$treated_pay)] <- 1
  agg$dx <-tmp
  
  
  #Selectors for Dx
  tmp <-matrix(0,ref_all$i$nstates,ref_all$i$nstates)
  tmp[c(ref_all$s$STx,ref$s$Tx,ref$s$Tx2), ] <- 1
  tmp<-tmp*check
  sel$dx<- tmp - diag(diag(tmp))
  
  
  
  # Selector for ACF 
  tmp <- matrix(0,1,ref_all$i$nstates)
  tmp[1,c(ref_all$s$S, ref$s$Ia, ref$s$Is, ref$s$E)] <- 1
  agg$acf <-tmp
  
  tmp <-matrix(0,ref_all$i$nstates,ref_all$i$nstates)
  tmp[c(ref$s$notbtreated_pay,ref$s$treated_pay), ]<- 1
  tmp<-tmp*check
  sel$acf <-tmp - diag(diag(tmp))
  
  
  # Aggregator for patient-month of tx
  tmp <- matrix(0,2,ref_all$i$nstates)
  tmp[1,intersect(c(ref$s$notbtreated_pay,ref$s$treated_pay),c(ref$s$Tx, ref_all$s$STx))] <- 1
  tmp[2,intersect(ref$s$treated_pay,ref$s$Tx2)] <- 1
  agg$pmo <-tmp
  
  #calib parameters
  xi <-list()
  
  xnames <-c('beta','beta_mdr','bRRhi','fast','kappa',
             'careseek_lo','careseek_hi','cs2','RRcs_rs',
             'symp_del','imm','selfcure','muTB',
             'pu','pse_base','Dxpu','Txinitpu',
             'ntpcov','cross')
  
  xnums  <- rep(1,length(xnames))
  lim <- 0
  for (ii in 1:length(xnames)){
    inds <- lim + (1:xnums[ii])
    xi[[xnames[ii]]] <- inds
    lim <- inds[length(inds)]
  }
  
  
  xi$nx <- lim
  bds<-matrix(0,length(xnames),2)
  bds[xi$beta,]             <- c(2, 18)
  bds[xi$beta_mdr,]         <- c(2, 18)
  bds[xi$bRRhi,]            <- c(1, 20)
  bds[xi$fast,]             <- c(0.012, 4)
  bds[xi$kappa,]            <- c(0.2, 5)
  bds[xi$careseek_lo,]      <- c(0.5, 12)
  bds[xi$careseek_hi,]      <- c(0.5, 12)
  bds[xi$cs2,]              <- c(1, 5)
  bds[xi$RRcs_rs,]          <- c(0, 1)   #Factor reducing careseeking amongs resp symptomatic
  bds[xi$symp_del,]         <- c(1.2, 12)
  bds[xi$imm,]              <- c(0.25, 0.75)
  bds[xi$selfcure,]         <- c(0.1, 0.22)
  bds[xi$muTB,]             <- c(0.1, 0.22)
  bds[xi$pu,]               <- c(0, 1)   #pu_bds
  bds[xi$pse_base,]         <- c(0, 1)
  bds[xi$Dxpu,]             <- c(0.7, 1)
  bds[xi$Txinitpu,]         <- c(0.6, 0.9)
  bds[xi$ntpcov,]           <- c(0, 1)
  bds[xi$cross,]            <- c(0, 1)
  
  
  #__________________________________________________________________________
  #  Prepare fixed parameters
  #__________________________________________________________________________
  
  
  r<-list()
  p<-list()
  
  #Epi
  
  r$fast_react    <- 0       # Ragonnet 2018_Epidemics
  r$slow_react    <- (0.0012 + 0.0012)/2  # Ragonnet 2018_Epidemics
  r$slow          <- 1.9710             # Ragonnet 2018_Epidemics
  r$careseek      <-c(0,0)
  r$careseeking2  <- c(0, 0, 0)
  r$access        <- 0
  p$kappa         <- 1
  r$selfcure  <- 0
  r$relapse   <- c(0.032, 0.14, 0.0015)
  r$mort_TB   <- 1/6
  p$Fast      <- c(0.035, 0.084, 0.087)
  p$imm       <- 0.5
  p$crossg    <- 0.3
  
  # Diagnosis stage
  p$strd_pr <-0
  if(strcmp(data_raw[location,'stdard_pr'],'good')){
    p$strd_pr<-1}
  else{
    p$strd_pr<-0.5
  }
  
  p$pu <- 0
  p$pse<- 0
  r$Dx <- 52
  p$Dx <- 0.83*c(1,p$strd_pr,p$strd_pr)
  r$actC  <-1
  p$Tx_init <- c(0.88, 0.75, 0.75)
  p$smear_sens<-0.8 #Smear test speccificity (Swai F 2011, BMC resreach notes)
  p$smear_spec<-0.94#Smear test speccificity (Swai F 2011, BMC resreach notes)
  p$xpert_sens<-0.9 #Smear test speccificity (Swai F 2011, BMC resreach notes)
  p$xpert_spec<-0.99#Smear test speccificity (Swai F 2011, BMC resreach notes)
  p$xray_sens<-0.9 #Xray test speccificity (Swai F 2011, BMC resreach notes)
  p$xray_spec<-0.5 #Xray test speccificity (Swai F 2011, BMC resreach notes)
  
  p$smear<-c(0.8 , 0, 0)   #prob of microscopy by sector
  p$xray <-c(0,1,1)          #prob of Xray as Dx upfront
  p$xpert_upf<-c(0.2,0,0)  #prob of xpert upfront
  p$xpert_fup<-c(0.1,0,0)  #prob of xpert follow-up aas confirmaion of smear
  
  # Treatment stage
  r$Tx <- 2
  p$pdef <- (1-data_raw[location,'tsr_fl'])*c(1, 1.75, 1.75)
  r$default <- r$Tx*p$pdef/(1-p$pdef)
  p$cure      <- c(1, 1, 1)
  
  #MDR
  r$MDR_acqu <- 0.05 # MDR acquisition
  p$SL_trans <-c(0.8, 0, 0) # Transfer from FL
  r$Tx2      <-0.5
  p$newmol   <- 0 #Fraction recieving new molecuiles for their SL treatment
  p$pdef2 <- (1-data_raw[location,'tsr_sl'])*c(1, 1.75, 1.75) 
  r$default2 <- r$Tx2*p$pdef2/(1-p$pdef2)
  p$cure2    <- c(0.5, 0, 0)
  
  
  # Demographic terms
  lex        <- 67    # Life expectancy
  r$mort     <- 1/lex #[(0.0253+	1.0000e-03)/2	(0.0065+	0.0933)/2];%[ 1/lex 1/(lex-15) 1/(lex-65)];    % Non-disease-related mortality
  p$growth   <- 0     #pars{location,'growth'};
  p$frac_pop <- c(0.0910+0.1824	, 0.6643+0.0623) #	Population fraction <15 2016;
  p$hi       <- data_raw[location,'size_hirisk'] # size of the vulnerable population
  r$wloss    <- 0.5 # Rate of 2 years fro untreated TBG to move from normal BMI to malnourished
  p$resp_symptomatic<-0.045 # Fraction Respiratoty symnptomatic (From kenya survey (eligible only by symptoms)(TB survey)
  #p$slum<-strcmp(data_raw[location,'setting'],'slum')
  p$popN<-data_raw[location,'popN']
  
  #intervs
  r$acf_asym<-c(0, 0)
  r$acf_sym <-c(0, 0)
  p$acfhi <-0.5
  p$xray_acf<-0
  p$xpert_acf<-1
  p$acf_k<-0 # Cumulative losses
  
  
  p$cfy_all<-0  #Switch parametyer to do contact trace
  p$cfy    <-0  #Switch parametyer to do contact trace
  r$nutri  <-0
  p$prevComm<-1 #TB prevalence in the community for CFY
  p$pteff  <-0.63
  p$ptreach<-0 # 0.1 for children, 1 for all population
  r$outpt <- 12/9 
  p$ptdef <- 1-0.75 
  r$ptdefault <- r$outpt*p$ptdef/(1-p$ptdef)
  r$ptrate<-0
  p$hcscreen<-1 
  
  
  # Unit costs (USD temporary costs )
  
  p$u_smear   <-10.87               #;cost.smear=rnu(u_smear,nruns);
  p$u_xpert   <-32.24               #cost.xpert=rnu(u_xpert,nruns);
  p$u_dst     <-20.03               #cost.dst=rnu(u_dst,nruns);
  p$u_xray    <-15               #cost.dst=rnu(u_dst,nruns);
  p$u_flmo    <- 33.91              #cost.flmoadu=rnu(u_flmoadu,nruns);
  p$u_slmo    <- 134.61             #cost.slmoadu=rnu(u_slmoadu,nruns);
  p$u_ipt     <-52.20               #cost.ipt=rnu(u_ipt,nruns);
  p$u_acf     <- ((30.69+1.93)) + 3 #cost.acf=rnu(u_acf,nruns);
  
  
  
  
  prm<-list(p=p, r=r, bds=t(bds),xnames=xnames )
  ref$i_all<-ref_all$i
  ref$s_all<-ref_all$s
  ref$xi<-xi
  
  return(list(prm=prm, ref=ref, lhd=lhd, sel=sel, agg=agg, gps=gps, datapoints=datapoints))
  
  
}
