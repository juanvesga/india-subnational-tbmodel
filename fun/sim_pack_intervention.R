sim_pack_intervention <-
  function (x0,
            prm,
            ref,
            sel,
            agg,
            gps) {
    
    itv <- prm$itv
    r <- prm$r
    p <- prm$p
    i <- ref$i
    s <- ref$s
    i_all = ref$i_all
    s_all <- ref$s_all
    xi <- ref$xi
    
    
    tmp <- allocate_parameters(x0, r, p, xi)
    r <- tmp$r
    p <- tmp$p
    
    
    #--- Baseline
    p0 <- p
    r0 <- r
    
    M0 <- make_model(p0, r0, i, i_all, s, gps)
    Mset <- list()
    Mset[[1]] <- M0
    
    
    #--- Package of interventions
    p1 <- p
    r1 <- r
    
    #FL TSR
    p1$pdef[1] <- itv$pdef
    r1$default <-  r1$Tx * p1$pdef / (1 - p1$pdef)
    
    #SL TSR
    p1$pdef2[1] <- itv$pdef2
    r1$default2 <- r1$Tx2 * p1$pdef2 / (1 - p1$pdef2)
    
    #Xprt DST
    p1$xpert_upf[1] <- itv$xpert_uf
    p1$xpert_fup[1] <- itv$xpert_fu
    
    #PSE
    p1$pse          <- itv$pse
    p1$xpert_upf[3] <- itv$xpert_pe
    p1$xpert_fup[3] <- itv$xpert_pe
    p1$pdef[3]      <- itv$pdef_pe
    p1$pdef2[3]     <- itv$pdef_pe
    
    
    #PT
    r1$ptrate    <- itv$ptreach * itv$hcscreen
    p1$ptdef     <- itv$ptdef
    r1$ptdefault <- r$outpt * p1$ptdef / (1 - p1$ptdef)
    
    # ACF
    r1$acf_sym[2]  <- itv$acfhi
    r1$acf_asym[2] <- itv$acfhi
    p1$acf_k       <- itv$acf_loss
    p1$xray_acf    <- itv$xray_acf
    p1$xpert_acf   <- itv$xpert_acf
    
    M1 <- make_model(p1, r1, i, i_all, s, gps)
    Mset[[2]] <- M1
    
    
    # Run model pre-2017
    obj       <-  function(x)
      get_objective(x, prm, ref, sel, agg, gps, c(), FALSE)
    
      
      runs <- return_output(obj, x0, ref)
      init <- unlist(runs$sfin)
    
    #run models
    fx <- goveqs_basis
    
    times  <- seq(2017, 2025, by = 1)           # time scale
    intT   <- 3
    t.interv   <- c(times[2], times[2] + intT)
    allsol <- c()
    
    parms0 <-
      list(
        agg = agg,
        sel = sel,
        s = s,
        p = p0,
        r = r0,
        i = i,
        i_all = i_all,
        M = M0
      )
    
    for (mi in 1:length(Mset)) {
      Mfinal <- Mset[[mi]]
      fx_scale <-
        function(t, state, Mfinal)
          scale_up(t, state, Mfinal, parms0, t.interv, fx)
      out <- as.data.frame(lsoda(
        y = as.numeric(init),
        times = times,
        func = fx_scale,
        parms = Mfinal
      ))   #
      soln <- out[, 2:ncol(out)]
      t   <- out[, 1]
      
      if (mi == 1) {
        allsol <- soln
      } else{
        allsol <-  abind(x = allsol, y = soln, along = 3)
        
      }
    }
    
    
    # Get Model outputs
    tmp <- allsol[, 1:i$nstates, ]
    tmp <- apply(tmp, c(1, 3), sum)
    pop <- (tmp[1:nrow(tmp) - 1, ] + tmp[2:nrow(tmp), ]) / 2
    
    inc_n <- diff(allsol[, i$aux$inc[1], ])
    inc_rate <- 1e5 * (inc_n / pop)
    inc_cum <- colSums(inc_n)
    inc_av_pr <- (1 - inc_cum[2] / inc_cum[1])
    inc_av_n <- (inc_cum[1] - inc_cum[2])* p$popN
    
    mort_n <- diff(allsol[, i$aux$mort, ])
    mort_rate <- 1e5 * (mort_n / pop)
    mort_cum <- colSums(mort_n)
    mort_av_pr <- (1 - mort_cum[2] / mort_cum[1])
    mort_av_n <- (mort_cum[1] - mort_cum[2])* p$popN
    
   
    # Count Units
    
    #Dx (semar, xpert, xray)
    smear_n <- diff(allsol[, i$aux$dx[1], ]) * p$popN
    xpert_n <- diff(allsol[, i$aux$dx[2], ]) * p$popN
    xray_n <- diff(allsol[, i$aux$dx[2], ]) * p$popN
    acf_n <- diff(allsol[, i$aux$acf[1], ]) * p$popN
    fl_pmo <- diff(allsol[, i$aux$pmo[1], ]) * p$popN * 12
    sl_pmo <- diff(allsol[, i$aux$pmo[2], ]) * p$popN * 12
    
    
  
    tmp <- colSums(smear_n)
    incr_sm <- (tmp[2] - tmp[1]) * p$u_smear
    
    tmp <- colSums(xpert_n)
    incr_xp <- (tmp[2] - tmp[1]) * p$u_xpert
    
    tmp <- colSums(xray_n)
    incr_xr <- (tmp[2] - tmp[1]) * p$u_xray
    
    
    tmp <- colSums(acf_n)
    incr_acf <- (tmp[2] - tmp[1]) * p$u_acf
    
    
    tmp <- colSums(fl_pmo)
    incr_fl <- (tmp[2] - tmp[1]) * p$u_flmo
    
    tmp <- colSums(sl_pmo)
    incr_sl <- (tmp[2] - tmp[1]) * p$u_slmo
    
    res <- list()
    
    res$inc_base <- as.numeric(inc_rate[,1])
    res$inc_itv <- as.numeric(inc_rate[,2])
    res$mort_base <- as.numeric(mort_rate[,1])
    res$mort_itv <- as.numeric(mort_rate[,2])
    res$inc_av_pr <- as.numeric(inc_av_pr)
    res$mort_av_pr <- as.numeric(mort_av_pr)
    res$inc_av_n <- as.numeric(inc_av_n)
    res$mort_av_n <- as.numeric(mort_av_n)
    res$ic_fl  <- as.numeric(incr_fl)
    res$ic_sl  <- as.numeric(incr_sl)
    res$ic_sm  <- as.numeric(incr_sm)
    res$ic_xp  <- as.numeric(incr_xp)
    res$ic_xr  <- as.numeric(incr_xr)
    res$icr_all <- as.numeric(incr_fl + incr_sl + incr_sm + incr_xp + incr_xr + incr_acf)
    
    
    
    return(res)
  }
