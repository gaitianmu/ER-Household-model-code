Findzero <- function(pars_estimate,file_name, init_settings) {
  idx_to_date <- function(d) {
    d <- as.numeric(d)
    mydate <- paste(rep(c("Jan/20", "Feb/20", "Mar/20", "Apr/20", "May/20", "Jun/20", "Jul/20", "Aug/20", "Sep/20","Oct/20","Nov/20","Dec/20","Jan/21", "Feb/21", "Mar/21", "Apr/21", "May/21", "Jun/21", "Jul/21", "Aug/21", "Sep/21","Oct/21","Nov/21","Dec/21"),c(21, 29, 31, 30, 31, 30, 31, 31, 30,31,30,31,31, 28, 31, 30, 31, 30, 31, 31, 30,31,30,31)),
                    c(1:21, 1:29, 1:31, 1:30, 1:31, 1:30, 1:31, 1:31, 1:30,1:31,1:30,1:31,1:31, 1:28, 1:31, 1:30, 1:31, 1:30, 1:31, 1:31, 1:30,1:31,1:30,1:31))
    return(mydate[d])
  }
  ##
  init_settings$days_to_fit <- 1:700
  init_settings$Di<-18
  init_settings$Dp<-0.8*init_settings$Dp
  init_settings$De<-0.8*init_settings$De
  init_settings$alpha<-1
  Dq_vec <- init_settings$Dq
  Di <- init_settings$Di
  Dp <- init_settings$Dp
  De <- init_settings$De
  alpha <- init_settings$alpha
  ##
  est_mat <- apply(pars_estimate, 1, function(x) SEIRsimu(pars = x, init_settings = init_settings, num_periods = 2)[, c("E", "P", "I", "A","ES", "PS", "IS", "AS","AAS","EL", "PL", "IL", "AL","AAL", "Onset_expect", "SOnset_expect", "LOnset_expect","expect_ESiH","ST","STC","STJ","BL","LTC","LTJ")])
  est_mat_I <-round(est_mat[9801:10500, ],5)+ round(est_mat[10501:11200, ],5)+ round(est_mat[11201:11900, ],5) 
  #est_mat_EPIA <-  round(est_mat[201:400, ],5) + round(est_mat[401:600, ],5) + round(est_mat[601:800, ],5)+round(est_mat[1001:1200, ],5)+round(est_mat[1201:1400, ],5)+round(est_mat[1401:1600, ],5)+round(est_mat[1601:1800, ],5)+round(est_mat[2001:2200, ],5)+round(est_mat[2201:2400, ],5)+round(est_mat[2401:2600, ],5)+round(est_mat[2601:2800, ],5)
  est_mat_EPIA <- round(est_mat[1:700, ],5) + round(est_mat[701:1400, ],5) + round(est_mat[1401:2100, ],5) + round(est_mat[2101:2800, ],5)+round(est_mat[2801:3500, ],5)+round(est_mat[3501:4200, ],5)+round(est_mat[4201:4900, ],5)+round(est_mat[4901:5600, ],5)+round(est_mat[5601:6300, ],5)+round(est_mat[6301:7000, ],5)+round(est_mat[7001:7700, ],5)+round(est_mat[7701:8400, ],5)+round(est_mat[8401:9100, ],5)+round(est_mat[9101:9800, ],5)
  
  
  init_settings$days_to_fit <- 1:80
  
  NSS_mat <- matrix(0, 45,4000)
  NAS_mat <- matrix(0, 45, 4000)
  NA_mat <- matrix(0, 45, 4000)
  NSL_mat <- matrix(0, 45,4000)
  NASL_mat <- matrix(0, 45, 4000)
  NAL_mat <- matrix(0, 45, 4000)
  
  for  (d in 1:45){
    est_mat_rr<-apply(pars_estimate, 1, function(x) SEIRsimuR(pars = x, init_settings = init_settings, num_periods = 2,t=d)[, c("saf","ssf","laf","lsf")])
    e_saf <-  est_mat_rr[1:80,]
    e_ssf <- est_mat_rr[81:160,]
    e_laf <- est_mat_rr[161:240,]
    e_lsf <- est_mat_rr[241:320,]
  
    Dq = Dq_vec[d%/%5.00001+1]
    for (k in (round(d+De)+1):(round(d+De+Dp))) {
      
      NSS_mat[d,]<-NSS_mat[d,]+e_saf[k,]
    }
    
    for (j in (round(d+De+Dp)+1):(round(d+De+Dp+(1/Di+1/Dq)^(-1)))) {
      
      NSS_mat[d,]<-NSS_mat[d,]+e_ssf[j,]
    }
    for (p in (round(d+De)+1):(round(d+De+Dp+(1/Di+1/Dq)^(-1)))) {
      
      NAS_mat[d,]<-NAS_mat[d,]+e_saf[p,]
    }
    for (v in (round(d+De)+1):(round(d+De+Dp+Di))) {
      
      NA_mat[d,]<-NA_mat[d,]+e_saf[v,]
    }
    
    for (k in (round(d+De)+1):(round(d+De+Dp))) {
      
      NSL_mat[d,]<-NSL_mat[d,]+e_laf[k,]
    }
    
    for (j in (round(d+De+Dp)+1):(round(d+De+Dp+(1/Di+1/Dq)^(-1)))) {
      
      NSL_mat[d,]<-NSL_mat[d,]+e_lsf[j,]
    }
    for (p in (round(d+De)+1):(round(d+De+Dp+(1/Di+1/Dq)^(-1)))) {
      
      NASL_mat[d,]<-NASL_mat[d,]+e_laf[p,]
    }
    
    for (v in (round(d+De)+1):(round(d+De+Dp+Di))) {
      
      NAL_mat[d,]<-NAL_mat[d,]+e_laf[v,]
    }
    
    
    
    
    
  }
  

  init_settings$days_to_fit <- 1:45
  e_fix<-apply(pars_estimate, 1, function(x) SEIRsimu(pars = x, init_settings = init_settings, num_periods = 2)[, c("ps","pl","r","rh","b")])
  e_ps <- e_fix[1:45,]
  e_pl <-  e_fix[46:90,]
  r <- e_fix[91:135,]
  rh <- e_fix[136:180,]
  b <-e_fix[181:225,]
  
  est_Rs<-matrix(0, 45, 4000)
  for  (d in 1:45){
  
    est_Rs[d,]<- e_ps[d,]*(r[d,]* NSS_mat[d,]+(rh[d,]-r[d,])*NAS_mat[d,]+(1-rh[d,])*NA_mat[d,])
    
  }
  
 
  
  est_Rl<-matrix(0, 45, 4000)
  
  
  
  for  (d in 1:45){
   
    est_Rl[d,]<- e_pl[d,]*(r[d,]* NSL_mat[d,]+(rh[d,]-r[d,])*NASL_mat[d,]+(1-rh[d,])*NAL_mat[d,])
    
  }
  
  est_Rhh<-est_Rl+est_Rs
  
  
 
  est_Rcc<-matrix(0, 45, 4000)
  for  (d in 1:45){
   
    Dq = Dq_vec[d%/%5.00001+1]
    
    est_Rcc[d,]<-(e_ps[d,]+e_pl[d,])*( r[d,]*b[d,]*(Dp*alpha+(1/Di+1/Dq)^(-1)) + (rh[d,]-r[d,])*alpha*b[d,]*(Dp+(1/Di+1/Dq)^(-1)) + (1-rh[d,])*alpha*b[d,]*(Dp+Di) ) +  (1-e_ps[d,]-e_pl[d,])*(  r[d,]*b[d,]*(Dp*alpha+(1/Di+1/Dq)^(-1)) +(1-r[d,])*alpha*b[d,]*(Dp+Di)  ) 
    
  }
  print(round(apply(est_Rcc, 1, mean), 3))
  est_Rtt=est_Rhh+est_Rcc
  print("######################### R #####################")
 
   Rs_mean <- round(apply(est_Rs, 1, mean), 3)
  Rs_up <- round(apply(est_Rs, 1, function(x) quantile(x, 0.975)), 3)
  Rs_low <- round(apply(est_Rs, 1, function(x) quantile(x, 0.025)), 3)
  
  Rl_mean <- round(apply(est_Rl, 1, mean), 3)
  Rl_up <- round(apply(est_Rl, 1, function(x) quantile(x, 0.975)), 3)
  Rl_low <- round(apply(est_Rl, 1, function(x) quantile(x, 0.025)), 3)
  
  Rhh_mean <- round(apply(est_Rhh, 1, mean), 3)
  Rhh_up <- round(apply(est_Rhh, 1, function(x) quantile(x, 0.975)), 3)
  Rhh_low <- round(apply(est_Rhh, 1, function(x) quantile(x, 0.025)), 3)
  
  Rcc_mean <- round(apply(est_Rcc, 1, mean), 3)
  Rcc_up <- round(apply(est_Rcc, 1, function(x) quantile(x, 0.975)), 3)
  Rcc_low <- round(apply(est_Rcc, 1, function(x) quantile(x, 0.025)), 3)
  
  Rt_mean <- round(apply(est_Rtt, 1, mean), 3)
  Rt_up <- round(apply(est_Rtt, 1, function(x) quantile(x, 0.975)), 3)
  Rt_low <- round(apply(est_Rtt, 1, function(x) quantile(x, 0.025)), 3)
  
  
  print("------------ Rhs----------")
  print(Rs_mean)
  print(Rs_up)
  print(Rs_low)
  
  print("------------ Rhl----------")
  print(Rl_mean)
  print(Rl_up)
  print(Rl_low)
  
  print("------------ Rh----------")
  print(Rhh_mean)
  print(Rhh_up)
  print(Rhh_low)
  
  print("------------ Rc----------")
  print(Rcc_mean)
  print(Rcc_up)
  print(Rcc_low)
  
  print("------------ R----------")
  
  print(Rt_mean)
  print(Rt_up)
  print(Rt_low)
  
  RR<-cbind(Rs_mean,Rs_up,Rs_low,Rl_mean,Rl_up,Rl_low,Rhh_mean,Rhh_up,Rhh_low,Rcc_mean,Rcc_up,Rcc_low,Rt_mean,Rt_up,Rt_low)
  write.csv(RR, paste0("../RT/Rt_all_",file_name,".csv"), quote = F, row.names = F)
  

  print(est_mat_I[,2])
  print(est_mat_EPIA[,2])
  ## I = 0
  #est <- apply(est_mat_I, 2, function(x) which(x == 0)[1])
  est <- apply(est_mat_I, 2,function(x) tail(which(x != 0), 1) + 1)
  #est <- apply(est_mat_I, 2,function(x) tail(which(x >= 2.1), 1) + 1)
  #print(est_mat_I)
  lowdate <- idx_to_date(round(quantile(est, 0.025), 0))
  #meandate <- idx_to_date(round(mean(est), 0))
  meandate <- idx_to_date(round(quantile(est, 0.5), 0))
  upperdate <- idx_to_date(round(quantile(est, 0.975), 0))
  intdate1 <- paste(meandate, " (", lowdate, " to ", upperdate,")", sep = "")
  rm(est)
  # E + P + I+ A = 0
  est <- apply(est_mat_EPIA, 2, function(x) tail(which(x!= 0), 1) + 1)
  lowdate <- idx_to_date(round(quantile(est, 0.025), 0))#0.025
  #meandate <- idx_to_date(round(mean(est), 0))
  meandate <- idx_to_date(round(quantile(est, 0.5), 0))
  upperdate <- idx_to_date(round(quantile(est, 0.975), 0))#0.975
  intdate2 <- paste(meandate, " (", lowdate, " to ", upperdate,")", sep = "")
  rm(est)
  
  estp <- apply(est_mat_EPIA, 2, which.max)
  #estp <- apply(est_mat_EPIA, 2, function(x) which.max)
  estact<- as.matrix(apply(est_mat_EPIA, 2, function(x) max(x)))
  upestact<-apply(estact, 2, function(x) quantile(x, 0.975))
  lowestact<-apply(estact, 2, function(x) quantile(x, 0.025))
  meanestact<-apply( estact, 2, mean)
  print("---------EPIA PEAK START------------")
  #print(estact)
  print(lowestact)
  print(meanestact)
  print(upestact)
  
  print("---------EPIA PEAK END------------")
  
  ptime <- 1:45
  mydate <- c(paste("Jan", 11:31), paste("Feb", 1:29), paste("Mar", 1:31), paste("Apr", 1:30), paste("May", 1:31), paste("Jun", 1:14))
  
  pdf(paste0("../output/RFigure_", file_name, ".pdf"), width = 16, height = 5)
  par(mfrow = c(1, 3), mar = c(4, 5, 3, 2))
  
  plot(ptime, Rt_mean, ylim = c(0, 1.1*max(Rt_up[1:45])), xlab = "", ylab = "", type = "p", col = "white", pch = 16, xaxt="n", cex = 0.5)
  mtext("Date (2020)", side = 1, line  = 3, cex = 1.01)
  mtext("Effective Reproductive Number", side = 2, line = 3, cex = 1.01)
  axis(1, at = c(1,15,45), labels = mydate[c(1,15,45)])
  abline(v = c(15), lty = 3, lwd = 3, col = "red")
  polygon(c(ptime[1:45], rev(ptime[1:45])), c(Rt_up[1:45], rev(Rt_low[1:45])), col = "#CEE6F4FF", border = NA)
  lines(ptime,Rt_mean,col = "#0072B5FF", lwd = 2)
  points(ptime[1:45], Rt_mean, col = "#0072B5FF", pch = 16, cex = 1.5)
  text(par()$usr[1] - (par()$usr[2] -par()$usr[1]) * 0.12, par()$usr[4] + (par()$usr[4] - par()$usr[3]) * 0.06, labels = "R", xpd = T, cex = 2)
  #########RH########
  plot(ptime, Rhh_mean, ylim = c(0, 1.1*max(Rhh_up[1:45])), xlab = "", ylab = "", type = "p", col = "white", pch = 16, xaxt="n", cex = 0.5)
  mtext("Date (2020)", side = 1, line  = 3, cex = 1.01)
  mtext("Effective Household Reproductive Number", side = 2, line = 3, cex = 1.01)
  axis(1, at = c(1,15,45), labels = mydate[c(1,15,45)])
  abline(v = c(15), lty = 3, lwd = 3, col = "red")
  polygon(c(ptime[1:45], rev(ptime[1:45])), c(Rhh_up[1:45], rev(Rhh_low[1:45])), col ="#F39B7FB2", border = NA)
  lines(ptime,Rhh_mean,col = "#BC3C29FF", lwd = 2)
  points(ptime[1:45], Rhh_mean, col = "#BC3C29FF", pch = 16, cex = 1.5)
  text(par()$usr[1] - (par()$usr[2] -par()$usr[1]) * 0.12, par()$usr[4] + (par()$usr[4] - par()$usr[3]) * 0.06, labels = expression("R"["H"]), xpd = T, cex = 2)
  ############RC#########
  plot(ptime, Rcc_mean, ylim = c(0, 1.1*max(Rcc_up[1:45])), xlab = "", ylab = "", type = "p", col = "white", pch = 16, xaxt="n", cex = 0.5)
  mtext("Date (2020)", side = 1, line  = 3, cex = 1.01)
  mtext("Effective Community Reproductive Number", side = 2, line = 3, cex = 1.01)
  axis(1, at = c(1,15,45), labels = mydate[c(1,15,45)])
  abline(v = c(15), lty = 3, lwd = 3, col = "red")
  polygon(c(ptime[1:45], rev(ptime[1:45])), c(Rcc_up[1:45], rev(Rcc_low[1:45])), col ="#FFDC91FF", border = NA)
  #polygon(c(ptime[1:45], rev(ptime[1:45])), c(Rcc_up[1:45], rev(Rcc_low[1:45])), col ="#CAFF70", border = NA)
  lines(ptime,Rcc_mean,col = "#E18727FF", lwd = 2)
  #lines(ptime,Rcc_mean,col = "#20854EFF", lwd = 2)
  points(ptime[1:45], Rcc_mean, col ="#E18727FF", pch = 16, cex = 1.5)
  #points(ptime[1:45], Rcc_mean, col ="#20854EFF", pch = 16, cex = 1.5)
  text(par()$usr[1] - (par()$usr[2] -par()$usr[1]) * 0.12, par()$usr[4] + (par()$usr[4] - par()$usr[3]) * 0.06, labels = expression("R"["C"]), xpd = T, cex = 2)
  
  dev.off()
  
  
  print("---------EPIA 2.24 start------------")
  est_end<-est_mat_EPIA[45,]
  #print(est_end)
  esend<-as.matrix(est_end)
  #print(esend)
  upend<-apply(esend, 2, function(x) quantile(x, 0.975))
  lowend<-apply(esend, 2, function(x) quantile(x, 0.025))
  meanend<-apply( esend, 2, mean)
  print( upend)
  print( meanend)
  print( lowend)
  print("---------EPIA 2.24 end------------")
  #est <- apply(est_mat_EPIA, 2, function(x) tail(which(x>= 2.1), 1) + 1)
  lowdate <- idx_to_date(round(quantile(estp, 0.025), 0))
  #meandate <- idx_to_date(round(mean(est), 0))
  meandate <- idx_to_date(round(quantile(estp, 0.5), 0))
  upperdate <- idx_to_date(round(quantile(estp, 0.975), 0))
  intdate3 <- paste(meandate, " (", lowdate, " to ", upperdate,")", sep = "")
  
  
  #
  intdate <- c(intdate1, intdate2,intdate3)
  names(intdate) <- c("I=0", "E+P+I+A=0","Active Infection Peak")
  return(intdate)
  
   
}







