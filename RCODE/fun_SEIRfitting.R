library(cairoDevice)

default_pars_density <- function(pars) {
  d_vec <- rep(NA,8)
   d_vec[1] <- dunif(pars[1],0,0.55, log = T)
   d_vec[2] <- dunif(pars[2], 0,1, log = T)
   r12 <- pars[3]
   d_vec[3] <- dunif(r12, 0,0.85, log = T)
   d_vec[4]<-dunif(pars[4], 0,1, log = T) 
   d_vec[5] <- dunif(pars[5],0,1, log = T)
   d_vec[6] <- dunif(pars[6],0,1, log = T)
   d_vec[7] <- dunif(pars[7],0,1, log = T)
   d_vec[8] <- dunif(pars[8], 0,1, log = T)
 
  return(sum(d_vec))
}

default_pars_sampler <- function(n = 1) {
  s_vec <- matrix(NA, n,8)
    s_vec[, 1] <- runif(n, 0,0.55) 
    s_vec[, 2] <- runif(n, 0,1) 
    r12<-runif(n, 0,0.85) 
    s_vec[, 3] <- r12
    s_vec[, 4] <-runif(n, 0,1) 
    s_vec[, 5] <- runif(n,0,1)
    s_vec[,6] <- runif(n,0,1)
    s_vec[, 7] <- runif(n,0,1)
    s_vec[, 8] <- runif(n,0,1)
    
  return(s_vec)
}


SEIRfitting=function(init_sets_list, 
                     randomize_startValue=F, 
                     startValue=NA, 
                     output_ret=T, 
                     run_id=0, 
                     skip_MCMC=F, 
                     panel_B_R_ylim=4,
                     plot_combined_fig=T,
                     pars_density=default_pars_density,
                     pars_sampler=default_pars_sampler,
                     pars_name=c("b12", "b3","r12","rh","bs1","bs2","bl1","bl2"),#
                     calc_clearance=T,
                     n_burn_in=16000,#16000,200000
                     n_iterations=200000) {
  if (randomize_startValue & !is.na(startValue)) {
    print("startValue will be ignored since randomize_startValue is set to TRUE!")
  } else if (!randomize_startValue & is.na(startValue)) {
    print("Please specify a startValue since you have set randomize_startValue to FALSE! Exiting!")
    q(save="no")
  }
  
  onset_obs <-init_sets_list$daily_new_case
  sonset_obs <- init_sets_list$Sdaily_new_case     #########small
  lonset_obs <- init_sets_list$Ldaily_new_case     #########large
  all_obs=init_sets_list$daily_new_case_all 
  init_states <- init_sets_list$init_states
  n_pars = length(pars_name)
  n_stage = length(init_sets_list$stage_intervals)
  
  ################################################################################################################
  ## likelihood function
  loglh_func <- function(pars){
    Aypred <- SEIRpred(pars, init_settings = init_sets_list)
   
    ypred <- Aypred[, "Onset_expect"]
    sypred <- Aypred[, "SOnset_expect"]
    lypred <- Aypred[, "LOnset_expect"]
    
    suppressWarnings(p <-dpois(onset_obs, ypred)*dpois(lonset_obs,lypred)*dpois(sonset_obs,sypred))
    if(any(p == 0) || any(is.nan(p))){
      logL <- -Inf
    }else{
      
      logL <- sum(log(p))
    }
    return(logL)
    print(logL)
  }
  
  pars_prior <- createPrior(density = pars_density, sampler = pars_sampler, 
                            lower = init_sets_list$par_lower, upper = init_sets_list$par_upper)
  
  if (!skip_MCMC) {
    bayesSEIR <- createBayesianSetup(loglh_func, prior = pars_prior)
  
    if (randomize_startValue) {  
      startValue=pars_sampler()
      while (is.infinite(loglh_func(startValue))) {
        startValue=pars_sampler()
      }
    }
    
   
    mh_settings = list(startValue = startValue,optimize = T, 
                       adapt = T, DRlevels = 2, iterations = n_iterations, thin = 10)
    mh_out <- runMCMC(bayesianSetup = bayesSEIR, sampler = "Metropolis", settings = mh_settings)
    mcmc_pars_estimate <- getSample(mh_out,start = n_burn_in+2, thin = 1) 
    mpe <- getSample(mh_out,parametersOnly =F,  start = n_burn_in+2, thin = 1) 
    mpelike<-mpe[,10]
    #write.csv(mpelike,file = "mpelike.csv",row.names = F)
    #write.csv(mpe,file = "ampe.csv",row.names = F)
    lax<-max(mpelike)
    mpe=data.frame(mpe)
    lamean<-mean(mpelike)
    print(lax)
    print(lamean)
    
    mcmc_pars_estimate <- round(mcmc_pars_estimate, 3)
    
    
    colnames(mcmc_pars_estimate) <- pars_name 
  
    if (output_ret) {
      write.table(mcmc_pars_estimate, paste0("/home/gaixin/output/pars_est_run_",run_id,".txt"), quote = F, row.names = F, sep = "\t")
    }
  } else {
    mcmc_pars_estimate = read.table(paste0("/home/gaixin/output/pars_est_run_",run_id,".txt"), header = T)
    pars_name = names(mcmc_pars_estimate)
  }
  
  summary_string=paste0(paste(pars_name, collapse = ","), "\n")
 
  par_str=list()
  for (i_par in 1:n_pars) {
    par_str[[i_par]]=paste0(round(mean(mcmc_pars_estimate[,i_par]),2), " (",
                            round(quantile(mcmc_pars_estimate[,i_par],0.025),2)," - " , 
                            round(quantile(mcmc_pars_estimate[,i_par],0.975),2), ")")
  }

  summary_string = paste0(summary_string, paste(par_str,collapse = ", "),"\n\n")
  
  mcmc_pars_estimate1<-mcmc_pars_estimate
  mcmc_pars_estimate1[,4]<-mcmc_pars_estimate1[,3]+(1-mcmc_pars_estimate1[,3])* mcmc_pars_estimate1[,4]
  mcmc_pars_estimate1[,6]<-mcmc_pars_estimate1[,5]* mcmc_pars_estimate1[,6]
  mcmc_pars_estimate1[,8]<-mcmc_pars_estimate1[,7]* mcmc_pars_estimate1[,8]
  for (i_par in 1:n_pars) {
    par_str[[i_par]]=paste0(round(mean(mcmc_pars_estimate1[,i_par]),2), " (",
                            round(quantile(mcmc_pars_estimate1[,i_par],0.025),2)," - " , 
                            round(quantile(mcmc_pars_estimate1[,i_par],0.975),2), ")")
  }
  summary_string = paste0(summary_string, paste(par_str,collapse = ", "),"\n\n")
  
  
  if (calc_clearance) {
    clearance_date = Findzero(mcmc_pars_estimate,file_name = run_id, init_sets_list)
    
    summary_string = paste0(summary_string, paste(names(clearance_date), collapse = ", "))
    
    summary_string = paste0(summary_string, "\n", paste(clearance_date, collapse = ", "), "\n")
  }
  
 
  
  write_file(summary_string, paste0("/home/gaixin/output/summary_run_",run_id,".txt"))
  
  cairo_pdf(paste0("/home/gaixin/output/par_cor_run_",run_id,".pdf"),width=10,height=10)
  correlationPlot_modified(mcmc_pars_estimate, scaleCorText = F)
  dev.off()
  
  png(paste0("/home/gaixin/output/par_hist_run_",run_id,".png"))
  par(mfrow = c(4, 3))
  for(i in 1:n_pars) {
    hist(mcmc_pars_estimate[, i], xlab = pars_name[i], main = "", col = "red")
    rm(i)
  }
  dev.off()
  
  png(paste0("/home/gaixin/output/par_traj_run_",run_id,".png"), width=1000, height=500)
  par(mfrow = c(4, 3))
  for(i in 1:n_pars) {
    plot(1:nrow(mcmc_pars_estimate), mcmc_pars_estimate[, i], ylab = pars_name[i], xlab = "iter", main = "", type = "l")
    rm(i)
  }
  dev.off()
  
 if (plot_combined_fig) {
    SEIRplot(pars_estimate = mcmc_pars_estimate, file_name = run_id, init_settings = init_sets_list, panel_B_R_ylim = panel_B_R_ylim)
  }
  
 
  
}
