SEIRpred <- function(pars, 
                     init_settings) {
  tmp_ret=init_settings$var_trans_fun(pars)
  b_vec=tmp_ret[[1]]
  r_vec=tmp_ret[[2]]
  rh_vec=tmp_ret[[3]]
  
  bs_vec=tmp_ret[[4]]
  bl_vec=tmp_ret[[5]]
  stage_intervals=init_settings$stage_intervals
  n_stage=length(stage_intervals)
  ##
  Di <- init_settings$Di
  Dp <- init_settings$Dp
  De <- init_settings$De
  Dq_vec <- init_settings$Dq
  Dqs_vec <- init_settings$Dqs
  alpha <- init_settings$alpha
  Dh <- init_settings$Dh
  MS<-init_settings$MS
  ML<-init_settings$ML
  N <- init_settings$N
  flowN_vec <- init_settings$flowN
  init_states <- init_settings$init_states
  days_to_fit <- init_settings$days_to_fit
  update_func <- function(stage_pars, states_old) {
   
    b <- stage_pars[1]
    r <- stage_pars[2]
    rh <- stage_pars[3]
    Dq <- stage_pars[4]
    Dqs <- stage_pars[5]
    n <- stage_pars[6]    
    bs<-stage_pars[7]
    bl<-stage_pars[8]

    S <- states_old[1]
    E <- states_old[2]
    P <- states_old[3]
    I <- states_old[4]
    A <- states_old[5]
    AA<-states_old[6]
    H <- states_old[7]
    R <- states_old[8]
    SS <- states_old[9]
    ST<-states_old[10]
    ES <- states_old[11]
    PS <- states_old[12]
    IS <- states_old[13]
    AS <- states_old[14]
    AAS <- states_old[15]
    HS <- states_old[16]
    RS <- states_old[17]
    SL <- states_old[18]
    LT <- states_old[19]
    EL <- states_old[20]
    PL <- states_old[21]
    IL <- states_old[22]
    AL <- states_old[23]
    AAL <- states_old[24]
    HL <- states_old[25]
    RL <- states_old[26]
    
    
    
    ## S
    ## meaning S->EA, S->EI,S->, S->S
    pS_vec <- c(b * (alpha * (P+PS+PL) +alpha * (A+AS+AL)+alpha * (AA+AAS+AAL)) / N, b * (I+IS+IL) / N, n /N, 1 - b *(alpha * (P+PS+PL) +alpha *(A+AS+AL)+alpha * (AA+AAS+AAL)+ (I+IS+IL)) /N - n / N)
    sample_S <- pS_vec*S
    ## E
    ## meaning E->P, E->, E->E
    pE_vec <- c(1 / De, n / N, 1 - 1 / De - n / N)
    sample_E <- pE_vec*E
    ## P
    ## meaning P->I, P->A, P->, P->P
    pP_vec <- c(r / Dp, (1 - r) / Dp, n/N, 1 - 1 / Dp - n/N)
    sample_P <- pP_vec*P
    ## I
    ## meaning I->H, I->R, I->I
    pI_vec <- c(1 / Dq, 1 / Di, 1 - 1 / Dq - 1 / Di)
    sample_I<- pI_vec*I

    ## A
    ## meaning A->R, A->, A->A
    pA_vec <- c(1 / Di, n /N, 1 - 1 / Di - n /N)
    sample_A <-pA_vec*A
    pAA_vec <- c(0)
    sample_AA <-pA_vec*AA
    
    ## H
    ## meaning H->R, H->H
    pH_vec <- c(1 / Dh, 1 - 1 / Dh)
    sample_H <- pH_vec*H
    ## R
    ## meaning R->, R->R
    pR_vec <- c(n /N, 1 - n /N)
    sample_R <- pR_vec*R
    
    ###########SAMll####################################################################
    #SS#
    ##SS->ASP,SS->ISP,SS->,SS->SS
    pSS_vec <- c(MS*b * (alpha * (P+PS+PL) +alpha * (A+AS+AL) +alpha * (AA+AAS+AAL)) / N, MS*b *(I+IS+IL) /N, n /N, 1 - MS*b * (alpha * (P+PS+PL) +alpha * (A+AS+AL) +alpha * (AA+AAS+AAL)) / N-MS*b *(I+IS+IL) /N- n / N)  #(SL+EL+PL+ALe+ALl+ILe+ILl)

    sample_SS <-SS*pSS_vec
    ##ST->ASP,SL->ISP,SS->ASH,SS->ISH,SL->,SL->SL
    pST_vec <- c(bs * (alpha * PS +alpha * AS+alpha * AAS)/(ST+ES+PS+AS+AAS+IS),bs*IS/(ST+ES+PS+AS+AAS+IS), n /N, 1 - bs* (alpha * PS +alpha * AS+alpha * AAS+ IS ) /(ST+ES+PS+AS+AAS+IS) - n / N)  #(SL+EL+PL+ALe+ALl+ILe+ILl)

    sample_ST <-ST*pST_vec
    ##ES#
    ## meaning ES->PS, E->, ES->ES
    pES_vec <-  c(1 / De, n / N, 1 - 1 / De - n / N)

    sample_ES<-ES*pES_vec
   
    ####PS###
    ## meaning PS->IS, PS->AS,PS->AAS, PL->, PL->PL
    pPS_vec <- c(r / Dp, (1 - rh) / Dp,(rh-r)/Dp, n/N, 1 - 1 / Dp - n/N)
    sample_PS <-pPS_vec *PS
   
    ## IS
    ## meaning IS->HS, IS->RS, IS->IS
    pIS_vec <- c(1 / Dqs, 1 / Di, 1 - 1 / Dqs - 1 / Di )
    sample_IS <-pIS_vec*IS
    ## AS
    ## meaning AS->R, AS->, AS->AS
    pAS_vec <- c(1 / Di, n /N, 1 - 1 / Di - n /N )
    sample_AS<-pAS_vec*AS
    ## AAS
    ## meaning AAS->HS, AAS->RS, AIS->AIS
    pAAS_vec <- c(1 / Dqs, 1 / Di, 1 - 1 / Dqs - 1 / Di)
    sample_AAS <-pAAS_vec*AAS

    ## HS
    ## meaning HS->RS, HS->HS
    pHS_vec <- c(1 / Dh, 1 - 1 / Dh)
    sample_HS <-pHS_vec*HS
    ## RS
    ## meaning RS->, RS->RS
    pRS_vec <- c(n /N, 1 - n /N)
    sample_RS <-pRS_vec*RS
    ###########large#####################################################################################
    #SL#
    ##SL->ALP,SL->ILP,SL->,SL->SL
    pSL_vec <- c(ML*b * (alpha * (P+PS+PL) +alpha * (A+AS+AL) +alpha * (AA+AAS+AAL)) / N, ML*b *(I+IS+IL) /N, n /N, 1 - ML*b * (alpha * (P+PS+PL) +alpha * (A+AS+AL) +alpha * (AA+AAS+AAL)) / N-ML*b *(I+IS+IL) /N- n / N)  
    sample_SL <-SL*pSL_vec
    ##LT->ALP,SL->ILP,SL->ALH,SL->ILH,SL->,SL->SL
    pLT_vec <- c(bl *(alpha * PL +alpha * AL+alpha * AAL)/(LT+EL+PL+AL+AAL+IL),bl*IL/(LT+EL+PL+AL+AAL+IL), n /N, 1 - bl* (alpha * PL +alpha * AL+alpha * AAL+ IL ) /(LT+EL+PL+AL+AAL+IL) - n / N)  
    sample_LT <-LT*pLT_vec
    ##EL#
    ## meaning EL->PL, E->, EL->EL
    pEL_vec <-  c(1 / De, n / N, 1 - 1 / De - n / N)

    sample_EL<-EL*pEL_vec
    
    ####PL###
    ## meaning PL->IL, PL->AL,PL->AAL, PL->, PL->PL
    pPL_vec <- c(r / Dp, (1 - rh) / Dp,(rh-r)/Dp, n/N, 1 - 1 / Dp - n/N)
    sample_PL <-pPL_vec *PL
    
    ## IL
    ## meaning IL->HL, IL->RL, IL->IL
    pIL_vec <- c(1 / Dqs, 1 / Di, 1 - 1 / Dqs - 1 / Di )
    sample_IL <-pIL_vec*IL
    ## AL
    ## meaning AL->R, AL->, AL->AL
    pAL_vec <- c(1 / Di, n /N, 1 - 1 / Di - n /N )
    sample_AL<-pAL_vec*AL
    ## AAL
    ## meaning AAL->HL, AAL->RL, AAL->AAL
    pAAL_vec <- c(1 / Dqs, 1 / Di, 1 - 1 / Dqs - 1 / Di)
    sample_AAL <-pAAL_vec*AAL
    ## HL
    ## meaning HL->RL, HL->HL
    pHL_vec <- c(1 / Dh, 1 - 1 / Dh)
    sample_HL <-pHL_vec*HL
    ## RL
    ## meaning RL->, RL->RL
    pRL_vec <- c(n /N, 1 - n /N)
    sample_RL <-pRL_vec *RL
    
    
    #################################################################################

    S_new <- sample_S[4] + 0.06*n
    E_new <- sample_E[3] + sample_S[1] +sample_S[2]
    P_new <- sample_P[4] + sample_E[1]
    I_new <- sample_I[3] + sample_P[1]
    A_new <- sample_A[3] + sample_P[2]
    AA_new <- sample_AA[1]
    H_new <- sample_H[2] + sample_I[1]
    R_new <- sample_R[2] + sample_I[2] + sample_A[1] + sample_H[1]
    Onset_expect <- sample_P[1]
    expect_H<-sample_I[1]
    expect_A<-sample_P[2]
    expect_Ea<-sample_S[1]
    expect_Ei<-sample_S[2]
    ########LARGR###########################################################
    SL_new <- sample_SL[4] + 0.48*n
    LT_new <- (sample_SL[1]+sample_SL[2])*(ML-1)/ML+sample_LT[4]
    
    EL_new <- sample_EL[3] + sample_LT[1] +sample_LT[2]+ (sample_SL[1] +sample_SL[2])/ML
    PL_new <- sample_PL[5] + sample_EL[1]
    
    IL_new <- sample_IL[3] + sample_PL[1]
    AL_new <- sample_AL[3] + sample_PL[2]
    
    AAL_new <- sample_AAL[3] + sample_PL[3]
    
    
    HL_new <- sample_HL[2] + sample_IL[1] + sample_AAL[1]
    RL_new <- sample_RL[2] + sample_IL[2] + sample_AL[1]+ sample_AAL[2] + sample_HL[1]
    
    
    LOnset_expect <- sample_PL[1]+ sample_PL[3]
    
    expect_HL<-sample_IL[1]+sample_AAL[1]
    expect_AL<-sample_PL[2]
    
    expect_ELaP<-sample_SL[1]/ML
    expect_ELiP<-sample_SL[2]/ML
    expect_ELaH<-sample_LT[1]
    expect_ELiH<-sample_LT[2]
    
    #######SMALL############################################################################

    SS_new <- sample_SS[4] + 0.48*n
    ST_new <- (sample_SS[1]+sample_SS[2])*(MS-1)/MS+sample_ST[4]
    
    ES_new <- sample_ES[3] + sample_ST[1] +sample_ST[2]+ (sample_SS[1] +sample_SS[2])/MS
    PS_new <- sample_PS[5] + sample_ES[1]
    
    IS_new <- sample_IS[3] + sample_PS[1]
    AS_new <- sample_AS[3] + sample_PS[2]
    
    AAS_new <- sample_AAS[3] + sample_PS[3]
    
    
    HS_new <- sample_HS[2] + sample_IS[1] + sample_AAS[1]
    RS_new <- sample_RS[2] + sample_IS[2] + sample_AS[1]+ sample_AAS[2] + sample_HS[1]
    
    
    SOnset_expect <- sample_PS[1]+ sample_PS[3]
    
    expect_HS<-sample_IS[1]+sample_AAS[1]
    expect_AS<-sample_PS[2]
    
    expect_ESaP<-sample_SS[1]/MS
    expect_ESiP<-sample_SS[2]/MS
    expect_ESaH<-sample_ST[1]
    expect_ESiH<-sample_ST[2]
    
    ##
    
    
    return(c(S_new, E_new, P_new, I_new, A_new, AA_new, H_new, R_new,SS_new,ST_new,ES_new, PS_new, IS_new, AS_new,AAS_new, HS_new, RS_new,SL_new,LT_new, EL_new,PL_new, IL_new, AL_new,AAL_new, HL_new, RL_new, Onset_expect,expect_H,expect_A,expect_Ea,expect_Ei,SOnset_expect,expect_HS,expect_AS,expect_ESaP,expect_ESiP,expect_ESaH,expect_ESiH, LOnset_expect,expect_HL,expect_AL,expect_ELaP, expect_ELiP,expect_ELaH,expect_ELiH))
  }
  
  states_mat <- matrix(0, length(days_to_fit), length(init_states) + 20)

  states_mat[, 1] <- days_to_fit
  
  colnames(states_mat) <- c("time", "S", "E", "P", "I", "A","AA", "H", "R", "SS","ST", "ES", "PS", "IS", "AS","AAS", "HS", "RS","SL","LT", "EL", "PL", "IL", "AL","AAL", "HL", "RL", "Onset_expect","expect_H","expect_A","expect_Ea","expect_Ei", "SOnset_expect","expect_HS","expect_AS", "expect_ESaP","expect_ESiP","expect_ESaH","expect_ESiH","LOnset_expect","expect_HL","expect_AL","expect_ELaP", "expect_ELiP","expect_ELaH","expect_ELiH")
  
  stage_start <- c(1,15)               
  stage_end <- c(14,length(days_to_fit))
  
  myold_states <- init_states
  
  for (i_stage in 1:n_stage) {
    
    for (d in stage_intervals[[i_stage]][["start"]]:stage_intervals[[i_stage]][["end"]]) {
      
      dd<-d
      if(dd >= 45) { 
        dd <- 45
      }
      
      stage_pars_setings <- c(b = b_vec[i_stage], r = r_vec[1],rh = rh_vec[i_stage],Dq = Dq_vec[dd%/%5.00001+1],Dqs = Dqs_vec[dd%/%5.00001+1], n = flowN_vec[i_stage],bs=bs_vec[i_stage],bl=bl_vec[i_stage])
     
      states_mat[d, -1] <- update_func(stage_pars = stage_pars_setings, states_old = myold_states)
      myold_states <- states_mat[d, -1]
    }
  }
  
  states_mat[,"Onset_expect"]<-states_mat[,"Onset_expect"]*qd
  states_mat[,"SOnset_expect"]<-states_mat[,"SOnset_expect"]*qd
  states_mat[,"LOnset_expect"]<-states_mat[,"LOnset_expect"]*qd
  write.csv(states_mat[,"Onset_expect"],file = "/home/gaixin/data/of.csv",row.names = F)
  write.csv(states_mat[,"SOnset_expect"],file = "/home/gaixin/data/Sf.csv",row.names = F)
  write.csv(states_mat[,"LOnset_expect"],file = "/home/gaixin/data/Lf.csv",row.names = F)
  return(states_mat)
 
}

