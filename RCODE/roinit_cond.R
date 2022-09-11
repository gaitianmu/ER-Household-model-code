generate_init_condi <- function(
                                Di = 10,
                                
                                Dp = 2.3,
                                De = 2.9,#2.9
                                Dq = c(16,12,9,7,6,5,5,3,2),#c(16,12,9,7,6,5,5,4,2)
                                Dqs= c(16,12,9,7,6,5,5,3,2),
                                
                                alpha = 0.55,#0.55
                                Dh = 30,
                                MS=2.5,
                                ML=5.1,
                                N = 4000000,
                                NS = 1840000,
                                NL = 1920000,#NL = 1920000,
                                flowN = c(0,0)
) {
  
  stopifnot(Di>=0  & Dp>=0 & De>=0 & all(Dq>=0)& all(Dqs>=0) & alpha>=0 & alpha<=1 & Dh>=0 & N>=0  & NS>=0 & NL>=0 & all(flowN>=0))
  
  
  R0 <- 0
  RS0<-3
  RL0<-3
  H0 <-0
  HS0<-3
  HL0<-3
  realData_all <- read.csv("/home/gaixin/data/roCovid19CasesYC.csv", row.names = 1)  
  realData <- realData_all 
  daily_new_case <- realData[8:52, 2]
  Sdaily_new_case <- realData[8:52, 3]
  Ldaily_new_case <- realData[8:52, 4]
  daily_new_case_all <- realData[8:52, 1]
  
  US<-1   #1,2,3
  UL<-9   #6,9,12
  
  E0<-US
  ES0<-UL
  
  EL0<-UL
 
  P0 <-US  
  PS0<-UL
  
  PL0<-UL
  
  I0 <-US
  IS0 <-UL
  IL0 <-UL
 
  A0<-0
  
  AA0<-0
  AS0<-US
  AAS0<-0
  AL0<-US
  AAL0<-0
  
  S0 <- 240000
  SS0 <- 1840000
  SL0 <- 1920000
  ST0<-(RS0+HS0+ES0+PS0+IS0+AS0+AAS0)*1.5        #51
  LT0<-(RL0+HL0+EL0+PL0+IL0+AL0+AAL0)*4.1        #139
  
  init_states <- round(c(S = S0, E = E0, P = P0, I = I0, A=A0,AA=AA0, H = H0, R = R0,SS = SS0,ST=ST0, ES = ES0, PS = PS0, IS = IS0,  AS=AS0, AAS = AAS0, HS = HS0, RS = RS0,SL = SL0,LT=LT0, EL = EL0, PL = PL0, IL = IL0,AL=AL0, AAL = AAL0, HL = HL0, RL = RL0), 0)
  print(I0)
 
  transform_var_main_5stage=function(pars) {
    b_vec <- pars[1:2]
    b_vec <- c(b_vec[1], b_vec[2])
    ##
    r12 <- pars[3]
    
    r_vec <- c(r12)
   
    rh<-r12+(1-r12)*pars[4]
    rh_vec<-c(rh,rh)
    bs1_vec<-pars[5]
   bs2_vec<-bs1_vec*pars[6]
    
    
    bs_vec<-c(bs1_vec,bs2_vec)
    bl1_vec<-pars[7]
    
    bl2_vec<-bl1_vec*pars[8]
    
    
    bl_vec<-c(bl1_vec,bl2_vec)
   
    return(list(b_vec, r_vec,rh_vec,bs_vec,bl_vec))
  }
  
  return(list(Di=Di,
             
              Dp=Dp,
              De=De,
              Dq=Dq,
              Dqs=Dqs,
              alpha=alpha,
              Dh=Dh,
              MS=MS,
              ML=ML,
              N=N,
              NS=NS,
              NL=NL,
              flowN=flowN,
              daily_new_case = daily_new_case, 
              Sdaily_new_case = Sdaily_new_case,
              Ldaily_new_case = Ldaily_new_case,
              daily_new_case_all = daily_new_case_all, 
              init_states = init_states,
              days_to_fit=1:45,
              stage_intervals=list(
                c(start=1, end=14),
                
                c(start=15, end=45)
                
              ),
              var_trans_fun=transform_var_main_5stage,
              par_lower = c(b12 = 0, b3 = 0,r12 = 0,rh=0,bs1=0,bs2=0,bl1=0,bl2=0),
              par_upper = c(b12 = 1, b3 =1, r12 = 1,rh=1, bs1=1,bs2=1,bl1=1,bl2=1)))
 
}


get_init_sets_list = generate_init_condi


