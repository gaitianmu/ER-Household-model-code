##  Please set code_root variable properly. 
## For more information, please read the "Model-code-note"
code_root="/home/gaixin/"

library(BayesianTools)
library("corrplot")
library(readr)
library(cairoDevice)
##
source(paste0(code_root, "RCODE/fun_SEIRpred.R"))
source(paste0(code_root, "RCODE/fun_SEIRsimu_ode.R"))
source(paste0(code_root, "RCODE/Rfun_SEIRsimu_ode.R"))
source(paste0(code_root, "RCODE/fun_SEIRfitting.R"))
source(paste0(code_root, "RCODE/M2init_cond.R"))
source(paste0(code_root, "RCODE/correlationPlot_modified.R"))
source(paste0(code_root, "RCODE/fun_SEIRplot.R"))
source(paste0(code_root, "RCODE/fun_Findzero.R"))

##

init_sets_list=get_init_sets_list(alpha = 0.55)
probbase<- read.csv("/home/gaixin/data/probability-14.csv", row.names = 1)
prob<-probbase[1:52,]
qd<-prob[8:52]
set.seed(66)
SEIRfitting(init_sets_list, randomize_startValue = T,
            run_id = "S17-ode", output_ret = T, skip_MCMC=F)

