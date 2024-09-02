# 1) Setting
library(rstan)
library(Rcpp)
library(RcppEigen)
library(StanHeaders)
library(BH)
library(inline)
library(plyr)
library(R.matlab)
library(loo)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

### Input variables

ChainNum  <- 1; # example, [1,2,3,4]
root <- paste("C:\\Users\\mh-hwang\\Dropbox\\i4Replication\\_DNCE",sep="") # example 

###

# setwd("/home/samba/Office/Users/minho/server/_Garcia/_Model_Fitting")
# setwd("/Volumes/samba/Office/Users/minho/server/_Garcia/_Model_Fitting")
setwd(paste(root,"\\git_hub\\stan",sep=""))


# args <- commandArgs(trailingOnly = TRUE);
# if (length(args)==0) {
#   stop("At least one argument must be supplied (input file).n", call.=FALSE)
# } else {
#   ChainNum  <- args[1];
# }
# show(ChainNum)

loadData <- paste(getwd(),"/Data/stanBehav_LE_choices.mat", sep = "");
allParticipantBehavior <- readMat(loadData);
Total_behavCol <-  allParticipantBehavior$stanAccum.choice;
Total_subjCol  <-  allParticipantBehavior$stanAccum.subject;
Total_ntrialCol<-  allParticipantBehavior$stanAccum.ntrial;
Total_conCol   <-  allParticipantBehavior$stanAccum.con;
Total_cfCol    <-  allParticipantBehavior$stanAccum.cf;
Total_outCol   <-  allParticipantBehavior$stanAccum.out;
Total_cfoutCol <-  allParticipantBehavior$stanAccum.cfout;
Total_nsubCol  <-  allParticipantBehavior$stanAccum.nsub;
Total_nconCol  <-  allParticipantBehavior$stanAccum.ncon;
Total_sessionCol  <-  allParticipantBehavior$stanAccum.session;
Total_expCol  <-  allParticipantBehavior$stanAccum.exp.num;


selected_exp = c(1, 2, 3, 4, 5, 6.1, 6.2, 7.1, 7.2, 8.1, 8.2, 9.1, 9.2);
for(exp_num in 1:length(selected_exp)){
  show(selected_exp[exp_num])
  tmp = Total_nsubCol[Total_expCol==selected_exp[exp_num]];
  NumSubject = tmp[1];
  
  Group_behavCol =Total_behavCol[Total_expCol==selected_exp[exp_num]];
  Group_subjCol  =Total_subjCol[Total_expCol==selected_exp[exp_num]];
  Group_ntrialCol=Total_ntrialCol[Total_expCol==selected_exp[exp_num]];
  Group_conCol   =Total_conCol[Total_expCol==selected_exp[exp_num]];
  Group_cfCol    =Total_cfCol[Total_expCol==selected_exp[exp_num]];
  Group_outCol   =Total_outCol[Total_expCol==selected_exp[exp_num]];
  Group_cfoutCol =Total_cfoutCol[Total_expCol==selected_exp[exp_num]];
  Group_nsubCol  =Total_nsubCol[Total_expCol==selected_exp[exp_num]];
  Group_nconCol  =Total_nconCol[Total_expCol==selected_exp[exp_num]];
  Group_sessionCol=Total_sessionCol[Total_expCol==selected_exp[exp_num]];
  Group_expCol   =Total_expCol[Total_expCol==selected_exp[exp_num]];
  
  NumObservation = length(Group_behavCol);
  RA_data=list("NumSubject"=NumSubject, "NumObservation"=NumObservation, "Group_subjCol" = Group_subjCol,
               "Group_behavCol"=as.vector(Group_behavCol), "Group_ntrialCol"=Group_ntrialCol,
               "Group_conCol"=Group_conCol, "Group_cfCol"=Group_cfCol,
               "Group_outCol"=Group_outCol,"Group_cfoutCol"=Group_cfoutCol,"Group_nsubCol"=Group_nsubCol,
               "Group_nconCol"=Group_nconCol, "Group_sessionCol"=Group_sessionCol, "Group_expCol"=Group_expCol,
               "Current_exp" = selected_exp[exp_num]);
  
  RA_fit=stan(file = paste(getwd(),"/Estim_RW_Model.stan",sep=""), data=RA_data, verbose=FALSE, warmup=2000, iter=5000, chains=1, init="random", control=list(max_treedepth=15, adapt_delta=0.99, stepsize=0.001))
  fname = paste("modelfitting_LE_fit_chain",ChainNum,"_exp_",selected_exp[exp_num],".rdata",sep="");
  save(RA_fit,file=fname)
}