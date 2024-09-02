rm(list=ls())
ls()
library(rstan)
library(Rcpp)
library(RcppEigen)
library(StanHeaders)
library(BH)
library(inline)
library(plyr)
library(R.matlab)
library(loo)
library(shinystan)

# setwd("/home/samba/Office/Users/minho/server/_Garcia/_Model_Fitting")
# setwd("/Volumes/samba/Office/Users/minho/server/_Garcia/_Model_Fitting")
setwd("Z:/Office/Users/minho/server/_Garcia/_Model_Fitting")

selected_exp = c(1, 2, 3, 4, 5, 6.1, 6.2, 7.1, 7.2, 8.1, 8.2, 9.1, 9.2);
for(exp_num in 1:length(selected_exp)){
  
  models <- paste("modelfitting_ES_fit_exp_",selected_exp[exp_num],sep="")
  
  # for(model_index in 1:length(models) ){
  #   chain_list = list.files(pattern=models[model_index])
  #   print(chain_list)
  #   fname = paste(models[model_index],".rdata",sep="")
  #   print(fname)
  #   if (length(chain_list) != 4){
  #     print(paste(fname,"has ", length(chain_list),"chains",sep=" "))
  #     next
  #   }
  # }
  chain_list<-c()
  for(chain_index in 1:4 ){
    chain_list[chain_index] <- paste("modelfitting_ES_fit_chain",chain_index,"_exp_",selected_exp[exp_num],".rdata",sep="")
  }
  fname = paste("modelfitting_LE_fit_exp_",selected_exp[exp_num],".rdata",sep="")
  
  for (ii in 1:4) {
    load(chain_list[ii])
    if(ii == 1){
      f1 <- RA_fit
    } else if(ii==2){
      f2 <- RA_fit
    }else if(ii==3){
      f3 <- RA_fit
    } else if(ii==4){
      f4 <- RA_fit
    }
  }
  RA_combined = sflist2stanfit(list(f1, f2, f3, f4))
  save(RA_combined,file=fname)
  
  printRhat<-summary(RA_combined)$summary[,"Rhat"]
  print(printRhat[1:10])
  printNeff<-summary(RA_combined)$summary[,"n_eff"]
  print(printNeff[1:10])
  
  printSummary<-summary(RA_combined)
  print(printSummary[1])
  
  sum(get_sampler_params(RA_combined, inc_warmup = FALSE)[[1]][,"divergent__"])
  sum(get_sampler_params(RA_combined, inc_warmup = FALSE)[[2]][,"divergent__"])
  sum(get_sampler_params(RA_combined, inc_warmup = FALSE)[[3]][,"divergent__"])
  sum(get_sampler_params(RA_combined, inc_warmup = FALSE)[[4]][,"divergent__"])
  
  png(file=paste("Fig.ES_fit_exp_",selected_exp[exp_num],".trace.png",sep=""), width=600, height=350,unit="px",bg="transparent")
  traceplot(RA_combined,pars=c("mu_p","s_p"));
  dev.off()
  
  Sys.sleep(0.1)
  
  png(file=paste("Fig.ES_fit_exp_",selected_exp[exp_num],".pairs.png",sep=""), width=1000, height=1000,unit="px",bg="transparent")
  pairs(RA_combined,pars=c("mu_p","s_p","lp__"))
  dev.off()
  
  Sys.sleep(0.1)
  
  png(file=paste("Fig.ES_fit_exp_",selected_exp[exp_num],".dense.png",sep=""), width=600, height=350,unit="px",bg="transparent")
  stan_dens(RA_combined, pars=c('mu_p','s_p'), separate_chains = TRUE)
  dev.off()
  
  Sys.sleep(0.1)
  
  models <- fname
  load(models)
  mylog_lik <- extract_log_lik(RA_combined,parameter_name = "log_lik")
  ll1 = sum(mylog_lik)
  loo_power1 <- loo(mylog_lik)
  
  parameterExtract<-extract(RA_combined,pars=c("beta1","midPoints1"))
  rawPosteriorMidPoints1 <- parameterExtract$midPoints1
  medianPosteriorBeta1  <- apply(parameterExtract$beta1,2,median)
  mylog_lik <- extract_log_lik(RA_combined,parameter_name = "log_lik")
  myloo <- loo_power1$looic;
  fname = paste(getwd(),"/new_posterior_ES_models_5000permute_exp_",selected_exp[exp_num],".mat",sep="");
  writeMat(fname, rawPosteriorMidPoints1=rawPosteriorMidPoints1, medianPosteriorBeta1=medianPosteriorBeta1,log_lik=mylog_lik, loo_looic=myloo);
}
