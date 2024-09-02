// RW model
data {
  int<lower=1> NumObservation;
  int<lower=1> NumSubject;
  real Current_exp;
  int<lower=1,upper=2> Group_behavCol[NumObservation]; //choice
  int<lower=1> Group_subjCol[NumObservation];
  int<lower=1> Group_ntrialCol[NumObservation];
  int<lower=1> Group_conCol[NumObservation];
  int Group_cfCol[NumObservation];
  int<lower=0,upper=1> Group_outCol[NumObservation];
  int<lower=0,upper=1> Group_cfoutCol[NumObservation];
  int<lower=1> Group_nsubCol[NumObservation];
  int<lower=1> Group_nconCol[NumObservation];
  int Group_sessionCol[NumObservation];
  real Group_expCol[NumObservation];
}


parameters {
  vector[2] mu_p; 		      //mean of parameters
  vector<lower=0>[2] s_p; 	//variance of parameters
  vector[NumSubject] alpha1_raw;
  vector[NumSubject] beta1_raw;
}

transformed parameters { //transformations on parameters (as name suggests)
vector[NumSubject] alpha1;      
vector[NumSubject] beta1;      

for (subjIndex in 1:NumSubject) {
  alpha1[subjIndex]= Phi_approx(mu_p[1] + s_p[1] * alpha1_raw[subjIndex]);
  beta1[subjIndex] = 50*Phi_approx(mu_p[2] + s_p[2] * beta1_raw[subjIndex]);
} 
}

model {
  int initIndex;
  int tmpCF[NumObservation];
  int xCF;
  int tmpNcon[NumObservation];
  int nCon;
  int tmpNtrial[NumObservation];
  int nTrial;
  int index;
  // hyperparameters
  mu_p ~ normal(0,10);						//top level priors- v.wide
  s_p  ~ cauchy(0,2.5);
  
  alpha1_raw ~ normal(0,1.0);
  beta1_raw ~ normal(0,1.0);

  tmpCF= Group_cfCol;
  xCF = tmpCF[1];
  
  tmpNcon = Group_nconCol;
  nCon = tmpNcon[1];
  tmpNtrial = Group_ntrialCol;
  nTrial = tmpNtrial[1];
  
  initIndex = 1;
  for (subjIndex in 1:NumSubject) {
    real Qvalue1[nTrial+1,nCon];
    real Qvalue2[nTrial+1,nCon];
    for (conIndex in 1:nCon) {
      Qvalue1[1,conIndex] = 0.5; 
      Qvalue2[1,conIndex] = 0.5; 
    }
    index = 1;
    for (trialIndex in initIndex:(initIndex+nTrial-1)) { 
      int BehavCol;
      int tmpCho;
      int tmpCR;
      int tmpUR;
      int tmpCon;
      int ConArray[nCon];
      real tmpPEc;
      real tmpPEu;
      real tmpQ1;
      real tmpQ2;
      for (conIndex in 1:nCon) {
        ConArray[conIndex] = conIndex;
      }
      tmpCon = Group_conCol[trialIndex];
      tmpQ1  = Qvalue1[index,tmpCon];
      tmpQ2  = Qvalue2[index,tmpCon];
      tmpCho = Group_behavCol[trialIndex];
      tmpCR  = Group_outCol[trialIndex];
      tmpUR  = Group_cfoutCol[trialIndex];
      if (xCF==1){
        if (tmpCho==1){
          tmpPEc = tmpCR - tmpQ1;
          tmpPEu = tmpUR - tmpQ2;
          tmpQ1 = tmpQ1 + alpha1[subjIndex] * tmpPEc;
          tmpQ2 = tmpQ2 + alpha1[subjIndex] * tmpPEu;
        } else {
          tmpPEu = tmpUR - tmpQ1;
          tmpPEc = tmpCR - tmpQ2;
          tmpQ1 = tmpQ1 + alpha1[subjIndex] * tmpPEu;
          tmpQ2 = tmpQ2 + alpha1[subjIndex] * tmpPEc;
        }
      } else{
        if (tmpCho==1){
          tmpPEc = tmpCR - tmpQ1;
          tmpPEu = 0;
          tmpQ1 = tmpQ1 + alpha1[subjIndex] * tmpPEc;
          tmpQ2 = tmpQ2 + alpha1[subjIndex] * tmpPEu;
        } else {
          tmpPEu = 0;
          tmpPEc = tmpCR - tmpQ2;
          tmpQ1 = tmpQ1 + alpha1[subjIndex] * tmpPEu;
          tmpQ2 = tmpQ2 + alpha1[subjIndex] * tmpPEc;
        }
      }
      for (conIndex in 1:nCon) {
        if (conIndex==tmpCon){
          Qvalue1[index+1,conIndex] = tmpQ1;
          Qvalue2[index+1,conIndex] = tmpQ2;
        }else{
          Qvalue1[index+1,conIndex] = Qvalue1[index,conIndex];
          Qvalue2[index+1,conIndex] = Qvalue2[index,conIndex];
        }
      }
      
      BehavCol = Group_behavCol[trialIndex]==1;
      target += bernoulli_logit_lpmf(BehavCol | beta1[subjIndex] * (Qvalue1[index,tmpCon] - Qvalue2[index,tmpCon])); //functions as softmax
      index = index+1;
    }
    initIndex = initIndex+nTrial;
  }
}

generated quantities {
  real log_lik[NumSubject];
  int initIndex;
  int tmpCF[NumObservation];
  int xCF;
  int tmpNcon[NumObservation];
  int nCon;
  int tmpNtrial[NumObservation];
  int nTrial;
  int index;

  tmpCF= Group_cfCol;
  xCF = tmpCF[1];
  
  tmpNcon = Group_nconCol;
  nCon = tmpNcon[1];
  tmpNtrial = Group_ntrialCol;
  nTrial = tmpNtrial[1];
  
  initIndex = 1;
  for (subjIndex in 1:NumSubject) {
    real Qvalue1[nTrial+1,nCon];
    real Qvalue2[nTrial+1,nCon];
    log_lik[subjIndex] = 0;
    for (conIndex in 1:nCon) {
      Qvalue1[1,conIndex] = 0.5; 
      Qvalue2[1,conIndex] = 0.5; 
    }
    index = 1;
    for (trialIndex in initIndex:(initIndex+nTrial-1)) { 
      int BehavCol;
      int tmpCho;
      int tmpCR;
      int tmpUR;
      int tmpCon;
      int ConArray[nCon];
      real tmpPEc;
      real tmpPEu;
      real tmpQ1;
      real tmpQ2;
      for (conIndex in 1:nCon) {
        ConArray[conIndex] = conIndex;
      }
      tmpCon = Group_conCol[trialIndex];
      tmpQ1  = Qvalue1[index,tmpCon];
      tmpQ2  = Qvalue2[index,tmpCon];
      tmpCho = Group_behavCol[trialIndex];
      tmpCR  = Group_outCol[trialIndex];
      tmpUR  = Group_cfoutCol[trialIndex];
      if (xCF==1){
        if (tmpCho==1){
          tmpPEc = tmpCR - tmpQ1;
          tmpPEu = tmpUR - tmpQ2;
          tmpQ1 = tmpQ1 + alpha1[subjIndex] * tmpPEc;
          tmpQ2 = tmpQ2 + alpha1[subjIndex] * tmpPEu;
        } else {
          tmpPEu = tmpUR - tmpQ1;
          tmpPEc = tmpCR - tmpQ2;
          tmpQ1 = tmpQ1 + alpha1[subjIndex] * tmpPEu;
          tmpQ2 = tmpQ2 + alpha1[subjIndex] * tmpPEc;
        }
      } else{
        if (tmpCho==1){
          tmpPEc = tmpCR - tmpQ1;
          tmpPEu = 0;
          tmpQ1 = tmpQ1 + alpha1[subjIndex] * tmpPEc;
          tmpQ2 = tmpQ2 + alpha1[subjIndex] * tmpPEu;
        } else {
          tmpPEu = 0;
          tmpPEc = tmpCR - tmpQ2;
          tmpQ1 = tmpQ1 + alpha1[subjIndex] * tmpPEu;
          tmpQ2 = tmpQ2 + alpha1[subjIndex] * tmpPEc;
        }
      }
      for (conIndex in 1:nCon) {
        if (conIndex==tmpCon){
          Qvalue1[index+1,conIndex] = tmpQ1;
          Qvalue2[index+1,conIndex] = tmpQ2;
        }else{
          Qvalue1[index+1,conIndex] = Qvalue1[index,conIndex];
          Qvalue2[index+1,conIndex] = Qvalue2[index,conIndex];
        }
      }  
      BehavCol = Group_behavCol[trialIndex]==1;
      log_lik[subjIndex] = log_lik[subjIndex] + bernoulli_logit_lpmf(BehavCol | beta1 * (Qvalue1[index,tmpCon] - Qvalue2[index,tmpCon])); //functions as softmax
      index = index+1;
    }
    initIndex = initIndex+nTrial;
  }
}

