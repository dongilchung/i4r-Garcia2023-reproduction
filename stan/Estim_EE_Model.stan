// EE model
data {
  int<lower=1> NumObservation;
  int<lower=1> NumSubject;
  int<lower=1> NumP1;
  real ListP1[NumP1];
  real Current_exp;
  real Group_p1Col[NumObservation];
  real Group_p2Col[NumObservation];
  int<lower=1,upper=2> Group_behavCol[NumObservation]; //choice
  int<lower=1> Group_subjCol[NumObservation];
  int<lower=1> Group_ntrialCol[NumObservation];
  int<lower=1> Group_nsubCol[NumObservation];
  int Group_sessionCol[NumObservation];
  real Group_expCol[NumObservation];
}

parameters {
  vector[NumP1+1] mu_p; 		      //mean of parameters
  vector<lower=0>[NumP1+1] s_p; 	//variance of parameters
  vector[NumSubject] beta1_raw;
  matrix[NumSubject,NumP1] midPoints1_raw;
}

transformed parameters { //transformations on parameters (as name suggests)
vector[NumSubject] beta1;
matrix[NumSubject,NumP1] midPoints1;

for (subjIndex in 1:NumSubject) {
  beta1[subjIndex] = 2.5*10^6*Phi_approx(mu_p[1] + s_p[1] * beta1_raw[subjIndex]);
  for (optIndex in 1:NumP1){
    midPoints1[subjIndex,optIndex] = Phi_approx(mu_p[1+optIndex] + s_p[1+optIndex] * midPoints1_raw[subjIndex,optIndex]);
  }
} 
}
model {
  int initIndex;
  int tmpNtrial[NumObservation];
  int nTrial;
  int index;
  // hyperparameters
  mu_p ~ normal(0,10);						//top level priors- v.wide
  s_p  ~ cauchy(0,2.5);
  
  beta1_raw ~ normal(0,1.0);
  for (optIndex in 1:NumP1){
    midPoints1_raw[:,optIndex] ~ normal(0,1.0);
  }

  tmpNtrial = Group_ntrialCol;
  nTrial = tmpNtrial[1];
  
  initIndex = 1;
  for (subjIndex in 1:NumSubject) {
    index = 1;
    for (trialIndex in initIndex:(initIndex+nTrial-1)) { 
      int BehavCol;
      real tmpUtil1;
      real tmpUtil2;
      for (optIndex in 1:NumP1){
        if(ListP1[optIndex]==Group_p1Col[trialIndex]){
          tmpUtil1 = midPoints1[subjIndex,optIndex];
        }
      }
      tmpUtil2 = Group_p2Col[trialIndex];
      
      BehavCol = Group_behavCol[trialIndex]==1;
      target += bernoulli_logit_lpmf(BehavCol | beta1[subjIndex] * (tmpUtil1 - tmpUtil2)); //functions as softmax
      index = index+1;
    }
    initIndex = initIndex+nTrial;
  }
}
generated quantities {
  real log_lik[NumSubject];
  int initIndex;
  int tmpNtrial[NumObservation];
  int nTrial;
  int index;
  
  { //local section, this saves time and space
  tmpNtrial = Group_ntrialCol;
  nTrial = tmpNtrial[1];
  
  initIndex = 1;
  for (subjIndex in 1:NumSubject) {
    log_lik[subjIndex] = 0;
    index = 1;
    for (trialIndex in initIndex:(initIndex+nTrial-1)) { 
      int BehavCol;
      real tmpUtil1;
      real tmpUtil2;
      for (optIndex in 1:NumP1){
        if(ListP1[optIndex]==Group_p1Col[trialIndex]){
          tmpUtil1 = midPoints1[subjIndex,optIndex];
        }
      }
      tmpUtil2 = Group_p2Col[trialIndex];
      
      BehavCol = Group_behavCol[trialIndex]==1;
      log_lik[subjIndex] = log_lik[subjIndex] + bernoulli_logit_lpmf(BehavCol | beta1[subjIndex] * (tmpUtil1 - tmpUtil2)); //functions as softmax
      index = index+1;
    }
    initIndex = initIndex+nTrial;
  }
  }
}