data {
  int<lower=0> M;  // Number of candidates
  int<lower=0> K;  // Number of raters
  int<lower=0> D;  // Number of wins or losses in sample
  int<lower=0> U;  // Number of undecided head-to-heads in sample
  vector[D] yd;     // Observed win-functions for decided matches
  int rd[D];     // Observed raters for decided matches (index)
  int ru[U];     // Observed raters for undecided matches (index)
  int ad[D];     // candidate for decided matches (index)
  int au[U];     // cadidate for undecided matches (index)
  int bd[D];     // opponent for decided matches (index)
  int bu[U];     // opponent for undecided matches (index)
  vector[D] futd; // follow-up times for the decided
  vector[U] futu; // follow-up times for the undecided
}

parameters {
  real lambda[M];
  //real<lower=0> lambda[M]; // severities
  //vector[K-1] eta_raw; // rater offsets
  real eta[K]; // rater offsets
}

transformed parameters {
  
  matrix[M,K] rater_lambda;
  
  for (k in 1:K){
  for (m in 1:M){
    rater_lambda[m,k] = exp(lambda[m] + eta[k]);
  }  
  }
  
}

model {
  // Prior distributions:
  for (m in 1:M){
    //lambda ~ gamma(1,1);
    lambda ~ normal(0,10);
  }
  
  //for (k in 1:(K-1)){
  //  eta_raw ~ normal(0,1);
  //}
  
  eta ~ normal(0, inv(sqrt(1 - inv(K)))); 

  // log-likelihood
  for (u in 1:U) {
    target += -futu[u]*(rater_lambda[au[u],ru[u]]+rater_lambda[bu[u],ru[u]]);
  }
  
  for (d in 1:D) {
    target += log(yd[d]*rater_lambda[ad[d],rd[d]] + (1-yd[d])*rater_lambda[bd[d],rd[d]]) - log(rater_lambda[ad[d],rd[d]] + rater_lambda[bd[d],rd[d]]) + log(1-exp(-futd[d]*(rater_lambda[ad[d],rd[d]]+rater_lambda[bd[d],rd[d]])));
  }
}
