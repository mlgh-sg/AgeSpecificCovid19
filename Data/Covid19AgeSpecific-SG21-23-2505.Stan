functions {
  vector fmax_vector(vector a, real b) {
    vector[size(a)] result;
    for (i in 1:size(a)) {
      result[i] = fmax(a[i], b);
    }
    return result;
  }
}

data {
  int<lower=1> Seed;
  int<lower=1> TVaccine;
  int<lower=1> TOmicron;
  int<lower=1> N0; // number of initial days for which to estimate infections
  int<lower=1> N; // days of observed data for age group a. each entry must be <= N2
  int<lower=1> N2; // days of observed data + # of days to forecast
  int<lower=1> A; // number of age bands
  int<lower=1> W; // number of weeks
  array[N2+Seed] int<lower=1> WeekId;
  int<lower=1> SI_CUT; // number of days in serial interval to consider
  //           data
  vector<lower=1>[A] Pop; // Population in each age group
  vector<lower=0, upper=1>[A] popByAge; // proportion of age bracket in population
  array[N2, A] int casesByAge; //reported cases
  array[N2, A] int hospsByAge; //reported hospitalization
  //           priors
  matrix<lower=0>[A,A] cntct_mean;
  vector<lower=0>[A] ICR;
  vector<lower=0>[A] IHR;
  row_vector[SI_CUT] rev_serial_interval_Omicron; // fixed pre-calculated serial interval of Omicron using empirical data from Neil in reverse order
  row_vector[SI_CUT] rev_serial_interval_Delta; // fixed pre-calculated serial interval of Delta using empirical data from Neil in reverse order
  row_vector[N2] rev_ihr_daysSinceInfection; // Reversed ihr distribution
  row_vector[N2] rev_icr_daysSinceInfection; // Reversed icr distribution
}
 
transformed data{
  vector[A] ones_vector_A = rep_vector(1.,A);
  real<lower=0> avg_cntct = (popByAge') * (cntct_mean * ones_vector_A);
}
 
parameters {
  real<lower=0> kappa;
  real R0;
  real<lower=0> phi; // overdispersion parameter for likelihood model
  real<lower=0> phi2; // overdispersion parameter for likelihood model
  real<lower=0> phi3;
  vector<lower=0>[A] icr_noise;
  vector<lower=0>[A] ihr_noise;
  row_vector[2] log_relsusceptibility_age_reduced;
  matrix<lower=0>[W,A] weekly_effect;
  real<lower=0, upper=1> weekly_rho;
  real<lower=0, upper=1> weekly_rho1;
  real<lower=0> weekly_sd;
  vector<lower=0>[A] tau;
  matrix<lower=0>[N0,A] e_cases_N0;  // Initial values of E_infeByAge
}
 
transformed parameters {
  // expected new infections by calendar day, and age under self-renewal model
  // and a container to store the precomputed infections by age
  matrix<lower=0>[(N2+Seed),A] E_infeByAge;
  matrix<lower=0>[A, (N2+Seed)] suscep;
  // initial values of E_infeByAge
  E_infeByAge[1:N0,:] = e_cases_N0;
  vector[A] zero_vector_A = rep_vector(0.,A);

  // Time-varying contact matrix
  array[W] matrix<lower=0>[A, A] ConMatTran;
  vector<lower=1>[A] Pop1 = Pop*((1-0.96) + 0.96*phi3);
 
  for (t in 1:W){
    ConMatTran[t][1, 1] = cntct_mean[1, 1];
    ConMatTran[t][1, 2:A] = (cntct_mean[1, 2:A] .* weekly_effect[t, 2:A]);
    ConMatTran[t][2:A, 1] = cntct_mean[2:A, 1] .* (weekly_effect[t, 2:A])';
    ConMatTran[t][2:A, 2:A] = cntct_mean[2:A, 2:A] .* (weekly_effect[t, 2:A]'*weekly_effect[t, 2:A]);
  }
  
  row_vector<lower=0>[A] tmp_row_vector_A;
  vector<lower=0>[A] tmp_vector_A;
  matrix<lower=0>[A, A] tmp_matrix;
  // Reproduction number by day, age
  vector<lower=0>[N2+Seed] Rt; // Rt in general
  matrix<lower=0>[(N2+Seed),A] RtByAge; // Rt for each age band
  // probability of infection given contact
  real<lower=0> rho0;
  // define probability of infection given contact
  rho0 = exp(R0) ./ avg_cntct;
  // transform log_relsusceptibility_age
  row_vector[A] log_relsusceptibility_age;
  log_relsusceptibility_age[1] = log_relsusceptibility_age_reduced[1];
  log_relsusceptibility_age[2:5] = rep_row_vector(0., 4);
  log_relsusceptibility_age[6] = log_relsusceptibility_age_reduced[2]*0.5;
  log_relsusceptibility_age[7] = log_relsusceptibility_age_reduced[2];
  // t = 1
  suscep[:,1] = rep_vector(1, 7);
  tmp_matrix =  rep_matrix(exp(log_relsusceptibility_age), A).*ConMatTran[1];
  RtByAge[1,:] = (tmp_matrix * ones_vector_A)'*rho0;
  //t = 2
  suscep[:,2] = fmax_vector(1 - (E_infeByAge[1,:])'./ Pop, 0);
  tmp_matrix = rep_matrix(suscep[:,2], A) .* rep_matrix(exp(log_relsusceptibility_age), A).*ConMatTran[1];
  RtByAge[2,:] = (tmp_matrix * ones_vector_A)'*rho0;
  tmp_vector_A = (E_infeByAge[1,:])';
  for (t in 3:N0){
    tmp_vector_A = tmp_vector_A + (E_infeByAge[(t-1),:])';
    suscep[:,t] = fmax_vector(1 - tmp_vector_A ./ Pop, 0);
    tmp_matrix = rep_matrix(suscep[:,3], A) .* rep_matrix(exp(log_relsusceptibility_age), A).*ConMatTran[WeekId[t]];
    RtByAge[t,:] = (tmp_matrix * ones_vector_A)'*rho0;
  }
  // calculate expected infections by age under self-renewal model after first N0 days
  // and adjusted for saturation
  for (t in (N0+1):(TVaccine+Seed)){
    tmp_row_vector_A = rev_serial_interval_Delta[max(1,(SI_CUT-t+2)):SI_CUT]*E_infeByAge[max(1,(t-SI_CUT)):(t-1),:]*rho0;
    E_infeByAge[t,:] = tmp_row_vector_A * ConMatTran[WeekId[t]];
    tmp_vector_A = tmp_vector_A + (E_infeByAge[(t-1),:])';
    suscep[:,t] = fmax_vector(1 - tmp_vector_A ./ Pop, 0);
    E_infeByAge[t,:] .*= (suscep[:,t])';
    E_infeByAge[t,:] .*= exp(log_relsusceptibility_age);
    tmp_matrix = rep_matrix(suscep[:,t], A).*rep_matrix(exp(log_relsusceptibility_age), A) .*ConMatTran[WeekId[t]];
    RtByAge[t,:] = (tmp_matrix * ones_vector_A)'*rho0;
  }
  for (t in (TVaccine+1+Seed):(TOmicron+Seed)){
    tmp_row_vector_A = rev_serial_interval_Delta[max(1,(SI_CUT-t+2)):SI_CUT]*E_infeByAge[max(1,(t-SI_CUT)):(t-1),:]*rho0;
    E_infeByAge[t,:] = tmp_row_vector_A * ConMatTran[WeekId[t]];
    tmp_vector_A = tmp_vector_A + (E_infeByAge[(t-1),:])';
    suscep[:,t] = fmax_vector(1 - tmp_vector_A ./ Pop1, 0);
    E_infeByAge[t,:] .*= (suscep[:,t])';
    E_infeByAge[t,:] .*= exp(log_relsusceptibility_age);
    tmp_matrix = rep_matrix(suscep[:,t], A).*rep_matrix(exp(log_relsusceptibility_age), A) .*ConMatTran[WeekId[t]];
    RtByAge[t,:] = (tmp_matrix * ones_vector_A)'*rho0;
  }
  for (t in (TOmicron+1+Seed):(N2+Seed)){
    tmp_row_vector_A = rev_serial_interval_Omicron[max(1,(SI_CUT-t+2)):SI_CUT]*E_infeByAge[max(1,(t-SI_CUT)):(t-1),:]*rho0;
    E_infeByAge[t,:] = tmp_row_vector_A * ConMatTran[WeekId[t]];
    tmp_vector_A = tmp_vector_A + (E_infeByAge[(t-1),:])';
    suscep[:,t] = fmax_vector(1 - tmp_vector_A ./ Pop1, 0);
    E_infeByAge[t,:] .*= (suscep[:,t])';
    E_infeByAge[t,:] .*= exp(log_relsusceptibility_age);
    tmp_matrix = rep_matrix(suscep[:,t], A).*rep_matrix(exp(log_relsusceptibility_age), A) .*ConMatTran[WeekId[t]];
    RtByAge[t,:] = (tmp_matrix * ones_vector_A)'*rho0;
  }
  Rt = RtByAge * popByAge;
  // calculate expected cases by age
  // expected cases by calendar day (1st dim) age (2nd dim) , under self-renewal model
  matrix<lower=0>[N2,A] E_caseByAge;
  vector<lower=0>[N2] E_case;
  E_caseByAge[1,:] = (ICR)' .*  E_infeByAge[(1+Seed),:];
  for (t in 2:N2)
  {
    E_caseByAge[t,:] = rev_icr_daysSinceInfection[(N2-(t-1)+1):N2] * E_infeByAge[(1+Seed):(t-1+Seed),:] .* (ICR .* icr_noise)';
  }
    // calculate expected cases
  E_case = E_caseByAge * ones_vector_A;
   // calculate expected deaths by age
   // expected deaths by calendar day (1st dim) age (2nd dim) , under self-renewal model
  matrix<lower=0>[N2,A] E_hospsByAge = rep_matrix( 0., N2, A );
  vector<lower=0>[N2] E_hosps = rep_vector(0., N2);
  E_hospsByAge[1,:] = (IHR)' .* E_infeByAge[(1+Seed),:];
  for (t in 2:N2)
  {
    E_hospsByAge[t,:] = rev_ihr_daysSinceInfection[(N2-(t-1)+1):N2] * E_infeByAge[(1+Seed):(t-1+Seed),:] .* (IHR .* ihr_noise)';
  }
  E_hosps = E_hospsByAge * ones_vector_A;
}
model {
  // priors
  tau ~ exponential([0.45326160,0.03696352,0.05222578,0.17913391,0.41611017,0.40742085,2.11974178]');
  for (a in 1:A){
    e_cases_N0[:,a] ~ exponential(1/tau[a]);
  }
  weekly_sd ~ normal(0,0.2);
  weekly_rho ~ normal(0.8, 0.05);
  weekly_rho1 ~ normal(0.1, 0.05);
  weekly_effect[1,] ~ normal(1, 0.01);
  weekly_effect[2,] ~ normal(1,weekly_sd*sqrt(1-pow(weekly_rho,2)-pow(weekly_rho1,2)-2*pow(weekly_rho,2)*weekly_rho1/(1-weekly_rho1)));
  for (a in 1:A){
    weekly_effect[3:W,a] ~ normal(weekly_effect[2:(W-1),a]*weekly_rho + weekly_effect[1:(W-2),a]*weekly_rho1,
                              weekly_sd *sqrt(1-pow(weekly_rho,2)-pow(weekly_rho1,2)-2*pow(weekly_rho,2)*weekly_rho1/(1-weekly_rho1)));
  }
  phi ~ normal(0,5);
  phi2 ~ normal(0,2);
  phi3 ~ normal(0.8, 0.05);
  ihr_noise ~ normal(1,0.1);
  icr_noise ~ normal(1,0.1);
  log_relsusceptibility_age_reduced[1] ~ normal(-1.07, 0.22);//citation: Zhang et al Science
  log_relsusceptibility_age_reduced[2] ~ normal(0.38, 0.16);//citation: Zhang et al Science
  kappa ~ normal(0.05,0.01);
  R0 ~ normal(1, kappa);
  // likelihood
  for (a in 1:A){
    casesByAge[:,a] ~ neg_binomial(E_caseByAge[:,a]/phi, inv(phi));
    hospsByAge[:,a] ~ neg_binomial(E_hospsByAge[:,a]/phi2, inv(phi2));
  }
}
generated quantities {
  matrix<lower=0>[N2,A] neg_binomial_case_samples = rep_matrix( 0., N2, A );
  matrix<lower=0>[N2,A] neg_binomial_hosp_samples = rep_matrix( 0., N2, A );
  for (a in 1:A){
    for (t in 1:N2){
        neg_binomial_case_samples[t, a] = neg_binomial_rng(E_caseByAge[t,a]/phi, inv(phi));
        neg_binomial_hosp_samples[t, a] = neg_binomial_rng(E_hospsByAge[t,a]/phi2, inv(phi2));
    }
  }
}
