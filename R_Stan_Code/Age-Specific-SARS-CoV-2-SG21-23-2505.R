library(readxl)
library(EnvStats)
library(cmdstanr)
library(ggplot2)
library(patchwork)
library(multipanelfigure)
library(tidyr)
library(grid)
options(mc.cores=parallel::detectCores())

################### Data
set.seed(10)
N0 = 14 # // number of initial days for which to estimate infections
N = N2 = 658 # // days of observed data. Each entry must be <= N2
A = 7 # // number of age bands

W = 110
WeekId = 1
for (t in 1:W){
  WeekId = c(WeekId, rep(t,7))
}
WeekId = WeekId[-1]

################## Important Date
### from 2021 Sept to 2023 April  (-182 days)
Seed = 112
TVaccine = 365 # Time as total population have completed their vaccination regimen
# and the time as reducing suspetibe size of people
TOmicron = 427

# Serial Interval
serial_interval <- readRDS("~/Desktop/Covid19-Age/Covid19Data/serial-interval.rds")
SI = serial_interval$fit
SI_CUT = length(SI) #// number of days in serial interval to consider data

#serial_interval_Delta <- readRDS("~/Desktop/Covid19-Age/Covid19Data/serial_interval_Delta.rds")
#SI_Delta = serial_interval_Delta

#serial_interval_Omicron <- readRDS("~/Desktop/Covid19-Age/Covid19Data/serial_interval_Omicron.rds")
#SI_Omicron = serial_interval_Omicron

# Proportion of age bracket in population
# Age-specific population from https://www.statista.com/statistics/624913/singapore-population-by-age-group/
PopTotal = c(175.27,201.97,204.46,209.58,229.07,277.89,323.63,300.04,310.17,302.52,302.27,300.67,293.79, 259.04,197.89,260.92)*1000
Pop = c(581.7, 716.54, 623.67, 612.69, 603.04, 552.83, 458.81)*1000
popByAge = Pop / sum(Pop)
epidemicStart = 1

# Actual deaths
CovidDeathSG <- read_excel("~/Desktop/Covid19-Age/Covid19Data/CovidDeathSG212223.xlsx")
CovidDeathSG <- CovidDeathSG[188:845, ]

NewCol <- apply(CovidDeathSG[,2:3], 1, sum) + CovidDeathSG[,4]/2
NewCol1 <- CovidDeathSG[,4]/2 + CovidDeathSG[,5]
NewCol2 <- apply(CovidDeathSG[,10:11], 1, sum)
deathsByAge <- cbind(NewCol, NewCol1, CovidDeathSG[,6:9], NewCol2)
deathsByAge <- round(deathsByAge)
deaths <- apply(deathsByAge, 1, sum)
DeathsByAge <- apply(data.matrix(deathsByAge), 2, cumsum)

###############
CovidLocalSG <- read_excel("~/Desktop/Covid19-Age/Covid19Data/CovidLocalSG212223.xlsx")
CovidImportSG <- read_excel("~/Desktop/Covid19-Age/Covid19Data/CovidImportSG212223.xlsx")
CovidLocalSG <- CovidLocalSG[188:845, ]
CovidImportSG <- CovidImportSG[188:845, ]

NewColLocal <- apply(CovidLocalSG[,2:3], 1, sum) + CovidLocalSG[,4]/2
NewColLocal1 <- CovidLocalSG[,4]/2 + CovidLocalSG[,5]
NewColLocal2 <- apply(CovidLocalSG[,10:11], 1, sum)
LocalByAge <- cbind(NewColLocal,NewColLocal1, CovidLocalSG[,6:9], NewColLocal2)

NewColImport <- apply(CovidImportSG[,2:3], 1, sum) + CovidImportSG[,4]/2
NewColImport1 <- CovidImportSG[,4]/2 + CovidImportSG[,5]
NewColImport2 <- apply(CovidImportSG[,10:11], 1, sum)
importByAge <- cbind(NewColImport,NewColImport1,CovidImportSG[,6:9], NewColImport2)

casesByAge <- LocalByAge+importByAge
casesByAge <- round(casesByAge)
initialcase <- casesByAge[1:N0,]
cases <- apply(casesByAge, 1, sum)
CasesByAge = apply(data.matrix(casesByAge), 2, cumsum)

#####################################################
CovidActiveHospSG <- read_excel("~/Desktop/Covid19-Age/Covid19Data/CovidActiveHospSG212223.xlsx")
CovidNewHospSG <- read_excel("~/Desktop/Covid19-Age/Covid19Data/CovidNewHospSG212223.xlsx")
CovidActiveHospSG <- CovidActiveHospSG[188:845, ]
CovidNewHospSG <- CovidNewHospSG[188:845, ]
NewColActiveHospSG <- apply(CovidActiveHospSG[,2:3], 1, sum) + CovidActiveHospSG[,4]/2
NewColActiveHospSG1 <- CovidActiveHospSG[,4]/2 + CovidActiveHospSG[,5]
NewColActiveHospSG2 <- apply(CovidActiveHospSG[,10:11], 1, sum)
ActiveHospByAge <- cbind(NewColActiveHospSG,NewColActiveHospSG1, CovidActiveHospSG[,6:9], NewColActiveHospSG2)

NewColNewHospSG <- apply(CovidNewHospSG[,2:3], 1, sum) + CovidNewHospSG[,4]/2
NewColNewHospSG1 <- CovidNewHospSG[,4]/2 + CovidNewHospSG[,5]
NewColNewHospSG2 <- apply(CovidNewHospSG[,10:11], 1, sum)
NewHospByAge <- cbind(NewColNewHospSG,NewColNewHospSG1, CovidNewHospSG[,6:9], NewColNewHospSG2)

hospsByAge <- ActiveHospByAge+NewHospByAge
hospsByAge <- round(hospsByAge)
hosps <- apply(hospsByAge, 1, sum)
HospsByAge = apply(data.matrix(hospsByAge), 2, cumsum)

################
#Estimates of the severity of coronavirus disease 2019:a model-based analysis
CFR = rep(0, A)
CFR1 = c(2.6e-05, 1.48e-04*204.46/(204.46+209.58))
CFR[1] = sum(CFR1*c(377.24,204.46)/(377.24+204.46))
CFR2 = c(1.48e-04*209.58/(204.46+209.58), 6e-04)
CFR[2] = sum(CFR2*c(209.58,506.96)/(209.58+506.96))

CFR[3:6] = c(0.146e-02, 0.3e-02, 1.3e-02, 4e-02)
WeightLast = c(318.89, 139.92)
WeightLast = WeightLast/sum(WeightLast)
CFR[A] = sum(c(8.6e-02, 13.4e-02)*WeightLast)

# priors - Baseline Contact Matrix
ConMat = read_excel("~/Desktop/Covid19-Age/Covid19Data/ConMat.xlsx")
ProbTotal = rep(0, length(PopTotal))
Group = list(a = 1:3, b = 4:6, c = 7:8, d = 9:10, e = 11:12, f = 13:14, g = 15:16)
for (a in 1:A){
  ProbTotal[unlist(Group[a])] = PopTotal[unlist(Group[a])]/sum(PopTotal[unlist(Group[a])])
}
ProbMat = outer(ProbTotal, ProbTotal)
Cntct_Mean = ProbMat*ConMat   # // mean of prior contact rates between age groups
cntct_mean = array(rep(0, A*A), dim=c(A, A))
for (i in 1:A){
  for (j in 1:A){
    cntct_mean[i,j] = sum(Cntct_Mean[unlist(Group[i]), unlist(Group[j])])
  }
}

# priors - Lockdown Contact Matrix
ConMatLock = read_excel("~/Desktop/Covid19-Age/Covid19Data/ConMat_Lockdown.xlsx")
ProbTotal1 = rep(0, length(PopTotal)-1)
Group = list(a = 1:3, b = 4:6, c = 7:8, d = 9:10, e = 11:12, f = 13:14, g = 15)
for (a in 1:A){
  ProbTotal1[unlist(Group[a])] = PopTotal[unlist(Group[a])]/sum(PopTotal[unlist(Group[a])])
}
ProbMat1 = outer(ProbTotal1, ProbTotal1)
Cntct_Mean1 = ProbMat1*ConMatLock
cntct_mean_Lockdown = array(rep(0, A*A), dim=c(A, A))
for (i in 1:A){
  for (j in 1:A){
    cntct_mean_Lockdown[i,j] = sum(Cntct_Mean1[unlist(Group[i]), unlist(Group[j])])
  }
}

### Infection-to-death distribution
# Age-specific IFR and ICR
# Age-specific IFR is from ifr-by-age-prior_Levin_201117 in github of covid19model. 
IFR_mean = c(7.17380079770443e-06,1.31216329755262e-05,2.40008420078084e-05,4.39000546466559e-05,8.02977955900133e-05,
             0.000146873093210161,0.000268646281212568,0.00049138218536973,0.000898789502680301,0.00164398025721312,
             0.00300701205432415,0.00550014019012451,0.0100603329539299, 0.0184014030694962,0.0336580936908722,
             0.061564181804657,0.11260733795166, 0.302116652443295)
# Weighted IFR
IFR = rep(0, A)
for (a in 1:(A-1)){
  IFR[a] = sum(IFR_mean[unlist(Group[a])]*ProbTotal[unlist(Group[a])])
}
WeightLast = c(197.89, 121, 76.17, 63.75)
WeightLast = WeightLast/sum(WeightLast)
IFR[A] = sum(WeightLast*IFR_mean[15:18])
ICR = IFR/CFR

############# Calculation of IHR
# Monitoring the proportion of the population infected by
#SARS-CoV-2 using age-stratified hospitalisation and
#serological data: a modelling study
IHR = rep(0, A)
IHR = c(0.4e-02, 0.4e-02, 0.7e-02, 1e-02, 3e-02, 7e-02, 17.6e-02)

mean1 = 5.1; cv1 = 0.86; # infection to onset
mean2 = 5.83; cv2 = 0.88 # onset to case
mean3 = 3.88; cv3 = 0.44
mean4 = 4.31; cv4 = 0.685 # SI_Delta
mean5 = 2.64; cv5 = 0.399 # SI_Omicron
x1 = rgammaAlt(1e7,mean1,cv1) # infection-to-onset distribution
x2 = rgammaAlt(1e7,mean2,cv2) # onset-to-case distribution
x3 = rgammaAlt(1e7,mean3,cv3) # onset-to-hospitalization distribution
x4 = rgammaAlt(1e7,mean4,cv4) # SI_Delta
x5 = rgammaAlt(1e7,mean5,cv5) # SI_Omicron

ecdf.saved = ecdf(x1+x2)
ecdf.saved1 = ecdf(x1+x3)
ecdf.saved2 = ecdf(x4)
ecdf.saved3 = ecdf(x5)

### Estimate SI of Delta and Omicron
SI_Delta = rep(0, SI_CUT)
rev_serial_interval_Delta = rep(0, SI_CUT)
convolution = function(u) (ecdf.saved2(u))
SI_Delta[1] = (convolution(1.5) - convolution(0))
for(t in 2:100) {
  SI_Delta[t] = (convolution(t+0.5) - convolution(t-0.5))
}
rev_serial_interval_Delta = rev(SI_Delta)

SI_Omicron = rep(0, SI_CUT)
rev_serial_interval_Omicron = rep(0, SI_CUT)
convolution = function(u) (ecdf.saved3(u))
SI_Omicron[1] = (convolution(1.5) - convolution(0))
for(t in 2:100) {
  SI_Omicron[t] = (convolution(t+0.5) - convolution(t-0.5))
}
rev_serial_interval_Omicron = rev(SI_Omicron)

### Estimate distributions for infection-to-case by composing two distributions:
# (1) time from infection to symptom onset, and
# (2) time from symptom onset to case.
# from literature "Statistical Deconvolution for Inference of Infection Time Series"
icr_daysSinceInfection = rep(0, N2)
rev_icr_daysSinceInfection = rep(0, N2)

convolution = function(u) (ecdf.saved(u))
icr_daysSinceInfection[1] = (convolution(1.5) - convolution(0))
for(t in 2:N2) {
  icr_daysSinceInfection[t] = (convolution(t+0.5) - convolution(t-0.5))
}
rev_icr_daysSinceInfection = rev(icr_daysSinceInfection)

### Estimate distributions for infection-to-hospitalization by composing two distributions:
# (1) time from infection to symptom onset, and
# (2) time from symptom onset to hospitalization
ihr_daysSinceInfection = rep(0, N2)
rev_ihr_daysSinceInfection = rep(0, N2)

convolution = function(u) (ecdf.saved1(u))
ihr_daysSinceInfection[1] = (convolution(1.5) - convolution(0))
for(t in 2:N2) {
  ihr_daysSinceInfection[t] = (convolution(t+0.5) - convolution(t-0.5))
}
rev_ihr_daysSinceInfection = rev(ihr_daysSinceInfection)

### Cmdstan runs 4 Markov chains each with 1500 iterations, 500 iterations are discarded as "warm-up".
### Specify control parameters
control_params <- list(adapt_delta = 0.9, max_treedepth = 14)
start.time <- Sys.time()
file <- "/Users/jingyan/Desktop/Covid19-Age/Covid19-Age_specific/Singapore/testing_0823.stan"
mod <- cmdstan_model(file)
data_list <-
  list(
    Seed=Seed,
    TVaccine=TVaccine,
    TOmicron=TOmicron,
    N0=N0,
    N=N,
    N2=N2,
    A=A,
    W=W,
    WeekId=WeekId,
    SI_CUT=SI_CUT,
    Pop=Pop,
    popByAge=popByAge,
    cntct_mean=cntct_mean,
    ICR=ICR,
    IHR=IHR,
    rev_serial_interval_Delta=rev_serial_interval_Delta,
    rev_serial_interval_Omicron=rev_serial_interval_Omicron,
    casesByAge = as.matrix(casesByAge),
    hospsByAge = as.matrix(hospsByAge),
    rev_ihr_daysSinceInfection=rev_ihr_daysSinceInfection,
    rev_icr_daysSinceInfection=rev_icr_daysSinceInfection)

m_hier <- mod$sample(
  seed = 100,
  data=data_list,
  iter_warmup=500,
  iter_sampling=500,
  chains=4,
  parallel_chains = 4,
  max_treedepth = 14,
  adapt_delta = 0.9,
  refresh = 50)
end.time <- Sys.time()
time.taken <- end.time - start.time
print(time.taken)
m_hier$save_object(file = "cases_hospitalizations_ihr0.3_phi22_TimeVaryInfe_58_fit_1203.RDS")

m_hier = cases_hospitalizations_ihr0.8_phi21_3_TimeVaryInfe_60_fit_1202
E_infeByAge_samples <- m_hier$draws("E_infeByAge")
E_hospsByAge_samples <- m_hier$draws("E_hospsByAge")
E_caseByAge_samples <- m_hier$draws("E_caseByAge")
RtByAge_samples <- m_hier$draws("RtByAge")
E_case_samples <- m_hier$draws("E_case")
E_hosp_samples <- m_hier$draws("E_hosps")
Rt_samples <- m_hier$draws("Rt")
Suscep_samples <- m_hier$draws("suscep")
rho0_samples <- m_hier$draws("rho0")
log_relsusceptibility_sampels <- m_hier$draws("log_relsusceptibility_age")
ConMatTran_samples <- m_hier$draws("ConMatTran")
e_cases_N0_samples <- m_hier$draws("e_cases_N0")
weekly_effect_samples <- m_hier$draws("weekly_effect")
TimeVaryInfeByAge_samples <- m_hier$draws("TimeVaryInfeByAge")

############################################# Calculation of proportions of infections
### Parameters
rho0 = mean(rho0_samples)
SuscepByAge = matrix(rep(0, A*(N2+Seed)), nrow = N2+Seed)
log_relsusceptibility = rep(0, A)
for (a in 1:A){
  log_relsusceptibility[a] = mean(log_relsusceptibility_sampels[,,a])
}
#rho0*exp(log_relsusceptibility)
ConMatTran <- array(data = rep(0,W*A*A), dim = c(W, A, A))
for (w in 1:W){
  for (a in 1:A){
    for (b in 1:A){
      ConMatTran[w, a, b] = mean(ConMatTran_samples[,,((a-1)*A+(b-1))*A + w])
    }
  }
}

TimeVaryInfeByAge <- array(data = rep(0,(N2+Seed)*A*A), dim = c((N2+Seed), A, A))
for (t in 1:(N2+Seed)){
  for (a in 1:A){
    for (b in 1:A){
      TimeVaryInfeByAge[t, a, b] = mean(TimeVaryInfeByAge_samples[,,((a-1)*A+(b-1))*A + t])
    }
  }
}

e_cases_N0 = matrix(rep(0, A*N0), nrow = N0)
for (t in 1:N0){
  for (a in 1:A){
    e_cases_N0[t,a] = mean(e_cases_N0_samples[,,(a-1)*N0 + t])
  }
}

for (t in 1:(N2+Seed)){
  for (a in 1:A){
    SuscepByAge[t,a] = mean(Suscep_samples[,,(a-1)*N2 + t])
  }
}

### Proportions of infections from different age groups
Proportion_infection = array(data = rep(0,(N2+Seed)*A*A), dim = c((N2+Seed), A, A))
for (a in 1:A){
  Proportion_infection[1:N0,a,] = matrix(rep(e_cases_N0[,a]/A, A), nrow = N0)
}

for (t in (N0+1):(TOmicron+Seed)){
  for (a in 1:A){
    tmp_row_vector_A = t(Proportion_infection[max(1,(t-SI_CUT)):(t-1),a,]) %*% as.matrix(rev_serial_interval_Delta[max(1,(SI_CUT-t+2)):SI_CUT], nrow = A)
    Proportion_infection[t,a,] = ConMatTran[WeekId[t],,] %*% tmp_row_vector_A *SuscepByAge[t,]*rho0*exp(log_relsusceptibility)
  }
}

for (t in (TOmicron+1+Seed):(N2+Seed)){
  for (a in 1:A){
    tmp_row_vector_A = t(Proportion_infection[max(1,(t-SI_CUT)):(t-1),a,]) %*% as.matrix(rev_serial_interval_Omicron[max(1,(SI_CUT-t+2)):SI_CUT], nrow = A)
    Proportion_infection[t,a,] = ConMatTran[WeekId[t],,] %*% tmp_row_vector_A *SuscepByAge[t,]*rho0*exp(log_relsusceptibility)
  }
}

############################################# Bayes Plot from

Neg_Binomial_case = m_hier$draws("neg_binomial_case_samples")
Neg_Binomial_hosp = m_hier$draws("neg_binomial_hosp_samples")

RtByAge = matrix(rep(0, A*N2), nrow = N2)
Rt = rep(0, N2)
CaseByAge = matrix(rep(0, A*N2), nrow = N2)
Case = rep(0, N2)
HospByAge = matrix(rep(0, A*N2), nrow = N2)
Hosp = rep(0, N2)

MeanDeathByAge = matrix(rep(0, A*N2), nrow = N2)
InfeByAge = matrix(rep(0, A*(N2+Seed)), nrow = N2+Seed)
MeanHospByAge = matrix(rep(0, A*N2), nrow = N2)

RtUp = matrix(rep(0, A*N2), nrow = N2)
RtDown = matrix(rep(0, A*N2), nrow = N2)
RtUp1 = rep(0, N2)
RtDown1 = rep(0, N2)

CaseUp = matrix(rep(0, A*N2), nrow = N2)
CaseDown = matrix(rep(0, A*N2), nrow = N2)
#DeathUp = matrix(rep(0, A*N2), nrow = N2)
#DeathDown = matrix(rep(0, A*N2), nrow = N2)
HospUp = matrix(rep(0, A*N2), nrow = N2)
HospDown = matrix(rep(0, A*N2), nrow = N2)

#RtUp_NB = matrix(rep(0, A*N2), nrow = N2)
#RtDown_NB = matrix(rep(0, A*N2), nrow = N2)
CaseUp_NB = matrix(rep(0, A*N2), nrow = N2)
CaseDown_NB = matrix(rep(0, A*N2), nrow = N2)
CaseUp_NB1 = rep(0, N2)
CaseDown_NB1 = rep(0, N2)
HospUp_NB = matrix(rep(0, A*N2), nrow = N2)
HospDown_NB = matrix(rep(0, A*N2), nrow = N2)

weekly_effect = matrix(rep(0, A*W), nrow = W)
#DeathUp_NB = matrix(rep(0, A*N2), nrow = N2)
#DeathDown_NB = matrix(rep(0, A*N2), nrow = N2)

##############################################
for (t in 1:W){
  for (a in 1:A){
    weekly_effect[t, a] = mean(weekly_effect_samples[,,(a-1)*W + t])
  }
}

for (t in 1:(N2+Seed)){
  for (a in 1:A){
    InfeByAge[t,a] = mean(E_infeByAge_samples[,,(a-1)*(N+Seed) + t])
  }
}

for (t in 1:N2){
  for (a in 1:A){
    ### InfectionByAge
    #CaseByAge[t,a] = mean(Normal_case[,,(a-1)*N + t])
    #MeanDeathByAge[t,a] = mean(E_deathsByAge_samples[,c(1,4),(a-1)*N + t])
    #RtByAge[t,a] = mean(Normal_rt[,,(a-1)*N + t])
    #MeanHospByAge[t,a] = mean(E_hospsByAge_samples[,,(a-1)*N + t])
    
    CaseByAge[t,a] = mean(E_caseByAge_samples[,,(a-1)*N + t])
    RtByAge[t,a] = mean(RtByAge_samples[,,(a-1)*(N+Seed) + (t+Seed)])
    HospByAge[t,a] = mean(E_hospsByAge_samples[,,(a-1)*N + t])
    
    #Credible Intervals
    CaseDown[t,a] = quantile(E_caseByAge_samples[,,(a-1)*N + t], probs = c(0.025, 0.975))[1]
    CaseUp[t,a] = quantile(E_caseByAge_samples[,,(a-1)*N + t], probs = c(0.025, 0.975))[2]
    HospDown[t,a] = quantile(E_hospsByAge_samples[,,(a-1)*N + t], probs = c(0.025, 0.975))[1]
    HospUp[t,a] = quantile(E_hospsByAge_samples[,,(a-1)*N + t], probs = c(0.025, 0.975))[2]
    
    #DeathDown[t,a] = quantile(E_deathsByAge_samples[,c(1,4),(a-1)*N + t], probs = c(0.025, 0.975))[1]
    #DeathUp[t,a] = quantile(E_deathsByAge_samples[,c(1,4),(a-1)*N + t], probs = c(0.025, 0.975))[2]
    #RtDown[t,a] = quantile(RtByAge_samples[,,(a-1)*N + t], probs = c(0.025, 0.975))[1]
    #RtUp[t,a] = quantile(RtByAge_samples[,,(a-1)*N + t], probs = c(0.025, 0.975))[2]
    RtDown[t,a] = quantile(RtByAge_samples[,,(a-1)*(N+Seed) + (t+Seed)], probs = c(0.025, 0.975))[1]
    RtUp[t,a] = quantile(RtByAge_samples[,,(a-1)*(N+Seed) + (t+Seed)], probs = c(0.025, 0.975))[2]
    
    #CaseDown_NB[t,a] = quantile(Neg_Binomial[1001:2000,t,a], probs = c(0.025, 0.975))[1]
    #CaseUp_NB[t,a] = quantile(Neg_Binomial[1001:2000,t,a], probs = c(0.025, 0.975))[2]
    CaseDown_NB[t,a] = quantile(Neg_Binomial_case[,,(a-1)*N + t], probs = c(0.025, 0.975))[1]
    CaseUp_NB[t,a] = quantile(Neg_Binomial_case[,,(a-1)*N + t], probs = c(0.025, 0.975))[2]
    HospDown_NB[t,a] = quantile(Neg_Binomial_hosp[,,(a-1)*N+t], probs = c(0.025, 0.975))[1]
    HospUp_NB[t,a] = quantile(Neg_Binomial_hosp[,,(a-1)*N+t], probs = c(0.025, 0.975))[2]
  }
}

for (t in 1:N2){
  Case[t] = mean(E_case_samples[,,t])
  Rt[t] = mean(Rt_samples[,,(t+Seed)])
  Hosp[t] = mean(E_hosp_samples[,,t])
  RtDown1[t] = quantile(Rt_samples[,,(t+Seed)], probs = c(0.025, 0.975))[1]
  RtUp1[t] = quantile(Rt_samples[,,(t+Seed)], probs = c(0.025, 0.975))[2]
}

CaseDown_NB1 = apply(CaseDown_NB, 1, sum)
CaseUp_NB1 = apply(CaseUp_NB, 1, sum)
HospDown_NB1 = apply(HospDown_NB, 1, sum)
HospUp_NB1 = apply(HospUp_NB, 1, sum)

###############################
### Plot of Case fit
time_series_data1_case <- expand.grid(
  time = seq(as.Date("2021-07-07"), by = "days", length.out = N)
)

Age = 1
# Plot of Age-specific case fit
time_series_data1_case$values <- as.vector(CaseByAge[c(2,2:N),Age])
time_series_data1_case$lower <- as.vector(CaseDown[1:N,Age])
time_series_data1_case$upper <- as.vector(CaseUp[1:N,Age])
time_series_data1_case$lower_NB <- as.vector(CaseDown_NB[1:N,Age])
time_series_data1_case$upper_NB <- as.vector(CaseUp_NB[1:N,Age])
time_series_data1_case$observe <- as.vector(data.matrix(casesByAge[1:N,Age]))

# Plot of case fit across all age groups
time_series_data1_case$values <- as.vector(Case[c(2,2:N)])
time_series_data1_case$lower_NB <- as.vector(CaseDown_NB1[1:N])
time_series_data1_case$upper_NB <- as.vector(CaseUp_NB1[1:N])
time_series_data1_case$observe <- as.vector(cases[1:N])

ggplot(time_series_data1_case, aes(x = time, y = values)) +
  geom_line(color = "#0C7BDC", linewidth = 0.8, alpha = 0.6) +
  geom_col(aes(x = time, y = observe), fill = "#FFC20A", alpha=0.6, size = 0.5) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill=  "#0C7BDC", alpha = 0.4) +
  labs(
    title = NULL,
    x = NULL,
    y = "Daily number of cases"
  ) +
  scale_x_continuous(limits = c(as.Date("2021-07-07"), as.Date("2023-04-27")), breaks = seq(as.Date("2021-07-07"), as.Date("2023-04-27"), 90)) +
  scale_y_continuous(limits = c(0, 8000), breaks = seq(0, 8000, 500)) +
  theme(
    axis.text = element_text(size = 10, face = "bold"),
    axis.title = element_text(size = 13, face = "bold"),
    plot.title = element_text(size = 16, face = "bold"),
    plot.background = element_rect(fill = "white", color = "white")  # Set entire plot background to white
  )

# Create a custom labeling function for months and years
custom_labels <- function(x) {
  # Extract the first letter of the month (e.g., "J", "F", "M", ...) and year
  months <- toupper(substr(format(x, "%b"), 1, 1))  # Get first letter of the month
  years <- format(x, "%Y")  # Get the year (e.g., "2021", "2022", etc.)
  
  # Initialize a vector to store the labels
  labels <- sapply(x, function(date) {
    month_letter <- toupper(substr(format(date, "%b"), 1, 1))  # First letter of the month
    year <- format(date, "%Y")  # Year
    
    # For January and July, return the month letter with year below it
    if (format(date, "%m") == "01" || format(date, "%m") == "07") {
      return(paste(month_letter, "\n", year))  # Display "J" or "J" with the year below for January and July
    } else {
      return(month_letter)  # For other months, show only the first letter of the month
    }
  })
  
  return(labels)
}

ggplot(time_series_data1_case, aes(x = time, y = values)) +
  geom_line(color = "#0C7BDC", linewidth = 0.8, alpha = 0.6) +
  geom_col(aes(x = time, y = observe), fill = "#FFC20A" , alpha=0.6, size = 0.5) +
  geom_ribbon(aes(ymin = lower_NB, ymax = upper_NB), fill = "#0C7BDC",alpha = 0.4) +
  labs(
    title = NULL,
    x = NULL,
    y = "Daily number of cases"
  ) +
  scale_x_date(
    limits = c(as.Date("2021-07-07"), as.Date("2023-04-27")),  # Set the x-axis limits
    breaks = seq(as.Date("2021-07-07"), as.Date("2023-04-27"), by = "1 month"),  # Breaks every month
    labels = custom_labels
  ) +
  scale_y_continuous(limits = c(0, 22000), breaks = seq(0, 22000, 2000)) +
  theme(
    axis.text = element_text(size = 28, face = "bold"),  # Axis text styling
    axis.title = element_text(size = 28, face = "bold"),  # Axis title styling
    plot.title = element_text(size = 30, face = "bold"),  # Plot title styling
    plot.background = element_rect(fill = "white", color = "white"),  # White plot background
    panel.grid.major = element_line(color = "grey", size = 0.5),  # Major grid lines in white
    panel.grid.minor = element_line(color = "grey", size = 0.25),  # Minor grid lines in white
    panel.background = element_rect(fill = "white", color = "white"),  # White background for the plot area
    axis.text.x = element_text(angle = 45, hjust = 1),  # Adjust labels' angle, horizontal and vertical alignment
    axis.title.x = element_text(size = 16, face = "bold"),
    #axis.text.x = element_blank(),  # Removes x-axis labels
    #axis.title.x = element_blank()  # Removes x-axis title
  ) 

###############################
### Plot of Number of infectious individuals
time_series_data1_infe <- expand.grid(
  time = seq(as.Date("2021-07-07"), by = "days", length.out = N)
)

Age = 1
# Plot of Age-specific number of infectious individuals
time_series_data1_infe$values <- as.vector(InfeByAge[(Seed+1):(Seed+N),Age])

# Plot of number of infectious individuals across all age groups
#time_series_data1_infe$values <- apply(InfeByAge[(Seed+1):(Seed+N),], 1, sum)

ggplot(time_series_data1_infe, aes(x = time, y = values)) +
  geom_line(color = "#088F8F", linewidth = 0.8, alpha = 0.6) +
  #  geom_col(aes(x = time, y = observe), fill = "lightcoral", alpha=0.6, size = 0.5) +
  labs(
    title = NULL,
    x = NULL,
    y = "Daily number of infectious individuals"
  ) +
  scale_x_continuous(limits = c(as.Date("2021-07-07"), as.Date("2023-04-27")), breaks = seq(as.Date("2021-07-07"), as.Date("2023-04-27"), 90)) +
  scale_y_continuous(limits = c(0, 46000), breaks = seq(0, 46000, 2000)) +
  theme(
    axis.text = element_text(size = 10, face = "bold"),
    axis.title = element_text(size = 13, face = "bold"),
    plot.title = element_text(size = 16, face = "bold"),
    plot.background = element_rect(fill = "white", color = "white")  # Set entire plot background to white
  )

ggplot(time_series_data1_infe, aes(x = time, y = values)) +
  geom_line(color = "#088F8F", linewidth = 0.8, alpha = 0.6) +
  #  geom_col(aes(x = time, y = observe), fill = "lightcoral" , alpha=0.4, size = 0.5) +
  labs(
    title = NULL,
    x = NULL,
    y = "Daily number of infectious individuals"
  ) +
  scale_x_continuous(limits = c(as.Date("2021-07-07"), as.Date("2023-04-27")), breaks = seq(as.Date("2021-07-07"), as.Date("2023-04-27"), 180)) +
  scale_y_continuous(limits = c(0, 5000), breaks = seq(0, 5000, 500)) +
  theme(
    axis.text = element_text(size = 28, face = "bold"),  # Axis text styling
    axis.title = element_text(size = 28, face = "bold"),  # Axis title styling
    plot.title = element_text(size = 30, face = "bold"),  # Plot title styling
    plot.background = element_rect(fill = "white", color = "white"),  # White plot background
    panel.grid.major = element_line(color = "grey", size = 0.5),  # Major grid lines in white
    panel.grid.minor = element_line(color = "grey", size = 0.25),  # Minor grid lines in white
    panel.background = element_rect(fill = "white", color = "white"),  # White background for the plot area
    #axis.text.x = element_text(angle = 0, hjust = -0.8, vjust = 1),  # Adjust labels' angle, horizontal and vertical alignment
    #axis.title.x = element_text(size = 16, face = "bold"),
    axis.text.x = element_blank(),  # Removes x-axis labels
    axis.title.x = element_blank()  # Removes x-axis title
  )

###############################
### Plot of Hosp fit
time_series_data1_hosp <- expand.grid(
  time = seq(as.Date("2021-07-07"), by = "days", length.out = N)
)

Age = 1
# Plot of Age-specific hosp fit
time_series_data1_hosp$values <- as.vector(HospByAge[c(2,2:N),Age])
time_series_data1_hosp$lower <- as.vector(HospDown[1:N,Age])
time_series_data1_hosp$upper <- as.vector(HospUp[1:N,Age])
time_series_data1_hosp$lower_NB <- as.vector(HospDown_NB[1:N,Age])
time_series_data1_hosp$upper_NB <- as.vector(HospUp_NB[1:N,Age])
time_series_data1_hosp$observe <- as.vector(data.matrix(hospsByAge[1:N,Age]))

# Plot of hosp fit across all age groups
time_series_data1_hosp$values <- as.vector(Hosp[c(2,2:N)])
time_series_data1_hosp$lower_NB <- as.vector(HospDown_NB1[1:N])
time_series_data1_hosp$upper_NB <- as.vector(HospUp_NB1[1:N])
time_series_data1_hosp$observe <- as.vector(hosps[1:N])

ggplot(time_series_data1_hosp, aes(x = time, y = values)) +
  geom_line(color = "#0C7BDC", linewidth = 0.8, alpha = 0.6) +
  geom_col(aes(x = time, y = observe), fill = "#FFC20A", alpha=0.6, size = 0.5) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill=  "#0C7BDC", alpha = 0.4) +
  labs(
    title = NULL,
    x = NULL,
    y = "Daily number of hospitalizations"
  ) +
  scale_x_continuous(limits = c(as.Date("2021-07-07"), as.Date("2023-04-27")), breaks = seq(as.Date("2021-07-07"), as.Date("2023-04-27"), 90)) +
  scale_y_continuous(limits = c(0, 300), breaks = seq(0, 300, 50)) +
  theme(
    axis.text = element_text(size = 28, face = "bold"),
    axis.title = element_text(size = 28, face = "bold"),
    plot.title = element_text(size = 30, face = "bold"),
    plot.background = element_rect(fill = "white", color = "white")  # Set entire plot background to white
  )

# Create a custom labeling function for months and years
custom_labels <- function(x) {
  # Extract the first letter of the month (e.g., "J", "F", "M", ...) and year
  months <- toupper(substr(format(x, "%b"), 1, 1))  # Get first letter of the month
  years <- format(x, "%Y")  # Get the year (e.g., "2021", "2022", etc.)
  
  # Initialize a vector to store the labels
  labels <- sapply(x, function(date) {
    month_letter <- toupper(substr(format(date, "%b"), 1, 1))  # First letter of the month
    year <- format(date, "%Y")  # Year
    
    # For January and July, return the month letter with year below it
    if (format(date, "%m") == "01" || format(date, "%m") == "07") {
      return(paste(month_letter, "\n", year))  # Display "J" or "J" with the year below for January and July
    } else {
      return(month_letter)  # For other months, show only the first letter of the month
    }
  })
  
  return(labels)
}

ggplot(time_series_data1_hosp, aes(x = time, y = values)) +
  geom_line(color = "#0C7BDC", linewidth = 0.8, alpha = 0.6) +
  geom_col(aes(x = time, y = observe), fill = "#FFC20A" , alpha=0.6, size = 0.5) +
  geom_ribbon(aes(ymin = lower_NB, ymax = upper_NB), fill = "#0C7BDC",alpha = 0.4) +
  labs(
    title = NULL,
    x = NULL,
    y = "Daily number of hospitalizations"
  ) +
  scale_x_date(
    limits = c(as.Date("2021-07-07"), as.Date("2023-04-27")),  # Set the x-axis limits
    breaks = seq(as.Date("2021-07-07"), as.Date("2023-04-27"), by = "1 month"),  # Breaks every month
    labels = custom_labels  # Custom label function: First letter in uppercase
  ) +
  scale_y_continuous(limits = c(0, 4000), breaks = seq(0, 4000, 400)) +
  theme(
    axis.text = element_text(size = 28, face = "bold"),  # Axis text styling
    axis.title = element_text(size = 28, face = "bold"),  # Axis title styling
    plot.title = element_text(size = 30, face = "bold"),  # Plot title styling
    plot.background = element_rect(fill = "white", color = "white"),  # White plot background
    panel.grid.major = element_line(color = "grey", size = 0.5),  # Major grid lines in white
    panel.grid.minor = element_line(color = "grey", size = 0.25),  # Minor grid lines in white
    panel.background = element_rect(fill = "white", color = "white"),  # White background for the plot area
    axis.text.x = element_text(angle = 45, hjust = 1),  # Adjust labels' angle, horizontal and vertical alignment
    axis.title.x = element_text(size = 16, face = "bold"),
    #axis.text.x = element_blank(),  # Removes x-axis labels
    #axis.title.x = element_blank()  # Removes x-axis title
  ) 

###############################
### Plot of Rt estimation
time_series_data1_rt <- expand.grid(
  time = seq(as.Date("2021-07-07"), by = "days", length.out = N)
)

Age = 1
# Plot of Age-specific Rt estimation
time_series_data1_rt$values <- as.vector(RtByAge[1:N, Age])
time_series_data1_rt$lower <- as.vector(RtDown[1:N,Age])
time_series_data1_rt$upper <- as.vector(RtUp[1:N,Age])
#time_series_data1_rt$lower_NB <- as.vector(RtDown_NB[1:N,Age])
#time_series_data1_rt$upper_NB <- as.vector(RtUp_NB[1:N,Age])

# Plot of Rt estimation across all age groups
time_series_data1_rt$values <- as.vector(Rt[1:N])
time_series_data1_rt$lower <- as.vector(RtDown1[1:N])
time_series_data1_rt$upper <- as.vector(RtUp1[1:N])

ggplot(time_series_data1_rt, aes(x = time, y = values)) +
  geom_line(color = "#0C7BDC", linewidth = 0.8, alpha = 0.6) +
  geom_hline(yintercept = 1, color = "black") +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "#0C7BDC", alpha = 0.4) +
  labs(
    title = NULL,
    x = NULL,
    y = "Rt"
  ) +
  scale_x_date(
    limits = c(as.Date("2021-07-07"), as.Date("2023-04-27")),  # Set the x-axis limits
    breaks = seq(as.Date("2021-07-07"), as.Date("2023-04-27"), by = "1 month"),  # Breaks every month
    labels = function(x) toupper(substr(format(x, "%B"), 1, 1))  # Custom label function: First letter in uppercase
  ) +
  scale_y_continuous(limits = c(0, 2.5), breaks = seq(0, 2.5, 0.5)) +
  theme(
    axis.text = element_text(size = 14, face = "bold"),  # Axis text styling
    axis.title = element_text(size = 15, face = "bold"),  # Axis title styling
    plot.title = element_text(size = 18, face = "bold"),  # Plot title styling
    plot.background = element_rect(fill = "white", color = "white"),  # White plot background
    panel.grid.major = element_line(color = "grey", size = 0.5),  # Major grid lines in white
    panel.grid.minor = element_line(color = "grey", size = 0.25),  # Minor grid lines in white
    panel.background = element_rect(fill = "white", color = "white"),  # White background for the plot area
    axis.text.x = element_text(angle = 0, hjust = -0.8, vjust = 1),  # Adjust labels' angle, horizontal and vertical alignment
    axis.title.x = element_text(size = 16, face = "bold"),
    #axis.text.x = element_blank(),  # Removes x-axis labels
    #axis.title.x = element_blank()  # Removes x-axis title
  )

# Create a custom labeling function for months and years
custom_labels <- function(x) {
  # Extract the first letter of the month (e.g., "J", "F", "M", ...) and year
  months <- toupper(substr(format(x, "%b"), 1, 1))  # Get first letter of the month
  years <- format(x, "%Y")  # Get the year (e.g., "2021", "2022", etc.)
  
  # Initialize a vector to store the labels
  labels <- sapply(x, function(date) {
    month_letter <- toupper(substr(format(date, "%b"), 1, 1))  # First letter of the month
    year <- format(date, "%Y")  # Year
    
    # For January and July, return the month letter with year below it
    if (format(date, "%m") == "01" || format(date, "%m") == "07") {
      return(paste(month_letter, "\n", year))  # Display "J" or "J" with the year below for January and July
    } else {
      return(month_letter)  # For other months, show only the first letter of the month
    }
  })
  
  return(labels)
}

ggplot(time_series_data1_rt, aes(x = time, y = values)) +
  geom_line(color = "#0C7BDC", linewidth = 0.8, alpha = 0.6) +
  geom_hline(yintercept = 1, color = "black") +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "#0C7BDC", alpha = 0.4) +
  labs(
    title = NULL,
    x = NULL,
    y = "Rt"
  ) +
  scale_x_date(
    limits = c(as.Date("2021-07-07"), as.Date("2023-04-27")),  # Set the x-axis limits
    breaks = seq(as.Date("2021-07-07"), as.Date("2023-04-27"), by = "1 month"),  # Breaks every month
    labels = custom_labels  # Custom label function: First letter in uppercase
  ) +
  scale_y_continuous(limits = c(0, 3), breaks = seq(0, 3, 0.5)) +
  theme(
    axis.text = element_text(size = 28, face = "bold"),  # Axis text styling
    axis.title = element_text(size = 28, face = "bold"),  # Axis title styling
    plot.title = element_text(size = 30, face = "bold"),  # Plot title styling
    plot.background = element_rect(fill = "white", color = "white"),  # White plot background
    panel.grid.major = element_line(color = "grey", size = 0.5),  # Major grid lines in white
    panel.grid.minor = element_line(color = "grey", size = 0.25),  # Minor grid lines in white
    panel.background = element_rect(fill = "white", color = "white"),  # White background for the plot area
    axis.text.x = element_text(angle = 45, hjust = 1),  # Adjust labels' angle, horizontal and vertical alignment
    axis.title.x = element_text(size = 16, face = "bold"),
    #axis.text.x = element_blank(),  # Removes x-axis labels
    #axis.title.x = element_blank()  # Removes x-axis title
  )

###############################
### Plot of weekly_effect
time_series_data1_weekly_effect <- expand.grid(
  time = seq(as.Date("2021-07-07"), by = "weeks", length.out = W-16),
  Age = 1:7
)

time_series_data1_weekly_effect$values <- as.vector(weekly_effect[17:W, 1:7])

custom_colors <- c(
  "1" = "#0C7BDC",  # Color for Age group 1
  "2" = "#FFC20A",  # Color for Age group 2
  "3" = "#40B0A6",  # Color for Age group 3
  "4" = "#E66100",  # Color for Age group 4
  "5" = "#5D3A9B",  # Color for Age group 5
  "6" = "#994F00",  # Color for Age group 6
  "7" = "#D35FB7"   # Color for Age group 7
)

ggplot(time_series_data1_weekly_effect, aes(x = time, y = values, color = as.factor(Age))) +
  geom_line(linewidth = 0.8, alpha = 0.6) +
  labs(
    title = NULL,
    x = NULL,
    y = "Weekly effects in time-varying contact matrix",
    color = "Age Group"  # Label the legend for the color
  ) +
  scale_x_date(
    limits = c(as.Date("2021-07-07"), as.Date("2023-04-19")),
    breaks = seq(as.Date("2021-07-07"), as.Date("2023-04-19"), by = "3 months"),
    date_labels = "%b %Y"  # Format date labels to show month and year
  ) +
  scale_y_continuous(limits = c(0, 4.5), breaks = seq(0, 4.5, 0.5)) +
  scale_color_manual(values = custom_colors) +  # Apply the custom color palette
  theme_minimal() +
  theme(
    axis.text = element_text(size = 14, face = "bold"),  # Axis text styling
    axis.title = element_text(size = 15, face = "bold"),  # Axis title styling
    plot.background = element_rect(fill = "white", color = "white"),  # White plot background
    legend.title = element_text(size = 14),  # Legend title styling
    legend.text = element_text(size = 12)  # Legend text styling
  )

###############################
### Plot of proportion of infection, Singapore average 2021-2023
# Example data for the bar chart
# TimeVaryInfeByAge

SumCase = rep(0, A)
LOWERCI = rep(0, A)
UPPERCI = rep(0, A)

SumCase = apply(TimeVaryInfeByAge[(Seed+TVaccine):(Seed+N),7,], 2, mean)
SdCase = apply(TimeVaryInfeByAge[(Seed+TVaccine):(Seed+N),7,], 2, sd)

for (a in 1:A){
  UPPERCI = SumCase + 1.96*SdCase/sqrt(dim(TimeVaryInfeByAge[(Seed+TVaccine):(Seed+N),7,])[1])
  LOWERCI = SumCase - 1.96*SdCase/sqrt(dim(TimeVaryInfeByAge[(Seed+TVaccine):(Seed+N),7,])[1])
}
#714

# Example data for the bar chart
data <- data.frame(
  Category = c("0-14", "15-29", "30-39", "40-49", "50-59", "60-69", "70+"),
  Proportions_of_infections_after_intervention_of_vaccine = SumCase / sum(SumCase),
  Share_in_population = popByAge, 
  LowerCI = LOWERCI / sum(SumCase),
  UpperCI = UPPERCI / sum(SumCase)
)

# Reshape data for side-by-side bars using pivot_longer
data_long <- pivot_longer(data, 
                          cols = c(Proportions_of_infections_after_intervention_of_vaccine, Share_in_population), 
                          names_to = "Bar_Type", 
                          values_to = "Value")

# Create the bar chart with two groups: one with blue and one with white background
p <- ggplot(data_long, aes(x = Category, y = Value, fill = Bar_Type)) +
  # Add white bars for Share_in_population
  geom_bar(data = subset(data_long, Bar_Type == "Share_in_population"),
           stat = "identity", 
           color = "black", fill = "white", size = 1, 
           position = position_dodge(width = 0.8),
           show.legend = FALSE) +  # Black border for white bars
  # Add blue bars for Proportions_of_infections_before_intervention_of_vaccine
  geom_bar(data = subset(data_long, Bar_Type == "Proportions_of_infections_after_intervention_of_vaccine"),
           stat = "identity", 
           color = "transparent", fill = "#56B4D3", size = 1,
           position = position_dodge(width = 0.8),
           show.legend = TRUE) +  # Light ocean blue color
  geom_bar(data = subset(data_long, Bar_Type == "Share_in_population"),
           stat = "identity", 
           color = "black", fill = "transparent", size = 1,
           position = position_dodge(width = 0.8),
           show.legend = FALSE) +  # Black border for white bars
  geom_errorbar(aes(ymin = LowerCI, ymax = UpperCI), width = 0.4, color = "black") +  # Error bars for CI
  # Labels and title
  labs(title = "Estimated sources of SARS-CoV-2 cases, Singapore average from 7 July 2021 until 30 April 2023", 
       x = "Age group", y = "Proportions of infections for age group [70+] caused by (%)") + 
  scale_fill_manual(values = c("Proportions_of_infections_after_intervention_of_vaccine" = "#56B4D3", "Share_in_population" = "white")) + 
  # Theme adjustments
  theme(
    axis.text = element_text(size = 18, face = "bold"),   # Enlarge axis tick labels
    axis.title = element_text(size = 18, face = "bold"),  # Enlarge axis labels
    plot.title = element_text(size = 0, face = "bold"),  # Enlarge plot title
    plot.margin = margin(10, 10, 50, 10),  # Increase bottom margin for annotation space
    
    # Set background color to white
    plot.background = element_rect(fill = "white", color = "white"),  # Set entire plot background to white
    panel.background = element_rect(fill = "white", color = "white"), # Set the panel background to white
    
    # Customize the grid lines (keep them visible, but change their appearance)
    panel.grid.major = element_line(color = "grey80", size = 0.5),  # Major grid lines with a lighter grey color and thinner size
    panel.grid.minor = element_line(color = "grey90", size = 0.25), # Minor grid lines with an even lighter grey color
    
    legend.position = "bottom",  # Move the legend to the bottom
    legend.title = element_text(size = 18, face = "bold"),  # Customize legend title size and boldness
    legend.text = element_text(size = 18, face = "bold", color = "black")  # Customize legend text size and color
  )
print(p)
