# This script computes the following operating characteristics for the augmented simulation results
# outputted by the "interim_pval_fn" function in the Simulate_FreqGS.R script: type I error, power,
# the number of people in the trial (n), the number of people in the experimental (nE) and control (nC) 
# arms, the number of survivors in the experimental (yE) and control (yC) arms, the number of deaths 
# in the trial, and the overall number of deaths in the fixed cohort population (pop), assuming the 
# recommended treatment at the end of each trial is given to the remaining cohort population
# Both Pocock (aggressive early stopping) and OBF (conservative early stopping) stopping boundaries
# are considered in conjunction with interim analyses K = 3, 5, and 10
# pop = the maximum sample size from the 12 group sequential designs for a given treatment effect (OBF with K=3,5,10,
# Pocock with K=3,5,10, and frequentist and Bayesian trial designs)

setwd("C:/Users/jenni/OneDrive/Documents/Research/Mortality Comparison_RAR vs. Freq Design/Frequentist GS Designs")
source("Final Code/GSdesign nmax and stop bounds/Simulate Frequentist Design/Simulate_FreqGS.R")

###################################################################################################
########################## Function for Alternative Scenario Simulations ##########################
###################################################################################################

Freq_ComputeStats_ALT = function(altSims, alt.RR, null.RR, bounds, pop){
   
  #convert z-score upper stopping boundaries to two-sided p-values
  p.bounds = (1-pnorm(bounds))*2 
  
  foo = sapply(altSims, function(w){
      
      #determine when each trial is stopped
      stop.at = which.max(w[,5] <= c(p.bounds[-length(p.bounds)],1))

      #obtain desired statistics when trial stopped
      nE = w[stop.at,1]; nC = w[stop.at,3]; yE = w[stop.at,2]; yC = w[stop.at,4]; n = nE+nC; deaths.trial = n-(yE+yC)
      
      #determine number of deaths in fixed population
      Remaining = pop - n
      if((sum(w[,5] < p.bounds) > 0) & w[stop.at,6] == 1){ 
        #when treatment is deemed efficacious -> apply treatment to remaining pop
        deaths.cohort = rbinom(1, Remaining, prob=(1-alt.RR)) + deaths.trial
      }else{ 
        #when treatment is deemed harmful or null is accepted -> apply control to remaining pop
        deaths.cohort = rbinom(1, Remaining, prob=(1-null.RR)) + deaths.trial
      }
  return(c(n, nE, nC, yE, yC, deaths.trial, deaths.cohort))})
  
  #compute power
  Power = mean(sapply(altSims,function(x) max(x[,5] < p.bounds) == 1 & x[which.max(x[,5] <= c(p.bounds[-length(p.bounds)],1)),6] == 1))
 
  #Prepare results for output
  rownames(foo) = c("Trial N","Treatment N","Control N","Treatment Survivors",
                    "Control Survivors","Trial Deaths","Cohort Deaths")
  colnames(foo) = NULL
  foo = list(foo,Power); names(foo) = c("Stats","Power")
  return(foo)
}

###################################################################################################
########################## Function for Null Scenario Simulations #################################
###################################################################################################

Freq_ComputeStats_NULL = function(nullSims, null.RR, bounds, pop){
  
  #convert z-score upper stopping boundaries to two-sided p-values
  p.bounds = (1-pnorm(bounds))*2 
  
  foo = sapply(nullSims, function(w){
    
    #determine when each trial is stopped
    stop.at = which.max(w[,5] <= c(p.bounds[-length(p.bounds)],1))
    
    #obtain desired statistics when trial stopped
    nE = w[stop.at,1]; nC = w[stop.at,3]; yE = w[stop.at,2]; yC = w[stop.at,4]; n = nE+nC; deaths.trial = n-(yE+yC)
    
    #determine number of deaths in fixed population
    Remaining = pop - n
    
    #In the null scenario, the death rate among the remaining population will be the same, 
    #regardless of which treatment, if any, is erroneously deemed efficacious
    deaths.cohort = rbinom(1, Remaining, prob=(1-null.RR)) + deaths.trial
    return(c(n, nE, nC, yE, yC, deaths.trial, deaths.cohort))})
  
  #compute type I error
  TIE = mean(sapply(nullSims,function(x) max(x[,5] < p.bounds) == 1))
  
  #Prepare results for output
  rownames(foo) = c("Trial N","Treatment N","Control N","Treatment Survivors",
                    "Control Survivors","Trial Deaths","Cohort Deaths")
  colnames(foo) = NULL
  foo = list(foo,TIE); names(foo) = c("Stats","Type I Error")
  return(foo)
}


###################################################################################################
########################################## Example ################################################
###################################################################################################

#Obtain the augmented simulation results using the "simulate_FreqGS" and "interim_pval_fn" functions,
#and then obtain operating characteristics using the "Freq_ComputeStats_NULL" and "Freq_ComputeStats_ALT"
#functions

ntrials = 10
set.seed(1994)
Example_alt = lapply(1:ntrials, function(trial) simulate.gs.trial(response.probs = c(0.37,0.12), nmax = 153, k = 10))
Example_null = lapply(1:ntrials, function(trial) simulate.gs.trial(response.probs = c(0.12,0.12), nmax = 153, k = 10))
Example_pval_alt = interim_pval_fn(sim = Example_alt, interim.test = "prop.test")
Example_pval_null = interim_pval_fn(sim = Example_null, interim.test = "prop.test")

Example_res_alt = Freq_ComputeStats_ALT(altSims = Example_pval_alt, alt.RR = 0.37, null.RR = 0.12, 
                                        bounds = rep(2.55013,10), pop = 179)
Example_res_null = Freq_ComputeStats_NULL(nullSims = Example_pval_null, null.RR = 0.12, 
                                     bounds = rep(2.55013,10), pop = 179)
