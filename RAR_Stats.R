# The functions below take in null or alternative simulation results and output the number of people
# in each trial (n), the number of people in the experimental (nE) and control (nC) arms, the number 
# of survivors in the experimental (yE) and control (yC) arms, the number of deaths in the trial, and
# the overall number of deaths in the fixed cohort population (pop), assuming the recommended treatment 
# at the end of each trial is given to the remaining cohort population
# The "RAR_ComputeStats_ALT" and "RAR_ComputeStats_Null" functions compute power and type I error, respectively 
# pop = the maximum sample size from the 12 group sequential designs for a given treatment effect (OBF with K=3,5,10,
# Pocock with K=3,5,10, and frequentist and Bayesian trial designs)

###################################################################################################
########################## Function for Alternative Scenario Simulations ##########################
###################################################################################################

RAR_ComputeStats_ALT = function(altSims, alt.RR, null.RR, bounds, pop){

  foo = sapply(altSims, function(w){
    
    #determine when each trial is stopped
    stop.at = which.max((w[,1] > c(bounds[-length(bounds)],0)) | (w[,1] < (1-bounds)))
    
    #obtain desired statistics when trial stopped
    nE = w[stop.at,2]; nC = w[stop.at,4]; yE = w[stop.at,3]; yC = w[stop.at,5]; n = nE+nC; deaths.trial = n-(yE+yC)
    
    #determine number of deaths in the remaining cohort population
    Remaining = pop - n
    if(sum(w[,1] > bounds) > 0){ 
      #when treatment is deemed efficacious -> apply treatment to remaining pop
      deaths.cohort = rbinom(1, Remaining, prob=(1-alt.RR)) + deaths.trial
    }else{ 
      #when treatment is deemed harmful or null is accepted -> apply control to remaining pop
      deaths.cohort = rbinom(1, Remaining, prob=(1-null.RR)) + deaths.trial 
    }
  return(c(n, nE, nC, yE, yC, deaths.trial, deaths.cohort))})
   
  #compute power (proportion of trials under the alternative in which treatment deemed efficacious)
  Power = mean(sapply(altSims,function(x) max(x[,1] > bounds) == 1))
  
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

RAR_ComputeStats_NULL = function(nullSims, null.RR, bounds, pop){
  
  foo = sapply(nullSims, function(w){
    
    #determine when each trial is stopped
    stop.at = which.max((w[,1] > c(bounds[-length(bounds)],0)) | (w[,1] < (1-bounds)))
    
    #obtain desired statistics when trial stopped
    nE = w[stop.at,2]; nC = w[stop.at,4]; yE = w[stop.at,3]; yC = w[stop.at,5]; n = nE+nC; deaths.trial = n-(yE+yC)
    
    #determine number of deaths in the remaining cohort population
    Remaining = pop - n
    
    #In the null scenario, the death rate among the remaining population will be the same, regardless
    #of which treatment, if any, is erroneously deemed efficacious
    deaths.cohort = rbinom(1, Remaining, prob=(1-null.RR)) + deaths.trial 
      
  return(c(n, nE, nC, yE, yC, deaths.trial, deaths.cohort))})
  
  #compute type I error (proportion of trials under the null in which treatment or control 
  #erroneously deemed efficacious)
  TIE = mean(sapply(nullSims,function(x) max(x[,1] > bounds | x[,1] < 1-bounds) == 1))
  
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

set.seed(1994)
Ex.RAR.nmax = find.nmax(nmax=153,ntrials=100,null.RR=c(0.12,0.12),alt.RR=c(0.37,0.12),k=10,bound.type = "Pocock")
Ex_Null = RAR_ComputeStats_NULL(nullSims = Ex.RAR.nmax[[5]], null.RR = 0.12, bounds = Ex.RAR.nmax[[4]], pop = 179)
Ex_Alt = RAR_ComputeStats_ALT(altSims = Ex.RAR.nmax[[6]], null.RR = 0.12, alt.RR = 0.37, bounds = Ex.RAR.nmax[[4]], pop = 179)
