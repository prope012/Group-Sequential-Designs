# This script simulates a classical frequentist group sequential (gs) clinical trial using 1:1 allocation. 
# A trial is stopped early if the two-sided p-value from the chosen two-sample binomial proportion test 
# is smaller than the two-sided p-value corresponding to the z-score stopping boundary at a given interim analysis. 
# The following metrics are computed:
# trial sample size, number of participants assigned to the treatment arm, number of participants 
# assigned to the control arm, number of responders in the treatment arm, number of responders in 
# the control arm, total number of non-responders in the actual trial, total number of non-responders 
# in the potential study sample, a binary indicator for whether the trial was stopped early for 
# efficacy, and a binary indicator for whether the trial was stopped early for harm or futility 
# as appropriate.


###################################################################################################
######### Function to Simulate a Frequentist GS Trial Using Potential Outcomes Framework ##########
###################################################################################################

### The function below has the following inputs:
###   - potential.outcomes = a dataframe with 2 columns containing the potential outcomes on the 
###     treatment (column 1) and control (column 2) for each individual in the potential study sample 
###     (generated using the function "Run.Trial.Get.Stats.Freq" below)
###   - nmax = the maximum sample size for the trial (determined using the gsDesign package)
###   - j = the number of interim analyses for the trial
###   - rand.prob = the randomization probability for the treatment arm
###   - max.deviation = alpha parameter in the mass-weighted urn design randomization scheme

### The function below has the following outputs:
###   - a jx4 matrix containing nE,nC,yE,yC at each interim analysis where nE and nC are the number of 
###     trial participants on the treatment and control, respectively, and yE and yC are the number of 
###     responders on the treatment and control, respectively.

### Purpose of function:
###   - This function simulates a 2-arm gs clinical trial with no early stopping (i.e. stopping 
###     boundaries are not yet applied) using the mass-weighted urn randomization scheme and 1:1 allocation
###   - Rather than generating outcomes for each participant within the trial, the outcomes for each 
###     participant are set equal to their potential outcome for the arm to which they were assigned

simulate.freq.trial.PO <- function(potential.outcomes, nmax, j, rand.prob=0.5, max.deviation = 3){
  
  #generate group sizes
  ns = round(seq(nmax/j,nmax,by=nmax/j)) 
  
  #Initialize Data
  n = c(0,0); y = c(0,0);
  
  #Storage object for number assigned to each arm, and number of survivors on each arm
  stats = matrix(NA,nrow=length(ns),ncol=4)
  colnames(stats) = c("nE","yE","nC","yC")
  rownames(stats) = 1:length(ns)
  
  #Generate data using MWUD and true response probabilities
  for(group in 1:nrow(stats)){
    
    #Number of new enrollees during current group
    n.new = c(0,ns)[group+1]-c(0,ns)[group]
    
    #Grab potential outcomes for current group
    start.row = c(0,ns)[group]+1
    end.row = c(0,ns)[group+1]
    p.outcomes = potential.outcomes[start.row:end.row,]
    
    #Generate assignments for current interval assuming a fixed, equal target allocation (1:1)
    assignments = NULL; rand.prob.temp = rand.prob
    for(i in 1:n.new){
      assignments[i] = rbinom(1,1,min(max(rand.prob.temp,0),1))
      rand.prob.temp = rand.prob.temp+(rand.prob-1*(assignments[i]==1))/max.deviation
    }
    nE.new = sum(assignments)
    nC.new = n.new - nE.new
    
    #Determine outcomes assuming consistency (Y_i = Y^{Zi}, where subject i was assigned treatment Z)
    outcomes = NULL
    for(i in 1:n.new){
      if(assignments[i]==1){outcomes[i] = p.outcomes[i,"potential_Ye"]} 
      else{outcomes[i] = p.outcomes[i,"potential_Yc"]}
    }
    yE.new = sum(assignments==1 & outcomes==1)
    yC.new = sum(assignments==0 & outcomes==1)
    
    #Update dataset
    n = n+c(nE.new,nC.new)
    y = y+c(yE.new,yC.new)
    
    #Store Relevant Statistics
    stats[group,] = c(n[1],y[1],n[2],y[2])  
  }
  return(stats)
}


###################################################################################################
##### Function to Perform 2-Sample Binomial Proportion Test Using Trial Data (Interim Analyses)####
#################### and obtain desired statistics from the trial #################################
###################################################################################################

### The function below has the following inputs:
###   - trial = the trial stats from the "simulate.freq.trial.PO" function
###   - interim.test = the 2-sample proportion test to be computed at each interim analysis (Fisher's exact test, Chi-square test, or the unconditional exact test)
###   - potential.outcomes = a dataframe containing the potential outcomes for each potential participant on the treatment and control (generated using the "Run.Trial.Get.Stats.Freq" function below)
###   - u.bound = the upper stopping boundary for the trial (found using the "Run.Trial.Get.Stats.Freq" function below)
###   - l.bound = the lower stopping boundary for the trial (found using the "Run.Trial.Get.Stats.Freq" function below)
###   - symmetric.bounds = TRUE or FALSE; indicates whether symmetric or asymmetric boundaries are being used
###   - ns = the sample size at each interim analysis (generally equal to round(seq(nmax/k,nmax,by=nmax/k))) 
###   - j = the number of interim analyses
###   - cohort.size = the size of the potential study sample for the target effect size
###   - response.probs = c(treatment arm response rate, control arm response rate)

### The function below has the following outputs:
###   - n = number of trial participants after stopping boundaries have been applied
###   - nE and nC = number of participants respectively assigned to the treatment and control arms after stopping boundaries have been applied 
###   - yE and yC = number of trial responders respectively on the treatment and control arms after stopping boundaries have been applied
###   - deaths.trial = total number of non-responders within the trial after stopping boundaries have been applied
###   - deaths.cohort = total number of non-responders within the potential study sample assuming the individuals in 
###                     the potential study sample but not in the trial itself receive the treatment recommended by the trial
###   - efficacy = a binary indicator for whether the trial stopped early due to treatment efficacy
###   - harm = a binary indicator for whether the trial stopped early due to treatment harm (for symmetric boundaries only)
###   - futility = a binary indicator for whether the trial stopped early due to treatment futility (for asymmetric boundaries only)

### Purpose of function: 
###   - This function computes the desired 2-sample binomial proportion test at each interim analysis using
###     the trial data outputted by the "simulate.freq.trial.PO" function
###   - This function then determines if the Z-score from the 2-sample binomial proportion test exceeded a stopping boundary
###   - Within-trial statistics (nE, nC, yE, yC, non-responders within trial) are reported at the time the trial was stopped
###   - The number of non-responders within the potential study sample is found.

library(gsDesign)
library(exact2x2)

prop_test_get_stats_fn = function(trial, interim.test=c("fisher","chisquare","uncondExact"), potential.outcomes,
                                  l.bound, u.bound, symmetric.bounds = TRUE, ns, j, cohort.size, response.probs){
  
  if(interim.test == "chisquare"){ #obtain p-value and test sign using prop.test
    pval = apply(trial,1,function(v){
      res = prop.test(x=c(v[2],v[4]), n=c(v[1],v[3]),alternative = "two.sided")
      tsign = sign(res$estimate[1] - res$estimate[2])
      pvalue = ifelse(is.nan(res$p.value)==TRUE,1,res$p.value) #prop.test returns "NaN" when there are 0 successes in both groups
      return(c(pvalue,tsign))})}
  
  if(interim.test == "uncondExact"){ #obtain p-value and test sign using unconditional exact test
    pval = apply(trial,1,function(v){
      res = uncondExact2x2(x1=v[4],n1=v[3],x2=v[2],n2=v[1],parmtype = "difference",method="FisherAdj",midp = FALSE)
      pvalue = res$p.value
      tsign = sign(res$estimate)
      return(c(pvalue,tsign))})}
  
  if(interim.test == "fisher"){ #obtain p-value and test sign using fisher's exact test
    pval = apply(trial,1,function(v){
      res = fisher.test(x = matrix(c(v[2],v[1]-v[2],v[4],v[3]-v[4]),ncol=2), alternative = "two.sided")
      pvalue = res$p.value
      tsign = sign(log(res$estimate))
      return(c(pvalue,tsign))})} 
  
  trial = cbind(trial,t(data.frame(pval)))
  colnames(trial)[5:6] = c("pvalue","tsign")
  
  # convert z-score upper stopping boundary to two-sided p-values for symmetric boundaries
  if(symmetric.bounds == TRUE){p.bounds = (1-pnorm(u.bound))*2}
  # convert z-score stopping boundaries to one-sided p-values for asymmetric boundaries
  if(symmetric.bounds == FALSE){
    p.u.bounds <- 1-pnorm(u.bound)
    p.l.bounds <- ifelse(sign(l.bound) == 1, 1-pnorm(l.bound), pnorm(l.bound))
  }
  
  # compute desired stats from trial
  if(symmetric.bounds == TRUE){stop.at = which.max(trial[,"pvalue"] <=  c(p.bounds[-length(p.bounds)],1))}
  if(symmetric.bounds == FALSE){
    
    # obtain 1-sided p-value
    trial[,"pvalue"] <- trial[,"pvalue"]/2
    
    # determine when boundaries were crossed
    stop.high <- (trial[,"pvalue"] <= c(p.u.bounds[-length(p.u.bounds)],1)) & (trial[,"tsign"] == 1)
    stop.low <- ifelse(sign(l.bound) == 1, ((trial[,"pvalue"] > p.l.bounds) & (trial[,"tsign"] == 1)) | (trial[,"tsign"] == -1), (trial[,"pvalue"] <= p.l.bounds) & (trial[,"tsign"] == -1))
    stop.at = which.max(stop.high | stop.low)
  }
  nE = trial[stop.at,"nE"]; nC = trial[stop.at,"nC"]; yE = trial[stop.at,"yE"]; yC = trial[stop.at,"yC"]; 
  n = nE+nC; deaths.trial = n-(yE+yC); maxN = max(ns)
  
  # determine whether trial stopped due to efficacy, harm, or futility
  if(symmetric.bounds == TRUE){
    efficacy = harm = 0 #null accepted
    if(sum(trial[,"pvalue"] <= p.bounds)>0 & trial[stop.at,"tsign"]==1) {efficacy=1; harm=0} #treatment efficacious
    if(sum(trial[,"pvalue"] <= p.bounds)>0 & trial[stop.at,"tsign"]==-1)  {efficacy=0; harm=1} #treatment harmful
  }
  if(symmetric.bounds == FALSE){
    cross.upper <- min(which(stop.high==TRUE))
    cross.lower <- min(which(stop.low==TRUE))
    efficacy = futility = 0
    if(cross.upper < cross.lower){efficacy = 1}else{futility = 1}
  }
  
  # determine number of deaths in the potential cohort, taking into account the appropriate hypothesis
  deaths.cohort = deaths.trial + (cohort.size-n)*(1-response.probs[2]-efficacy*(response.probs[1]-response.probs[2]))
  
  #return appropriate information
  if(symmetric.bounds == T){return(c(n, nE, nC, yE, yC, deaths.trial, deaths.cohort, efficacy, harm))}
  if(symmetric.bounds == F){return(c(n, nE, nC, yE, yC, deaths.trial, deaths.cohort, efficacy, futility))}
}

### Note: Power is equal to the proportion of trials where efficacy = 1 under alternative (for symmetric and asymmetric boundaries)
### Note: TIE is equal to the proportion of trials where efficacy = 1 or harm = 1 under null (for symmetric boundaries)
### Note: TIE is equal to the proportion of trials where efficacy = 1 under null (for asymmetric boundaries)


###################################################################################################
############ Function to Run a Trial and Get Stats Under Potential Outcomes Framework #############
###################################################################################################

### The function below has the following inputs:
###   - response.probs = c(treatment arm response rate, control arm response rate)
###   - cohort.size = the size of the potential study sample for the target effect size
###   - j = the number of interim analyses
###   - bound.type = "OF","Pocock","asymmetric-OF","asymmetric-PK"; specifies whether to use symmetric OBF-like or Pocock-like stopping boundaries, an 
###                   OBF-like upper boundary with a binding lower futility boundary, or a Pocock-like upper boundary with a binding lower futility boundary###   - interim.test = the 2-sample proportion test to be computed at each interim analysis (Fisher's exact test, Chi-square test, or the unconditional exact test)
###   - interim.test = the 2-sample proportion test to be computed at each interim analysis (Fisher's exact test, Chi-square test, or the unconditional exact test)
###   - symmetric.bounds = TRUE or FALSE; indicates whether symmetric or asymmetric boundaries are being used

### This function has the following outputs:
###   - n = number of trial participants after stopping boundaries have been applied
###   - nE and nC = number of participants respectively assigned to the treatment and control arms after stopping boundaries have been applied 
###   - yE and yC = number of trial responders respectively on the treatment and control arms after stopping boundaries have been applied
###   - deaths.trial = total number of non-responders within the trial after stopping boundaries have been applied
###   - deaths.cohort = total number of non-responders within the potential study sample assuming the individuals in 
###                     the potential study sample but not in the trial itself receive the treatment recommended by the trial
###   - efficacy = a binary indicator for whether the trial stopped early due to treatment efficacy
###   - harm = a binary indicator for whether the trial stopped early due to treatment harm (for symmetric boundaries only)
###   - futility = a binary indicator for whether the trial stopped early due to treatment futility (for asymmetric boundaries only)

### Purpose of Function:
###   - This function generates the potential outcomes for each individual in the potential study 
###     sample using the response.probs vector & cohort.size
###   - This function then runs the "simulate.freq.trial.PO" and "prop_test_get_stats_fn" functions to get
###     all desired sim stats

Run.Trial.Get.Stats.Freq = function(response.probs, bound.probs, cohort.size, j, bound.type, symmetric.bounds = TRUE,
                                    interim.test=c("fisher","chisquare","uncondExact")){
  
  #Generate Potential Outcomes
  potential.outcomes =  data.frame(potential_Ye = rbinom(cohort.size, 1, prob = response.probs[1]),
                                   potential_Yc = rbinom(cohort.size, 1, prob = response.probs[2]))

  # determine group size, max N, and stopping boundaries needed for 90% power, one-sided alpha = 0.025
  if(symmetric.bounds == TRUE){
    fixed.ss = 2*(power.prop.test(n = NULL, p1 = bound.probs[1], p2 = bound.probs[2], sig.level = 0.05, power = 0.90, alternative = "two.sided"))$n
    GSdesign = gsDesign(k=j, test.type=2, alpha = 0.025, beta = 0.10, sfu = bound.type, n.fix = fixed.ss)
    ns = ceiling(GSdesign$n.I)
    u.bound = GSdesign$upper$bound
    l.bound = GSdesign$lower$bound
  }
  if(symmetric.bounds == FALSE & bound.type == "asymmetric-OF"){
  
    ### function to compute type I error rate for gs design
    t1e = function(c,j){
      t=1:j/j; b=c*sqrt(1/t); a = c*(2-sqrt(1/t)) 
      sum(gsProbability(k=j,n.I=t,a=a,b=b,theta=0)$upper$prob)
    } 
    
    ### function to get multiplicative constant a that controls one-sided type I error at 0.025
    get.constant = function(alpha=0.025) uniroot(function(x) alpha-t1e(c=x,j=j),interval=c(0,4))$root
    a <- get.constant()
    u.bound <- sqrt(j/1:j)*a
    l.bound <- a*(2-sqrt(j/1:j))
    
    ### function to compute type II error rate for same design that controls type I error rate at 0.025
    t2e = function(theta,alpha=0.025,j){
      c=get.constant(alpha); t=1:j/j; b=c*sqrt(1/t); a = c*(2-sqrt(1/t)) 
      sum(gsProbability(k=j,n.I=t,a=a,b=b,theta=theta)$lower$prob)
    } 
    
    ### function to get sample size inflation factor (R = n.gs/n.fix)
    get.inf.factor = function(beta=0.10,alpha=0.025,j) (uniroot(function(x) beta-t2e(theta=x,alpha=alpha,j=j),interval=c(0,10))$root/(qnorm(1-alpha)+qnorm(1-beta)))^2
    R = get.inf.factor(j=j)
    fixed.ss = 2*(power.prop.test(n = NULL, p1 = bound.probs[1], p2 = bound.probs[2], sig.level = 0.05, power = 0.90, alternative = "two.sided"))$n
    gs.ss <- fixed.ss*R
    ns <- round(seq(gs.ss/j, gs.ss, gs.ss/j))
  }
  
  if(symmetric.bounds == FALSE & bound.type == "asymmetric-PK"){
    
    ### function to compute type I error rate for gs design
    t1e = function(c,j){
      t=1:j/j; b=rep(c,j); a = c*(2-sqrt(1/t)) 
      sum(gsProbability(k=j,n.I=t,a=a,b=b,theta=0)$upper$prob)
    } 
    
    ### function to get multiplicative constant a that controls one-sided type I error at 0.025
    get.constant = function(alpha=0.025) uniroot(function(x) alpha-t1e(c=x, j=j),interval=c(0,4))$root
    a <- get.constant()
    u.bound <- rep(a,j)
    l.bound <- a*(2-sqrt(j/1:j))
    
    ### function to compute type II error rate for same design that controls type I error rate at 0.025
    t2e = function(theta,alpha=0.025,j){
      c=get.constant(alpha); t=1:j/j; b=rep(c,j); a = c*(2-sqrt(1/t)) 
      sum(gsProbability(k=j,n.I=t,a=a,b=b,theta=theta)$lower$prob)
    } 
    
    ### function to get sample size inflation factor (R = n.gs/n.fix)
    get.inf.factor = function(beta=0.10,alpha=0.025,j) (uniroot(function(x) beta-t2e(theta=x,alpha=alpha,j=j),interval=c(0,10))$root/(qnorm(1-alpha)+qnorm(1-beta)))^2
    R = get.inf.factor(j=j)
    fixed.ss = 2*(power.prop.test(n = NULL, p1 = bound.probs[1], p2 = bound.probs[2], sig.level = 0.05, power = 0.90, alternative = "two.sided"))$n
    gs.ss <- fixed.ss*R
    ns <- round(seq(gs.ss/j, gs.ss, gs.ss/j))
  }
  
  #Simulate Trial
  trial = simulate.freq.trial.PO(potential.outcomes = potential.outcomes, nmax = max(ns), j=j)
  
  #Compute Stats using user-specified 2-sample binomial proportion test
  stats = prop_test_get_stats_fn(trial=trial, interim.test = interim.test, potential.outcomes = potential.outcomes,
                                 l.bound = l.bound, u.bound = u.bound, symmetric.bounds = symmetric.bounds, ns=ns, j=j, 
                                 cohort.size=cohort.size, response.probs=response.probs)
  
  #Label stats as appropriate
  if(symmetric.bounds == TRUE){
  names(stats) =  c("Trial N","Treatment N","Control N","Treatment Survivors",
                    "Control Survivors","Trial Deaths","Cohort Deaths","Efficacy","Harm")
  }
  if(symmetric.bounds == FALSE){
    names(stats) =  c("Trial N","Treatment N","Control N","Treatment Survivors",
                      "Control Survivors","Trial Deaths","Cohort Deaths","Efficacy","Futility")
  }
  return(stats)
}

###################################################################################################
############################## Example Using the Above Functions ##################################
###################################################################################################

### simulate 1000 trials and get stats - ALT
Ex.Freq.alt = data.frame(do.call(rbind, lapply(1:1000, function(x){
  set.seed(x)
  Run.Trial.Get.Stats.Freq(response.probs = c(0.37,0.12), bound.probs = c(0.37,0.12), cohort.size=447, j=10, bound.type = "OF",
                           symmetric.bounds = TRUE, interim.test = "uncondExact")})))

### simulate 1000 trials and get stats - NULL
Ex.Freq.null = data.frame(do.call(rbind, lapply(1:1000, function(x){
  set.seed(x)
  Run.Trial.Get.Stats.Freq(response.probs = c(0.12,0.12), bound.probs = c(0.37,0.12), cohort.size=447, j=10, bound.type = "OF",
                           symmetric.bounds = TRUE, interim.test = "uncondExact")})))