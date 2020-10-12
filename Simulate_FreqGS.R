# This script simulates a frequentist group sequential design employing 1:1 randomization within each block
# 3 different proportion tests (Prop.test, Barnard's test, fisher's exact test) are considered
# and used to test the null that the proportion of successes in the treatment and control groups
# are the same
# The function "interim_pval_fn" augments the simulation results from the "simulate.gs.trial"
# function with the p-values and the sign of the difference in proportions (+1 = the treatment group 
# had a greater success rate, -1 = the control group had a greater success rate)for each interim 
# analysis proportion test
# ns = sequential group sizes
# nmax = maximum sample size for each trial
# k = the number of interim analyses
# response.probs = (treatment arm response rate, control arm response rate)

###################################################################################################
###################### Simulate frequentist group sequential trial ################################
###################################################################################################

simulate.gs.trial <- function(response.probs, nmax, k){
  
  #generate group sizes
  ns = round(seq(nmax/k,nmax,by=nmax/k)) 
  
  #Initialize Data
  n = c(0,0); y = c(0,0);
  
  #Storage object for number assigned to each arm, and number of survivors on each arm
  stats = matrix(NA,nrow=length(ns),ncol=4)
  colnames(stats) = c("nE","yE","nC","yC")
  rownames(stats) = 1:length(ns)
  
  #Generate data under 1:1 randomization scheme and true response probabilities
  for(group in 1:nrow(stats)){
    #Number of new enrollees during current group
    n.new = c(0,ns)[group+1]-c(0,ns)[group]
    
    #Generate assignments and outcomes for current interval
    if(n.new %% 2 == 0){nE.new = n.new/2; nC.new = n.new/2} else{
      a = rbinom(1,1,0.5)
      nE.new = ceiling(n.new/2)*a + floor(n.new/2)*(1-a)
      nC.new = n.new - nE.new}
    
    yE.new = rbinom(1,nE.new,response.probs[1])
    yC.new = rbinom(1,nC.new,response.probs[2])
    
    #Update dataset
    n = n+c(nE.new,nC.new)
    y = y+c(yE.new,yC.new)
    
    #Store Relevant Statistics
    stats[group,] = c(n[1],y[1],n[2],y[2])  
  }
  
  return(stats)
}

###################################################################################################
###################### Compute proportion tests for each interim analyses #########################
###################################################################################################

library(Barnard)
interim_pval_fn = function(sim, interim.test){
  results = lapply(sim, function(w){
    
    if(interim.test == "prop.test"){ #obtain p-value and test sign using prop.test
      pval = apply(w,1,function(v){
        res = prop.test(x=c(v[2],v[4]), n=c(v[1],v[3]),alternative = "two.sided")
        tsign = sign(res$estimate[1] - res$estimate[2])
        pvalue = res$p.value
        pvalue = ifelse(is.nan(pvalue) == TRUE, 1, pvalue) #prop.test returns "NaN" when there are 0 successes in both groups
        return(c(pvalue,tsign))})}
    
    if(interim.test == "barnard"){ #obtain p-value and test sign using barnard's test
      pval = apply(w,1,function(v){
        if(v[2]+v[4] == 0){pvalue = 1; tsign = 1} else{ #barnard's test returns a warning message when there are 0 successes in both groups
          res = barnard.test(n1=v[2], n2=v[4], n3=(v[1]-v[2]), n4=(v[3]-v[4]))
          pvalue = res$p.value[2]
          tsign = (sign(as.numeric(res$statistic)))*(-1)}
          return(c(pvalue,tsign))})}
    
    if(interim.test == "fisher"){ #obtain p-value and test sign using fisher's exact test
      pval = apply(w,1,function(v){
        res = fisher.test(x = matrix(c(v[2],v[1]-v[2],v[4],v[3]-v[4]),ncol=2), alternative = "two.sided")
        pvalue = res$p.value
        tsign = sign(log(res$estimate))
        return(c(pvalue,tsign))})} 
    
    foo = cbind(w,t(data.frame(pval)))
    colnames(foo)[5:6] = c("pvalue","tsign")
    return(foo)})
  return(results)
}

###################################################################################################
########################################## Example ################################################
###################################################################################################

set.seed(1994)
Ex.Freq.Sims = lapply(1:100, function(trial) simulate.gs.trial(response.probs = c(0.27,0.12),nmax = 368, k = 10))
Ex.Freq.Res = interim_pval_fn(sim = Ex.Freq.Sims, interim.test = "fisher")
