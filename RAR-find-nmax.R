# This script determines the maximum sample size (nmax) for a clinical trial implementing a 
# Bayesian response adaptive randomization (RAR) group-sequential (gs) design using Pocock-like 
# (aggressive early stopping) and OBF-like (conservative early stopping) symmetric stopping boundaries for
# a specific number of interim analyses K. If symmetric boundaries are undesirable, the user may implement 
# a binding lower futility boundary.
#   For symmetric boundaries:
#     1) Pocock-like boundaries: we will use the same posterior probability (PP) stopping boundary throughout 
#        the trial that controls two-sided type I error at 0.05.
#     2) OBF-like boundaries: we will use c_k = pnorm(sqrt(K/k)*a), where c_k is the value of the kth stopping 
#        boundary and a is an appropriate constant that leads to a two-sided type I error rate of 0.05; 
#        PP boundaries start high and get progressively smaller.
#   For asymmetric boundaries: 
#     1) Pocock-like upper boundary and binding lower futility boundary: we will set the upper boundary equal to
#        u_k = pnorm(a) and lower boundary equal to l_k = pnorm(2a-sqrt(K/k)*a) where a is an appropriate constant that 
#        leads to a one-sided type I error rate of 0.025.
#     2) OBF-like upper boundary and binding lower futility boundary: we will set the upper boundary equal to
#        u_k = pnorm(sqrt(K/k)*a) and lower boundary equal to l_k = pnorm(2a-sqrt(K/k)*a) where a is an appropriate constant that 
#        leads to a one-sided type I error rate of 0.025.
# nmax is further determined to obtain 90% power
# PP = Pr(Experimental Treatment > Control|data)
# Code uses a logistic regression probability model and mass-weighted urn randomization scheme

###################################################################################################
####################### Function to Employ Logistic Regression Probability Model ##################
###################################################################################################

### The function below has the following inputs:
###   - n = c(number of participants on treatment arm, number of participants on control arm)
###   - y = c(number of treatment responders, number of control responders)
###   - X = design matrix for the logistic regression model (default: rbind(c(1,0.5),c(1,-0.5)))
###   - df = degrees of freedom of the prior t-distribution for the intercept and slope 
###   - m0 = location parameter of the prior t-distribution for the intercept and slope 
###   - s0 = scale parameter of the prior t-distribution for the intercept and slope 
###   - nsamps = number of posterior samples using the Gibbs sampler
###   - warmup = number of posterior samples to be used as a burn-in period (not included in final posterior draws)

### The function below has the following outputs:
###   - a 2 x nsamps dataframe of posterior samples for beta0 and beta1

### Purpose of function:
###   - This function executes the logistic regression probability model and returns beta posterior samples

library(BayesLogit)
library(locfit)
X = rbind(c(1,0.5),c(1,-0.5)) #establish design matrix for tlr model

# This function returns beta posterior samples
tlr.sampler = function(n=NULL,y,X,df=rep(7,ncol(X)),m0=rep(0,ncol(X)),s0=rep(2.5,ncol(X)),nsamps=2000,warmup=100){
  
  #All trials are of size 1
  if(is.null(n)) n = rep(1,length(y))
  
  # Pseudo-outcome vector
  kappa = cbind(y-0.5*n)
  
  # Storage object
  beta.samps = matrix(NA,nrow=ncol(X),ncol=nsamps)
  
  # Initial values
  beta = rep(0,ncol(X))
  
  # Gibbs Sampler 
  for(samp in 1:(nsamps+warmup)){
    
    # Sample beta
    omega = rpg.devroye(ncol(X),n,as.vector(X%*%cbind(beta))) # generate auxillary PG parameters, i.e. omegas
    tau = rgamma(ncol(X),shape=(df+1)/2,rate=(df*s0^2+(beta-m0)^2)/2) # generate auxillary Gamma parameters, i.e. taus
    V = chol2inv(chol(t(X)%*%(X*omega) + diag(tau))) # Posterior variance matrix
    m = V%*%(t(X)%*%kappa+diag(tau)%*%cbind(m0)) # Posterior mean vector
    beta = as.vector(m + t(chol(V))%*%rnorm(ncol(X))) # generate regression parameters, i.e. beta
    
    # Save posterior sample 
    if(samp>warmup){
      beta.samps[,samp-warmup] = beta
    } 
  }
  
  # Return posterior samples
  return(beta.samps) 
}

###################################################################################################
############### Function to Simulate a Bayesian RAR Group Sequential Trial ########################
###################################################################################################

### The function below has the following inputs:
###   - response.probs = c(treatment arm response rate, control arm response rate)
###   - ns = the sample size at each interim analysis (generally equal to round(seq(nmax/k,nmax,by=nmax/k))) 
###   - max.ar = the maximum randomization probability to be used by the "Restrict" RAR design; default is 0.75
###   - rand.type = "Coin" for weighted coin flipping, "Urn" for mass-weighted urn design
###   - max.deviation = alpha parameter in the mass-weighted urn design randomization scheme
###   - model = "ibb" for independent beta-binomial probability model vs. "tlr" for logistic regression probabiltiy model
###   - PP.mod = specify which RAR adaptive modification to use ("restrict","exponent","neyman","optimal","DABCD-ney","DABCD-opt","ERADE-ney","ERADE-opt");
###              "exponent" refers to the modification that uses the tuning parameter
###   - pi.star = the prior mean for the ibb model
###   - pess = the prior effective sample size for the ibb model
###   - df = degrees of freedom of the prior t-distribution for the intercept and slope 
###   - m0 = location parameter of the prior t-distribution for the intercept and slope 
###   - s0 = scale parameter of the prior t-distribution for the intercept and slope 
###   - nsamps = number of posterior samples using the Gibbs sampler
###   - warmup = number of posterior samples to be used as a burn-in period (not included in final posterior draws)

### The function below has the following outputs: stats[group,] = c(post.prob,n[1],y[1],n[2],y[2])  
###   - a k x 5 matrix containing the posterior probability that the treatment is better than the control,
###     nE, nC, yE, and yC at each interim analysis where nE and nC are the number of 
###     trial participants on the treatment and control, respectively, and yE and yC are the number of 
###     responders on the treatment and control, respectively.

### Purpose of function:
###   - This function simulates a 2-arm Bayesian RAR gs clinical trial with no early stopping (i.e. stopping 
###     boundaries are not yet applied) using the mass-weighted urn randomization scheme

simulate.trial <- function(response.probs, ns, max.ar=0.75, rand.type=c("Coin","Urn","Block"), max.deviation=3, model=c("tlr","ibb"),
                           PP.mod = c("restrict","exponent","neyman","optimal","DABCD-ney","DABCD-opt","ERADE-ney","ERADE-opt"),
                           pi.star=0.5, pess=2, df=c(7,7), m0=c(log(0.12/0.88),0), s0=c(2.5,2.5), nsamps=2000, warmup=100){
  
  #Initialize Data
  n = c(0,0); y = c(0,0); rand.prob = 0.5
  
  #Storage object for posterior probabilities
  stats = matrix(NA,nrow=length(ns),ncol=5)
  colnames(stats) = c("PP","nE","yE","nC","yC")
  rownames(stats) = 1:length(ns)
  
  #Generate Data under Block AR scheme and true response probabilities
  for(group in 1:nrow(stats)){
    #Number of new enrollees during current group
    n.new = c(0,ns)[group+1]-c(0,ns)[group]
    
    #Generate assignments and outcomes for current interval
    if(rand.type == "Coin"){
      nE.new = rbinom(1,n.new,rand.prob)
      yE.new = rbinom(1,nE.new,response.probs[1])
      nC.new = n.new - nE.new
      yC.new = rbinom(1,nC.new,response.probs[2])
    }
    if(rand.type == "Urn"){
      z = NULL; rand.prob.temp = rand.prob
      for(i in 1:n.new){
        z[i] = rbinom(1,1,min(max(rand.prob.temp,0),1))
        rand.prob.temp = rand.prob.temp+(rand.prob-1*(z[i]==1))/max.deviation
      }
      nE.new = sum(z); nC.new = n.new - nE.new
      yE.new = rbinom(1,nE.new,response.probs[1])
      yC.new = rbinom(1,nC.new,response.probs[2])
    }
    if(rand.type == "Block"){
      u = rbinom(1,1,n.new*rand.prob - floor(n.new*rand.prob))
      nE.new = u*ceiling(n.new*rand.prob)+(1-u)*floor(n.new*rand.prob)
      nC.new = n.new-nE.new
      yE.new = rbinom(1,nE.new,response.probs[1])
      yC.new = rbinom(1,nC.new,response.probs[2])
    }
    #Update dataset
    n = n+c(nE.new,nC.new)
    y = y+c(yE.new,yC.new)
    
    #Update posterior probability that E is superior to C and randomization probability for next group
    if(model=="ibb") post.prob = get.post.prob(n=n,y=y,pi.star=pi.star,pess=pess)
    if(model=="tlr"){
      gibbs.draws = tlr.sampler(n=n,y=y,X=rbind(c(1,0.5),c(1,-0.5)),df=df,m0=m0,s0=s0,nsamps=nsamps,warmup=warmup)
      post.prob = rowMeans(gibbs.draws>0)[2]
    }
    
    #Determine updated rand prob by stabilizing the current post prob via user-specified modification
    C = 0.5*(group)/(length(ns)-1)
    if(PP.mod == "restrict") rand.prob = min(max.ar,max(1-max.ar,post.prob))
    if(PP.mod == "exponent") rand.prob = (post.prob^C)/(post.prob^C + (1-post.prob)^C)
    if(PP.mod == "neyman"){
      pE = expit(gibbs.draws[1,] + 0.5*gibbs.draws[2,]); qE = 1-pE
      pC = expit(gibbs.draws[1,] - 0.5*gibbs.draws[2,]); qC = 1-pC
      frac = sqrt(pE*qE)/(sqrt(pE*qE)+sqrt(pC*qC)) 
      frac = ifelse((pE==1 & pC==1)==1, 0.5, frac) # when pE & pC = 1, then 0/0 = NaN; set value to 0.5
      rand.prob = mean(frac) #posterior mean for Neyman allocation
      }
    if(PP.mod == "optimal"){
      pE = expit(gibbs.draws[1,] + 0.5*gibbs.draws[2,])
      pC = expit(gibbs.draws[1,] - 0.5*gibbs.draws[2,])
      rand.prob = mean(sqrt(pE)/(sqrt(pE)+sqrt(pC))) # posterior mean for optimal allocation
    }
    if(PP.mod == "DABCD-ney"){
      pE = expit(gibbs.draws[1,] + 0.5*gibbs.draws[2,]); qE = 1-pE
      pC = expit(gibbs.draws[1,] - 0.5*gibbs.draws[2,]); qC = 1-pC
      frac = sqrt(pE*qE)/(sqrt(pE*qE)+sqrt(pC*qC)) # Neyman allocation
      frac = ifelse((pE==1 & pC==1)>=1, 0.5, frac) # when pE & pC = 1, then 0/0 = NaN; set value to 0.5
      target = mean(frac) # posterior mean for Neyman allocation  
      realized.prop = n[1]/sum(n) # realized proportion of participants on the treatment arm
      alloc.fn.num = target*(target/realized.prop)^2 # numerator of Hu & Zhang allocation function
      alloc.fn.den = target*(target/realized.prop)^2 + (1-target)*((1-target)/(1-realized.prop))^2 # denominator of Hu & Zhang allocation function
      rand.prob = alloc.fn.num/alloc.fn.den 
    }
    if(PP.mod == "DABCD-opt"){
      pE = expit(gibbs.draws[1,] + 0.5*gibbs.draws[2,])
      pC = expit(gibbs.draws[1,] - 0.5*gibbs.draws[2,])
      target = mean(sqrt(pE)/(sqrt(pE)+sqrt(pC))) # posterior mean for optimal allocation
      realized.prop = n[1]/sum(n) # realized proportion of participants on the treatment arm
      alloc.fn.num = target*(target/realized.prop)^2 # numerator of Hu & Zhang allocation function
      alloc.fn.den = target*(target/realized.prop)^2 + (1-target)*((1-target)/(1-realized.prop))^2 # denominator of Hu & Zhang allocation function
      rand.prob = alloc.fn.num/alloc.fn.den 
    }    
    if(PP.mod == "ERADE-ney"){
      pE = expit(gibbs.draws[1,] + 0.5*gibbs.draws[2,]); qE = 1-pE
      pC = expit(gibbs.draws[1,] - 0.5*gibbs.draws[2,]); qC = 1-pC
      frac = sqrt(pE*qE)/(sqrt(pE*qE)+sqrt(pC*qC)) # Neyman allocation
      frac = ifelse((pE==1 & pC==1)>=1, 0.5, frac) # when pE & pC = 1, then 0/0 = NaN; set value to 0.5
      target = mean(frac) # posterior mean for Neyman allocation  
      realized.prop = n[1]/sum(n) # realized proportion of participants on the treatment arm
      if(realized.prop > target){rand.prob = (2/3)*target}
      if(realized.prop == target){rand.prob = target}
      if(realized.prop < target){rand.prob = 1-(2/3)*(1-target)}
    }
    if(PP.mod == "ERADE-opt"){
      pE = expit(gibbs.draws[1,] + 0.5*gibbs.draws[2,]); qE = 1-pE
      pC = expit(gibbs.draws[1,] - 0.5*gibbs.draws[2,]); qC = 1-pC
      target = mean(sqrt(pE)/(sqrt(pE)+sqrt(pC))) # posterior mean for optimal allocation
      realized.prop = n[1]/sum(n) # realized proportion of participants on the treatment arm
      if(realized.prop > target){rand.prob = (2/3)*target}
      if(realized.prop == target){rand.prob = target}
      if(realized.prop < target){rand.prob = 1-(2/3)*(1-target)}
    }
    
    #Store Relevant Statistics
    stats[group,] = c(post.prob,n[1],y[1],n[2],y[2])  
  }
  
  return(stats)
}

###################################################################################################
################################# Function to Compute Power #######################################
###################################################################################################

### This function has the following inputs:
###   - nmax = maximum trial sample size
###   - ntrials = number of simulated trials
###   - null.RR = response.probs vector under the null hypothesis
###   - alt.RR = response.probs vector under the alternative hypothesis
###   - k = number of interim analyses
###   - bound.type = "OBF","Pocock","asymmetric-OBF","asymmetric-PK"; specifies whether to use symmetric OBF-like or Pocock-like stopping boundaries, an 
###                   OBF-like upper boundary with a binding lower futility boundary, or a Pocock-like upper boundary with a binding lower futility boundary
###   - PP.mod = specify which RAR adaptive modification to use ("restrict","exponent","neyman","optimal","DABCD-ney","DABCD-opt","ERADE-ney","ERADE-opt");
###              "exponent" refers to the modification that uses the tuning parameter

### This function has the following outputs:
###   - t1e = type I error of design
###   - pow = power of design
###   - stop.bounds = stopping boundaries calibrating the 2-sided type I error rate at 0.05 (for symmetric boundaries)
###                   or 1-sided type I error rate at 0.025 (for asymmetric boundaries)
###   - null.sim = trial results for the null scenario simulations
###   - alt.sim = trial results for the alternative scenario simulations

### Purpose of function:
###   - This function finds the stopping boundaries that calibrate type I error at 0.05 (symmetric boundaries) or 0.025 (asymmetric boundaries)
###   - The power associated with the identified stopping boundaries is then found

find.pow = function(nmax,ntrials,null.RR,alt.RR,k,bound.type,PP.mod){
  
  #group size at current nmax
  ns = round(seq(nmax/k,nmax,by=nmax/k)) 
  
  #simulate null and alternative scenarios at current nmax
  null.sim = lapply(1:ntrials,function(trial) simulate.trial(response.probs = null.RR, ns = ns, rand.type = "Urn", model="tlr", PP.mod = PP.mod))
  alt.sim = lapply(1:ntrials,function(trial) simulate.trial(response.probs = alt.RR, ns = ns, rand.type = "Urn", model="tlr", PP.mod = PP.mod))
  
  #determine stopping boundaries such that type I error = 0.05
  if(bound.type == "Pocock"){
    options = seq(0.9700,1,by=0.0001)
    foo = rep(0,length(options));
    for(i in 1:length(options)){
      bounds = c(rep(options[i],length(ns)))
      foo[i] = mean(sapply(null.sim,function(x) max(x[,1] > bounds | x[,1] < 1-bounds) == 1))
    }
    loc = which.min(abs(foo-0.05)) #determine closest type I error value to 0.05
    t1e = foo[loc] #store type I error
    stop.bounds = c(rep(options[loc],length(ns))) #store corresponding stopping boundary
    
    #compute power at the determined stopping boundary
    pow = mean(sapply(alt.sim,function(x) max(x[,1] > stop.bounds) == 1))
  }
  
  if(bound.type == "OBF"){
    bounds = list()
    a = seq(1.8,2.7,by=0.0001)
    foo = rep(0,length(a))
    for(i in 1:length(a)){
      bounds[[i]] = pnorm(sqrt(k/1:k)*a[i])
      foo[i] = mean(sapply(null.sim,function(x) max(x[,1] > bounds[[i]] | x[,1] < 1-bounds[[i]]) == 1))
    }
    loc = which.min(abs(foo-0.05)) #determine closest type I error value to 0.05
    t1e = foo[loc] #store type I error
    stop.bounds = bounds[[loc]] #store corresponding stopping boundary
    
    #compute power at the determined stopping boundary
    pow = mean(sapply(alt.sim,function(x) max(x[,1] > stop.bounds) == 1))
  }
  
  if(bound.type=="asymmetric-OBF" | bound.type=="asymmetric-PK" ){
    u.bound = l.bound = list()
    a = seq(1.8,2.7,by=0.0001)
    foo = rep(0,length(a))
    for(i in 1:length(a)){
      if(bound.type=="asymmetric-OBF"){u.bound[[i]] <- pnorm(sqrt(k/1:k)*a[i])} else{u.bound[[i]] <- rep(pnorm(a[i]),k)}
      l.bound[[i]] <- pnorm(2*a[i]-sqrt(k/1:k)*a[i])
      pass.upper.time <- as.numeric(sapply(null.sim,function(x){if(sum(x[,1] > u.bound[[i]]) > 0){which.max(x[,1] > u.bound[[i]])} else{0}})) # determine the first time the upper boundary was crossed
      pass.lower.time <- as.numeric(sapply(null.sim,function(x){if(sum(x[,1] < l.bound[[i]]) > 0){which.max(x[,1] < l.bound[[i]])} else{0}})) # determine the first time the lower boundary was crossed
      foo[i] <- mean((pass.upper.time > 0) & ifelse(pass.lower.time == 0, 1, ifelse(pass.upper.time < pass.lower.time,1,0))) # determine percentage of times the upper boundary was crossed first
    }
    loc = which.min(abs(foo-0.025)) #determine closest one-sided type I error value to 0.025
    t1e = foo[loc] #store type I error
    stop.bounds = list(l.bound[[loc]],u.bound[[loc]]) #store corresponding stopping boundary
    
    #compute power at the determined stopping boundary
    pass.upper.time <- as.numeric(sapply(alt.sim,function(x){if(sum(x[,1] > stop.bounds[[2]]) > 0){which.max(x[,1] > stop.bounds[[2]])} else{0}})) # determine the first time the upper boundary was crossed
    pass.lower.time <- as.numeric(sapply(alt.sim,function(x){if(sum(x[,1] < stop.bounds[[1]]) > 0){which.max(x[,1] < stop.bounds[[1]])} else{0}})) # determine the first time the lower boundary was crossed
    pow <- mean((pass.upper.time > 0) & ifelse(pass.lower.time == 0, 1, ifelse(pass.upper.time < pass.lower.time,1,0))) # determine percentage of times the upper boundary was crossed first
  }
  
  #return type I error, power, stopping boundaries, null simulations, alternative simulations
  return(c(list(t1e,pow,stop.bounds,null.sim,alt.sim)))
}

###################################################################################################
################# Function to Determine Maximum Sample Size & Stopping Boundaries #################
###################################################################################################

### This function has the following inputs:
###   - nmax = maximum trial sample size
###   - ntrials = number of simulations
###   - null.RR = response.probs vector under the null hypothesis
###   - alt.RR = response.probs vector under the alternative hypothesis
###   - k = number of interim analyses
###   - bound.type = "OBF","Pocock","asymmetric-OBF","asymmetric-PK"; specifies whether to use symmetric OBF-like or Pocock-like stopping boundaries, an 
###                   OBF-like upper boundary with a binding lower futility boundary, or a Pocock-like upper boundary with a binding lower futility boundary
###   - PP.mod = specify which RAR adaptive modification to use ("restrict","exponent","neyman","optimal");
###              "exponent" refers to the modification that uses the tuning parameter

### This function has the following outputs:
###   - nmax = maximum sample size of trial that achieves 90% power
###   - t1e = type I error of design
###   - pow = power of design
###   - stop.bounds = stopping boundaries calibrating the 2-sided type I error rate at 0.05 (for symmetric boundaries)
###                   or 1-sided type I error rate at 0.025 (for asymmetric boundaries)###   - null.sim = trial results for the null scenario simulations
###   - alt.sim = trial results for the alternative scenario simulations

### Purpose of function:
###   - This function first executes the "find.pow" function using the user-specified nmax
###   - if the power of the design is less than 0.89 or greater than 0.90, the "find.pow" function is re-performed
###     using a new estimate of nmax
###   - The above step is repeated until 0.89<power<0.91

find.nmax = function(nmax,ntrials,null.RR,alt.RR,k,bound.type,PP.mod){
  Pow.res = find.pow(nmax=nmax,ntrials=ntrials,null.RR=null.RR,alt.RR=alt.RR,k=k,bound.type = bound.type, PP.mod = PP.mod)
  while(Pow.res[[2]] < 0.89 | Pow.res[[2]] > 0.91){
    R = ((qnorm(0.975)+qnorm(0.90))^2)/((qnorm(0.975)+qnorm(Pow.res[[2]]))^2)
    nmax = ceiling(nmax*R)
    Pow.res = find.pow(nmax=nmax,ntrials=ntrials,null.RR=null.RR,alt.RR=alt.RR,k=k,bound.type = bound.type, PP.mod = PP.mod)
  }  
  return(c(list(nmax,Pow.res[[1]],Pow.res[[2]],Pow.res[[3]],Pow.res[[4]],Pow.res[[5]])))
}

###################################################################################################
###################################### Example ####################################################
###################################################################################################

### Determine nmax for a Bayesian RAR gs trial with a target effect size of 25% and K=10,5,3
### We recommend starting with nmax equal to the maximum sample size from the gsDesign package for the analogous classical frequentist design

library(gsDesign)
fixed.ss = 2*(power.prop.test(n = NULL, p1 = 0.12, p2 = 0.37, sig.level = 0.05, power = 0.90, alternative = "two.sided"))$n
GSdesign = gsDesign(k=10, test.type=2, alpha = 0.025, beta = 0.10, sfu = "Pocock", n.fix = fixed.ss)
max(ceiling(GSdesign$n.I))

set.seed(1994)
RAR.Pocock.10.TE25 = find.nmax(nmax=153,ntrials=10000,null.RR=c(0.12,0.12),alt.RR=c(0.37,0.12),k=10,bound.type = "Pocock")
RAR.Pocock.5.TE25 = find.nmax(nmax=146,ntrials=10000,null.RR=c(0.12,0.12),alt.RR=c(0.37,0.12),k=5,bound.type = "Pocock")
RAR.Pocock.3.TE25 = find.nmax(nmax=139,ntrials=10000,null.RR=c(0.12,0.12),alt.RR=c(0.37,0.12),k=3,bound.type = "Pocock")