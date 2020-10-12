# This script determines the maximum sample size (nmax) for a clinical trial implementing a 
# Bayesian response adaptive randomization (RAR) group-sequential (gs) design using Pocock-like 
# (aggressive early stopping) and OBF-like (conservative early stopping) stopping boundaries for
# a specific number of interim analyses K
#     1) For Pocock: we will use the same posterior probability (PP) stopping boundary throughout 
#        the trial that controls type I error at 0.05
#     2) For OBF: we will use c_k = pnorm(sqrt(K/k)*a), where c_k is the value of the kth stopping 
#        boundary and a is an appropriate constant that leads to a type I error rate of 0.05; 
#        PP boundaries start high and get progressively smaller
# nmax is further determined to obtain a power of 90%
# The function "find.nmax" outputs the following: nmax, type I error, power, stopping boundaries,
# null scenario simulations, and alternative scenario simulations
# PP = Pr(Experimental Treatment > Control|data); two-sided symmetric stopping boundaries
# Code uses a logistic regression probability model and mass-weighted urn randomization scheme
# The PP will be stabilized before using as a randomization probability via 1 of the 2 following modifications:
#     1) Restrict: 1-max.ar < rand.prob < max.ar
#     2) Exponent: rand.prob = (PP^c)/(PP^c + (1-PP)^c), with c = 0.5*(k-1)/(K-1); k = current group to be enrolled
#                               and K = the total number of groups (Per Thall and Wathen paper)

###################################################################################################
####################### Function to Employ Logistic Regression Probability Model ##################
###################################################################################################

library(BayesLogit)
n = c(8,7); y = c(1,0); X = rbind(c(1,0.5),c(1,-0.5))#establish design matrix for tlr model

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
####################### Function to Simulate a Bayesian RAR Trial #################################
###################################################################################################

# This function returns posterior probabilities, number assigned to each arm, number of survivors on each arm
simulate.trial <- function(response.probs, ns, max.ar=0.75, rand.type=c("Coin","Urn","Block"), max.deviation=3, model=c("tlr","ibb"),
                  PP.mod = c("restrict","exponent"),pi.star=0.5, pess=2, df=c(7,7), m0=c(log(0.12/0.88),0), s0=c(2.5,2.5), nsamps=2000, warmup=100){
  
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
      #z = sample(c(rep(1,nE.new),rep(0,nC.new))) #Group Randomization Schedule
      yE.new = rbinom(1,nE.new,response.probs[1])
      yC.new = rbinom(1,nC.new,response.probs[2])
    }
    #Update dataset
    n = n+c(nE.new,nC.new)
    y = y+c(yE.new,yC.new)
    
    #Update posterior probability that E is superior to C and randomization probability for next group
    if(model=="ibb") post.prob = get.post.prob(n=n,y=y,pi.star=pi.star,pess=pess)
    if(model=="tlr") post.prob = rowMeans(tlr.sampler(n=n,y=y,X=rbind(c(1,0.5),c(1,-0.5)),df=df,m0=m0,s0=s0,nsamps=nsamps,warmup=warmup)>0)[2]
    
    #Determine updated rand prob by stabilizing the current post prob via user-specified modification
    C = 0.5*(group)/(length(ns)-1)
    if(PP.mod == "restict") rand.prob = min(max.ar,max(1-max.ar,post.prob))
    if(PP.mod == "exponent") rand.prob = (post.prob^C)/(post.prob^C + (1-post.prob)^C)
    
    #Store Relevant Statistics
    stats[group,] = c(post.prob,n[1],y[1],n[2],y[2])  
  }
  
  return(stats)
}

###################################################################################################
################################# Function to Compute Power #######################################
###################################################################################################

# This function computes power using a PP stopping boundary calibrated to control type I error at 0.05
# alternative & null simulations performed using a maximum sample size and interim look number specified by the user
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
  }
  
  if(bound.type == "OBF"){
    bounds = list()
    options = seq(1.8,2.3,by=0.0001)
    foo = rep(0,length(options))
    for(i in 1:length(options)){
      bounds[[i]] = pnorm(sqrt(k/1:k)*options[i])
      foo[i] = mean(sapply(null.sim,function(x) max(x[,1] > bounds[[i]] | x[,1] < 1-bounds[[i]]) == 1))
    }
    loc = which.min(abs(foo-0.05)) #determine closest type I error value to 0.05
    t1e = foo[loc] #store type I error
    stop.bounds = bounds[[loc]] #store corresponding stopping boundary
  }
  
  #compute power at the determined stopping boundary
  pow = mean(sapply(alt.sim,function(x) max(x[,1] > stop.bounds) == 1))
  
  #return type I error, power, stopping boundaries, null simulations, alternative simulations
  return(c(list(t1e,pow,stop.bounds,null.sim,alt.sim)))
}

###################################################################################################
########################## Function to Determine Maximum Sample Size ##############################
###################################################################################################

# This function determines nmax so that power and type I error are 0.90 and 0.05, respectively
find.nmax = function(nmax,ntrials,null.RR,alt.RR,k,bound.type,PP.mod){
  Pow.res = find.pow(nmax=nmax,ntrials=ntrials,null.RR=null.RR,alt.RR=alt.RR,k=k,bound.type = bound.type, PP.mod = PP.mod)
  while(Pow.res[[2]] < 0.89 | Pow.res[[2]] > 0.91){
    R = ((qnorm(0.975)+qnorm(0.90))^2)/((qnorm(0.975)+qnorm(Pow.res[[2]]))^2)
    nmax = ceiling(nmax*R)
    Pow.res = find.pow(nmax=nmax,ntrials=ntrials,null.RR=null.RR,alt.RR=alt.RR,k=k,bound.type = bound.type, PP.mod = PP.mod)
  }  
  return(c(list(nmax,Pow.res[[1]],Pow.res[[2]],Pow.res[[3]],Pow.res[[4]],Pow.res[[5]])))
}

