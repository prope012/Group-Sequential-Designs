
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
X = rbind(c(1,0.5),c(1,-0.5))#establish design matrix for tlr model

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
########### Function to Simulate a Bayesian RAR Trial Using Potential Outcomes Framework ##########
###################################################################################################

### The function below has the following inputs:
###   - potential.outcomes = a dataframe with 2 columns containing the potential outcomes on the 
###     treatment (column 1) and control (column 2) for each individual in the potential study sample 
###     (generated using the function "Run.Trial.Get.Stats" below)
###   - ns = the sample size at each interim analysis (generally equal to round(seq(nmax/k,nmax,by=nmax/k))) 
###   - max.ar = the maximum randomization probability to be used by the "Restrict" RAR design; default is 0.75
###   - rand.type = "Coin" for weighted coin flipping, "Urn" for mass-weighted urn design
###   - max.deviation = alpha parameter in the mass-weighted urn design randomization scheme
###   - model = "ibb" for independent beta-binomial probability model vs. "tlr" for logistic regression probabiltiy model
###   - PP.mod = specify which RAR adaptive modification to use ("restrict","exponent","neyman","optimal");
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
###   - Rather than generating outcomes for each participant within the trial, the outcomes for each 
###     participant are set equal to his or her potential outcome for the arm to which he or she was assigned

simulate.RAR.trial.PO <- function(potential.outcomes, ns, max.ar=0.75, rand.type=c("Coin","Urn","Block"), max.deviation=3, model=c("tlr","ibb"),
                           PP.mod = c("restrict","exponent","neyman","optimal"),pi.star=0.5, pess=2, df=c(7,7), m0=c(log(0.12/0.88),0), s0=c(2.5,2.5), nsamps=2000, warmup=100){
  
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
    
    #Grab potential outcomes for current group
    start.row = c(1,ns)[group]
    end.row = c(1,ns)[group+1]
    p.outcomes = potential.outcomes[start.row:end.row,]
    
    #Generate assignments for current interval
    if(rand.type == "Coin"){
      assignments = rbinom(n.new,1,rand.prob) #flip coin to determine assignments for current group
      nE.new = sum(assignments)
      nC.new = n.new - nE.new
    }
    if(rand.type == "Urn"){
      assignments = NULL; rand.prob.temp = rand.prob
      for(i in 1:n.new){
        assignments[i] = rbinom(1,1,min(max(rand.prob.temp,0),1))
        rand.prob.temp = rand.prob.temp+(rand.prob-1*(assignments[i]==1))/max.deviation
      }
      nE.new = sum(assignments)
      nC.new = n.new - nE.new
    }
    if(rand.type == "Block"){
      u = rbinom(1,1,n.new*rand.prob - floor(n.new*rand.prob))
      nE.new = u*ceiling(n.new*rand.prob)+(1-u)*floor(n.new*rand.prob)
      nC.new = n.new-nE.new
      assignments = sample(c(rep(1,nE.new),rep(0,nC.new))) #random permutation to determine assignments for current group
    }
    
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
      rand.prob = mean(sqrt(pE*qE)/(sqrt(pE*qE)+sqrt(pC*qC))) #posterior mean for Neyman allocation
    }
    if(PP.mod == "optimal"){
      pE = expit(gibbs.draws[1,] + 0.5*gibbs.draws[2,])
      pC = expit(gibbs.draws[1,] - 0.5*gibbs.draws[2,])
      rand.prob = mean(sqrt(pE)/(sqrt(pE)+sqrt(pC))) #posterior mean for optimal allocation
    }
    
    #Store Relevant Statistics
    stats[group,] = c(post.prob,n[1],y[1],n[2],y[2])  
  }
  
  return(stats)
}


###################################################################################################
############ Function to Obtain Desired Statistics for a Bayesian RAR Trial #######################
############################ Using Potential Outcomes Framework ###################################
###################################################################################################

### The function below has the following inputs:
###   - w = the trial stats from the "simulate.RAR.trial.PO" function
###   - potential.outcomes = a dataframe with 2 columns containing the potential outcomes on the 
###     treatment (column 1) and control (column 2) for each individual in the potential study sample 
###     (generated using the function "Run.Trial.Get.Stats" below)
###   - bounds = the efficacy stopping boundaries found using the "RAR_findnmax_functions.R" script
###   - cohort.size = the size of the potential study sample for the target effect size
###   - hypothesis = "alt" to run simulations under the target alternative vs. "null" to run simulations under the null scenario

### The function below has the following outputs:
###   - n = number of trial participants after efficacy & harm stopping boundaries have been applied
###   - nE and nC = number of participants respectively assigned to the treatment and control arms 
###     after efficacy & harm stopping boundaries have been applied 
###   - yE and yC = number of trial responders respectively on the treatment and control arms after 
###     efficacy & harm stopping boundaries have been applied
###   - deaths.trial = total number of non-responders within the trial after efficacy & harm stopping 
###     boundaries have been applied
###   - deaths.cohort = total number of non-responders within the potential study sample assuming 
###     the individuals in the potential study sample but not in the trial itself receive the treatment 
###     recommended by the trial
###   - efficacy = a binary indicator for whether the trial stopped early due to treatment efficacy
###   - harm = a binary indicator for whether the trial stopped early due to treatment harm

### Purpose of function:
###   - This function applies the appropriate stopping boundaries to the trial data from the "simulate.RAR.trial.PO" 
###     function and determines whether a trial was stopped early for efficacy or harm
###   - Within-trial statistics (nE, nC, yE, yC, non-responders within trial) are reported at the time the trial was stopped
###   - The treatment recommended by the trial is applied to the remaining potential study sample

RAR_ComputeStat = function(w, potential.outcomes, bounds, cohort.size, hypothesis=c("null","alt")){
    
    #determine when trial is stopped
    stop.at = which.max((w[,1] > c(bounds[-length(bounds)],0)) | (w[,1] < (1-bounds)))
    
    #obtain desired statistics when trial stopped
    nE = w[stop.at,2]; nC = w[stop.at,4]; yE = w[stop.at,3]; yC = w[stop.at,5]; n = nE+nC; deaths.trial = n-(yE+yC)
    
    #determine whether trial stopped due to efficacy, harm, or neither
    efficacy = harm = 0 #null accepted
    if(sum(w[,1] > bounds)>0) {efficacy=1; harm=0} #treatment efficacious
    if(sum(w[,1] < 1-bounds)>0) {efficacy=0; harm=1} #treatment harmful
  
    #determine number of deaths in the potential cohort, taking into account the appropriate hypothesis
    if(hypothesis=="alt" & n!=cohort.size){
      if(efficacy==1){deaths.cohort = sum(potential.outcomes[(n+1):cohort.size,"potential_Ye"]==0) + deaths.trial}
      if(efficacy==0){deaths.cohort = sum(potential.outcomes[(n+1):cohort.size,"potential_Yc"]==0) + deaths.trial}
    }
    #determine number of deaths in the potential cohort, taking into account the appropriate hypothesis
    if(hypothesis=="null" & n!=cohort.size){
      deaths.cohort = sum(potential.outcomes[(n+1):cohort.size,"potential_Yc"]==0) + deaths.trial
    }
    if(n==cohort.size){deaths.cohort=deaths.trial}
    return(c(n, nE, nC, yE, yC, deaths.trial, deaths.cohort, efficacy, harm))
}

### Note: Power is equal to the proportion of trials where efficacy = 1 under alternative
### Note: TIE is equal to the proportion of trials where efficacy = 1 or harm = 1 under null


###################################################################################################
############ Function to Run a Trial and Get Stats Under Potential Outcomes Framework #############
###################################################################################################

### The function below has the following inputs:
###   - response.probs = c(treatment arm response rate, control arm response rate)
###   - cohort.size = the size of the potential study sample for the target effect size
###   - nmax = maximum trial sample size
###   - k = the number of interim analyses
###   - PP.mod = specify which RAR adaptive modification to use ("restrict","exponent","neyman","optimal");
###              "exponent" refers to the modification that uses the tuning parameter###   - bounds = the efficacy stopping boundaries found using the "RAR_findnmax_functions.R" script
###   - hypothesis = "alt" to run simulations under the target alternative vs. "null" to run simulations under the null scenario

### This function has the following outputs:
###   - n = number of trial participants after efficacy & harm stopping boundaries have been applied
###   - nE and nC = number of participants respectively assigned to the treatment and control arms 
###     after efficacy & harm stopping boundaries have been applied 
###   - yE and yC = number of trial responders respectively on the treatment and control arms after 
###     efficacy & harm stopping boundaries have been applied
###   - deaths.trial = total number of non-responders within the trial after efficacy & harm stopping 
###     boundaries have been applied
###   - deaths.cohort = total number of non-responders within the potential study sample assuming 
###     the individuals in the potential study sample but not in the trial itself receive the treatment 
###     recommended by the trial
###   - efficacy = a binary indicator for whether the trial stopped early due to treatment efficacy
###   - harm = a binary indicator for whether the trial stopped early due to treatment harm

### Purpose of Function:
###   - This function generates the potential outcomes for each individual in the potential study 
###     sample using the response.probs vector & cohort.size
###   - This function then runs the "simulate.RAR.trial.PO" and "RAR_ComputeStat" functions to get
###     all desired sim stats

Run.Trial.Get.Stats = function(response.probs, cohort.size, nmax, k, PP.mod, bounds, 
                               hypothesis = c("null","alt")){

  #Get Group Size
  ns = round(seq(nmax/k,nmax,by=nmax/k))
  
  #Generate Potential Outcomes
  potential.outcomes =  data.frame(potential_Ye = rbinom(cohort.size, 1, prob = response.probs[1]),
                                   potential_Yc = rbinom(cohort.size, 1, prob = response.probs[2]))
  if(hypothesis=="null") potential.outcomes[,1] = potential.outcomes[,2]
  
  #Simulate Trial
  trial = simulate.RAR.trial.PO(potential.outcomes = potential.outcomes, ns = ns, rand.type="Urn", 
                        max.deviation=3, model="tlr",PP.mod = PP.mod)
  
  #Compute Stats
  stats = RAR_ComputeStat(w = trial, potential.outcomes = potential.outcomes, bounds = bounds, 
                          cohort.size = cohort.size, hypothesis = hypothesis)
  
  names(stats) =  c("Trial N","Treatment N","Control N","Treatment Survivors",
                       "Control Survivors","Trial Deaths","Cohort Deaths","Efficacy","Harm")
  return(stats)
}

###################################################################################################
############################## Example Using the Above Functions ##################################
###################################################################################################

### Determine nmax and boundaries needed to achieve 90% power, TIE = 0.05 under desired scenario
Ex.RAR.nmax = find.nmax(nmax=153,ntrials=100,null.RR=c(0.12,0.12),alt.RR=c(0.37,0.12),k=10,bound.type = "Pocock",PP.mod="restrict")

### simulate 100 trials and get stats - ALT
Ex.RAR.alt = data.frame(do.call(rbind, lapply(1:100, function(x){
             set.seed(x)
             Run.Trial.Get.Stats(response.probs = c(0.37,0.12), cohort.size = 200, nmax = Ex.RAR.nmax[[1]],
                                 k=10, PP.mod="restrict", bounds = Ex.RAR.nmax[[4]], hypothesis = "alt")})))

### simulate 100 trials and get stats - NULL
Ex.RAR.null = data.frame(do.call(rbind, lapply(1:100, function(x){
              set.seed(x)
              Run.Trial.Get.Stats(response.probs = c(0.37,0.12), cohort.size = 200, nmax = Ex.RAR.nmax[[1]],
                      k=10, PP.mod="restrict", bounds = Ex.RAR.nmax[[4]], hypothesis = "null")})))
