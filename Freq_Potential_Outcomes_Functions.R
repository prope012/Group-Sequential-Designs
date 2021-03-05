
###################################################################################################
######### Function to Simulate a Frequentist GS Trial Using Potential Outcomes Framework ##########
###################################################################################################

### The function below has the following inputs:
###   - potential.outcomes = a dataframe with 2 columns containing the potential outcomes on the 
###     treatment (column 1) and control (column 2) for each individual in the potential study sample 
###     (generated using the function "Run.Trial.Get.Stats.Freq" below)
###   - nmax = the maximum sample size for the trial (determined using the gsDesign package)
###   - k = the number of interim analyses for the trial
###   - rand.prob = the randomization probability for the treatment arm

### The function below has the following outputs:
###   - a kx4 matrix containing nE,nC,yE,yC at each interim analysis where nE and nC are the number of 
###     trial participants on the treatment and control, respectively, and yE and yC are the number of 
###     responders on the treatment and control, respectively.

### Purpose of function:
###   - This function simulates a 2-arm gs clinical trial with no early stopping (i.e. stopping 
###     boundaries are not yet applied) using the mass-weighted urn randomization scheme and 1:1 allocation
###   - Rather than generating outcomes for each participant within the trial, the outcomes for each 
###     participant are set equal to his or her potential outcome for the arm to which he or she was assigned

simulate.freq.trial.PO <- function(potential.outcomes, nmax, k, rand.prob=0.5){
  
  #generate group sizes
  ns = round(seq(nmax/k,nmax,by=nmax/k)) 
  
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
    start.row = c(1,ns)[group]
    end.row = c(1,ns)[group+1]
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
###   - potential.outcomes = a dataframe containing the potential outcomes for each potential participant on the treatment and control(generated using the "Run.Trial.Get.Stats.Freq" function below)
###   - bounds = the efficacy stopping boundaries for the trial (found using the "Run.Trial.Get.Stats.Freq" function below)
###   - ns = the sample size at each interim analysis (generally equal to round(seq(nmax/k,nmax,by=nmax/k))) 
###   - k = the number of interim analyses
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
###   - This function computes the desired 2-sample binomial proportion test at each interim analysis using
###     the trial data outputted by the "simulate.freq.trial.PO" function
###   - This function then determines if the Z-score from the 2-sample binomial proportion test
###     exceeded an efficacy or harm stopping boundary
###   - Within-trial statistics (nE, nC, yE, yC, non-responders within trial) are reported at the time the trial was stopped
###   - The treatment recommended by the trial is applied to the remaining potential study sample 

library(gsDesign)
library(exact2x2)

prop_test_get_stats_fn = function(trial, interim.test=c("fisher","chisquare","uncondExact"),potential.outcomes,
                         bounds, ns, k, cohort.size, hypothesis=c("alt","null")){
  
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
    
    #convert z-score upper stopping boundaries to two-sided p-values
    p.bounds = (1-pnorm(bounds))*2 
    
    #compute desired stats from trial
    stop.at = which.max(trial[,"pvalue"] <=  c(p.bounds[-length(p.bounds)],1))
    nE = trial[stop.at,"nE"]; nC = trial[stop.at,"nC"]; yE = trial[stop.at,"yE"]; yC = trial[stop.at,"yC"]; 
    n = nE+nC; deaths.trial = n-(yE+yC); maxN = max(ns)
    
    #determine whether trial stopped due to efficacy, harm, or neither
    efficacy = harm = 0 #null accepted
    if(sum(trial[,"pvalue"] <= p.bounds)>0 & trial[stop.at,"tsign"]==1) {efficacy=1; harm=0} #treatment efficacious
    if(sum(trial[,"pvalue"] <= p.bounds)>0 & trial[stop.at,"tsign"]==-1)  {efficacy=0; harm=1} #treatment harmful
      
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
###   - k = the number of interim analyses
###   - bound.type = "OF" or "Pocock"; specifies whether to use an OBF or Pocock stopping boundary
###   - interim.test = the 2-sample proportion test to be computed at each interim analysis (Fisher's exact test, Chi-square test, or the unconditional exact test)
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
###   - This function then runs the "simulate.freq.trial.PO" and "prop_test_get_stats_fn" functions to get
###     all desired sim stats

Run.Trial.Get.Stats.Freq = function(response.probs, cohort.size,k,bound.type=c("OF,Pocock"),
                                    interim.test=c("fisher","chisquare","uncondExact"),hypothesis=c("alt","null")){
  
  #Generate Potential Outcomes
  potential.outcomes =  data.frame(potential_Ye = rbinom(cohort.size, 1, prob = response.probs[1]),
                                   potential_Yc = rbinom(cohort.size, 1, prob = response.probs[2]))
  if(hypothesis=="null") potential.outcomes[,1] = potential.outcomes[,2]
  
  #determine group size, max N, and stopping boundaries needed for 90% power, one-sided alpha = 0.025
  fixed.ss = 2*(power.prop.test(n = NULL, p1 = response.probs[1], p2 = response.probs[2], sig.level = 0.05, power = 0.90, alternative = "two.sided"))$n
  GSdesign = gsDesign(k=k, test.type=2, alpha = 0.025, beta = 0.10, sfu = bound.type, n.fix = fixed.ss)
  ns = ceiling(GSdesign$n.I)
  bounds = GSdesign$upper$bound
  
  #Simulate Trial
  trial = simulate.freq.trial.PO(potential.outcomes = potential.outcomes, nmax = max(ns), k=k)

  #Compute Stats using user-specified 2-sample binomial proportion test
  stats = prop_test_get_stats_fn(trial=trial, interim.test = interim.test, potential.outcomes = potential.outcomes,
                                 bounds = bounds, ns=ns, k=k, cohort.size=cohort.size, hypothesis=hypothesis)
  
  names(stats) =  c("Trial N","Treatment N","Control N","Treatment Survivors",
                    "Control Survivors","Trial Deaths","Cohort Deaths","Efficacy","Harm")
  return(stats)
}

###################################################################################################
############################## Example Using the Above Functions ##################################
###################################################################################################

### simulate 100 trials and get stats - ALT
Ex.Freq.alt = data.frame(do.call(rbind, lapply(1:100, function(x){
  set.seed(x)
  Run.Trial.Get.Stats.Freq(response.probs = c(0.37,0.12),cohort.size=447,k=10,bound.type = "OF",
                           interim.test = "uncondExact",hypothesis = "alt")})))

### simulate 100 trials and get stats - NULL
Ex.Freq.null = data.frame(do.call(rbind, lapply(1:100, function(x){
  set.seed(x)
  Run.Trial.Get.Stats.Freq(response.probs = c(0.37,0.12),cohort.size=447,k=10,bound.type = "OF",
                           interim.test = "uncondExact",hypothesis = "null")})))
