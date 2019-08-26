# install.packages("dplyr")
library(truncdist)
library(rmutil)
library(rootSolve)
library(dplyr)

source("G:/My Drive/MDrive/dissertation/3/code/sample_size_calculation/first_method/r/single_pair_comparison/functions.R")
#############################
start.time <- Sys.time()
# constant
K = 3
LIST_OF_PIS=c(0.90,0.65,0.65)
LIST_OF_PIS=LIST_OF_PIS[order(LIST_OF_PIS,decreasing = T)]
  
LINKAGE_SAMPLE=1000

COVRAGE=0.9
# COVRAGE_BONFF=1-(1-COVRAGE)/(K-1)
# COVRAGE=COVRAGE_BONFF
POW = 0.8 #when POW is small, we need to decrease sample size lower limit, need to estimate the lower limit or there will not be a opposite sign

SAMPLE_SIZE_LLIMIT=1
SAMPLE_SIZE_ULIMIT=1000

CIL_MIN=0.1# can not be too small
CIL_MAX=2
CIL_STEP=0.01


###################
# check if pis are all the same
# see if the first two treatment are the same-if not the same then do computation, if the same, do computation with the first and last
# need all treatment arm for computation
if(LIST_OF_PIS[1]==LIST_OF_PIS[2]){
  piA=LIST_OF_PIS[1]
  piB=LIST_OF_PIS[K]
  piC=LIST_OF_PIS[setdiff(1:K,c(1,K))]
}else{
  piA=LIST_OF_PIS[1]
  piB=LIST_OF_PIS[2]
  piC=LIST_OF_PIS[setdiff(1:K,c(1,2))]
}

# generate truncated beta1 and beta0
set.seed(1991)
# generate pareto truncated at 1.5 with location 1 and scale 3, mean 1.5
beta1_sample=rtrunc(LINKAGE_SAMPLE, spec="pareto",a=1,b=1/max(c(piA,piB,piC)),1,3)
# beta1_sample=rep(1,LINKAGE_SAMPLE)
beta0_sample=rbeta(LINKAGE_SAMPLE,.5*2,2-.5*2)
# beta0_sample=rep(1,LINKAGE_SAMPLE)
# plot(hist(beta0))
# load functions

# for each pair of beta1 and beta0, solve for n then take average
sample_size_list_pair1=NULL

error_count_pair1=0
warn_count_pair1=0

sim_count=0

error_round_pair1=NULL
warn_round_pair1=NULL
error_mesg_pair1=NULL
warn_mesg_pair1=NULL

error_round_pair1=NULL
warn_round_pair1=NULL
error_mesg_pair1=NULL
warn_mesg_pair1=NULL

for(i in 1:LINKAGE_SAMPLE){
  print(i)
  # i=974 when some ill beta1 and beta0 generated so that we exaust ciL there is no solution, then we throw this draw and report this situation
  beta1=beta1_sample[i]
  beta0=beta0_sample[i]
  
  aA = piA*2
  bA = 2-aA
  # aA = 0.25*2
  # bA = 2-aA
  cA = cFunc(aA, bA, beta1)
  dA = dFunc(aA, bA, beta1)
  eA = eFunc(aA, bA, beta0, K)
  fA = fFunc(aA, bA, beta0, K)
  
  aB = piB*2
  bB = 2-aB
  # aB = 0.25*2
  # bB = 2-aB
  cB = cFunc(aB, bB, beta1)
  dB = dFunc(aB, bB, beta1)
  eB = eFunc(aB, bB, beta0, K)
  fB = fFunc(aB, bB, beta0, K)
  
  aC = piC*2
  bC = 2-aC
  # aC = 0.25*2
  # bC = 2-aC
  cC = cFunc(aC, bC, beta1)
  dC = dFunc(aC, bC, beta1)
  eC = eFunc(aC, bC, beta0, K)
  fC = fFunc(aC, bC, beta0, K)
  
  ciZ=qnorm(1-(1-COVRAGE)/2) # critical value of coverage
  
  for(ciL in seq(CIL_MIN,CIL_MAX,by=CIL_STEP))
  {
    # print(ciL)
    # get sample size solved
    # print(ciL)
    # ciL=0.26
    if(ciL/2==(max(c(piA,piB,piC))-min(c(piA,piB,piC)))){
      next
    } 
    tryCatch({
      # set.seed(i)
      # skip iteration and go to next iteration if ciL/2==(max(c(piA,piB,piC))-min(c(piA,piB,piC)))
      # sample_size_tmp_pair1 <- ceiling(uniroot(sample_size_equation, c(SAMPLE_SIZE_LLIMIT, SAMPLE_SIZE_ULIMIT))$root)
      sample_size_tmp_pair1 <- max(ceiling(uniroot.all(function(n) sample_size_equation(K,
                                                                                        piA, piB, piC,
                                                                                        beta1, beta0, 
                                                                                        aA, bA, cA, dA, eA, fA,
                                                                                        aB, bB, cB, dB, eB, fB,
                                                                                        n), c(SAMPLE_SIZE_LLIMIT, SAMPLE_SIZE_ULIMIT))))
      # All <- uniroot.all(sample_size_equation, c(0, 8))
      # calculate powerpower
      # 1-pnorm((ciL/2-(max(c(piA,piB)-min(piA,piB))))/sqrt(sigmaSqABDiffFunc(sample_size_tmp_pair1)))+pnorm((-ciL/2-(max(c(piA,piB)-min(piA,piB))))/sqrt(sigmaSqABDiffFunc(sample_size_tmp_pair1)))
      pow_pair1=1-pnorm((ciL/2-abs(muABDiffFunc(K,
                                                piA, piB, piC,
                                                beta1, beta0, 
                                                aA, bA, cA, dA, eA, fA,
                                                aB, bB, cB, dB, eB, fB,
                                                sample_size_tmp_pair1)))/sqrt(sigmaSqABDiffFunc(K,
                                                                                                piA, piB, piC,
                                                                                                beta1, beta0, 
                                                                                                aA, bA, cA, dA, eA, fA,
                                                                                                aB, bB, cB, dB, eB, fB,
                                                                                                sample_size_tmp_pair1)))
      if(pow_pair1<POW) break
    },
    error = function(c) {
      error_round_tmp_pair1=cbind(i,beta1,beta0,ciL)
      error_round_pair1<<-rbind(error_round_pair1,error_round_tmp_pair1)
      # next
    },
    warning = function(c) {
      warn_round_tmp_pair1=cbind(i,beta1,beta0,ciL)
      warn_round_pair1<<-rbind(warn_round_pair1,warn_round_tmp_pair1)
      # print(i)
      # print(ciL)
      # next
    },
    finally = {# posterior_sample_burn<<-window(posterior_sample,start=BURNING, end=MCMC_SAMPLE)
      # posterior_sample_cmb<<-do.call(rbind, posterior_sample_burn)
    }
    )
    
  }
  # catch the first time power is lower than target
  tryCatch({
    ciL=ciL-CIL_STEP
    sample_size_tmp_pair1 <- max(ceiling(uniroot.all(function(n) sample_size_equation(K,
                                                                                      piA, piB, piC,
                                                                                      beta1, beta0, 
                                                                                      aA, bA, cA, dA, eA, fA,
                                                                                      aB, bB, cB, dB, eB, fB,
                                                                                      n), c(SAMPLE_SIZE_LLIMIT, SAMPLE_SIZE_ULIMIT))))
    # calculate powerpower
    # 1-pnorm((ciL/2-(max(c(piA,piB)-min(piA,piB))))/sqrt(sigmaSqABDiffFunc(sample_size_tmp_pair1)))+pnorm((-ciL/2-(max(c(piA,piB)-min(piA,piB))))/sqrt(sigmaSqABDiffFunc(sample_size_tmp_pair1)))
    pow_pair1=1-pnorm((ciL/2-abs(muABDiffFunc(K,
                                              piA, piB, piC,
                                              beta1, beta0, 
                                              aA, bA, cA, dA, eA, fA,
                                              aB, bB, cB, dB, eB, fB,
                                              sample_size_tmp_pair1)))/sqrt(sigmaSqABDiffFunc(K,
                                                                                              piA, piB, piC,
                                                                                              beta1, beta0, 
                                                                                              aA, bA, cA, dA, eA, fA,
                                                                                              aB, bB, cB, dB, eB, fB,
                                                                                              sample_size_tmp_pair1)))
  },
  error = function(c) {
    error_mesg_pair1=rbind(error_mesg_pair1,c)
    # next
  },
  warning = function(c) {
    warn_mesg_pair1=rbind(error_mesg_pair1,c)
    # next
  },
  finally = {# posterior_sample_burn<<-window(posterior_sample,start=BURNING, end=MCMC_SAMPLE)
    # posterior_sample_cmb<<-do.call(rbind, posterior_sample_burn)
  }
  )
  
  cvrg=1-(1-round(2*pnorm(ciZ)-1,digits = 3))
  # sample
  
  sample_size_pair1=cbind(ss=sample_size_tmp_pair1,coverage=cvrg,power=pow_pair1,length_ci=ciL,
                          bt1=round(beta1,digits = 4),bt0=round(beta0,digits = 4))
  sample_size_list_pair1=rbind(sample_size_list_pair1,sample_size_pair1)
}

# generate summary sample size for different arms
index_error_pair1=unique(error_round_pair1[,1])
error_count_pair1=length(index_error_pair1)
index_warn_pair1=unique(warn_round_pair1[,1])
warn_count_pair1=length(index_warn_pair1)
index_error_warn_pair1=c(index_error_pair1,index_warn_pair1)
if(is.null(index_error_warn_pair1)){
  sample_size_list_pure_pair1=sample_size_list_pair1
}else{
  sample_size_list_pure_pair1=sample_size_list_pair1[-index_error_warn_pair1,]
}
# head(sample_size_list_pure_pair1)
sample_size_mean_pair1=ceiling(mean(sample_size_list_pure_pair1[,1]))
sample_size_median_pair1=ceiling(median(sample_size_list_pure_pair1[,1]))


# choose sample size

sample_size_mean=sample_size_mean_pair1
sample_size_median=sample_size_median_pair1
nrow(sample_size_list_pure_pair1)
index_error_pair1
index_warn_pair1
# error_mesg_pair1
# warn_mesg_pair1

head(sample_size_list_pure_pair1)

end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken

sample_size_mean
sample_size_median
