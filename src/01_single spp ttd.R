#' ---
#' title: "Single species TTD models"
#' author: "Dominic Henry"
#' date: "06/08/2020"
#' output: 
#'   html_document:
#'     toc: yes
#'     toc_float: yes
#'     toc_collapsed: yes
#'     code_folding: "show"
#'     toc_depth: 3
#'     number_sections: true
#' ---
#' 
## ----setup, include=FALSE-------------------------------------------------------------------------------------------------
knitr::opts_chunk$set(
  echo=TRUE,
  message=FALSE,
  warning=FALSE,
  eval=TRUE,
  fig.width=12,
  fig.height=9
  )

#' 
#' 
## -------------------------------------------------------------------------------------------------------------------------
library(nimble)    
library(AHMbook)
library(unmarked)
library(tidyverse)
library(boot)

#' 
#' # Simulate TTD data
## -------------------------------------------------------------------------------------------------------------------------
set.seed(123) # Needed to reproduce results of simulations 

data <- AHMbook::simOccttd(
  M = 250,            # Number of sites
  mean.psi = 0.4,     # Occupancy probability 
  mean.lambda = 0.3,  # Detection rate
  beta1 = 1,          # Slope of continuous covariate B on logit(psi)
  alpha1 = -1,        # Slope of continuous covariate A on log(lambda)
  Tmax = 10           # Maximum survey/search time
)

## Check all components of simulated data set
str(data)

#'  
#' # TTD model in unmarked
#' 
#' ## Compile data 
#' `NA` values in unmarked are typically used to indicate missing values as a result of sampling issues (i.e., lost data, equipment failure) but `simOccttd` indicates censored observations as `NA`. We therefore indicate a censored `y` in `occuTTD()` by setting `y = Tmax`
## -------------------------------------------------------------------------------------------------------------------------
um_ttd <- data$ttd
um_ttd[which(is.na(um_ttd))] <- data$Tmax

## Build unmarkedFrame
umf <- unmarkedFrameOccuTTD(
  y = um_ttd,           
  surveyLength = data$Tmax, 
  siteCovs = data.frame(
    cbind(covA = data$covA, 
          covB = data$covB)
    )
)

## Check data
head(umf)

#'  
#' ## Fit model 
#' Note that you can specify a dynamic occupancy model using the  `gammaformula` and `epsilonformula` arguments.
## -------------------------------------------------------------------------------------------------------------------------
um_fit <- occuTTD(psiformula = ~covB, detformula = ~covA, data = umf)
um_fit


#' 
#' ## Create summary table 
#' This is used to compare the results of the models to those of the simulated parameters. 
## -------------------------------------------------------------------------------------------------------------------------
summary <- tibble(parameter = c("int.psi", "beta1", "int.lambda", "alpha1"),
                  estimate = coef(um_fit)) %>%
  mutate(estimate_unmarked = c(inv.logit(.$estimate[1]), # Back transform to get estimate on probability scale
                               .$estimate[2],
                               exp(.$estimate[3]),       # Back transform to get estimate on probability scale    
                               .$estimate[4])) %>%
  mutate(truth = c(data$mean.psi, data$beta1, data$mean.lambda, data$alpha1)) %>% 
  select(parameter, truth, estimate_unmarked)

summary


#' 
#' ## Predict detection probabilites
## -------------------------------------------------------------------------------------------------------------------------
lam <- predict(um_fit, type='det', newdata=data.frame(covA=0.5, covB=0.1))$Predicted
survey_time <- 5
pexp(survey_time, lam)
head(getP(um_fit))


#' 
#' # TTD model in NIMBLE
#' 
#' ## Compile data
## -------------------------------------------------------------------------------------------------------------------------
nim_data <- list(ttd = data$ttd,   # TTD data
                 d = data$d,       # Censoring indicator (1 = censored / 0 = not censored)
                 covA = data$covA, # Values of covariate A
                 covB = data$covB, # Values of covariate A
                 nobs = data$M,    # Number of sites
                 Tmax = data$Tmax  # Maximum survey time
                 )

str(nim_data)


#'  
#'  
#' ## Define model 
## -------------------------------------------------------------------------------------------------------------------------
nim_model <- nimbleCode({
  
  # Priors
  int.psi ~ dunif(0, 1)               # Intercept occupancy on prob. scale
  beta1 ~ dnorm(0, 0.001)             # Slope coefficient in logit(occupancy)
  int.lambda ~ dgamma(0.0001, 0.0001) # Poisson rate parameter
  alpha1 ~ dnorm(0, 0.001)            # Slope coefficient in log(rate)
  
  # Likelihood
  for (i in 1:nobs){
   
    # Model for occurrence
    z[i] ~ dbern(psi[i])
    logit(psi[i]) <- logit(int.psi) + beta1 * covB[i]
    
    # Observation model
    # Exponential model for time to detection ignoring censoring
    ttd[i] ~ dexp(lambda[i])
    log(lambda[i]) <- log(int.lambda) + alpha1 * covA[i]
    # Model for censoring due to species absence and ttd>=Tmax
    d[i] ~ dbern(theta[i])
    theta[i] <- z[i] * step(ttd[i] - Tmax) + (1 - z[i])
  }
  # Derived quantities
  n.occ <- sum(z[1:nobs])              # Number of occupied sites among M
}
)


#'  
#' ## Model inits 
## -------------------------------------------------------------------------------------------------------------------------
## Initialize with z = 1 throughout all sites
zst <- rep(1, length(nim_data$ttd))
## All NA's due to censoring rather than non-occurrence
ttdst <-rep(nim_data$Tmax+1, data$M)
ttdst[nim_data$d == 0] <- NA 
## Function run at each iteration 
inits <- function(){list(z = zst, 
                         ttd = ttdst, 
                         int.psi = runif(1), 
                         int.lambda = runif(1))}


#' 
#' ## Monitored parameters and settings 
## -------------------------------------------------------------------------------------------------------------------------
## Parameters to estimate
params <- c("int.psi", "beta1", "int.lambda", "alpha1")

## MCMC settings
ni <- 12000 # Iterations
nt <- 10    # Thinning
nb <- 2000  # Burn-in
nc <- 3     # Chains


#'  
#' ## Fit NIMBLE model 
## ----message=TRUE---------------------------------------------------------------------------------------------------------
nim_fit <- nimbleMCMC(
  code = nim_model,      
  constants = nim_data ,    
  inits = inits,
  monitors = params,
  nburnin = nb,
  niter = ni,
  nchains = nc,
  thin = nt,
  samplesAsCodaMCMC = TRUE,
  summary = TRUE
)


#' 
#' ## Compare results 
## -------------------------------------------------------------------------------------------------------------------------
summary <- summary %>% 
  mutate(estimate_nimble = nim_fit$summary$all.chains[c(4,2,3,1),1])

summary

#' 
#' # Notes on NIMBLE
#' 
#' + NIMBLE [homepage](https://r-nimble.org/)
#' + nimbleEcology [package](https://r-nimble.org/nimbleecology-custom-nimble-distributions-for-ecologists)
#' + Installation [instructions](https://r-nimble.org/html_manual/cha-installing-nimble.html)
#' + [Converting](https://r-nimble.org/nimbleExamples/converting_to_nimble.html) from JAGS, OpenBUGS or WinBUGS
#' + [Parallelization](https://r-nimble.org/nimbleExamples/parallelizing_NIMBLE.html) with NIMBLE
#' + Ecological models from Applied Hierarchical Modelling in Ecology ([AHM book](https://www.elsevier.com/books/applied-hierarchical-modeling-in-ecology-analysis-of-distribution-abundance-and-species-richness-in-r-and-bugs/kery/978-0-12-801378-6)) which have been converted from JAGS to NIMBLE
#' + Ponisio, L.C., et al. 2020. One size does not fit all: Customizing MCMC methods for hierarchical models using NIMBLE. [Ecology & Evolution 10: 2385– 2416](https://onlinelibrary.wiley.com/doi/full/10.1002/ece3.6053)
#' 
#' Efficiency comparison:
#' 
#' ![](nimVsjags.jpg)
#' 
