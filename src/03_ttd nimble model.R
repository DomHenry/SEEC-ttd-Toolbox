#' ---
#' title: "Multi-species TTD NIMBLE model"
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
  fig.width=10,
  fig.height=7,
  rows.print= 20
  )

#' 
## -------------------------------------------------------------------------------------------------------------------------
library(tidyverse)
library(boot)
library(coda)
library(nimble)


#'  
#' # Import workspace
## -------------------------------------------------------------------------------------------------------------------------
load("data output/TTD_arrays.RData")

#' 
#' # Scale covariates 
## -------------------------------------------------------------------------------------------------------------------------
scaleCen <- function(array) {
  ref <- which(!is.na(array))
  array[ref] <- scale(array[ref])
  return(array)
}

temp_sc <- scaleCen(temp_arr); cloud_sc <- scaleCen(cloud_arr); wind_sc <- scaleCen(wind_arr)
min_sc <- scaleCen(mins_arr); jday_sc <- scaleCen(jday_arr)

site_covs_sc <- site_covs %>% 
  mutate_at(vars(map_ctn:ndvi),scale) %>%
  mutate_at(vars(map_ctn:ndvi),as.numeric)

#' 
#' # Define TTD NIMBLE model
#' 
#' Turn BUGS model code into an object for use in `nimbleModel()`
## -------------------------------------------------------------------------------------------------------------------------
ttd_model <- nimbleCode({
  
  # Prior for species-specific effects in Occ & Detec #
  # ************************************************** #    
  
  for(k in 1:n_spp){	                                 
    logLambda[k] ~ dnorm(mu.logLambda, tau.logLambda) 
    lpsi[k] ~ dnorm(mu.lpsi, tau.lpsi)                     
    beta1[k] ~ dnorm(mu.beta1, tau.beta1)                # ndvi
    beta2[k] ~ dnorm(mu.beta2, tau.beta2)                # map_ctn
    beta3[k] ~ dnorm(mu.beta3, tau.beta3)                # elev   
    beta4[k] ~ dnorm(mu.beta4, tau.beta4)                # tri
  }	
  
  # Hyperpriors for model of occupancy #
  # ********************************** #    
  
  # Dorazio, Gotelli, and Ellison (2011)
  
  mu.psi <- ilogit(mu.lpsi)
  mu.lpsi ~ dt(muT,tauT,kT)    
  
  muT <- 0
  tauT <- pow(sdT, -2)
  sdT <- 1.566
  kT <- 7.763
  
  x.lpsi ~ dt(0,1,1)
  sd.lpsi  <- abs(x.lpsi)  
  tau.lpsi <- pow(sd.lpsi, -2)
  
  mu.beta1 ~ dnorm(0, 0.1)
  x.beta1 ~ dt(0,1,1)
  sd.beta1 <- abs(x.beta1)  
  tau.beta1 <- pow(sd.beta1, -2)
  
  mu.beta2 ~ dnorm(0, 0.1)
  x.beta2 ~ dt(0,1,1)
  sd.beta2 <- abs(x.beta2)  
  tau.beta2 <- pow(sd.beta2, -2)
  
  mu.beta3 ~ dnorm(0, 0.1)
  x.beta3 ~ dt(0,1,1)
  sd.beta3 <- abs(x.beta3)  
  tau.beta3 <- pow(sd.beta3, -2)
  
  mu.beta4 ~ dnorm(0, 0.1)
  x.beta4 ~ dt(0,1,1)
  sd.beta4 <- abs(x.beta4)  
  tau.beta4 <- pow(sd.beta4, -2)
  
  # Hyperpriors for model of detection #
  # ********************************** #
  
  mu.logLambda <- log(lambdaP)
  lambdaP ~ dgamma(0.0001, 0.0001)
  
  x.logLambda ~ dt(0,1,1)
  sd.logLambda <- abs(x.logLambda) 
  tau.logLambda <- pow(sd.logLambda, -2)
  
  # Detection covariate priors #
  # ************************** #
  
  alpha1 ~ dnorm(0, 0.1)       # temp
  alpha2 ~ dnorm(0, 0.1)       # wind
  alpha3 ~ dnorm(0, 0.1)       # cloud
  alpha4 ~ dnorm(0, 0.1)       # mins
  alpha5 ~ dnorm(0, 0.1)       # mins (quad)
  alpha6 ~ dnorm(0, 0.1)       # jday
  
  # Occurrence sub-model #
  # ******************* #
  
  for(k in 1:n_spp){ 
    for (i in 1:n_pen) {
      ## Ecological model of true occurrence [Occ z at site i for spp k]
      z[i,k] ~ dbern(psi[i,k])                           
      logit(psi[i,k]) <- lpsi[k] + beta1[k]*ndvi[i] + beta2[k]*map_ctn[i] + 
                                   beta3[k]*elev[i] + beta4[k]*tri[i] 
      
    }
  }
  
  # Detection sub-model #
  # ******************* #
  for(k in 1:n_spp){
    for (i in 1:n_pen) {
      for (j in 1:transvec[i]){       # Loop over transects
        
        Y[i,j,k] ~ dexp(lambda[i,j,k])
        
        ## Linear model for log(rate)
        log(lambda[i,j,k]) <- logLambda[k] + alpha1*temp[i,j] + alpha2*wind[i,j] +
                              alpha3*cloud[i,j] + alpha4*mins[i,j] + 
                              alpha5*pow(mins[i,j],2) + alpha6*jday[i,j]
        
        ## Accommodation of z=0 and censoring
        d[i,j,k] ~ dbern(theta[i,j,k])
        theta[i,j,k] <- z[i,k] * step(Y[i,j,k] -  tmax) + (1 - z[i,k])
        
      }
    }
  }
  
  # Derived quantities #
  # ****************** #
  
  for(k in 1:n_spp){
    Nocc.fs[k] <- sum(z[1:n_pen,k])	        # Number of occupied sites
    lam.est.sp[k] <- exp(logLambda[k])
    p.sp[k] <- 1-exp(-lam.est.sp[k]*time)   # Mean species-specific detec prob after 700 seconds (time)
  }
  
  for (i in 1:n_pen) {                               
    Nsite[i] <- sum(z[i,1:n_spp])           # Number of occurring species at each site
  } 

})


#' 
#' **Prior choices based on**: Dorazio, R. M., Gotelli, N. J., & Ellison, A. M. (2011). Modern methods of estimating biodiversity from presence-absence surveys. In O. Grillo, & G. Venora (Eds.), Biodiversity loss in a changing planet (pp. 277–302). Rijeka, Croatia: IntechOpen.  
#'   
#' **Also see**: Banner, KM, Irvine, KM, Rodhouse, TJ. The use of Bayesian priors in Ecology: The good, the bad and the not great. Methods Ecol Evol. 2020; 11: 882– 889. https://doi.org/10.1111/2041-210X.13407  
#' 
#' # Monitored parameters 
## -------------------------------------------------------------------------------------------------------------------------
params <- c("mu.lpsi","mu.psi","sd.lpsi","mu.logLambda","sd.logLambda",
            "mu.beta1","sd.beta1", "mu.beta2","sd.beta2",
            "mu.beta3","sd.beta3","mu.beta4","sd.beta4",
            "alpha1","alpha2","alpha3","alpha4","alpha5","alpha6",
            "Nocc.fs",
            "p.sp", "lam.est.sp",
            "lpsi","logLambda",
            "Nsite",
            "beta1", "beta2", "beta3", "beta4") 


#'  
#' # Assign constants 
## -------------------------------------------------------------------------------------------------------------------------
occ_mod_consts <- list(n_spp = dim(Y)[3], 
                       transvec = pen_sum$n_survey,
                       n_pen = dim(Y)[1],
                       tmax = 665,
                       time = 700)


#' 
#' # Compile data 
## -------------------------------------------------------------------------------------------------------------------------
occ_mod_data <- list(Y = Y, 
                     d = d,
                     ndvi = site_covs_sc$ndvi, map_ctn = site_covs_sc$map_ctn, 
                     elev = site_covs_sc$elev, tri = site_covs_sc$tri_med,
                     temp = temp_sc, wind = wind_arr, cloud = cloud_arr, 
                     mins = min_sc, jday = jday_sc)


#'  
#' # Initial values
## -------------------------------------------------------------------------------------------------------------------------
zst <- array(NA, dim = c(dim(Y)[1],dim(Y)[3]))
zst[] <- 1

## TTD inits
ttdst <- Y
I <- which(is.na(Y))
ttdst[I] <- 800 # Make sure this is greater than Tmax (in our case 665)
ttdst[-I] <- NA

transvec <- pen_sum$n_survey
maxpen <- max(pen_sum$n_survey)

## Very important: Put a missing value (NA) for nodes of "ttdst" that you don't want 
## to initialize (i.e array values in which no transect was completed based on the  
## transvec vector)

for(i in 1:dim(ttdst)[3]){     # Species
  for(r in 1:dim(ttdst)[1]){   # Pentad
    if(transvec[r] < maxpen){
      ref <- (transvec[r]+1):maxpen
      ttdst[r,ref,i] <- NA
    }
  }
}


#'  
#' # Init function 
## -------------------------------------------------------------------------------------------------------------------------
occ_mod_inits <- list(z = zst,
                      Y = ttdst,
                      mu.psi = runif(1),
                      alpha1 = runif(1),
                      alpha2 = runif(1),
                      alpha3 = runif(1),
                      alpha4 = runif(1),
                      alpha5 = runif(1),
                      alpha6 = runif(1),
                      mu.beta1 = runif(1),
                      mu.beta2 = runif(1),
                      mu.beta3 = runif(1),
                      mu.beta4 = runif(1),
                      beta1 = rep(runif(1), occ_mod_consts$n_spp),
                      beta2 = rep(runif(1), occ_mod_consts$n_spp),
                      beta3 = rep(runif(1), occ_mod_consts$n_spp),
                      beta4 = rep(runif(1), occ_mod_consts$n_spp),
                      logLambda = rep(-6.5, occ_mod_consts$n_spp),
                      p.sp = rep(0.2, occ_mod_consts$n_spp), 
                      lam.est.sp = rep(0.0006, occ_mod_consts$n_spp),
                      mu.logLambda = -6.5,
                      sd.logLambda = 0.001,
                      lambdaP = 0.0001,
                      tau.logLambda = 0.0002)


#'  
#' # Build the NIMBLE model 
#' Processes BUGS model code and optional constants, data, and initial values
## -------------------------------------------------------------------------------------------------------------------------
s0 <- Sys.time()
occ_model <- nimbleModel(code = ttd_model, 
                         name = "Time-to-detection MSOM", 
                         constants = occ_mod_consts,
                         data = occ_mod_data, 
                         inits = occ_mod_inits,
                         calculate = TRUE) 
e0 <- Sys.time()
e0-s0

occ_model$initializeInfo()
occ_model$logProb_alpha1
occ_model$logLambda
occ_model$lam.est.sp
occ_model$mu.logLambda 
occ_model$sd.logLambda # Why is this NA?
occ_model$beta1

#'  
#' # Build an MCMC object for this model
#' Create an MCMC function from a NIMBLE model
## -------------------------------------------------------------------------------------------------------------------------
s1 <- Sys.time()
occ_model_MCMC <- buildMCMC(occ_model, monitors = params)
e1 <- Sys.time()
e1-s1


#'  
#' # Compile the model 
#' Generate and compile C++
## -------------------------------------------------------------------------------------------------------------------------
s2 <- Sys.time()
comp_occ_mod <- compileNimble(occ_model)
e2 <- Sys.time()
e2-s2

s3 <- Sys.time()
compile_occ_MCMC <- compileNimble(occ_model_MCMC, project = comp_occ_mod)
e3 <- Sys.time()
e3-s3

#' 
#' # MCMC settings
## -------------------------------------------------------------------------------------------------------------------------
## MCMC settings - TEST
ni <- 300 # Iterations
nt <- 5    # Thinning
nb <- 100  # Burn-in
nc <- 3     # Chains

## MCMC settings - FOR REAL REAL
# ni <- 60000  # Iterations
# nt <- 10     # Thinning
# nb <- 15000  # Burn-in
# nc <- 3      # Chains

## Number of posterior samples
floor((ni-nb)/nt)


#'  
#' # Run MCMC 
#' Run chains of an MCMC algorithm and return samples
## ----eval=FALSE-----------------------------------------------------------------------------------------------------------
## # 8 hours to run 60K with 15K burnin
## s4 <- Sys.time()
## runMCMC_samples <- runMCMC(compile_occ_MCMC,
##                            nburnin = nb,
##                            niter = ni,
##                            nchains = nc,
##                            thin = nt,
##                            summary = TRUE,
##                            samples = TRUE,
##                            samplesAsCodaMCMC = TRUE,
##                            progressBar = TRUE,
##                            WAIC = FALSE)
## e4 <- Sys.time()
## e4-s4

#' 
#' # Timings 
## ----eval=FALSE-----------------------------------------------------------------------------------------------------------
## e4-s4 + e3-s3 + e2-s2 + e1-s1 + e0-s0

#'  
#' # Write outputs to file 
## ----eval=FALSE-----------------------------------------------------------------------------------------------------------
## head(runMCMC_samples$summary$all.chains)
## 
## runMCMC_samples$summary$all.chains %>%
##   as_tibble %>%
##   mutate(param = rownames(runMCMC_samples$summary$all.chains)) %>%
##   select(param, everything()) %>%
##   write_csv("data output/ttd_nimble_mcmc_out_small.csv") # Used for testing
##   # write_csv("data output/ttd_nimble_mcmc_out.csv")
## 
## 
## saveRDS(runMCMC_samples, "data output/nimble_mcmc_samples_small.rds") # Used for testing
## # saveRDS(runMCMC_samples, "data output/nimble_mcmc_samples.rds")
## 

#' 
