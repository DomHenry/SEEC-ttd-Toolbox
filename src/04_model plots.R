#' ---
#' title: "Plotting TTD model outputs"
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
  fig.height=7
  )

#' 
#' 
## -------------------------------------------------------------------------------------------------------------------------

library(tidyverse)
library(gridExtra)
library(ggpubr)
library(boot)
library(janitor)
library(MCMCvis)

#'  
#' # Load data
## -------------------------------------------------------------------------------------------------------------------------
## Load workspace with array and input data as well as the MCMC object
load("data output/TTD_arrays.RData")
ttd_modfit <- readRDS("data output/nimble_mcmc_samples.rds")

#'  
#' # TTD data exploration
#' 
#' ## TTD frequency
## -------------------------------------------------------------------------------------------------------------------------
## Frequency plot of observed time-to-detection observations
ttd %>% 
  ggplot(aes(ttd))+
  geom_histogram(col = "black",alpha = 0.3, bins = 20) + 
  xlab("\nTime-to-detection (seconds)")+
  ylab("Frequency\n")+
  theme_classic()+
  theme(axis.text = element_text(size = 16, colour = "black"),
        axis.title = element_text(size = 18))

#'  
#' ## TTD species mean
## -------------------------------------------------------------------------------------------------------------------------
## Create a summary table with mean, sd, se, and number of counts
summary_ttd <- ttd %>% 
  group_by(species) %>% 
  summarise(mean_ttd = mean(ttd, na.rm = TRUE), 
            sd_ttd = sd(ttd),
            n_detect = length(ttd), 
            se_ttd = sd(ttd,na.rm = TRUE)/sqrt(length(ttd))) %>%
  arrange(mean_ttd)

summary_ttd

## In order to create plots that are more readable when you have a high 
## number of species you can split the data set and use these vectors
## to create separate plots. You can also plot all species together
## and adjust the axis label size accordingly.
odds <- seq(1,n_spp,2) 
even <- seq(2,n_spp,2)

summary_ttd %>% 
  filter(row_number() %in% even) %>% 
  ggplot(aes(x = forcats::fct_reorder(species, mean_ttd, .desc = TRUE), 
             y = mean_ttd)) +
  geom_bar(stat = "identity", width = 0.8) +
  coord_flip() +
  ylab("Mean time-to-detection (seconds)") +
  xlab("") +
  geom_errorbar(aes(ymin = mean_ttd, ymax = mean_ttd + sd_ttd),
                width = 0.2, position = position_dodge(0.9),
                col = "gray35") +
  scale_y_continuous(limits = c(0, 600)) +
  theme_classic() +
  theme(axis.text.x = element_text(size = 10, color = "black"),
        axis.text.y = element_text(size = 10))

#' 
#' ## TTD events
## -------------------------------------------------------------------------------------------------------------------------
## Plot the number of detections events for each species
summary_ttd %>% 
  filter(row_number() %in% odds) %>% 
  ggplot(aes(x= forcats::fct_reorder(species,n_detect,.desc = TRUE), 
             y= n_detect)) +
  geom_bar(stat="identity",width = 0.8)+
  coord_flip()+
  xlab("")+
  ylab("Detection events")+
  scale_y_continuous(limits = c(0, 250)) +
  theme_classic()+
  theme(axis.text.x = element_text(size = 10, color = "black"),
        axis.text.y = element_text(size = 10))

#' 
#' # MCMC diagnostics 
#' Check convergence
## -------------------------------------------------------------------------------------------------------------------------
MCMCvis::MCMCtrace(ttd_modfit$samples, 
          params = c("mu.psi","alpha1","mu.beta1","beta2\\[1\\]"),
          pdf = FALSE, 
          ind = TRUE,
          n.eff = TRUE,
          Rhat = TRUE,
          ISB = FALSE)


#' 
#' Write traceplots for all nodes to PDF 
## ----eval=FALSE-----------------------------------------------------------------------------------------------------------
## MCMCtrace(ttd_modfit$samples,
##           pdf = TRUE,
##           filename = "data output/nimble_traceplots.pdf",
##           ind = TRUE,
##           n.eff = TRUE,
##           Rhat = TRUE)
## 

#' 
#' 
#' # Create MCMC plotting data frame
## -------------------------------------------------------------------------------------------------------------------------
plot_data <- ttd_modfit$summary$all.chains %>% 
  as_tibble() %>% 
  clean_names() %>% 
  mutate(param = row.names(ttd_modfit$summary$all.chains)) %>% 
  select(param, mean, x95_percent_ci_low, x95_percent_ci_upp) %>% 
  rename(lower = x95_percent_ci_low, upper = x95_percent_ci_upp)

plot_data

#'  
#' # Detection
#' 
#' ## Detection coefficients estimates
## -------------------------------------------------------------------------------------------------------------------------
## Extract alpha parameters
alphas <- plot_data %>% 
  filter(str_detect(param,"alpha")) %>% 
  mutate(param = c("Temperature","Wind","Cloud","Time of day","Time of day (quadratic)","Julian Day"))

alphas

ggplot(alphas, aes(x = fct_reorder(param,mean), y = mean, ymin = lower, ymax = upper)) +
  geom_pointrange() + 
  scale_y_continuous(limits = c(-0.15,0.15))+
  coord_flip() + 
  geom_hline(yintercept = 0, linetype = "dotted") +
  ylab("Estimate (+- 95% CI)") + xlab("")+ 
  theme_classic()+
  theme(axis.title = element_text(size = 14),
        axis.text.x = element_text(size = 14, color = "black"),
        axis.text.y = element_text(size = 14, color = "black"))


#'  
#' ## Detection covariate predictions
## -------------------------------------------------------------------------------------------------------------------------
## Covariate arrays (from bird data prep)
covs <- list(temp_arr, wind_arr, cloud_arr, mins_arr, jday_arr)
## Range of values for prediction
ranges <- list(c(0, 40, 0.1), c(0, 6, 0.01), c(0, 100, 1), c(0, 500, 1), c(-1, 1, 0.01))
## Names of coefficients as per NIMBLE model definition
coeffs1 <- list("alpha1", "alpha2", "alpha3", "alpha4", "alpha6")
## Names of quadratic coefficients (if any)
coeffs2 <- list(NA, NA, NA, "alpha5", NA)
## Lamda rate intercept
intname <- "mu.logLambda"
## Axis labels
xlabs <- list("Temperature", "Wind", "Cloud", "Time of day (mins since dawn)", "Julian Day")
ylabs <- list(rep("Detection probability", 5))
## Axis limits
ylim <- list(c(0.40, 0.70))

## Helper function
source("src/functions/plot_detec_predic.R")

## Returns a list of plots
plots_out <- purrr::pmap(.l = list(cov = covs, covRange = ranges,
                       coeffname1 = coeffs1,coeffname2 = coeffs2,
                       intname = intname,xlab = xlabs,
                       ylab = ylabs, ylim = ylim),
                  .f = plot_detec_predic)

## Plot all figures together (otherwise write each to PDF or JPEG)
ggpubr::ggarrange(plots_out[[1]],
          plots_out[[2]] + rremove("y.text") + rremove("ylab"),
          plots_out[[3]],
          plots_out[[4]] + rremove("y.text") + rremove("ylab"),
          plots_out[[5]],
          ncol = 2, nrow = 3)


#'  
#' ## Lambda parameters
## -------------------------------------------------------------------------------------------------------------------------
## Extract species lambda estimates
lam_est_spp <- plot_data %>% 
  filter(str_detect(param, "lam.est.sp")) %>% 
  mutate(param = spp_vec) 

## Plot a subset of estimates
lam_est_spp %>% 
  filter(row_number() %in% odds) %>% 
  ggplot(aes(x = fct_reorder(param,mean), 
           y = mean, ymin = lower, ymax = upper)) +
  geom_pointrange(size = 1.0, fatten = 1.1) +
  scale_y_continuous(limits = c(0,0.0010))+
  coord_flip() + geom_hline(yintercept = 0, linetype = "dotted") +
  ylab(expression(lambda)) + xlab("")+ 
  theme_classic()+
  theme(axis.text.x = element_text(size = 10, color = "black"),
        axis.text.y = element_text(size = 10),
        axis.title.x.bottom = element_text(size = 22))



#'  
#' ## Species detection curves
## -------------------------------------------------------------------------------------------------------------------------
## Helper function
source("src/functions/plot_detec_probability.R")

## Select species of interest
spp_input <- c("Bokmakierie","Rufous-eared Warbler",
               "Cape Bunting","Yellow-bellied Eremomela")

## Or create plots for all species
spp_input <- spp_vec

## Define the maximum survey time for the predictions
tmax <-  1000

## Run plotting function
plots_out <- map2(.x = spp_input,
                 .y = tmax,
                 .f = plot_detec_probability)

## Plot on single page
ggarrange(plots_out[[2]] + rremove("xlab"),
          plots_out[[3]] + rremove("y.text") + rremove("ylab") + rremove("xlab"),
          plots_out[[1]], 
          plots_out[[4]] + rremove("y.text") + rremove("ylab"),
          ncol = 2, nrow = 2)


#'  
#' # Occupancy
#' 
#' ## Occupancy coefficient estimates
## -------------------------------------------------------------------------------------------------------------------------
## Extract beta parameters (random effects)
betas <- plot_data %>% 
  filter(str_detect(param,"mu.beta")) %>% 
  mutate(param = c("NDVI","Rainfall concentration","Elevation","TRI"))

betas

ggplot(betas, aes(x = fct_reorder(param,mean), y = mean, ymin = lower, ymax = upper)) +
  geom_pointrange() + 
  scale_y_continuous(limits = c(-1,1))+
  coord_flip() + 
  geom_hline(yintercept = 0, linetype = "dotted") +
  ylab("Estimate (+- 95% CI)") + xlab("")+ 
  theme_classic()+
  theme(axis.title = element_text(size = 14),
        axis.text.x = element_text(size = 14, color = "black"),
        axis.text.y = element_text(size = 14, color = "black"))


#' 
#' ## Predictions
## -------------------------------------------------------------------------------------------------------------------------
covs <- list(site_covs[["ndvi"]],site_covs[["map_ctn"]],site_covs[["elev"]],site_covs[["tri_med"]])
ranges <- list(c(0,0.5,0.001),c(30,50,1),c(300,1700,5),c(10,700,5))
coeffs1 <- list("mu.beta1","mu.beta2","mu.beta3","mu.beta4")
coeffs2 <- NA
intname <- "mu.psi"
xlabs <- list("NDVI","Rainfall concentration","Elevation","TRI")
ylabs <- list(rep("Occupancy probability",4))
ylim <- list(c(0,1))

source("src/functions/plot_occ_predic.R")

plots_out <- pmap(list(cov = covs, covRange = ranges,
                       coeffname1 = coeffs1,coeffname2 = coeffs2,
                       intname = intname,xlab = xlabs,
                       ylab = ylabs, ylim = ylim),
                  plot_occ_predic)

ggarrange(plots_out[[1]],
          plots_out[[2]] + rremove("y.text") + rremove("ylab"),
          plots_out[[3]],
          plots_out[[4]] + rremove("y.text") + rremove("ylab"),
          ncol = 2, nrow = 2)


#'  
#' 
#' ## Random effects 
## ----fig.height=12, fig.width=8-------------------------------------------------------------------------------------------
# Plot:: Occupancy random effects -----------------------------------------

mod <- rbind(as.matrix(ttd_modfit$samples$chain1),
             as.matrix(ttd_modfit$samples$chain2),
             as.matrix(ttd_modfit$samples$chain3))

att <- attributes(mod)$dimnames[2][[1]][]
coeffs <- c("beta1","beta2","beta3","beta4")
df_list <- list()

for(i in seq_along(coeffs)){
  ref <- which(grepl(coeffs[i],att) & !grepl(c("sd"),att) & !grepl(c("mu"),att)) # Exclude sd and mu
  df <- as_tibble(mod[,ref])
  
  mean  <- df %>% summarise_all(funs(mean)) %>% 
    gather(key = "coeff_name", value = "mean") %>% 
    select(mean)
  
  lower <- df %>% 
    summarise_all(funs(quantile(.,prob = c(0.025)))) %>% 
    gather(key = "coeff_name", value = "lower") %>% 
    select(lower)
  
  upper <- df %>% 
    summarise_all(funs(quantile(.,prob = c(0.975)))) %>% 
    gather(key = "coeff_name", value = "upper") %>% 
    select(upper)
  
  df_list[[i]] <- cbind(mean, lower, upper)
  
}

source("src/functions/plot_random_effects.R")

names(df_list) <-  coeffs 
coeffs <- list("beta1","beta2","beta3","beta4")
xlims <- list(c(-3,3),c(-2,2),c(-2,2),c(-3,3))
mainlabs <- list("NDVI","Rainfall concentration","Elevation","TRI")

par(mfrow = c(2,2), cex.lab = 1, cex.axis = 0.9)
pwalk(.l = list(pm = df_list, 
           coeff_name = coeffs, 
           main_lab = mainlabs, 
           xlims = xlims, 
           n = n_spp),
      .f = plot_random_effects)


#' 
#' ## Species probabilities
## -------------------------------------------------------------------------------------------------------------------------
## Extract species occupancy probability
psi_spp <- plot_data %>% 
  filter(str_detect(param, "lpsi")) %>% 
  filter(!str_detect(param,"sd|mu")) %>% 
  mutate(param = spp_vec) %>% 
  mutate_at(vars(mean:upper),inv.logit)

psi_spp %>%  
  filter(row_number() %in% odds) %>% 
  ggplot(aes(x = fct_reorder(param,mean), y = mean, ymin = lower, ymax = upper)) +
  geom_pointrange(size = 0.6) + theme_bw() + scale_y_continuous(limits = c(0,1))+
  coord_flip() + geom_hline(yintercept = 0, linetype = "dotted") +
  ylab("Occupancy probability") + xlab("")+
  theme(axis.text.y=element_text(size=13),
        axis.text.x=element_text(size=13),
        axis.title.x = element_text(size=13))

#' 
#' ## Number of sites occupied
## -------------------------------------------------------------------------------------------------------------------------
## Extract the number of sites occupied by each species
occ_fs <- plot_data %>% 
  filter(str_detect(param, "occ.fs")) %>% 
  mutate(param = spp_vec) 

ggplot(occ_fs, aes(x = fct_reorder(param,mean), y = mean, ymin = lower, ymax = upper)) +
  geom_pointrange() + theme_bw() + scale_y_continuous(limits = c(0,60))+
  coord_flip() + 
  ylab("Number of sites occupied by each species") + xlab("")+
  theme(axis.text.y=element_text(size=8),
        axis.text.x=element_text(size=16),
        axis.title=element_text(size=16))



#'  
#' ## Species richness
## -------------------------------------------------------------------------------------------------------------------------
## Extract species richness estimates
Nsite <- plot_data %>% 
  filter(str_detect(param,"Nsite")) %>% 
  mutate(param = as.factor(pen_vec))

ggplot(Nsite, aes(x = fct_reorder(param,mean), y = mean, ymin = lower, ymax = upper)) +
  geom_pointrange() + theme_bw() + scale_y_continuous(limits =  c(15,max(Nsite$upper)+2))+
  coord_flip() + geom_hline(yintercept = 0, linetype = "dotted") +
  ylab("Species richness") + xlab("Pentad")+
  theme(axis.text.y=element_text(size=8),
        axis.text.x=element_text(size=10),
        axis.title=element_text(size=17))



#'  
#' # Helper functions 
#' 
#' ## `plot_detec_predic`
#' 
## -------------------------------------------------------------------------------------------------------------------------
plot_detec_predic <- function(cov, covRange, coeffname1, coeffname2,
                              intname, xlab, ylab, ylim) {
  
  vec <- which(!is.na(cov))
  vec <- cov[vec]
  
  mean <- attr(scale(vec), "scaled:center") # mean 
  sd <- attr(scale(vec), "scaled:scale")  # standard deviation
  orig.pred <- seq(from = covRange[1], to = covRange[2], by = covRange[3])# Predictor (x-values) which are unscaled 
  sc.pred <- (orig.pred - mean)/sd 

  coef1samples <- c(ttd_modfit$samples[[1]][,coeffname1],
                    ttd_modfit$samples[[2]][,coeffname1],
                    ttd_modfit$samples[[3]][,coeffname1])
  
  intsamples <- c(ttd_modfit$samples[[1]][,intname],
                  ttd_modfit$samples[[2]][,intname],
                  ttd_modfit$samples[[3]][,intname])
  
  predictions <- array(dim = c(length(sc.pred), length(coef1samples)))
  dim(predictions)
  
  if (is.na(coeffname2)) {
    
    for(i in 1:length(sc.pred)){
      
      predictions[i,] <- 1/(1+ exp(-(exp(intsamples) + 
                                       coef1samples * sc.pred[i])))
    }  
  } else {
    
    coef2samples <- c(ttd_modfit$samples[[1]][,coeffname2],
                      ttd_modfit$samples[[2]][,coeffname2],
                      ttd_modfit$samples[[3]][,coeffname2])
    
    for(i in 1:length(sc.pred)){
      predictions[i,] <- 1/(1+ exp(-(exp(intsamples) + 
                                       coef1samples * sc.pred[i]+
                                       coef2samples * sc.pred[i]^2)))
    }  
    
  }
  
  LPB <-  apply(predictions, 1, quantile, probs = 0.025) # Lower bound
  UPB <-  apply(predictions, 1, quantile, probs = 0.975) # Upper bound
  y <- apply(predictions, 1, mean) # Upper bound
  ylim <- ylim
  
  plotdf <- tibble(xcol = orig.pred, ycol = y, lower = LPB, upper = UPB)
  
  finalplot <- ggplot(plotdf, aes(x = xcol, y = ycol)) +
    geom_line(colour="black", size = 1) +
    geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.2) +
    xlab(xlab) +
    ylab(ylab) +
    #scale_x_continuous(breaks = c(8:14))+
    scale_y_continuous(limits = ylim)+
    theme(axis.text.x=element_text(size=14, colour = "black"),
          axis.text.y=element_text(size=14, colour = "black"),
          axis.title = element_text(size = 16,margin = margin(t = 0, r = 20, b = 100, l = 20)),
          panel.grid = element_blank(),panel.background = element_blank(),
          panel.border = element_blank(),
          axis.line = element_line(size = 1),
          axis.ticks = element_line(color = "black",size = 1.2),
          axis.ticks.length = unit(0.2,"cm"))
  
  return(finalplot)
  
}

#' 
#' ## `plot_detec_probability`
#' 
## -------------------------------------------------------------------------------------------------------------------------
plot_detec_probability <- function(spp, tmax) {
  
  ref <- which(lam_est_spp$param %in% spp)
  mean.lambda <- as.numeric(lam_est_spp[ref,2])
  upper <- as.numeric(lam_est_spp[ref,4])
  lower <- as.numeric(lam_est_spp[ref,3])
  
  duration <- 1:tmax 
  p.pred <- 1 - exp(-mean.lambda*duration)
  upp <- 1-exp(-upper*duration)
  low <- 1-exp(-lower*duration)
    
  plotdf <- data_frame(xcol = duration, ycol = p.pred, lower = low, upper = upp)
    
  ggplot(plotdf, aes(x = xcol, y = ycol)) +
      geom_line(colour="black", size = 1) +
      geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.2) +
      labs(title = spp) +
      xlab("\nSurvey duration (seconds)") +
      ylab("Detection probability\n") +
      scale_y_continuous(limits = c(0,0.5))+
      theme(axis.text.x=element_text(size=12, colour = "black"),
            axis.text.y=element_text(size=12, colour = "black"),
            axis.title = element_text(size = 12,margin = margin(t = 0, r = 20, b = 100, l = 20)),
            panel.grid = element_blank(),panel.background = element_blank(),
            panel.border = element_blank(),
            axis.line = element_line(size = 1),
            axis.ticks = element_line(color = "black",size = 1.2),
            axis.ticks.length = unit(0.2,"cm"),
            plot.title = element_text(size = 15))
  
 
}

#' 
#' ## `plot_random_effects`
#' 
## -------------------------------------------------------------------------------------------------------------------------
plot_random_effects <- function(pm,coeff_name,main_lab,xlims,n){
  
    plot(pm[["mean"]] ,1:n, xlim = xlims, 
       xlab = "Parameter estimate", ylab = "", 
       main = main_lab, pch = 16)
  abline(v = 0, lwd = 2, col = "green")
  
  segments(pm[1:n,"lower"], 1:n, pm[1:n,"upper"], 1:n, col = "grey", lwd = 1.5)
  sig1 <- (pm[1:n,"lower"] * pm[1:n,"upper"]) > 0
  
  segments(pm[1:n,"lower"][sig1 == 1], (1:n)[sig1 == 1], pm[1:n,"upper"][sig1 == 1], (1:n)[sig1 == 1], col = "blue", lwd = 1.5)
  points(pm[["mean"]], 1:n, pch = 16)
  
  abline(v = plot_data[which(plot_data$param == paste0("mu.",coeff_name)),2], col = "red", lwd = 3)
  abline(v = plot_data[which(plot_data$param == paste0("mu.",coeff_name)),3:4], col = "red", lwd = 2, lty = 2)
  
  }

#' 
#' ## `plot_occ_predic`
#' 
## -------------------------------------------------------------------------------------------------------------------------
plot_occ_predic <- function(cov,covRange,coeffname1,coeffname2,intname,xlab,ylab,ylim) {
  
  vec <- which(!is.na(cov))
  vec <- cov[vec]
  
  mean <- attr(scale(vec), "scaled:center") # mean 
  sd <- attr(scale(vec), "scaled:scale")  # standard deviation
  orig.pred <- seq(from = covRange[1], to = covRange[2], by = covRange[3])# Predictor (x-values) which are unscaled 
  sc.pred <- (orig.pred - mean)/sd 
  
  coef1samples <- c(ttd_modfit$samples[[1]][,coeffname1],
                    ttd_modfit$samples[[2]][,coeffname1],
                    ttd_modfit$samples[[3]][,coeffname1])
  
  intsamples <- c(ttd_modfit$samples[[1]][,intname],
                  ttd_modfit$samples[[2]][,intname],
                  ttd_modfit$samples[[3]][,intname])
  
  predictions <- array(dim = c(length(sc.pred), length(coef1samples)))
  dim(predictions)
  
  if (is.na(coeffname2)) {
    
    for(i in 1:length(sc.pred)){
      
      predictions[i,] <- plogis(logit(intsamples) +
                                  coef1samples * sc.pred[i])
      
    }  
  } else {
    
    coef2samples <- c(ttd_modfit$samples[[1]][,coeffname2],
                      ttd_modfit$samples[[2]][,coeffname2],
                      ttd_modfit$samples[[3]][,coeffname2])
    
    for(i in 1:length(sc.pred)){
      predictions[i,] <- plogis(logit(intsamples) +
                                  coef1samples * sc.pred[i] + 
                                  coef2samples * sc.pred[i]^2)
    }  
    
  }
  
  LPB <-  apply(predictions, 1, quantile, probs = 0.025) # Lower bound
  UPB <-  apply(predictions, 1, quantile, probs = 0.975) # Upper bound
  y <- apply(predictions, 1, mean) 
  ylim <- ylim
  
  plotdf <- data_frame(xcol = orig.pred, ycol = y, lower = LPB, upper = UPB)
  
  finalplot <- ggplot(plotdf, aes(x = xcol, y = ycol)) +
    geom_line(colour="black", size = 2) +
    geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.2) +
    xlab(xlab) +
    ylab(ylab) +
    scale_y_continuous(limits = ylim)+
    theme(axis.text.x=element_text(size=16, colour = "black"),
          axis.text.y=element_text(size=16, colour = "black"),
          axis.title = element_text(size = 16,margin = margin(t = 0, r = 20, b = 100, l = 20)),
          panel.grid = element_blank(),panel.background = element_blank(),
          panel.border = element_blank(),
          axis.line = element_line(size = 1),
          axis.ticks = element_line(color = "black",size = 1.2),
          axis.ticks.length = unit(0.2,"cm"))
  
  return(finalplot)
  
}

#' 
