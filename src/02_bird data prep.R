#' ---
#' title: "Bird data preparation"
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
library(lubridate)


#' 
#' # Import data 
## -------------------------------------------------------------------------------------------------------------------------
## Point count data
birds <- read_csv("data input/bird_count_data.csv") 
birds

## Survey covariates
detec_covs <- read_csv("data input/survey_metadata.csv")
detec_covs

## Site covariates
site_covs <- read_csv("data input/site_covs.csv")

## map_ctn = Rainfall concentration
## elev = Mean elevation
## tri_med = Median Terrain Ruggedness Index
## ndvi = Mean NDVI over 2016 - 2018

site_covs


#' 
#' # Create TTD data frame 
## -------------------------------------------------------------------------------------------------------------------------

## Extract start times of each unique survey
start_times <- birds %>% 
  select (survey_id, datetime) %>% 
  group_by(survey_id) %>%  
  arrange(datetime) %>%
  filter(row_number()==1) %>% 
  rename(survey_start = datetime) %>% 
  ungroup

start_times

ttd <- left_join(birds, start_times, by = "survey_id") %>% 
  mutate(ttd = as.numeric(datetime - survey_start))  %>%  # ttd in seconds
  mutate(ttd = replace(ttd,ttd == 0, 10)) %>%             # Replace 0's with 10 secs
  filter(ttd <= 660) %>%                                  # Remove entries where ttd  > 11 mins (Tmax = 660s)
  group_by(survey_id, species) %>% 
  filter(row_number() == 1) %>%   # NB! Select the first occurrence of a species record in the survey - this is TTD!
  ungroup %>% 
  select(survey_id, pentad, species,ttd, datetime, survey_start)
  
ttd


#'  
#' ## Clean TTD data frame 
## -------------------------------------------------------------------------------------------------------------------------

## Remove rare, vagrant birds
rare_spp <- ttd %>% 
  count(species) %>% 
  filter(n < 20) %>% 
  pull(species)

rare_spp

ttd <- ttd %>% 
  filter(!species %in% rare_spp) %>% 
  filter(!is.na(species))


#' 
#' # Order detection covariates 
## -------------------------------------------------------------------------------------------------------------------------
## Match order of survey IDs in covariate and TTD data frames
detec_covs <- detec_covs %>% 
  filter(survey_id %in% ttd$survey_id)

#' 
#' # Create TTD variables and arrays 
## -------------------------------------------------------------------------------------------------------------------------

## Summary of the number of surveys in each site (pentad)
pen_sum <- ttd %>% 
  group_by(pentad) %>% 
  summarise(n_survey = length(unique(survey_id))) %>% 
  arrange(pentad) %>% 
  print(n = 20)


#' 
#' These are the values that are necessary to construct the arrays with the correct dimensions and extract the appropriate TTD and covariate values.
## -------------------------------------------------------------------------------------------------------------------------
## Create vector of pentad IDs
pen_vec <- pen_sum %>% 
  select(pentad) %>% 
  pull

## Create vector of species names
spp_vec <- ttd %>% 
  select(species) %>% 
  distinct(species) %>% 
  arrange(species) %>% 
  pull

## Number of species index
n_spp <- length(spp_vec)

## Number of pentads index
n_pen <- nrow(pen_sum)

## Maximum number of survey replicates across sites
n_counts <- max(pen_sum$n_survey)

## Create an empty 3D array for TTD observations
Y <- array(NA, dim = c(n_pen, n_counts,n_spp))

## Assign dimnames to array
dimnames(Y) <- list(pen_vec, NULL, spp_vec) 
## sites, survey replicates, species
dim(Y) 
dimnames(Y)

## Create detection covariate arrays
temp_arr <- wind_arr <- cloud_arr <- 
  mins_arr <- jday_arr <- array(NA, dim = c(n_pen, n_counts))

dimnames(temp_arr) <- dimnames(wind_arr) <- dimnames(cloud_arr) <- 
  dimnames(mins_arr) <- dimnames(jday_arr) <- list(pen_vec, NULL) 

#' 
#' # Populate the arrays 
## -------------------------------------------------------------------------------------------------------------------------

## Loop over each site
for(i in 1:n_pen){ 
  
  ## Extract TTD data for 1 site
  ttd_pen <- ttd %>% filter(pentad == pen_vec[i])
  
  ## Extract covariate data for 1 site
  detec_pen <- detec_covs %>% filter(pentad == pen_vec[i])
  
  ## Create a vector of survey IDs in the 1 site
  t_id <- as.character(unique(ttd_pen$survey_id))
 
  ## Loop over each survey within the site
  for(t in 1:length(t_id)){
    
    ## Extract TTD data for the survey
    ttd_surv <- ttd_pen %>% filter(survey_id  == t_id[t])
    
    ## Extract the species detected in the survey
    surv_spp <- as.character(ttd_surv$species)
   
     ## Assign those TTD value for the combination of site(i), 
    ## survey(t), and species(surv_spp) to the Y array
    Y[i,t,surv_spp] <- ttd_surv$ttd
    
  }
  
  ## Assign covariate data to arrays
  temp_arr[i,1:length(t_id)]  <- detec_pen$temperature
  wind_arr[i,1:length(t_id)]  <- detec_pen$wind
  cloud_arr[i,1:length(t_id)] <- detec_pen$cloud
  mins_arr[i,1:length(t_id)]  <- detec_pen$mins_after_dawn
  jday_arr[i,1:length(t_id)]  <- detec_pen$ju_rad
  
}


#'  
#' # Create indicator array
## -------------------------------------------------------------------------------------------------------------------------
d <- Y
d <- ifelse(d > 0, d <- 0, d <- 1)
ref <- which(is.na(d))
d[ref] <- 1

#' 
#' # Array checks 
## -------------------------------------------------------------------------------------------------------------------------
Y[1,1:10,1:4]

#'  
## -------------------------------------------------------------------------------------------------------------------------
d[1,1:10,1:4]

#'  
## -------------------------------------------------------------------------------------------------------------------------
temp_arr[1:10,1:5]

#'  
## -------------------------------------------------------------------------------------------------------------------------
mins_arr[1:10,1:5]

#' 
#' # Save workspace
## -------------------------------------------------------------------------------------------------------------------------
rm(list = c("t","ref","i","ttd_pen","detec_pen","t_id","ttd_surv", "surv_spp","start_times"))
save.image("data output/TTD_arrays.RData")

#'  
#'  
