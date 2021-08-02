# Author: Alex Paynter
# What: An R script which creates fake data, designed to be similar to eICE in 
#   a few important ways:  (1) variable number of home observations in each
#   subject and (2) higher variability in home spirometry data.

library(dplyr)
library(tibble)
library(magrittr)
library(stringr)
library(readr)
library(ggplot2) # optional - commented out plots

set.seed(6)




######################################
### Define a few helper functions: ###
######################################

# Each of these functions generates a subject-level parameter (usually latent).

# Random generation of home spirometry followup length.
# Our generation is a mixture distribution for three groups:
# 1. near full follow up (300-400 days), about 65%
# 2. partial followup (30-300 days), about 25%
# 3. almost no followup (<30 days), about 10%
r_home_t <- function(n_subj) {
    mix_vec <- runif(n = n_subj)
    dplyr::case_when(
        mix_vec < 0.65 ~ rnorm(n_subj, 355, 18),
        mix_vec < 0.9 ~ runif(n_subj, 30, 300),
        T ~ runif(n_subj, 0, 30)
    ) %>%
        floor
}

# A helper function for the multinomial distribution.
# Normalizes probabilities and draws one observation.
mult_helper <- function(n, p) {
    p <- p/sum(p) # normalize probability vector
    rmultinom(n = n, size = 1, prob = p) %>%
        apply(., MARGIN = 2, FUN = (function(x) which(x == 1)))
}

# random draw for the last clinic visit observed.
# In eICE (and obviously) this was strongly related to home followup length,
#  so one of the imputs is the home_fu generated with r_home_t().
r_max_visit <- function(n_subj, home_fu) {
    dplyr::case_when(
        # p at each cut is based on empirical counts. 
        home_fu < 100 ~ mult_helper(n_subj, p = c(8,4,2,0,8)),
        home_fu < 200 ~ mult_helper(n_subj, p = c(1,5,3,1,8)),
        home_fu < 300 ~ mult_helper(n_subj, p = c(0,1,2,4,4)),
        # in the case where people had full home followup virtually 100%
        #   also had full clinic followup:
        T ~ as.integer(5) # 5 meaning we assume 5 clinic visits.
    )
}







###############################
### subject level variables ###
###############################

n_subj <- 150 # number of subjects we want to simulate (~133 in eICE)

df_subj_level <- tibble(
        subject = paste0("s_", str_pad(1:n_subj, side = "left", pad = "0",
                                       width = ceiling(log10(n_subj+1))))
        ) %>%
    mutate(
        # helper functions above to signify days followed at home and
        #   last clinic visit:
        max_t_home = r_home_t(n_subj = n()),
        max_vis_clin = r_max_visit(n_subj = n(), home_fu = max_t_home),
        
        # probability of being observed each day from 0 to max_t_home
        freq_home = rbeta(n(), 2.81, 2.86)/2.5, 
        
        # subject-level true slopes (population truth = -2 ppFEV/yr):
        slope_clin = rnorm(n(), -2, 0.74),
        slope_home = rnorm(n(), -2, 0.74*2)
    )

# The subject-level intercept correlation obviously needs to be very high to 
#   be realistic, so we generate those from a multivariate normal:
df_ints <- MASS::mvrnorm(
    n = n_subj, 
    mu = c(78.7, 75.7), # True cross sectional difference: -3 ppFEV1
    Sigma = matrix(c(400, 350, 350, 400), ncol = 2)
) %>%
    as_tibble %>%
    set_names(c("intercept_clin", "intercept_home")) 
# add those intercepts to data:
df_subj_level %<>% bind_cols(., df_ints)
rm(df_ints) # cleanup







############################
### generating home data ###
############################

study_length <- 400 # in days - the maximum home followup permitted.
df_home <- df_subj_level %>%
    dplyr::select(subject, max_t_home, freq_home, slope_home, intercept_home) %>%
    # create a row for each day, each subject:
    slice(rep(1:n(), each = (study_length+1))) %>%
    # populated each subject/time row with a time marker:
    group_by(subject) %>%
    mutate(t_home = 0:(n()-1))

# We now use our subject-level characteristics to filter the dataframe down:
df_home %<>%
    filter(t_home <= max_t_home) %>%
    mutate(obs_rand_draw = runif(n())) %>%
    # keep rows randomly depending on that subject's observation frequency:
    # we also keep the time zero row (usually occurred in eICE):
    filter(obs_rand_draw <= freq_home | t_home == 0) %>%
    ungroup(.)

# now we can generate ppFEV over time for each subject:
df_home %<>%
    mutate(t_home_yr = t_home/365) %>%
    mutate(ppfev_home = intercept_home + slope_home*t_home_yr,
           # add residual noise on par with what we observed:
           ppfev_home = ppfev_home + rnorm(n(), 0, 7.9))

df_home %<>% # cleanup
    dplyr::select(subject, t_home, t_home_yr, ppfev_home)





##############################
### generating clinic data ###
##############################

max_visit <- 5 # by visit number - visit 5 was at about 1 year.
df_clin <- df_subj_level %>%
    dplyr::select(subject, max_vis_clin, slope_clin, intercept_clin) %>%
    # create a row for each visit, each subject:
    slice(rep(1:n(), each = max_visit)) %>%
    # populated each subject/time row with a time marker:
    group_by(subject) %>%
    mutate(visit_num = 1:n(),
           # timing for clinic visits was a +/-7d window, adherence to that timing
           #   was actually quite good.  And time zero is defined by the 
           #   first clinic visit, so we fix that firmly:
           t_clin = dplyr::case_when(
               visit_num %in% 1 ~ 0,
               T ~ round((visit_num-1)*91 + rnorm(n(), mean = 0, sd = 5))
           )
    )

# We now use our subject-level characteristics to filter the dataframe down:
df_clin %<>%
    filter(visit_num <= max_vis_clin) %>%
    mutate(obs_rand_draw = runif(n())) %>%
    # For the clinic data we will just put in a 5% probability of randomly
    #  missing followup visits (on top of the loss to followup above). 
    filter(obs_rand_draw <= 0.95 | visit_num == 1) %>%
    ungroup(.)

df_clin %<>%
    mutate(t_clin_yr = t_clin/365) %>%
    mutate(ppfev_clin = intercept_clin + slope_clin*t_clin_yr,
           # residual variation was substantially lower in clinic:
           ppfev_clin = ppfev_clin + rnorm(n(), 0, 4.3))

df_clin %<>% # cleanup
    dplyr::select(subject, visit_num, t_clin, t_clin_yr, ppfev_clin)


#####################
### save datasets ###
#####################
    
readr::write_csv(x = df_home, file = "home.csv")
readr::write_csv(x = df_clin, file = "clin.csv")


