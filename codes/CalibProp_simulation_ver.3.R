# some error occured when installing sl3... This code is for trying to fix that error
# this chunk is for removing all the packages installed
#
# ip <- as.data.frame(installed.packages())
# head(ip)
# # if you use MRO, make sure that no packages in this library will be removed
# ip <- subset(ip, !grepl("MRO", ip$LibPath))
# # we don't want to remove base or recommended packages either\
# ip <- ip[!(ip[,"Priority"] %in% c("base", "recommended")),]
# # determine the library where the packages are installed
# path.lib <- unique(ip$LibPath)
# # create a vector with all the names of the packages you want to remove
# pkgs.to.remove <- ip[,1]
# head(pkgs.to.remove)
# # remove the packages
# sapply(pkgs.to.remove, remove.packages, lib = path.lib)

rm(list = ls())
graphics.off()
#install/load relevant packages
list.of.packages <- c("tidyverse","caret","xgboost","bootstrap","mgcv","remotes", "ranger", "earth")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, type="binary")
lapply(list.of.packages, library, character.only = TRUE)


# usethis::edit_r_environ()
# input GITHUB_TOKEN="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"
remotes::install_github("tlverse/sl3")
library(sl3)

source("CalibProp_simulation_ver.2.3.R")

n_sample <- 50000
set.seed(342)
dat <- DGP(n_sample, overlap=1)

################################################################################
# some examples for estimating AIPW:

# given a dataset, we have two functions for estimating the ATE:

# 1. get_ATE_AIPW_with_pi(dat, outcome_reg_type, ps_est_type, calibration_method) is for those
#    ATE estimators first estimate pi, calibrate it and then transform it in to alpha; (We call this the first type AIPW estimator)
#    outcome_reg_type can takes value in c("logit", "GAM", "earth", "random forest", "xgboost", "oracle");
#    ps_est_type can also takes value in c("logit", "GAM", "earth", "random forest", "xgboost", "oracle");
#    calibration_method can takes value in c("deterministic truncation", "adaptive truncation", "isotonic regression", "none");

# 2. get_ATE_AIPW_with_alpha(dat, outcome_reg_type, alpha_est_type) is for those ATE estimators
#    uses directly estimated alpha; (We call this the second type AIPW estimator)
#    again, outcome_reg_type can takes value in c("logit", "GAM", "earth", "random forest", "xgboost", "oracle");
#    alpha_est_type takes value in c("random forest", "xgboost")

# Moreover, get_true_ATE() is for getting the true ATE value using simulation;
################################################################################

get_true_ATE()

get_ATE_AIPW_with_pi(dat, "xgboost", "xgboost", "none")
get_ATE_AIPW_with_pi(dat, "xgboost", "xgboost", "deterministic truncation")
get_ATE_AIPW_with_pi(dat, "xgboost", "xgboost", "adaptive truncation")
get_ATE_AIPW_with_pi(dat, "xgboost", "xgboost", "isotonic regression")

get_ATE_AIPW_with_pi(dat, "logit", "logit", "none")
get_ATE_AIPW_with_pi(dat, "logit", "logit", "deterministic truncation")
get_ATE_AIPW_with_pi(dat, "logit", "logit", "adaptive truncation")
get_ATE_AIPW_with_pi(dat, "logit", "logit", "isotonic regression")

get_ATE_AIPW_with_pi(dat, "GAM", "earth", "none")
get_ATE_AIPW_with_pi(dat, "GAM", "earth", "deterministic truncation")
get_ATE_AIPW_with_pi(dat, "GAM", "earth", "adaptive truncation")
get_ATE_AIPW_with_pi(dat, "GAM", "earth", "isotonic regression")

get_ATE_AIPW_with_pi(dat, "oracle", "oracle", "none")

################################################################################
# some examples for diagnostic Monte Carlo simulations:

# since the running time is too long, we're evaluating one AIPW estimator at a time here

# sim_one_ATE_est_with_pi(K_times, n_samples,  outcome_reg_type, ps_est_type, calibration_method)
# is for evaluating the a first type AIPW estimator, specified by outcome_reg_type, ps_est_type
# and calibration_method as before

# sim_one_ATE_est_with_alpha(K_times, n_samples,  outcome_reg_type, alpha_est_type)
# is for evaluating the a second type AIPW estimator, specified by outcome_reg_type, alpha_est_type as before

# In both functions, K_times determines the time for iteration, and n_samples determines
# how many data we generate in a single iteration
################################################################################


sim_one_ATE_est_with_pi(K_times=50, n_samples=10000, outcome_reg_type="earth", ps_est_type="earth", calibration_method="none")
sim_one_ATE_est_with_pi(K_times=50, n_samples=10000, outcome_reg_type="earth", ps_est_type="earth", calibration_method="deterministic truncation")
sim_one_ATE_est_with_pi(K_times=50, n_samples=10000, outcome_reg_type="earth", ps_est_type="earth", calibration_method="adaptive truncation")
sim_one_ATE_est_with_pi(K_times=50, n_samples=10000, outcome_reg_type="earth", ps_est_type="earth", calibration_method="isotonic regression")

sim_one_ATE_est_with_pi(K_times=100, n_samples=10000, outcome_reg_type="xgboost", ps_est_type="xgboost", calibration_method="none")
sim_one_ATE_est_with_pi(K_times=100, n_samples=10000, outcome_reg_type="xgboost", ps_est_type="xgboost", calibration_method="deterministic truncation")
sim_one_ATE_est_with_pi(K_times=100, n_samples=10000, outcome_reg_type="xgboost", ps_est_type="xgboost", calibration_method="adaptive truncation")
sim_one_ATE_est_with_pi(K_times=100, n_samples=10000, outcome_reg_type="xgboost", ps_est_type="xgboost", calibration_method="isotonic regression")

sim_one_ATE_est_with_pi(K_times=1000, n_samples=5000,  outcome_reg_type="logit", ps_est_type="xgboost", calibration_method="isotonic regression")

sim_one_ATE_est_with_pi(K_times=10, n_samples=5000,  outcome_reg_type="logit", ps_est_type="xgboost", calibration_method="isotonic regression")
sim_one_ATE_est_with_pi(K_times=10, n_samples=5000,  outcome_reg_type="xgboost", ps_est_type="xgboost", calibration_method="isotonic regression")

sim_one_ATE_est_with_pi(K_times=100, n_samples=10000,  outcome_reg_type="xgboost", ps_est_type="xgboost", calibration_method="isotonic regression")

sim_some_ATE_est_with_pi(K_times=5, n_samples=2000, outcome_reg_type="xgboost", ps_est_type="xgboost")
sim_some_ATE_est_with_pi(K_times=1000, n_samples=20000, outcome_reg_type="xgboost", ps_est_type="xgboost")

# quick example
my_sample_range <- seq(1000, 5000, by=2000)
sim_ATE_est_with_pi_n_varying(K_times=10, n_samples_range = my_sample_range, overlap=1, outcome_reg_type="xgboost", ps_est_type="xgboost")

# investigating in using isotonic regression for xgboost
my_sample_range <- c(500, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000,12000,14000,16000,18000,20000)
xgboost_evaluation_results <- sim_ATE_est_with_pi_n_varying(K_times=1000, n_samples_range = my_sample_range, overlap=1, outcome_reg_type="xgboost", ps_est_type="xgboost")
write.csv(xgboost_evaluation_results, "D:/InvProp_project/xgboost_evaluation_results.csv", row.names=FALSE)
xgb_plots <- make_plots_for_results(xgboost_evaluation_results)

my_sample_range2 <- c(500, 1000, 2000, 3000, 4000, 5000)
xgbooost_evaluation_results2 <- sim_ATE_est_with_pi_n_varying(K_times=1000, n_samples_range = my_sample_range2, overlap=1, outcome_reg_type="xgboost", ps_est_type="xgboost")
write.csv(xgbooost_evaluation_results2, "D:/InvProp_project/xgboost_evaluation_results2.csv", row.names=FALSE)
xgboost_plots2 <- make_plots_for_results(xgbooost_evaluation_results2)

my_sample_range3 <- c(10000,12000,14000,16000,18000,20000)
xgboost_evaluation_results3 <- sim_ATE_est_with_pi_n_varying(K_times=100, n_samples_range = my_sample_range3, overlap=1, outcome_reg_type="xgboost", ps_est_type="xgboost")
make_plots_for_results(xgboost_evaluation_results3)

# sim_one_ATE_est_with_pi(K_times=50, n_samples=10000,  outcome_reg_type="logit", ps_est_type="oracle", calibration_method="isotonic regression")
# 
# sim_one_ATE_est_with_alpha(K_times=20, n_samples=10000,  outcome_reg_type="xgboost", alpha_est_type="random forest")
# sim_one_ATE_est_with_alpha(K_times=20, n_samples=10000,  outcome_reg_type="xgboost", alpha_est_type="xgboost")

