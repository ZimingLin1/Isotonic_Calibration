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
list.of.packages <- c("data.table","tidyverse","caret","xgboost","mgcv","remotes", "ranger", "earth","glmnet")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(list.of.packages, library, character.only = TRUE)

# load datasets from ACIC 2017
# devtools::install_github("vdorie/aciccomp/2017")

library(sl3)
library(aciccomp2017)

source("CalibProp_Real_funcs.R")
############################################

# sbatch_index <- strtoi(commandArgs(trailingOnly=TRUE))
# sim_id <- readRDS("sim_ids.rds")[sbatch_index]
# 
# ACIC_2017_results <- sim_for_ACIC_2017_data(sim_id = sim_id, CATE_est_type = "xgb_ens", iters = 250)
# ATE_results <- ACIC_2017_results$evaluation_ATE
# CATE_results <- ACIC_2017_results$evaluation_CATE
# 
# ATE_dir <- paste0("results/ATE_xgb_", as.character(sim_id), ".csv")
# CATE_dir <- paste0("results/CATE_xgb_", as.character(sim_id), ".csv")
# 
# write.csv(ATE_results, ATE_dir, row.names=TRUE)
# write.csv(CATE_results, CATE_dir, row.names=TRUE)

sbatch_index <- strtoi(commandArgs(trailingOnly=TRUE))
sim_id <- readRDS("sim_ids.rds")[sbatch_index]

ATE_results <- sim_for_ACIC_2017_data2(sim_id = sim_id, iters = 250)

ATE_dir <- paste0("results/ATE_xgb_", as.character(sim_id), ".csv")
write.csv(ATE_results, ATE_dir, row.names=TRUE)

############################################


# ACIC_2017_results <- sim_for_ACIC_2017_data(sim_id = 17, CATE_est_type = "xgb_ens", iters = 1)

# load the saved CATE est.s
# sim_id <- 17
# true_CATE_list <- readRDS(paste0("results/CATE_est/", as.character(sim_id), "/", "true_CATE_list", ".rds"))
# calibration_types <- c("none", "deterministic truncation", "adaptive truncation", "isotonic regression")
# 
# for(calibration_method in calibration_types){
#   list_name <- paste0(calibration_method, "_list")
#   list_name <- gsub(" ", "_", list_name) # e.g. "isotonic_regression_list"
#   assign(list_name, readRDS(paste0("results/CATE_est/", as.character(sim_id), "/", list_name, ".rds")))
# }
# 
# anyNA(isotonic_regression_list, recursive = TRUE)
# 
# evaluation_CATE_all_calib_2017(250, true_CATE_list = true_CATE_list, none_list = none_list, deterministic_truncation_list = deterministic_truncation_list, adaptive_truncation_list = adaptive_truncation_list, isotonic_regression_list = isotonic_regression_list)
# 
# 
# for (sim_id in 17:17){
#   
#   true_CATE_list <- readRDS(paste0("results/CATE_est/", as.character(sim_id), "/", "true_CATE_list", ".rds"))
#   calibration_types <- c("none", "deterministic truncation", "adaptive truncation", "isotonic regression")
#   
#   for(calibration_method in calibration_types){
#     list_name <- paste0(calibration_method, "_list")
#     list_name <- gsub(" ", "_", list_name) # e.g. "isotonic_regression_list"
#     assign(list_name, readRDS(paste0("results/CATE_est/", as.character(sim_id), "/", list_name, ".rds")))
#   }
#   
#   CATE_results <- evaluation_CATE_all_calib_2017(250, true_CATE_list = true_CATE_list, none_list = none_list, deterministic_truncation_list = deterministic_truncation_list, adaptive_truncation_list = adaptive_truncation_list, isotonic_regression_list = isotonic_regression_list)
#   CATE_dir <- paste0("results/CATE_xgb_", as.character(sim_id), ".csv")
#   write.csv(CATE_results, CATE_dir, row.names=TRUE)
# }




# nsize <- 1000 # either 1000, 2500, 5000 or 10000
# params <- fread("acic2018_scaling_small/params.csv")
# ids <- params$ufid[params$size == nsize]
# iters <- seq_along(ids) # which is seq(nsize)
# 
# # simulation for instantiated dataset with size nsize, ATE case with xgb ensemble
# 
# true_ATE_list <- numeric(iters)
# 
# calibration_types <- c("none", "deterministic truncation", "adaptive truncation", "isotonic regression")
# df <- data.frame(matrix(ncol = length(calibration_types), nrow = iters))
# colnames(df) <- calibration_types
# 
# df_est <- df # a dataframe used to store estimated ATE
# df_std <- df # a dataframe used to store estimated std of ATE, for computing confidence intervals
# Multiplier_CI <-  qnorm(1-{1-0.95}/2)
# 
# # loop through the instantiated datasets of size n
# 
# true_CATE_list <- vector(mode='list', length=length(iters))
# for(calibration_method in calibration_types){
#   list_name <- paste0(calibration_method, "_list")
#   list_name <- gsub(" ", "_", list_name)
#   assign(list_name, vector(mode='list', length=length(iters)))
# }
# 
# 
# for(i in iters){
# 
#   id <- ids[i] # get simulation id
#   inst_ob <- fread(paste0("acic2018_scaling_small/factuals/", id, ".csv")) # i.e. the instantiated observation file, contains the treatment and outcome,
#   inst_cov <- fread(paste0("acic2018_scaling_small/x.csv")) # i.e. instantiated covariate file, full covariate information for the instantiated dataset
#   inst_cov <- inst_cov[match(inst_ob$sample_id, inst_cov$sample_id)] # subset to ids that appear in data, since the label files contan fewer samples than the covariate file
#   
#   inst_cf <- fread(paste0("counterfactuals/", id, "_cf.csv")) # i.e. counterfactuals for the instantiated dataset
#   inst_cf[, c("ITE")] <- inst_cf$y1 - inst_cf$y0
#   inst_cf <- inst_cf[, -c("y0", "y1")]
#   
#   inst_list_f <- list(inst_ob, inst_cov)
#   inst_dat <- inst_list_f %>% reduce(full_join, by="sample_id")
#   
#   # so up until now, we have a instantiated dataset and its corresponding covariate file and observation file
# 
#   # get those four DR estimated value of ATE with different calibration methods
# 
#   mu <- get_mu(inst_dat)
#   pi <- get_pi_est(inst_dat)
# 
#   for(calibration_method in calibration_types){
# 
#     alpha <- get_alpha_star(inst_dat, pi, calibration_method)
#     reg_dif <- mu$diff
#     riesz_prod <- alpha * (inst_dat$y - mu$mu_xa)
#     pseudo_outcome_est <- reg_dif + riesz_prod
#     inst_dat_pso <- inst_dat
#     inst_dat_pso[,c('pseudo_outcome_est')] <- pseudo_outcome_est
#     inst_dat_pso[,c('CATE_est')] <- get_CATE_pred(inst_dat_pso, regression_type = "xgb_ens")
#     
#     # record the estimated CATEs of this instantiated dataset in calib_method_list[[i]]
#     list_name <- paste0(calibration_method, "_list")
#     list_name <- gsub(" ", "_", list_name) #e.g. "isotonic_regression_list"
#     command_text <- paste0(list_name, "[[", as.character(i), "]] <- inst_dat_pso$CATE_est")
#     eval(parse(text = command_text))
# 
#   }
#   
#   # align the true ITE to estimated CATE according to the sample_id's
# 
#   inst_list_CATE <- list(inst_dat, inst_cf)
#   inst_dat_CATE <- inst_list_CATE %>% reduce(full_join, by="sample_id") 
#   true_CATE_list[[i]] <- inst_dat_CATE$ITE
#   
# }
# 
# # evaluations for CATE
# # draft
# for(calibration_method in calibration_types){
#   
#   list_name <- paste0(calibration_method, "_list")
#   list_name <- gsub(" ", "_", list_name) # list for estimated CATEs across all datasets with size n
#   
#   eval_name <- paste0(calibration_method, "_eval")
#   eval_name <- gsub(" ", "_", list_name)
#   
#   evaluation_CATE(true_CATE_list, get(list_name))
#   assign(eval_name, evaluation_CATE(true_CATE_list, get(list_name)))
# }
# 
# 
# 
# # evaluations for ATE
# evaluation_ATE(true_ATE_list, df_est$none, df_std$none)
# ATE_results_1000 <- evaluation_ATE_all_calib(true_ATE_list, df_est, df_std)
# write.csv(ATE_results_1000, "ATE_results_1000.csv", row.names=FALSE)
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# ############################################
# # for all evaluations
# 
# # test <- real_ATE_est_eval_n_varying()
# # ATE_eval_plots_test <- make_plots_for_ATE_results(test)
# # write.csv(test, "ATE_results_test.csv", row.names=FALSE)
# # ATE_eval_plots_test
# 
# # ATE_results <- real_ATE_est_eval_n_varying()
# # write.csv(ATE_results, "ATE_results.csv", row.names=FALSE)
# # ATE_eval_plots <- make_plots_for_ATE_results(ATE_results)
# # ATE_eval_plots
# # 
# # 
# # 
# 
# 
# ############################################
# 
# for(calibration_method in calibration_types){
#   
#   
#   # record the estimated CATEs of this instantiated dataset in calib_method_list[[i]]
#   list_name <- paste0(calibration_method, "_list")
#   list_name <- gsub(" ", "_", list_name) #e.g. "isotonic_regression_list"
#   command_text <- paste0(list_name, "[[", as.character(i), "]] <- inst_dat_pso$CATE_est")
#   
#   print(command_text)
#   eval(parse(text = command_text))
# }
# 
# 
# real_CATE_est_eval_test(1000)


# use dgp_2017(sim_id, i) to generate datasets; sim_id in 17:24, i in 1:250
# each datasets contains z, y, alpha (true CATE); 
# each sim_id correspond to a instantiated DGP
# to get covariates, use aciccomp2017::input_2017








# CATE_results <- real_CATE_est_eval_n_varying()
# CATE_eval <- CATE_results$CATE_eval
# CATE_ps <- CATE_results$ps_info
# 
# write.csv(CATE_eval, "results/CATE_evaluation_results.csv", row.names=FALSE)
# write.csv(CATE_ps, "CATE_ps_info.csv", row.names=FALSE)
# 
# CATE_eval <- read.csv("results/CATE_evaluation_results.csv", check.names = FALSE)
# 
# CATE_plots <- make_plots_for_CATE_results(CATE_eval)
# 
# pdf("CATE_plots.pdf")
# CATE_plots
# dev.off()
# 
# # for test
# sim_for_ACIC_2017_data(sim_id = 17, CATE_est_type = "xgb_ens", iters = 2)




############################################

# K_times <- 5000
# my_sample_range <- c(500, 1000, 1500, 2000, 2500, 3000, 3500, 4000, 4500, 5000)
# overlap_const <- 1
# 
# xgboost_evaluation_results_1 <- sim_ATE_est_with_pi_n_varying(K_times=K_times, n_samples_range = my_sample_range, overlap=overlap_const, outcome_reg_type="xgboost", ps_est_type="xgboost")
# xgb_plots_1 <- make_plots_for_results(xgboost_evaluation_results_1)
# 
# xgboost_evaluation_results_1
# xgb_plots_1
# 
# write.csv(xgboost_evaluation_results_1, "xgboost_evaluation_results_1.csv", row.names=FALSE)


###############

# dat_obs <- dgp_2017(17, 171)
# true_CATE_list[[k]] <- dat_obs[,c("alpha")]
# true_ATE[k] <- mean(dat_obs[,c("alpha")])
# dat <- list(dat_cov, dat_obs) %>% reduce(cbind)
# head(subset(dat, select = -c(alpha)))
# 
# 
# 
# calibration_method <- "isotonic regression"
# list_name <- paste0(calibration_method, "_list")
# list_name <- gsub(" ", "_", list_name) # e.g. "isotonic_regression_list"
# command_text <- paste0(list_name, "[[", as.character(17), "]] <- inst_dat_pso$CATE_est")
# 
# 
# # test with small iters (in 1:250)
# 
# sim_for_ACIC_2017_data(sim_id = 17, CATE_est_type = "xgb_ens", iters = 2)

