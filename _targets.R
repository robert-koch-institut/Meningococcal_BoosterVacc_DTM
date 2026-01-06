# _targets.R file
library(targets)
library(tarchetypes)
source("R/functions.R")
source("R/ode_funs.R")

# Load packages
library(tidyverse)
library(crew)
library(purrr)
library(deSolve)
library(mgcv)
library(hhh4contacts)
library(optimParallel)
library(parallel)
# Setup
theme_set(theme_bw())
options(dplyr.summarise.inform = FALSE)

#' For parallelization
#' Note that computation of some targets requires up to 10 cores
#' For 3 workers >30 cores are required
controller <- crew::crew_controller_local(
  name = "my_controller",
  workers = 3,
  seconds_idle = 20
)

tar_option_set(packages = c("tidyverse", 
                            "purrr", 
                            "deSolve", 
                            "mgcv",
                            "hhh4contacts",
                            "optimParallel",
                            "mvtnorm",
                            "parallel"),
               memory = "transient", 
               garbage_collection = TRUE,
               error = "trim",
               controller = controller)
# Helper functions and objects
make_csv = function(table, path) {
  write_csv2(x=table, file = path)
  path
}

make_png = function(fig, path, width=9, height=3) {
  ggsave(filename = path, plot = fig, width = width, height = height)
  path
}

vacc_scen_values = tibble( # Use all possible combinations of input settings.
  names = paste0("scen", 0:5),
  vacc_scen = 0:5,
  parms = c(1, 1, 2, 3, 2, 3)
)

# fit_scen_nb = tibble(
#   names = paste0("nb_", rep(c("lambda1", "lambda2"), each=2), "_",
#                  rep(c("age4", "age8"), times=2)),
#   mod = rep(c("lambda1", "lambda2"), each=2),
#   age_cat = rep(c("4", "8"), times=2))

# Analysis pipeline
list(
  # Load and preprocess data
  tar_target(pop_dat_destatis_file,
             "./data/12411-0005.csv", format = "file"),
  #' Create contact matrix
  tar_target(cont_mat_main,
             make_cont_mat(filepath_destatis = pop_dat_destatis_file)),
  
  ############
  # Simulation
  ############
  # Population size during simulation
  tar_target(pop_size_sim_dest,
             read_destatis_20_49(pop_dat_destatis_file = pop_dat_destatis_file)),
  #' Parameters main simulation
  #' Population data
  tar_target(pop_size_sim_dest_tibble,
             tibble(year = rep(2020:2049, 86),
                    age = rep(0:85, each = 30),
                    n = pop_size_sim_dest %>% unlist())
  ),
  #' Model parameters
  tar_target(parms_sim_main_01_24_35,
             list(
               list(r_m = c(2, 2, 2),
                    vacc_eff_mmCc = c(0.85, 0, 0),
                    vacc_eff_pmACWYc = c(0.16, 0.16, 0),
                    waning_mmCc = c(rep(1/4, 12), rep(1/10, 74)),
                    waning_pmACWYc = c(rep(1/4, 12), rep(1/10, 74)),
                    cont_matrix = cont_mat_main,
                    pop_size = pop_size_sim_dest),
               # MenC Booster
               list(r_m = c(2, 2, 2),
                    vacc_eff_mmCc = c(0.85, 0, 0),
                    vacc_eff_pmACWYc = c(0.16, 0.16, 0),
                    waning_mmCc = c(rep(1/4, 12), rep(1/10, 74)),
                    waning_pmACWYc = c(rep(1/4, 12), rep(1/10, 74)),
                    cont_matrix = cont_mat_main,
                    pop_size = pop_size_sim_dest),
               # MenACWY Booster
               list(r_m = c(2, 2, 2),
                    vacc_eff_mmCc = c(0.85, 0, 0),
                    vacc_eff_pmACWYc = c(0.16, 0.16, 0),
                    waning_mmCc = c(rep(1/4, 12), rep(1/10, 74)),
                    waning_pmACWYc = c(rep(1/4, 12), rep(1/10, 74)),
                    cont_matrix = cont_mat_main,
                    pop_size = pop_size_sim_dest))),
  #' Load objects from fitting
  #' Summary File from fitting model 2B
  tar_target(fit_file,
             "./data/fit_nb_lambda2_age8", format = "file"),
  tar_target(fit_nb_lambda2_age8,
             readRDS(fit_file)),
  #' Integrated model 2B
  tar_target(int_file,
             "./data/int_nb_lambda2_age8", format = "file"),
  tar_target(int_nb_lambda2_age8,
             readRDS(int_file)),
  #' Compartments over time during calibration
  tar_target(comp_file,
             "./data/comp_nb_lambda2_age8", format = "file"),
  tar_target(comp_nb_lambda2_age8,
             readRDS(comp_file)),
  #' Parameter samples from MLE
  tar_target(parms_file,
             "./data/parameter_samples_nb_lambda2_age8", format = "file"),
  tar_target(parameter_samples_nb_lambda2_age8,
             readRDS(parms_file)),
  #' Integrated model 2B for each parameter sample
  tar_target(int_parms_file,
             "./data/int_samples_nb_lambda2_age8_small", format = "file"),
  tar_target(int_samples_nb_lambda2_age8,
             readRDS(int_parms_file)),
  #' Case-carrier ratio data
  tar_target(ccr_file,
             "./data/ccr_main"),
  tar_target(ccr_main,
             readRDS(ccr_file)),
  #' Vaccine effectiveness against IMD
  tar_target(ve_d_file,
             "./data/ve_d_main"),
  tar_target(ve_d_main,
             readRDS(ve_d_file)),
  #' CFR data
  tar_target(cfr_file,
             "./data/cfr"),
  tar_target(cfr,
             readRDS(cfr_file)),
  #' Sequelae probability
  tar_target(seq_prob_file,
             "./data/seq_prob"),
  tar_target(seq_prob,
             readRDS(seq_prob_file)),
  # Simulate Scenarios
  # Simulation with NegBinom Fit
  tar_map(
    values = vacc_scen_values,
    names = "names", # Select columns from `values` for target names.
    tar_target(sim_nb_lambda2_age8,
               forward_sim_lambda2_8agecat(fit_nb_lambda2_age8$par[1:24],
                                           int_nb_lambda2_age8,
                                           vacc_scen = vacc_scen,
                                           grid_len=1/12,
                                           parms=parms_sim_main_01_24_35[[parms]],
                                           rootfun=rootfunc)),
    tar_target(sim_samples_nb_lambda2_age8,
               mclapply(1:length(parameter_samples_nb_lambda2_age8),
                        function(x) {
                          forward_sim_lambda2_8agecat(parameter_samples_nb_lambda2_age8[[x]][1:24],
                                                      int_samples_nb_lambda2_age8[[x]],
                                                      vacc_scen = vacc_scen,
                                                      grid_len=1/12,
                                                      parms=parms_sim_main_01_24_35[[parms]],
                                                      rootfun=rootfunc)
                        },
                        mc.cores = 10)),
    tar_target(comp_sim_nb_lambda2_age8,
               get_comp_time(integrate = sim_nb_lambda2_age8,
                             begin_year = 2020)),
    tar_target(comp_sim_samples_nb_lambda2_age8,
               mclapply(1:length(parameter_samples_nb_lambda2_age8),
                        function(x) {
                          get_comp_time(integrate = sim_samples_nb_lambda2_age8[[x]],
                                        begin_year = 2020)
                        },
                        mc.cores = 10)),
    tar_target(vacc_num_nb_lambda2_age8,
               get_vacc_time(sim_nb_lambda2_age8,
                             begin_year = 2020)),
    tar_target(vacc_num_samples_nb_lambda2_age8,
               mclapply(1:length(parameter_samples_nb_lambda2_age8),
                        function(x) {
                          get_vacc_time(integrate = sim_samples_nb_lambda2_age8[[x]],
                                        begin_year = 2020)
                        },
                        mc.cores = 10)),
    tar_target(carr_inc_sim_nb_lambda2_age8,
               get_exp_new_inc_carr(sim_nb_lambda2_age8, years = 1:30)),
    tar_target(carr_inc_sim_samples_nb_lambda2_age8,
               mclapply(1:length(parameter_samples_nb_lambda2_age8),
                        function(x) {
                          get_exp_new_inc_carr(integrate = sim_samples_nb_lambda2_age8[[x]],
                                               years = 1:30)
                        },
                        mc.cores = 10)),
    tar_target(imd_sim_nb_lambda2_age8,
               get_exp_new_imd_vacc_b(carr_inc = carr_inc_sim_nb_lambda2_age8,
                                      ccr = ccr_main,
                                      ve_d = ve_d_main)),
    tar_target(imd_sim_samples_nb_lambda2_age8,
               mclapply(1:length(parameter_samples_nb_lambda2_age8),
                        function(x) {
                          get_exp_new_imd_vacc_b(carr_inc = carr_inc_sim_samples_nb_lambda2_age8[[x]],
                                                 ccr = ccr_main,
                                                 ve_d = ve_d_main)
                        },
                        mc.cores = 10))
  ),
  #'  #' Simulation results
  #'  #' Tables
  tar_target(sim_res_lambda2_age8,
             sim_res(imd_sim_scen1 = imd_sim_nb_lambda2_age8_scen1,
                     imd_sim_scen2 = imd_sim_nb_lambda2_age8_scen2,
                     imd_sim_scen3 = imd_sim_nb_lambda2_age8_scen3,
                     imd_sim_scen4 = imd_sim_nb_lambda2_age8_scen4,
                     imd_sim_scen5 = imd_sim_nb_lambda2_age8_scen5,
                     imd_sim_samples_scen1 = imd_sim_samples_nb_lambda2_age8_scen1,
                     imd_sim_samples_scen2 = imd_sim_samples_nb_lambda2_age8_scen2,
                     imd_sim_samples_scen3 = imd_sim_samples_nb_lambda2_age8_scen3,
                     imd_sim_samples_scen4 = imd_sim_samples_nb_lambda2_age8_scen4,
                     imd_sim_samples_scen5 = imd_sim_samples_nb_lambda2_age8_scen5,
                     vacc_num_scen1 = vacc_num_nb_lambda2_age8_scen1,
                     vacc_num_scen2 = vacc_num_nb_lambda2_age8_scen2,
                     vacc_num_scen3 = vacc_num_nb_lambda2_age8_scen3,
                     vacc_num_scen4 = vacc_num_nb_lambda2_age8_scen4,
                     vacc_num_scen5 = vacc_num_nb_lambda2_age8_scen5,
                     vacc_num_samples_scen1 = vacc_num_samples_nb_lambda2_age8_scen1,
                     vacc_num_samples_scen2 = vacc_num_samples_nb_lambda2_age8_scen2,
                     vacc_num_samples_scen3 = vacc_num_samples_nb_lambda2_age8_scen3,
                     vacc_num_samples_scen4 = vacc_num_samples_nb_lambda2_age8_scen4,
                     vacc_num_samples_scen5 = vacc_num_samples_nb_lambda2_age8_scen5,
                     cfr = cfr,
                     seq_prob = seq_prob)),
   tar_target(sim_res_lambda2_age8_30y,
              sim_res(imd_sim_scen1 = imd_sim_nb_lambda2_age8_scen1,
                      imd_sim_scen2 = imd_sim_nb_lambda2_age8_scen2,
                      imd_sim_scen3 = imd_sim_nb_lambda2_age8_scen3,
                      imd_sim_scen4 = imd_sim_nb_lambda2_age8_scen4,
                      imd_sim_scen5 = imd_sim_nb_lambda2_age8_scen5,
                      imd_sim_samples_scen1 = imd_sim_samples_nb_lambda2_age8_scen1,
                      imd_sim_samples_scen2 = imd_sim_samples_nb_lambda2_age8_scen2,
                      imd_sim_samples_scen3 = imd_sim_samples_nb_lambda2_age8_scen3,
                      imd_sim_samples_scen4 = imd_sim_samples_nb_lambda2_age8_scen4,
                      imd_sim_samples_scen5 = imd_sim_samples_nb_lambda2_age8_scen5,
                      vacc_num_scen1 = vacc_num_nb_lambda2_age8_scen1,
                      vacc_num_scen2 = vacc_num_nb_lambda2_age8_scen2,
                      vacc_num_scen3 = vacc_num_nb_lambda2_age8_scen3,
                      vacc_num_scen4 = vacc_num_nb_lambda2_age8_scen4,
                      vacc_num_scen5 = vacc_num_nb_lambda2_age8_scen5,
                      vacc_num_samples_scen1 = vacc_num_samples_nb_lambda2_age8_scen1,
                      vacc_num_samples_scen2 = vacc_num_samples_nb_lambda2_age8_scen2,
                      vacc_num_samples_scen3 = vacc_num_samples_nb_lambda2_age8_scen3,
                      vacc_num_samples_scen4 = vacc_num_samples_nb_lambda2_age8_scen4,
                      vacc_num_samples_scen5 = vacc_num_samples_nb_lambda2_age8_scen5,
                      cfr = cfr,
                      seq_prob = seq_prob,
                      times = 1:30)),
   tar_target(sim_res_lambda2_age8_table,
              sim_res_lambda2_age8$prevented_cases %>%
                select(scen, type, prev = diff, prev_q025 = q025, prev_q975 = q975) %>%
                mutate(prev = format_averted(prev),
                       prev_q025 = format_averted(prev_q025),
                       prev_q975 = format_averted(prev_q975)) %>%
                left_join(sim_res_lambda2_age8$nnv %>%
                            select(scen, type, nnv, nnv_q025, nnv_q975) %>%
                            mutate(nnv = format_nnv(nnv),
                                   nnv_q025 = format_nnv(nnv_q025),
                                   nnv_q975 = format_nnv(nnv_q975)))),
   tar_target(save_sim_res_lambda2_age8_table,
              make_csv(sim_res_lambda2_age8_table, path = "./results/paper/table_1_effec_effic.csv"),
              format = "file"),

   tar_target(sim_res_lambda2_age8_wo_o_table,
              sim_res_lambda2_age8$prevented_cases %>%
                select(scen, type, prev = diff_wo_o, prev_q025 = wo_o_q025, prev_q975 = wo_o_q975) %>%
                mutate(prev = format_averted(prev),
                       prev_q025 = format_averted(prev_q025),
                       prev_q975 = format_averted(prev_q975)) %>%
                left_join(sim_res_lambda2_age8$nnv %>%
                            select(scen, type, nnv_wo_o, nnv_wo_o_q025, nnv_wo_o_q975) %>%
                            mutate(nnv = format_nnv(nnv_wo_o),
                                   nnv_q025 = format_nnv(nnv_wo_o_q025),
                                   nnv_q975 = format_nnv(nnv_wo_o_q975)) %>%
                            select(scen, type, nnv, nnv_q025, nnv_q975))),
   tar_target(save_sim_res_lambda2_age8_wo_o_table,
              make_csv(sim_res_lambda2_age8_wo_o_table, path = "./results/paper/supp_tab_4_effec_effic_wo_o.csv"), format = "file"),

  tar_target(write_sim_res_lambda2_age8_prev_by_sero,
             make_csv(sim_res_lambda2_age8$prevented_cases_by_sero %>%
                        mutate(prevented  = round(prevented, 1),
                               q025 = round(q025, 1),
                               q975 = round(q975, 1)) %>%
                        pivot_wider(id_cols = c(scen, sero),
                                    names_from = type, values_from = c("prevented", "q025", "q975")) %>%
                        select(scen, sero, contains("imd"), contains("sequelae"), contains("fatalities")),
                      path = "./results/paper/supp_tab_3_sim_res_lambda2_age8_tab_prev_by_sero.csv"),
             format = "file"),
  # 30 year simulation period
  tar_target(sim_res_lambda2_age8_table_30y,
             sim_res_lambda2_age8_30y$prevented_cases %>%
               select(scen, type, prev = diff, prev_q025 = q025, prev_q975 = q975) %>%
               mutate(prev = format_averted(prev),
                      prev_q025 = format_averted(prev_q025),
                      prev_q975 = format_averted(prev_q975)) %>%
               left_join(sim_res_lambda2_age8_30y$nnv %>%
                           select(scen, type, nnv, nnv_q025, nnv_q975) %>%
                           mutate(nnv = format_nnv(nnv),
                                  nnv_q025 = format_nnv(nnv_q025),
                                  nnv_q975 = format_nnv(nnv_q975)))),
  tar_target(save_sim_res_lambda2_age8_30y_table,
             make_csv(sim_res_lambda2_age8_table_30y, path = "./results/paper/supp_tab_5_sim_res_30y.csv"),
             format = "file"),
  tar_target(write_sim_res_lambda2_age8_30y_prev_by_sero,
             make_csv(sim_res_lambda2_age8_30y$prevented_cases_by_sero %>%
                        mutate(prevented  = round(prevented, 1),
                               q025 = round(q025, 1),
                               q975 = round(q975, 1)) %>%
                        pivot_wider(id_cols = c(scen, sero),
                                    names_from = type, values_from = c("prevented", "q025", "q975")) %>%
                        select(scen, sero, contains("imd"), contains("sequelae"), contains("fatalities")),
                      path = "./results/paper/supp_tab_6_sim_res_lambda2_age8_tab_prev_by_sero_30y.csv"),
             format = "file"),
  # Plot results
  tar_target(plot_sim_res_lambda2_8age_list,
             plots_sim_scen(imd_sim_nb_lambda2_age8_scen0,
                            imd_sim_nb_lambda2_age8_scen1,
                            imd_sim_nb_lambda2_age8_scen2,
                            imd_sim_nb_lambda2_age8_scen3,
                            imd_sim_nb_lambda2_age8_scen4,
                            imd_sim_nb_lambda2_age8_scen5,
                            imd_sim_samples_nb_lambda2_age8_scen0,
                            imd_sim_samples_nb_lambda2_age8_scen1,
                            imd_sim_samples_nb_lambda2_age8_scen2,
                            imd_sim_samples_nb_lambda2_age8_scen3,
                            imd_sim_samples_nb_lambda2_age8_scen4,
                            imd_sim_samples_nb_lambda2_age8_scen5,
                            times = 1:10)),

  # 30y
  tar_target(plot_sim_res_lambda2_8age_30y_list,
             plots_sim_scen(imd_sim_nb_lambda2_age8_scen0,
                            imd_sim_nb_lambda2_age8_scen1,
                            imd_sim_nb_lambda2_age8_scen2,
                            imd_sim_nb_lambda2_age8_scen3,
                            imd_sim_nb_lambda2_age8_scen4,
                            imd_sim_nb_lambda2_age8_scen5,
                            imd_sim_samples_nb_lambda2_age8_scen0,
                            imd_sim_samples_nb_lambda2_age8_scen1,
                            imd_sim_samples_nb_lambda2_age8_scen2,
                            imd_sim_samples_nb_lambda2_age8_scen3,
                            imd_sim_samples_nb_lambda2_age8_scen4,
                            imd_sim_samples_nb_lambda2_age8_scen5,
                            times = 1:30)),
  tar_target(plot_sim_res_lambda2_8age_30y,
             format_plot_sim_res_30y(plot_sim_res_lambda2_8age_30y_list)),
  tar_target(save_plot_sim_res_lambda2_8age_30y,
             make_png(plot_sim_res_lambda2_8age_30y,
                      path = "./results/paper/supp_fig_17_sim_res_lambda2_8age_30y.png", width = 12*0.75, height = 10*0.75), format = "file"),

  tar_target(plot_sim_res_lambda2_8age,
             format_plot_sim_res_fig_3(plot_sim_res_lambda2_8age_list,
                                       plot_sim_res_lambda2_8age_30y_list)),
  tar_target(save_plot_sim_res_lambda2_8age,
             make_png(plot_sim_res_lambda2_8age,
                      path = "./results/paper/fig_3_sim_res_lambda2_8age.png", width = 9, height = 10), format = "file"),
  tar_target(plot_comp_sim_scen_by,
             plot_comp_scen_by(comp_nb_lambda2_age8,
                               comp_sim_nb_lambda2_age8_scen1,
                               comp_sim_nb_lambda2_age8_scen2,
                               comp_sim_nb_lambda2_age8_scen3,
                               comp_sim_nb_lambda2_age8_scen4,
                               comp_sim_nb_lambda2_age8_scen5,
                               scens = paste0("scen_", 1:3),
                               years = c(2005, 2030),
                               birthyears = c(2020, 2010, 2000))),
  tar_target(plot_comp_sim_res_lambda2_8age,
             ggpubr::ggarrange(plot_comp_sim_scen_by$birthyear_2010 + ggtitle("Birthyear 2010") +
                                 scale_x_continuous(breaks = c(2005, 2015, 2025)) +
                                 scale_y_continuous(n.breaks = 4, labels = scales::percent) +
                                 theme(axis.text=element_text(size=10)),
                               plot_comp_sim_scen_by$birthyear_2000 + ggtitle("Birthyear 2000") +
                                 scale_x_continuous(breaks = c(2005, 2015, 2025)) +
                                 scale_y_continuous(n.breaks = 4, labels = scales::percent) +
                                 theme(axis.text=element_text(size=10)),
                               ncol = 2, common.legend = T, legend = "top", labels = "AUTO")),
  tar_target(save_plot_comp_sim_res_lambda2_8age,
             make_png(plot_comp_sim_res_lambda2_8age,
                      path = "./results/paper/fig_4_sim_res_lambda2_8age_carrier.png", width = 9, height = 4.5),
             format = "file"),
  # Carriage by age and simulation year
  tar_target(fig_carr_age_sim_main,
             plot_carr_age_year_sim(comp_sim_nb_lambda2_age8_scen1,
                                    comp_sim_nb_lambda2_age8_scen2,
                                    comp_sim_nb_lambda2_age8_scen3)),
  tar_target(save_plot_carr_age_sim_main,
             make_png(fig_carr_age_sim_main + theme(axis.text=element_text(size=10)),
                      path = "./results/paper/supp_fig_3_carr_age_sim_main.png", width = 9, height = 6),
             format = "file")
)
