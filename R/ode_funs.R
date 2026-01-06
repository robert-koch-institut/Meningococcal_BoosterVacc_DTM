#' calc_lambda_am_t Function to calculate lambda at time t per age (A different) and type 
#' of meningococcal serogroup (M different)
#' @param prob_inf_m probability of infection given contact to infectious per 
#' sero-group, vector of dim M (here 3, serogroup c, awy, o)
#' @param cont_matrix contact matrix, dimension A x A, now many people from age group a2 (columns)
#' does individual from age group a1 (rows) meet (e.g., per day) on average?
#' @param infec_am number of infectious per age and serogroup at time t, matrix dimension A x M
#' @param scale_cont age-specific scaling factor of infection probability per age-group, length A, 
#' Default NULL, i.e., no scaling
#' @param N_t number of individuals in population at time t, vector of length A
#' @return matrix with current force of infection per age group (A rows) 
#' serogroup (M columns)


calc_lambda_am_t = function(prob_inf_m, 
                            cont_matrix, 
                            infec_am,
                            scale_cont = NULL,
                            N_t){
  A = nrow(cont_matrix)
  M = length(prob_inf_m)
  stopifnot(dim(cont_matrix)==c(A, A))
  stopifnot(dim(infec_am)==c(A, M))
  stopifnot(A==length(N_t))
  if (is.null(scale_cont)) {
    out = matrix(rep(prob_inf_m, A), ncol = M, byrow=T) * 
      cont_matrix %*% (infec_am / N_t)
  } else {
    stopifnot(A==length(scale_cont))
    out = matrix(rep(prob_inf_m, A), ncol = M, byrow=T) * 
      ((matrix(scale_cont, nrow=A, ncol=A) * cont_matrix *
          matrix(scale_cont, nrow=A, ncol=A, byrow=T)) %*% (infec_am / N_t))
  }
  stopifnot(dim(out) == c(A, M))
  out
}


#' calc_lambda_am_t_2 Function to calculate lambda at time t per age (A different) and type 
#' of meningococcal serogroup (M different)
#' @param prob_inf_m age-specific probability of infection given contact to infectious per 
#' sero-group, matrix of dim AxM (86x3, agegroups x serogroup c, awy, o)
#' @param cont_matrix contact matrix, dimension A x A, now many people from age group a2 (columns)
#' does individual from age group a1 (rows) meet (e.g., per day) on average?
#' @param infec_am number of infectious per age and serogroup at time t, matrix dimension A x M
#' @param N_t number of individuals in population at time t, vector of length A
#' @return matrix with current force of infection per age group (A rows) 
#' serogroup (M columns)


calc_lambda_am_t_2 = function(prob_inf_a_m, 
                              cont_matrix, 
                              infec_am,
                              N_t){
  A = nrow(cont_matrix)
  M = dim(prob_inf_a_m)[2]
  stopifnot(dim(cont_matrix)==c(A, A))
  stopifnot(dim(infec_am)==c(A, M))
  stopifnot(A==length(N_t))
  stopifnot(M == dim(prob_inf_a_m)[2])
    out = prob_inf_a_m * 
      cont_matrix %*% (infec_am / N_t)
  stopifnot(dim(out) == c(A, M))
  out
}

#' calc_comp_dt_a Calculate derivative of compartments for age group a at time t
#' @param state current state of 20+15+2 compartments in age group a: 20 model 
#' compartments; 15 helper compartments for new incidence of 3 serotypes in 5 vaccine 
#' groups (1 unvaccinated, 2 vaccine-immune, 2 vaccine-waned). 2 helper compartments
#' for number of administered vaccines (C, ACWY). 
#' @param r_m vector of length 3, recovery rate for three sero-groups (c, awy, o)
#' @param vacc_eff_mmCc vaccine efficacy of mmCc vaccine (protection against
#' carrier). vector of length 3 for three sero-groups (c, awy, o)
#' @param vacc_eff_pmACWYc vaccine efficacy of pmACWYc vaccine. vector of length 3 
#' for three sero-groups (c, awy, o)
#' @param waning_mmCc real, waning rate of mmCc vaccine in age-group a
#' @param waning_pmACWYc real, waning rate of pmACWYc vaccine in age-group a
#' @param lambda_m_at force of infection for m sero-groups in age group a at time-point t.
#' Vector of length 3 (sero-groups c, awy, o)
#' @return vector of derivatives of 20 + 17 helper compartments in age group a at time t 

calc_comp_dt_a = function(state, 
                          r_m, 
                          vacc_eff_mmCc, 
                          vacc_eff_pmACWYc,
                          waning_mmCc,
                          waning_pmACWYc,
                          lambda_m_t){
  state = unname(state)
  lambda_m_t = unname(lambda_m_t)
  s_uv = state[1]
  s_vi_mmCc = state[2]
  s_vi_pmACWYc = state[3]
  s_vw_mmCc = state[4]
  s_vw_pmACWYc = state[5]
  car_uv_c = state[6]
  car_uv_awy = state[7]
  car_uv_o = state[8]
  car_vi_mmCc_c = state[9]
  car_vi_mmCc_awy = state[10]
  car_vi_mmCc_o = state[11]
  car_vi_pmACWYc_c = state[12]
  car_vi_pmACWYc_awy = state[13]
  car_vi_pmACWYc_o = state[14]
  car_vw_mmCc_c = state[15]
  car_vw_mmCc_awy = state[16]
  car_vw_mmCc_o = state[17]
  car_vw_pmACWYc_c = state[18]
  car_vw_pmACWYc_awy = state[19]
  car_vw_pmACWYc_o = state[20]
  # helper compartments carrier incidence
  car_inc_uv_c = state[21]
  car_inc_uv_awy = state[22]
  car_inc_uv_o = state[23]
  car_inc_vi_mmCc_c = state[24]
  car_inc_vi_mmCc_awy = state[25]
  car_inc_vi_mmCc_o = state[26]
  car_inc_vi_pmACWYc_c = state[27]
  car_inc_vi_pmACWYc_awy = state[28]
  car_inc_vi_pmACWYc_o = state[29]
  car_inc_vw_mmCc_c = state[30]
  car_inc_vw_mmCc_awy = state[31]
  car_inc_vw_mmCc_o = state[32]
  car_inc_vw_pmACWYc_c = state[33]
  car_inc_vw_pmACWYc_awy = state[34]
  car_inc_vw_pmACWYc_o = state[35]
  
  
  # Susceptible
  # Unvaccinated
  # "remove new infected, add recovered of all serotypes"
  d_s_uv = - s_uv * sum(lambda_m_t) + sum(r_m * c(car_uv_c, car_uv_awy, car_uv_o))
  # Vaccinated immune
  # mmCc
  # "remove new infected accounting for efficacy, add recovered of all serotypes"
  ds_vi_mmCc = - s_vi_mmCc * (sum(lambda_m_t * (1-vacc_eff_mmCc)) + waning_mmCc) +
    sum(r_m * c(car_vi_mmCc_c, car_vi_mmCc_awy, car_vi_mmCc_o))
  # pmACWY
  ds_vi_pmACWYc = - s_vi_pmACWYc * (sum(lambda_m_t * (1-vacc_eff_pmACWYc)) + waning_pmACWYc) +
    sum(r_m * c(car_vi_pmACWYc_c, car_vi_pmACWYc_awy, car_vi_pmACWYc_o))
  # Vaccinated waned
  # "remove new infected, add wained and recovered of all serotypes"
  # mmCc
  ds_vw_mmCc = - s_vw_mmCc * sum(lambda_m_t) + s_vi_mmCc * waning_mmCc +
    sum(r_m * c(car_vw_mmCc_c, car_vw_mmCc_awy, car_vw_mmCc_o))
  # pmACWY
  ds_vw_pmACWYc = - s_vw_pmACWYc * sum(lambda_m_t) + s_vi_pmACWYc * waning_pmACWYc +
    sum(r_m * c(car_vw_pmACWYc_c, car_vw_pmACWYc_awy, car_vw_pmACWYc_o))
  
  # Carrier
  # Unvaccinated
  # "remove recovered, add newly infected"
  dcar_uv_c = - car_uv_c * r_m[1] + s_uv * lambda_m_t[1]
  dcar_uv_awy = - car_uv_awy * r_m[2] + s_uv * lambda_m_t[2]
  dcar_uv_o = - car_uv_o * r_m[3] + s_uv * lambda_m_t[3]
  
  # Vaccinated immune
  # "remove recovered and wained, add newly infected (accounting for vacc efficacy" 
  # mmCc
  dcar_vi_mmCc_c = - car_vi_mmCc_c * (r_m[1] + waning_mmCc) + 
    s_vi_mmCc * lambda_m_t[1] * (1-vacc_eff_mmCc[1])
  dcar_vi_mmCc_awy = - car_vi_mmCc_awy * (r_m[2] + waning_mmCc) + 
    s_vi_mmCc * lambda_m_t[2] * (1-vacc_eff_mmCc[2])
  dcar_vi_mmCc_o = - car_vi_mmCc_o * (r_m[3] + waning_mmCc) + 
    s_vi_mmCc * lambda_m_t[3] * (1-vacc_eff_mmCc[3])
  # pmACWY
  dcar_vi_pmACWYc_c = - car_vi_pmACWYc_c * (r_m[1] + waning_pmACWYc) + 
    s_vi_pmACWYc * lambda_m_t[1] * (1-vacc_eff_pmACWYc[1])
  dcar_vi_pmACWYc_awy = - car_vi_pmACWYc_awy * (r_m[2] + waning_pmACWYc) + 
    s_vi_pmACWYc * lambda_m_t[2] * (1-vacc_eff_pmACWYc[2])
  dcar_vi_pmACWYc_o = - car_vi_pmACWYc_o * (r_m[3] + waning_pmACWYc) + 
    s_vi_pmACWYc * lambda_m_t[3] * (1-vacc_eff_pmACWYc[3])
  # Vaccinated waned
  # "remove recovered, add newly infected and waned" 
  # mmCc
  dcar_vw_mmCc_c = - car_vw_mmCc_c * r_m[1] + 
    s_vw_mmCc * lambda_m_t[1] + car_vi_mmCc_c * waning_mmCc
  dcar_vw_mmCc_awy = - car_vw_mmCc_awy * r_m[2] + 
    s_vw_mmCc * lambda_m_t[2] + car_vi_mmCc_awy * waning_mmCc
  dcar_vw_mmCc_o = - car_vw_mmCc_o * r_m[3] + 
    s_vw_mmCc * lambda_m_t[3] + car_vi_mmCc_o * waning_mmCc
  
  # pmACWY
  dcar_vw_pmACWYc_c = - car_vw_pmACWYc_c * r_m[1] + 
    s_vw_pmACWYc * lambda_m_t[1] + car_vi_pmACWYc_c * waning_pmACWYc
  dcar_vw_pmACWYc_awy = - car_vw_pmACWYc_awy * r_m[2] + 
    s_vw_pmACWYc * lambda_m_t[2] + car_vi_pmACWYc_awy * waning_pmACWYc
  dcar_vw_pmACWYc_o = - car_vw_pmACWYc_o * r_m[3] + 
    s_vw_pmACWYc * lambda_m_t[3] + car_vi_pmACWYc_o * waning_pmACWYc
  
  # helper compartments
  dcar_inc_uv_c = s_uv * lambda_m_t[1]
  dcar_inc_uv_awy = s_uv * lambda_m_t[2]
  dcar_inc_uv_o = s_uv * lambda_m_t[3]
  
  dcar_inc_vi_mmCc_c = s_vi_mmCc * lambda_m_t[1] * (1-vacc_eff_mmCc[1])
  dcar_inc_vi_mmCc_awy = s_vi_mmCc * lambda_m_t[2] * (1-vacc_eff_mmCc[2])
  dcar_inc_vi_mmCc_o = s_vi_mmCc * lambda_m_t[3] * (1-vacc_eff_mmCc[3])
  dcar_inc_vi_pmACWYc_c = s_vi_pmACWYc * lambda_m_t[1] * (1-vacc_eff_pmACWYc[1])
  dcar_inc_vi_pmACWYc_awy = s_vi_pmACWYc * lambda_m_t[2] * (1-vacc_eff_pmACWYc[2])
  dcar_inc_vi_pmACWYc_o = s_vi_pmACWYc * lambda_m_t[3] * (1-vacc_eff_pmACWYc[3])
  
  dcar_inc_vw_mmCc_c = s_vw_mmCc * lambda_m_t[1]
  dcar_inc_vw_mmCc_awy = s_vw_mmCc * lambda_m_t[2]
  dcar_inc_vw_mmCc_o = s_vw_mmCc * lambda_m_t[3]
  dcar_inc_vw_pmACWYc_c = s_vw_pmACWYc * lambda_m_t[1]
  dcar_inc_vw_pmACWYc_awy = s_vw_pmACWYc * lambda_m_t[2]
  dcar_inc_vw_pmACWYc_o = s_vw_pmACWYc * lambda_m_t[3]
  
  # Output, derivatives of 37 compartments
  c(# susceptible
    # unvaccinated
    d_s_uv=d_s_uv,
    # vaccinated immune
    ds_vi_mmCc=ds_vi_mmCc, 
    ds_vi_pmACWYc=ds_vi_pmACWYc,
    # vaccinated waned
    ds_vw_mmCc=ds_vw_mmCc, 
    ds_vw_pmACWYc = ds_vw_pmACWYc,
    # carrier 
    # unvaccinated
    dcar_uv_c=dcar_uv_c, 
    dcar_uv_awy=dcar_uv_awy, 
    dcar_uv_o=dcar_uv_o, 
    # vaccinated immune
    dcar_vi_mmCc_c=dcar_vi_mmCc_c, 
    dcar_vi_mmCc_awy=dcar_vi_mmCc_awy, 
    dcar_vi_mmCc_o=dcar_vi_mmCc_o, 
    dcar_vi_pmACWYc_c=dcar_vi_pmACWYc_c, 
    dcar_vi_pmACWYc_awy=dcar_vi_pmACWYc_awy, 
    dcar_vi_pmACWYc_o=dcar_vi_pmACWYc_o,
    # vaccinated waned
    dcar_vw_mmCc_c=dcar_vw_mmCc_c, 
    dcar_vw_mmCc_awy=dcar_vw_mmCc_awy, 
    dcar_vw_mmCc_o=dcar_vw_mmCc_o, 
    dcar_vw_pmACWYc_c=dcar_vw_pmACWYc_c, 
    dcar_vw_pmACWYc_awy=dcar_vw_pmACWYc_awy, 
    dcar_vw_pmACWYc_o=dcar_vw_pmACWYc_o,
    # helper compartments, carriage incidence
    dcar_inc_uv_c=dcar_inc_uv_c, 
    dcar_inc_uv_awy=dcar_inc_uv_awy, 
    dcar_inc_uv_o=dcar_inc_uv_o, 
    dcar_inc_vi_mmCc_c=dcar_inc_vi_mmCc_c, 
    dcar_inc_vi_mmCc_awy=dcar_inc_vi_mmCc_awy, 
    dcar_inc_vi_mmCc_o=dcar_inc_vi_mmCc_o, 
    dcar_inc_vi_pmACWYc_c=dcar_inc_vi_pmACWYc_c, 
    dcar_inc_vi_pmACWYc_awy=dcar_inc_vi_pmACWYc_awy, 
    dcar_inc_vi_pmACWYc_o=dcar_inc_vi_pmACWYc_o, 
    dcar_inc_vw_mmCc_c=dcar_inc_vw_mmCc_c, 
    dcar_inc_vw_mmCc_awy=dcar_inc_vw_mmCc_awy, 
    dcar_inc_vw_mmCc_o=dcar_inc_vw_mmCc_o, 
    dcar_inc_vw_pmACWYc_c=dcar_inc_vw_pmACWYc_c, 
    dcar_inc_vw_pmACWYc_awy=dcar_inc_vw_pmACWYc_awy, 
    dcar_inc_vw_pmACWYc_o=dcar_inc_vw_pmACWYc_o,
    # helper comp vacc c, acwy
    dvacc_mmCc = 0,
    dvacc_pmACWYc = 0)
}


#' calc_comp_dt Calculate derivative of compartments multiple age-groups at time t
#' calculates current number of infectious per sero-type and age-group and total
#' number of individuals per age-group, derives age-specific force of infection,
#' and the deriviatives of the compartments per age-group subsequently.
#' @param state current state of (20+17)*A compartments 
#' @param parms named list of (pre-specified) parameters of DTM model. Contains (at least):
#' "r_m": vector of length 3, recovery rate for three sero-groups (c, awy, o)
#' "vacc_eff_mmCc": vaccine efficacy of mmCc vaccine. vector of length 3 
#' for three sero-groups (c, awy, o)
#' "vacc_eff_pmACWYc": vaccine efficacy of pmACWYc vaccine. vector of length 3 
#' for three sero-groups (c, awy, o)
#' "waning_mmCc": vector of length A, waning rate of mmCc vaccine in age-group a
#' "waning_pmACWYc": vector of length A, waning rate of pmACWYc vaccine in age-group a
#' "cont_matrix": contact matrix between age groups, dimension A x A, 
#' how many people from age group a2 (columns)
#' does individual from age group a1 (rows) meet (e.g., per day) on average?
#' @param prob_inf_m infection probability for each sero-type,
#' vector of length 3 (sero-groups c, awy, o)
#' @param age_groups number of considered age groups (=A), integer
#' @return vector of derivatives of (20 + 17)*A compartments at time t 
calc_comp_dt = function(state_full, 
                        parms, 
                        prob_inf_m = c(0.25,.25,.25), 
                        age_groups=2,
                        scale_cont = NULL) {
  state_data = tibble(age_group = rep(1:age_groups, each=37),
                      comp_name = names(state_full),
                      n=state_full) %>%
    mutate(comp = ifelse(stringr::str_detect(comp_name, "s_"), "susc", "carrier"),
           comp = ifelse(stringr::str_detect(comp_name, "car_inc"), "car_inc", comp),
           comp = ifelse(stringr::str_detect(comp_name, "vacc_"), "vacc", comp),
           sero = ifelse(stringr::str_detect(comp_name, "_c"), "c", "awy"),
           sero = ifelse(stringr::str_detect(comp_name, "_o"), "o", sero),
           sero = ifelse(comp=="susc", NA, sero),
           sero = ifelse(stringr::str_detect(comp_name, "vacc_mmCc"), "mmCc", sero),
           sero = ifelse(stringr::str_detect(comp_name, "vacc_pmACWYc"), "pmACWYc", sero))
  infec_am = state_data %>% filter(comp=="carrier") %>%
    group_by(age_group, sero) %>%
    summarise(n_carrier=sum(n)) %>%
    ungroup() %>%
    pivot_wider(id_cols = age_group, names_from=sero, values_from = n_carrier) %>%
    select(-age_group) %>% 
    relocate(c, awy, o) %>% 
    as.matrix()
  N_t = state_data %>% filter(comp!="car_inc", comp != "vacc") %>% group_by(age_group) %>%
    summarise(n=sum(n)) %>% ungroup() %>% pull(n)
  lambda_am_t = calc_lambda_am_t(prob_inf_m = prob_inf_m, 
                                 infec_am = infec_am, 
                                 cont_matrix = parms$cont_matrix,
                                 scale_cont = scale_cont,
                                 N_t = N_t)
  dt_list = map(1:age_groups, function(x) calc_comp_dt_a(state=state_full[(1:37) + 37*(x-1)],
                                                         r_m = parms$r_m,
                                                         vacc_eff_mmCc = parms$vacc_eff_mmCc,
                                                         vacc_eff_pmACWYc = parms$vacc_eff_pmACWYc,
                                                         waning_mmCc = parms$waning_mmCc[x],
                                                         waning_pmACWYc = parms$waning_pmACWYc[x],
                                                         lambda_m_t = lambda_am_t[x,]))
  unlist(dt_list)
  
}


#' calc_comp_dt_2 Calculate derivative of compartments multiple age-groups at time t
#' calculates current number of infectious per sero-type and age-group and total
#' number of individuals per age-group, derives age-specific force of infection,
#' and the deriviatives of the compartments per age-group subsequently.
#' @param state current state of (20+17)*A compartments 
#' @param parms named list of (pre-specified) parameters of DTM model. Contains (at least):
#' "r_m": vector of length 3, recovery rate for three sero-groups (c, awy, o)
#' "vacc_eff_mmCc": vaccine efficacy of mmCc vaccine. vector of length 3 
#' for three sero-groups (c, awy, o)
#' "vacc_eff_pmACWYc": vaccine efficacy of pmACWYc vaccine. vector of length 3 
#' for three sero-groups (c, awy, o)
#' "waning_mmCc": vector of length A, waning rate of mmCc vaccine in age-group a
#' "waning_pmACWYc": vector of length A, waning rate of pmACWYc vaccine in age-group a
#' "cont_matrix": contact matrix between age groups, dimension A x A, 
#' how many people from age group a2 (columns)
#' does individual from age group a1 (rows) meet (e.g., per day) on average?
#' @param prob_inf_a_m age-specific infection probability for each sero-group,
#' matrix with dim AxM (age-group X sero-groups c, awy, o)
#' @param age_groups number of considered age groups (=A), integer
#' @return vector of derivatives of (20 + 17)*A compartments at time t 
calc_comp_dt_2 = function(state_full, 
                          parms, 
                          prob_inf_a_m = matrix(rep(c(0.25,.25,.25), 2),
                                                nrow = 2), 
                          age_groups=2) {
  state_data = tibble(age_group = rep(1:age_groups, each=37),
                      comp_name = names(state_full),
                      n=state_full) %>%
    mutate(comp = ifelse(stringr::str_detect(comp_name, "s_"), "susc", "carrier"),
           comp = ifelse(stringr::str_detect(comp_name, "car_inc"), "car_inc", comp),
           comp = ifelse(stringr::str_detect(comp_name, "vacc_"), "vacc", comp),
           sero = ifelse(stringr::str_detect(comp_name, "_c"), "c", "awy"),
           sero = ifelse(stringr::str_detect(comp_name, "_o"), "o", sero),
           sero = ifelse(comp=="susc", NA, sero),
           sero = ifelse(stringr::str_detect(comp_name, "vacc_mmCc"), "mmCc", sero),
           sero = ifelse(stringr::str_detect(comp_name, "vacc_pmACWYc"), "pmACWYc", sero))
  infec_am = state_data %>% filter(comp=="carrier") %>%
    group_by(age_group, sero) %>%
    summarise(n_carrier=sum(n)) %>%
    ungroup() %>%
    pivot_wider(id_cols = age_group, names_from=sero, values_from = n_carrier) %>%
    select(-age_group) %>% 
    relocate(c, awy, o) %>% 
    as.matrix()
  N_t = state_data %>% filter(comp!="car_inc", comp != "vacc") %>% group_by(age_group) %>%
    summarise(n=sum(n)) %>% ungroup() %>% pull(n)
  lambda_am_t = calc_lambda_am_t_2(prob_inf_a_m = prob_inf_a_m, 
                                   infec_am = infec_am, 
                                   cont_matrix = parms$cont_matrix,
                                   N_t = N_t)
  dt_list = map(1:age_groups, function(x) calc_comp_dt_a(state=state_full[(1:37) + 37*(x-1)],
                                                         r_m = parms$r_m,
                                                         vacc_eff_mmCc = parms$vacc_eff_mmCc,
                                                         vacc_eff_pmACWYc = parms$vacc_eff_pmACWYc,
                                                         waning_mmCc = parms$waning_mmCc[x],
                                                         waning_pmACWYc = parms$waning_pmACWYc[x],
                                                         lambda_m_t = lambda_am_t[x,]))
  unlist(dt_list)
  
}

#' get_exp_new_inc_carr get expected new incident carriers per year
#' sero-group and age from ODE integration

get_exp_new_inc_carr = function(integrate, years=1:9) {
  colnames(integrate)[-1] = paste0(colnames(integrate)[-1],"-", rep(0:85, each=37))
  new_carrier = integrate %>% as_tibble() %>% pivot_longer(cols = -time) %>%
    filter(!stringr::str_detect(name, "vacc_")) %>%
    separate(name, into = c("comp", "age"), sep = "-") %>%
    mutate(age=as.numeric(age),
           value = as.numeric(value),
           time = as.numeric(time)) %>%
    filter(grepl("car_inc", comp)) %>%
    filter(time %in% c((min(years)-1): max(years))) %>% 
    group_by(age, comp) %>% 
    arrange(time) %>% 
    mutate(new_carrier=c(NA, diff(value)),
           comp = str_replace(comp, "car_inc_", "")) %>% 
    filter(time %in% years)
  out = new_carrier %>% 
    mutate(sero_vacc_prot = 
             case_when(comp %in% c("uv_c", "vw_mmCc_c", "vw_pmACWYc_c") ~ "c_no",
                       comp %in% c("uv_awy", "vw_mmCc_awy", "vw_pmACWYc_awy") ~ "awy_no",
                       comp %in% c("uv_o", "vw_mmCc_o", "vw_pmACWYc_o") ~ "o_no",
                       comp %in% c("vi_mmCc_c") ~ "c_mmCc",
                       comp %in% c("vi_mmCc_awy") ~ "awy_mmCc",
                       comp %in% c("vi_mmCc_o") ~ "o_mmCc",                       
                       comp %in% c("vi_pmACWYc_c") ~ "c_pmACWYc",
                       comp %in% c("vi_pmACWYc_awy") ~ "awy_pmACWYc",
                       comp %in% c("vi_pmACWYc_o") ~ "o_pmACWYc")) %>%
    group_by(sero_vacc_prot, time, age) %>%
    summarise(new_carrier = sum(new_carrier)) %>%
    ungroup() %>%
    separate(sero_vacc_prot, sep = "_", into = c("sero", "vacc"))
  
  out
}

#' get_exp_new_imd get expected new invasive meningococcal cases per year
#' and sero-group from expected new carriers, case-carrier ratio and assumption on 
#' vaccine effectiveness against disease (given new carrier) accounting for MenB vaccination
#' @param carr_inc incident new carriers
#' @param ccr age- and serogroup-specific case carrier ratio
#' @param ve_d tibble of vaccine effectivenes against IMD from specific serogroups (c, awy, other)
#' and vaccines (mmCc, pmACWYc, no)
#' @param ve0_b vaccine effectiveness against MenB-IMD
#' @param waning_b inverse mean duration of protection for MenB vaccine in years (waning rate)
#' @param uptake_b uptake_b assumed vaccine uptake for the MenB vaccine among infants


get_exp_new_imd_vacc_b = function(carr_inc, 
                                  ccr, 
                                  ve_d = tibble(vacc = c("mmCc", "mmCc", "mmCc",
                                                         "no", "no", "no",
                                                         "pmACWYc", "pmACWYc", "pmACWYc"),
                                                sero = rep(c("c", "awy", "o"), 3),
                                                ve_d = c(0.8, 0, 0,
                                                         0, 0, 0,
                                                         0.5, 0.5, 0)),
                                  ve0_b = 0.75,
                                  waning_b = 1/10,
                                  uptake_b = 0.8) {
  prot_vac_b = function(ve0 = ve0_b, 
                       waning_rate = waning_b,
                       uptake = uptake_b) {
    int = function(r, ve, a, b) {
      ((ve*(exp(-a*r) - exp(-b*r)))/r)/(b-a)
    } 
    tibble(expand_grid(time=1:30, age = 0:85, sero = c("c", "awy", "o"))) %>%
      mutate(uptake = ifelse(age>0, uptake, 0),
             uptake = ifelse((time-age)<1, 0, uptake),
             a = pmax(age-1, 0),
             b = a+1,
             prot_b = int(r = waning_rate,
                        ve=ve0,
                        a = a,
                        b = b)*uptake*(sero=="o")) %>%
      select(time, age, sero, prot_b)
  }
  prot_b = prot_vac_b()
  dat = carr_inc %>% left_join(ccr) %>% left_join(ve_d) %>% left_join(prot_b)
  dat %>% mutate(exp_inv = new_carrier*ccr*(1-ve_d)*(1-prot_b)) %>%
    group_by(time, sero, age) %>%
    summarise(exp_inv = sum(exp_inv))
} 

#' rootfun indicator for step-changes in the ODE system
rootfunc <- function(t, population, params) {
  # The events in the stepChanges function are triggered every full year
  # i.e. when the time is an integer
  if(t==0) {
    return(1)
  } else {
    return(t %% 1)
  }
}


#' get_comp_time get compartment counts per vaccine group over time
#' sero-group and age from ODE integration

get_comp_time = function(integrate, begin_year = 2005) {
  colnames(integrate)[-1] = paste0(colnames(integrate)[-1],"-", rep(0:85, each=37))
  comp_num = integrate %>% as_tibble() %>% pivot_longer(cols = -time) %>%
    filter(!stringr::str_detect(name, "vacc_")) %>%
    separate(name, into = c("comp", "age"), sep = "-") %>%
    separate(comp, into = c("one", "two", "three", "four", "five"), sep = "_") %>%
    mutate(comp = one,
           comp = ifelse(one=="car" & two=="uv", paste0(one, "_", three), comp),
           comp = ifelse(one=="car" & two=="vi", paste0(one, "_", four), comp),
           comp = ifelse(one=="car" & two=="vw", paste0(one, "_", four), comp),
           vacc = ifelse(two == "uv", two, paste0(two, "_", three))) %>%
    filter(two!="inc") %>%
    mutate(age=as.numeric(age),
           value = as.numeric(value),
           time = as.numeric(time),
           year = begin_year + time,
           vacc_prot = ifelse(grepl("vi", vacc), "vi", "uv"),
           vacc_prot = ifelse(grepl("vw", vacc), "vw", vacc_prot),
           birthyear = floor(begin_year - age + time))
  comp_num
}


#' step_fun_a_sim Perform yearly step-changes for age-group a
#' @param state_a_t current state of (20+17) compartments in age group a "before ageing" 
#' @param stat_amin1_t current state of (20+17) compartments in age group a-1 "before ageing"
#' @param n_a_t population size in age group a after "ageing"
#' @param n_vac_amin1_t How many individuals have been vaccinated in age group a-1 during last year
#' per vaccine type
#' @param frac_boost_amin1_t What fraction of prev. vaccinated have been boostered in age-group
#' a-1
#' @return state after ageing
step_fun_a_sim = function(state_a_t,
                      state_amin1_t, 
                      n_a_t, 
                      n_vac_amin1_t=c(0, 0), 
                      frac_boost_amin1_t=c(0, 0),
                      boost_count_immune = TRUE) {
  
  # Vaccinate in age-group a-1
  # Primary
  state_amin_1_vacc = state_amin1_t[1:20]
  
  # Derive newly boostered from susceptible/carriers in waned compartments
  n_boost_vi_mmCc = state_amin_1_vacc[c("s_vi_mmCc", "s_vi_pmACWYc", 
                                        "car_vi_mmCc_c", "car_vi_mmCc_awy", "car_vi_mmCc_o",
                                        "car_vi_pmACWYc_c", "car_vi_pmACWYc_awy", "car_vi_pmACWYc_o")] *
    frac_boost_amin1_t[1]
  n_boost_vw_mmCc = state_amin_1_vacc[c("s_vw_mmCc", "s_vw_pmACWYc", 
                                        "car_vw_mmCc_c", "car_vw_mmCc_awy", "car_vw_mmCc_o",
                                        "car_vw_pmACWYc_c", "car_vw_pmACWYc_awy", "car_vw_pmACWYc_o")] *
    frac_boost_amin1_t[1]
  
  n_boost_vi_pmACWYc = state_amin_1_vacc[c("s_vi_mmCc", "s_vi_pmACWYc", 
                                           "car_vi_mmCc_c", "car_vi_mmCc_awy", "car_vi_mmCc_o",
                                           "car_vi_pmACWYc_c", "car_vi_pmACWYc_awy", "car_vi_pmACWYc_o")] *
    frac_boost_amin1_t[2]
  n_boost_vw_pmACWYc = state_amin_1_vacc[c("s_vw_mmCc", "s_vw_pmACWYc", 
                                           "car_vw_mmCc_c", "car_vw_mmCc_awy", "car_vw_mmCc_o",
                                           "car_vw_pmACWYc_c", "car_vw_pmACWYc_awy", "car_vw_pmACWYc_o")] *
    frac_boost_amin1_t[2]
  
  
  # Removed boostered from waned compartments
  state_amin_1_vacc[c("s_vw_mmCc", "s_vw_pmACWYc", 
                      "car_vw_mmCc_c", "car_vw_mmCc_awy", "car_vw_mmCc_o",
                      "car_vw_pmACWYc_c", "car_vw_pmACWYc_awy", "car_vw_pmACWYc_o")] = 
    state_amin_1_vacc[c("s_vw_mmCc", "s_vw_pmACWYc", 
                        "car_vw_mmCc_c", "car_vw_mmCc_awy", "car_vw_mmCc_o",
                        "car_vw_pmACWYc_c", "car_vw_pmACWYc_awy", "car_vw_pmACWYc_o")] -
    n_boost_vw_mmCc - n_boost_vw_pmACWYc
  
  # Add boostered to VI compartments
  state_amin_1_vacc[c("s_vi_mmCc", "car_vi_mmCc_c", "car_vi_mmCc_awy", "car_vi_mmCc_o")] =
    state_amin_1_vacc[c("s_vi_mmCc", "car_vi_mmCc_c", "car_vi_mmCc_awy", "car_vi_mmCc_o")] +
    c(sum(n_boost_vw_mmCc[c("s_vw_mmCc", "s_vw_pmACWYc")]),
      sum(n_boost_vw_mmCc[c("car_vw_mmCc_c", "car_vw_pmACWYc_c")]),
      sum(n_boost_vw_mmCc[c("car_vw_mmCc_awy", "car_vw_pmACWYc_awy")]),
      sum(n_boost_vw_mmCc[c("car_vw_mmCc_o", "car_vw_pmACWYc_o")])) -
    n_boost_vi_pmACWYc[c("s_vi_mmCc", "car_vi_mmCc_c", "car_vi_mmCc_awy", "car_vi_mmCc_o")] +
    n_boost_vi_mmCc[c("s_vi_pmACWYc", "car_vi_pmACWYc_c", "car_vi_pmACWYc_awy", "car_vi_pmACWYc_o")]
  
  
  state_amin_1_vacc[c("s_vi_pmACWYc", "car_vi_pmACWYc_c", "car_vi_pmACWYc_awy", "car_vi_pmACWYc_o")] =
    state_amin_1_vacc[c("s_vi_pmACWYc", "car_vi_pmACWYc_c", "car_vi_pmACWYc_awy", "car_vi_pmACWYc_o")] +
    c(sum(n_boost_vw_pmACWYc[c("s_vw_mmCc", "s_vw_pmACWYc")]),
      sum(n_boost_vw_pmACWYc[c("car_vw_mmCc_c", "car_vw_pmACWYc_c")]),
      sum(n_boost_vw_pmACWYc[c("car_vw_mmCc_awy", "car_vw_pmACWYc_awy")]),
      sum(n_boost_vw_pmACWYc[c("car_vw_mmCc_o", "car_vw_pmACWYc_o")])) -
    n_boost_vi_mmCc[c("s_vi_pmACWYc", "car_vi_pmACWYc_c", "car_vi_pmACWYc_awy", "car_vi_pmACWYc_o")] +
    n_boost_vi_pmACWYc[c("s_vi_mmCc", "car_vi_mmCc_c", "car_vi_mmCc_awy", "car_vi_mmCc_o")]
  
  # Primary
  # mmCc
  # Remove newly vaccinated from susceptible/carriers in UV compartments
  state_amin_1_vacc[c("s_uv", "car_uv_c", "car_uv_awy", "car_uv_o")] =
    state_amin_1_vacc[c("s_uv", "car_uv_c", "car_uv_awy", "car_uv_o")] - 
    (n_vac_amin1_t[1]*state_amin1_t[c("s_uv", "car_uv_c", "car_uv_awy", "car_uv_o")]/
       sum(state_amin1_t[c("s_uv", "car_uv_c", "car_uv_awy", "car_uv_o")]))
  
  # Add newly vaccinated to vaccinated immune compartment
  state_amin_1_vacc[c("s_vi_mmCc", "car_vi_mmCc_c", "car_vi_mmCc_awy", "car_vi_mmCc_o")] =
    state_amin_1_vacc[c("s_vi_mmCc", "car_vi_mmCc_c", "car_vi_mmCc_awy", "car_vi_mmCc_o")] + 
    (n_vac_amin1_t[1]*state_amin1_t[c("s_uv", "car_uv_c", "car_uv_awy", "car_uv_o")]/
       sum(state_amin1_t[c("s_uv", "car_uv_c", "car_uv_awy", "car_uv_o")]))
  # pmACWYc
  # Remove newly vaccinated from susceptible/carriers in UV compartments
  state_amin_1_vacc[c("s_uv", "car_uv_c", "car_uv_awy", "car_uv_o")] =
    state_amin_1_vacc[c("s_uv", "car_uv_c", "car_uv_awy", "car_uv_o")] - 
    (n_vac_amin1_t[2]*state_amin1_t[c("s_uv", "car_uv_c", "car_uv_awy", "car_uv_o")]/
       sum(state_amin1_t[c("s_uv", "car_uv_c", "car_uv_awy", "car_uv_o")])) 
  
  # Add newly vaccinated to vaccinated immune compartment
  state_amin_1_vacc[c("s_vi_pmACWYc", "car_vi_pmACWYc_c", "car_vi_pmACWYc_awy", "car_vi_pmACWYc_o")] =
    state_amin_1_vacc[c("s_vi_pmACWYc", "car_vi_pmACWYc_c", "car_vi_pmACWYc_awy", "car_vi_pmACWYc_o")] + 
    (n_vac_amin1_t[2]*state_amin1_t[c("s_uv", "car_uv_c", "car_uv_awy", "car_uv_o")]/
       sum(state_amin1_t[c("s_uv", "car_uv_c", "car_uv_awy", "car_uv_o")]))
  
  # ageing and adjust pop size
  n_boost_mmCc = ifelse(boost_count_immune, 
                        sum(n_boost_vi_mmCc) + sum(n_boost_vw_mmCc), 
                        sum(n_boost_vw_mmCc))
  
  n_boost_pmACWYc = ifelse(boost_count_immune, 
                           sum(n_boost_vi_pmACWYc) + sum(n_boost_vw_pmACWYc), 
                           sum(n_boost_vw_pmACWYc))
  
  
  state_a_new = c(state_amin_1_vacc/sum(state_amin_1_vacc)*n_a_t, 
                  state_a_t[21:35],
                  state_a_t[36] + n_vac_amin1_t[1] + n_boost_mmCc,
                  state_a_t[37] + n_vac_amin1_t[2] + n_boost_pmACWYc)
  names(state_a_new) = names(state_a_t)
  state_a_new
}

#' stepChanges_sim Perform stepChanges in the comparment counts depending on vaccination scenario
#' @param  description t time
#' @param y compartment vector
#' @param parms list of parameters of the DTM from external data
#' @param scen vaccination scenario

stepChanges_sim <- function(t, y, parms, scen = 0) {
  
  # print(paste("t = ", t))
  if(t>0) {
    vacc = array(0, dim = c(86, 30, 2))
    boost = array(0, dim = c(86, 30, 2))
    boost_count_array = array(0, dim = c(86, 30, 1))
    if (scen == 1) {
      vacc = array(c(rep(c(0, 0.8, rep(0, 84)), 30),
                     rep(c(rep(0, 86)), 30)),
                   dim = c(86, 30, 2))
      boost = array(0, dim = c(86, 30, 2))
    }
    if (scen == 2) {
      vacc = array(c(rep(c(0, 0.8, rep(0, 84)), 30),
                     rep(c(rep(0, 86)), 30)),
                   dim = c(86, 30, 2))
      boost = array(c(rep(c(rep(0, 12), 1, 0,  rep(0, 72)), 30),
                      rep(c(rep(0, 86)), 30)),
                    dim = c(86, 30, 2))
      boost_count_array = array(c(rep(c(rep(0, 12), 1, 0, rep(0, 72)), 30)),
                                dim = c(86, 30, 1))
    }
    if (scen == 3) {
      vacc = array(c(rep(c(0, 0.8, rep(0, 84)), 30),
                     rep(c(rep(0, 86)), 30)),
                   dim = c(86, 30, 2))
      boost = array(c(rep(c(rep(0, 86)), 30),
                      rep(c(rep(0, 12), 1, 0,  rep(0, 72)), 30)),
                    dim = c(86, 30, 2))
      boost_count_array = array(c(rep(c(rep(0, 12), 1, 0, rep(0, 72)), 30)),
                                dim = c(86, 30, 1))
    }
    if (scen == 4) {
      vacc = array(c(rep(c(0, 0.8, rep(0, 84)), 30),
                     rep(c(rep(0, 86)), 30)),
                   dim = c(86, 30, 2))
      boost = array(c(rep(c(rep(0, 15), 1, 0,  rep(0, 69)), 30),
                      rep(c(rep(0, 86)), 30)),
                    dim = c(86, 30, 2))
      boost_count_array = array(c(rep(c(rep(0, 15), 1, 0, rep(0, 69)), 30)),
                                dim = c(86, 30, 1))
    }
    if (scen == 5) {
      vacc = array(c(rep(c(0, 0.8, rep(0, 84)), 30),
                     rep(c(rep(0, 86)), 30)),
                   dim = c(86, 30, 2))
      boost = array(c(rep(c(rep(0, 86)), 30),
                      rep(c(rep(0, 15), 1, 0,  rep(0, 69)), 30)),
                    dim = c(86, 30, 2))
      boost_count_array = array(c(rep(c(rep(0, 15), 1, 0, rep(0, 69)), 30)),
                                dim = c(86, 30, 1))
    }
    
    pop_0 = c(parms$pop_size[[1]][t], rep(0, 19), y[21:37])
    names(pop_0) = names(y)[1:37]
    pop_1_84 = purrr::map(1:84, function(x) step_fun_a_sim(state_a_t = y[1:37 + x*37],
                                                           state_amin1_t = y[1:37 + (x-1)*37],
                                                           n_a_t = parms$pop_size[[x+1]][t], 
                                                           n_vac_amin1_t = vacc[x, t,] * 
                                                             parms$pop_size[[x]][t], 
                                                           frac_boost_amin1_t= boost[x,t,],
                                                           boost_count_immune = boost_count_array[x,t,])) %>%
      unlist()
    pop_85 = step_fun_a_sim(state_a_t = y[1:37 + 85*37],
                            state_amin1_t = y[1:37 + (85-1)*37],
                            n_a_t = parms$pop_size[[85+1]][t], 
                            n_vac_amin1_t = vacc[85, t,] * parms$pop_size[[85]][t] +
                              vacc[86, t,] * parms$pop_size[[86]][t], 
                            frac_boost_amin1_t=c(0, 0)) %>%
      unlist()
    pop_out = c(pop_0, pop_1_84, pop_85)
  } else {
    pop_out=y
  }
  return(pop_out)
}

#' get_vacc_time get compartment counts per vaccine group over time
#' sero-group and age from ODE integration
#' 
#' Result corrresponds to the number of newly vaccinated individuals
#' at the step change of the current year to the next year (e.g., 2020-2021) in age-group
#' age - 1 (i.e., aging and vaccination at the same time).
#' Interpretation: a number of 10.000 mmCc vaccinations in age = 2 for year = 2020
#' approximates 10.000 mmCc vaccinations administered to children of age 
#' month 12-23 (i.e., for "one-year-old") between Jan 1, 2020, and Dec 31, 2020.
#' 
get_vacc_time = function(integrate, begin_year = 2020) {
  colnames(integrate)[-1] = paste0(colnames(integrate)[-1],"-", rep(0:85, each=37))
  vacc_num = integrate %>% as_tibble() %>% pivot_longer(cols = -time) %>%
    filter(stringr::str_detect(name, "vacc_")) %>%
    separate(name, into = c("vacc", "age"), sep = "-") %>%
    mutate(vacc = str_replace(vacc, pattern = "vacc_", replacement = ""),
           age = as.numeric(age),
           value = as.numeric(value),
           time = as.numeric(time)) %>%
    filter((time %% 1) == 0) %>%
    group_by(vacc, age) %>%
    arrange(time) %>%
    mutate(value = c(NA, diff(value, 1))) %>%
    ungroup() %>% filter(time>0) %>%
    mutate(time = as.numeric(time),
           year = begin_year + time - 1)
  vacc_num
}