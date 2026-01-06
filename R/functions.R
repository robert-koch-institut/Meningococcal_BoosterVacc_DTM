#' read_destatis Function to read the destatis population data for specific years
#' @param filepath
#' @param from start year
#' @param to end year
#' @return list of length 86 with population size in year from-to for age-group
#' 0-(85+) years (86 age-groups)
read_destatis = function(filepath = "../data/12411-0005.csv", from=2010, to=2019){
  dat = read_delim(filepath, skip = 7, delim = ";", n_max = 86, col_select = -1, 
                   col_names = paste0("y_",2000:2022)) %>%
    mutate(age=0:85) %>% pivot_longer(cols=-age, names_to="year", 
                                      names_prefix = "y_", 
                                      names_transform = list(year=as.numeric),
                                      values_to = "n") %>%
    arrange(year, age) %>%
    filter(year>=from, year<=to)
  map(0:85, function(x) dat %>% filter(age==x) %>% arrange(year) %>% pull(n))
}

#' read_destatis Function to read the destatis population projection from 2020-2049 
#' for simulation period
read_destatis_20_49 = function(pop_dat_destatis_file) {
  dat = read_delim("./data/12421-0002_$F.csv", skip = 8, delim = ";", n_max = 315-9, col_select = -c(1:2))
  colnames(dat) = c("sex", "age", paste0("y_",2023:2071))
  pop_23_49 = dat %>% filter(sex == "Insgesamt") %>%
    filter(age != "Insgesamt") %>%
    select(-sex, -age) %>%
    mutate(age=0:100,
           age_group = cut(age, breaks = c(0:85, 101), right=FALSE, labels = c(0:85))) %>%
    select(-age) %>%
    rename(age=age_group) %>%
    pivot_longer(cols=-age, names_to="year",
                 names_prefix = "y_",
                 names_transform = list(year=as.numeric),
                 values_to = "n") %>%
    group_by(age, year) %>%
    summarise(n=sum(n)) %>%
    mutate(age=as.numeric(as.character(age))) %>%
    arrange(year, age) %>%
    filter(year>=2023, year<=2049) %>%
    mutate(n=n*100)

  pop_20_22 = read_destatis(pop_dat_destatis_file, from=2020, to=2022)
  pop_20_49 = tibble(age = rep(0:85, each=3),
                     year = rep(2020:2022, 86),
                     n=pop_20_22 %>% unlist()) %>%
    rbind(pop_23_49)

  dat = pop_20_49 %>% arrange(year, age)
  map(0:85, function(x) dat %>% filter(age==x) %>% arrange(year) %>% pull(n))
}

#' make_cont_mat Function to create a contact matrix one year age-groups from
#' 0-85+
#' @return matrix of dimension 86x86 with the daily number of contacts of an individual
#' with age i (row i) to individuals of age j (column j)
make_cont_mat = function(filepath_destatis){
  polymod_mat = hhh4contacts::contactmatrix(grouping = rep(1, 15))
  poly_long = tibble(ind_group = rep(rownames(polymod_mat), 15),
                     cont_group = rep(colnames(polymod_mat), each=15),
                     num_cont = polymod_mat %>% as.vector())
  pop_2010 = tibble(cont=0:85, 
                    n=read_destatis(filepath = filepath_destatis, 
                                    from=2010, 
                                    to=2010) %>% unlist()) %>%
    mutate(share_cont = n/sum(n),
           cont_group = cut(cont, breaks = c(seq(0,70, by=5), 86),
                            right=FALSE),
           cont_group = factor(cont_group, labels = rownames(polymod_mat)))
  pop_2010_group = pop_2010 %>% group_by(cont_group) %>%
    summarise(n=sum(n), share_cont = sum(share_cont))
  poly_long = poly_long %>% left_join(pop_2010_group %>% select(cont_group, share_cont))%>%
    mutate(rate=num_cont/share_cont)
  
  cont_long = expand_grid(ind=0:85,
                          cont=0:85) %>% 
    tibble() %>%
    mutate(ind_group = cut(ind, breaks = c(seq(0,70, by=5), 86),
                           right=FALSE, labels = rownames(polymod_mat)), 
           cont_group = cut(cont, breaks = c(seq(0,70, by=5), 86),
                            right=FALSE, labels = rownames(polymod_mat))
    )
  cont_long = cont_long %>% 
    left_join(pop_2010 %>% select(cont, share_cont)) %>%
    left_join(poly_long %>% select(ind_group, cont_group, rate))
  
  cont_long = cont_long %>%
    mutate(num_cont = rate*share_cont) %>%
    arrange(ind, cont)
  cont_mat = matrix(cont_long$num_cont, ncol=86, nrow=86, byrow=T, 
                    dimnames = list(ind=c(0:84, "85+"),
                                    cont=c(0:84, "85+")))
  cont_mat
}

#' sim_res Summarise simulation results in terms of prevented cases (IMD, sequelae, and death) 
#' and corresponding NNVs in total and by serogroup for each adolescent booster scenario (scen 2-5) compared to 
#' scenario 1 (no booster).
#' @param imd_sim_scen*1-5* Expected IMD counts from simulation for each scenario 
#' @param imd_sim_samples_scen*1-5* List of expected IMD counts from simulation 
#' for each sample of parameters governing simulation dynamics
#' @param vacc_num_scen*1-5* Vaccination counts
#' @param vacc_num_samples_scen*1-5* List of vaccination counts  
#' for each sample of parameters governing simulation dynamics from simulation 
#' @param cfr data.frame with assumed case fatality rate by age and serogroup
#' @param seq_prob probability of experiencing at least one sequelae among incident
#' IMD cases
#' @param times vector of years considered for summarizing simulation results (e.g.,
#' 1:10 corresponds to the first ten simulation years -> 2010-2019)
sim_res = function(imd_sim_scen1 = imd_sim_lambda2_age8_scen1,
                   imd_sim_scen2 = imd_sim_lambda2_age8_scen2,
                   imd_sim_scen3 = imd_sim_lambda2_age8_scen3,
                   imd_sim_scen4 = imd_sim_lambda2_age8_scen4,
                   imd_sim_scen5 = imd_sim_lambda2_age8_scen5,
                   imd_sim_samples_scen1 = imd_sim_samples_lambda2_age8_scen1,
                   imd_sim_samples_scen2 = imd_sim_samples_lambda2_age8_scen2,
                   imd_sim_samples_scen3 = imd_sim_samples_lambda2_age8_scen3,
                   imd_sim_samples_scen4 = imd_sim_samples_lambda2_age8_scen4,
                   imd_sim_samples_scen5 = imd_sim_samples_lambda2_age8_scen5,
                   vacc_num_scen1 = vacc_num_lambda2_age8_scen1,
                   vacc_num_scen2 = vacc_num_lambda2_age8_scen2,
                   vacc_num_scen3 = vacc_num_lambda2_age8_scen3,
                   vacc_num_scen4 = vacc_num_lambda2_age8_scen4,
                   vacc_num_scen5 = vacc_num_lambda2_age8_scen5,
                   vacc_num_samples_scen1 = vacc_num_samples_lambda2_age8_scen1,
                   vacc_num_samples_scen2 = vacc_num_samples_lambda2_age8_scen2,
                   vacc_num_samples_scen3 = vacc_num_samples_lambda2_age8_scen3,
                   vacc_num_samples_scen4 = vacc_num_samples_lambda2_age8_scen4,
                   vacc_num_samples_scen5 = vacc_num_samples_lambda2_age8_scen5,
                   cfr = cfr,
                   seq_prob = seq_prob,
                   times = 1:10) {
  
  est_prev = function(imd_sim_scen1,
                      imd_sim_scen2,
                      imd_sim_scen3,
                      imd_sim_scen4,
                      imd_sim_scen5,
                      cfr,
                      seq_prob,
                      times = 1:30) {
    prevented_imd = imd_sim_scen2 %>% mutate(scen="scen_2") %>%
      rbind(imd_sim_scen3 %>% mutate(scen="scen_3")) %>%
      rbind(imd_sim_scen4 %>% mutate(scen="scen_4")) %>%
      rbind(imd_sim_scen5 %>% mutate(scen="scen_5")) %>%
      filter(time %in% times) %>%
      mutate(exp = exp_inv) %>%
      left_join(imd_sim_scen1 %>%
                  rename(exp_1 = exp_inv)) %>%
      mutate(diff = exp_1 - exp) %>%
      mutate(type = "imd")
    
    prevented_fat = imd_sim_scen2 %>% mutate(scen="scen_2") %>%
      rbind(imd_sim_scen3 %>% mutate(scen="scen_3")) %>%
      rbind(imd_sim_scen4 %>% mutate(scen="scen_4")) %>%
      rbind(imd_sim_scen5 %>% mutate(scen="scen_5")) %>%
      filter(time %in% times) %>%
      mutate(exp = exp_inv) %>%
      left_join(imd_sim_scen1 %>%
                  rename(exp_1 = exp_inv)) %>%
      left_join(cfr %>% mutate(sero = factor(sero, levels = c("C", "AWY", "O"),
                                             labels = c("c", "awy", "o")))) %>%
      mutate(exp = exp*cfr,
             exp_1 = exp_1*cfr) %>%
      mutate(diff = exp_1 - exp) %>%
      mutate(type = "fatalities") %>%
      select(-cfr)
    
    prevented_seq = imd_sim_scen2 %>% mutate(scen="scen_2") %>%
      rbind(imd_sim_scen3 %>% mutate(scen="scen_3")) %>%
      rbind(imd_sim_scen4 %>% mutate(scen="scen_4")) %>%
      rbind(imd_sim_scen5 %>% mutate(scen="scen_5")) %>%
      filter(time %in% times) %>%
      mutate(exp = exp_inv) %>%
      left_join(imd_sim_scen1 %>%
                  rename(exp_1 = exp_inv)) %>%
      left_join(cfr %>% mutate(sero = factor(sero, levels = c("C", "AWY", "O"),
                                             labels = c("c", "awy", "o")))) %>%
      mutate(exp = exp*(cfr + (1-cfr)*seq_prob),
             exp_1 = exp_1*(cfr + (1-cfr)*seq_prob)) %>%
      mutate(diff = exp_1 - exp) %>%
      mutate(type = "sequelae") %>%
      select(-cfr)
    
    
    prevented_ov = prevented_imd %>% 
      rbind(prevented_seq) %>%
      rbind(prevented_fat) %>%
      group_by(scen, type) %>%
      summarise(exp = sum(exp),
                exp_1 = sum(exp_1),
                diff = sum(diff)) %>%
      mutate(type = factor(type, levels = c("imd", "sequelae", "fatalities"))) %>%
      arrange(scen, type) %>%
      left_join(prevented_imd %>% 
                  rbind(prevented_seq) %>%
                  rbind(prevented_fat) %>%
                  filter(sero!="o") %>%
                  group_by(scen, type) %>%
                  summarise(exp_wo_o = sum(exp),
                            exp_1_wo_o = sum(exp_1),
                            diff_wo_o = sum(diff)) %>%
                  mutate(type = factor(type, levels = c("imd", "sequelae", "fatalities"))) %>%
                  arrange(scen, type))
    prevented_ov
  }
  est_prev_by_sero = function(imd_sim_scen1,
                              imd_sim_scen2,
                              imd_sim_scen3,
                              imd_sim_scen4,
                              imd_sim_scen5,
                              cfr,
                              seq_prob,
                              times = 1:30) {
    prevented_imd = imd_sim_scen2 %>% mutate(scen="scen_2") %>%
      rbind(imd_sim_scen3 %>% mutate(scen="scen_3")) %>%
      rbind(imd_sim_scen4 %>% mutate(scen="scen_4")) %>%
      rbind(imd_sim_scen5 %>% mutate(scen="scen_5")) %>%
      filter(time %in% times) %>%
      mutate(exp = exp_inv) %>%
      left_join(imd_sim_scen1 %>%
                  rename(exp_1 = exp_inv)) %>%
      mutate(diff = exp_1 - exp) %>%
      mutate(type = "imd")
    
    prevented_fat = imd_sim_scen2 %>% mutate(scen="scen_2") %>%
      rbind(imd_sim_scen3 %>% mutate(scen="scen_3")) %>%
      rbind(imd_sim_scen4 %>% mutate(scen="scen_4")) %>%
      rbind(imd_sim_scen5 %>% mutate(scen="scen_5")) %>%
      filter(time %in% times) %>%
      mutate(exp = exp_inv) %>%
      left_join(imd_sim_scen1 %>%
                  rename(exp_1 = exp_inv)) %>%
      left_join(cfr %>% mutate(sero = factor(sero, levels = c("C", "AWY", "O"),
                                             labels = c("c", "awy", "o")))) %>%
      mutate(exp = exp*cfr,
             exp_1 = exp_1*cfr) %>%
      mutate(diff = exp_1 - exp) %>%
      mutate(type = "fatalities") %>%
      select(-cfr)
    
    prevented_seq = imd_sim_scen2 %>% mutate(scen="scen_2") %>%
      rbind(imd_sim_scen3 %>% mutate(scen="scen_3")) %>%
      rbind(imd_sim_scen4 %>% mutate(scen="scen_4")) %>%
      rbind(imd_sim_scen5 %>% mutate(scen="scen_5")) %>%
      filter(time %in% times) %>%
      mutate(exp = exp_inv) %>%
      left_join(imd_sim_scen1 %>%
                  rename(exp_1 = exp_inv)) %>%
      left_join(cfr %>% mutate(sero = factor(sero, levels = c("C", "AWY", "O"),
                                             labels = c("c", "awy", "o")))) %>%
      mutate(exp = exp*(cfr + (1-cfr)*seq_prob),
             exp_1 = exp_1*(cfr + (1-cfr)*seq_prob)) %>%
      mutate(diff = exp_1 - exp) %>%
      mutate(type = "sequelae") %>%
      select(-cfr)
    
    
    prevented_ov = prevented_imd %>% 
      rbind(prevented_seq) %>%
      rbind(prevented_fat) %>%
      group_by(scen, type, sero) %>%
      summarise(diff = sum(diff)) %>%
      mutate(type = factor(type, levels = c("imd", "sequelae", "fatalities")))
    prevented_ov
  }
  
  est_nnv = function(imd_sim_scen1,
                     imd_sim_scen2,
                     imd_sim_scen3,
                     imd_sim_scen4,
                     imd_sim_scen5,
                     vacc_num_scen1,
                     vacc_num_scen2,
                     vacc_num_scen3,
                     vacc_num_scen4,
                     vacc_num_scen5,
                     cfr,
                     seq_prob,
                     times = 1:30) {
    prevented = est_prev(imd_sim_scen1,
                         imd_sim_scen2,
                         imd_sim_scen3,
                         imd_sim_scen4,
                         imd_sim_scen5,
                         cfr,
                         seq_prob,
                         times = times)
    add_immu = vacc_num_scen2 %>% mutate(scen="scen_2") %>%
      rbind(vacc_num_scen3 %>% mutate(scen="scen_3")) %>%
      rbind(vacc_num_scen4 %>% mutate(scen="scen_4")) %>%
      rbind(vacc_num_scen5 %>% mutate(scen="scen_5")) %>%
      filter(time %in% times) %>%
      left_join(vacc_num_scen1 %>%
                  rename(value_1 = value)) %>%
      mutate(add_immu = value - value_1) %>%
      group_by(scen) %>%
      summarise(add_immu = sum(add_immu))
      
    nnv = prevented %>% left_join(add_immu) %>%
      mutate(nnv = add_immu/diff,
             nnv_wo_o = add_immu/diff_wo_o) %>%
      mutate(nnv = ifelse(nnv<0, Inf, nnv),
             nnv = ifelse(nnv<0, Inf, nnv))
    nnv
  }
  
  prevented = est_prev(imd_sim_scen1,
                       imd_sim_scen2,
                       imd_sim_scen3,
                       imd_sim_scen4,
                       imd_sim_scen5,
                       cfr,
                       seq_prob,
                       times = times)
  
  prevented_samples = do.call(rbind,
                              lapply(1:length(imd_sim_samples_scen1),
                                     function(i) {
                                       est_prev(imd_sim_samples_scen1[[i]],
                                                imd_sim_samples_scen2[[i]],
                                                imd_sim_samples_scen3[[i]],
                                                imd_sim_samples_scen4[[i]],
                                                imd_sim_samples_scen5[[i]],
                                                cfr = cfr,
                                                seq_prob = seq_prob,
                                                times = times
                                       ) %>%
                                         mutate(iter=i)
                                     }))
  prevented_quant = prevented_samples %>%
    group_by(scen, type) %>%
    summarise(q025 = quantile(diff, 0.025),
              q975 = quantile(diff, 0.975),
              wo_o_q025 = quantile(diff_wo_o, 0.025),
              wo_o_q975 = quantile(diff_wo_o, 0.975),
              exp_q025 = quantile(exp, 0.025),
              exp_q975 = quantile(exp, 0.975),
              exp_wo_o_q025 = quantile(exp_wo_o, 0.025),
              exp_wo_o_q975 = quantile(exp_wo_o, 0.975),
              exp_1_q025 = quantile(exp_1, 0.025),
              exp_1_q975 = quantile(exp_1, 0.975),
              exp_1_wo_o_q025 = quantile(exp_1_wo_o, 0.025),
              exp_1_wo_o_q975 = quantile(exp_1_wo_o, 0.975))
  
  prevented_cases = prevented %>% left_join(prevented_quant)
  
  # Prevented cases by serogroup
  prevented_by_sero = est_prev_by_sero(imd_sim_scen1,
                                       imd_sim_scen2,
                                       imd_sim_scen3,
                                       imd_sim_scen4,
                                       imd_sim_scen5,
                                       cfr,
                                       seq_prob, 
                                       times = times)
  prevented_by_sero_samples = do.call(rbind,
                                      lapply(1:length(imd_sim_samples_scen1),
                                             function(i) {
                                               est_prev_by_sero(imd_sim_samples_scen1[[i]],
                                                                imd_sim_samples_scen2[[i]],
                                                                imd_sim_samples_scen3[[i]],
                                                                imd_sim_samples_scen4[[i]],
                                                                imd_sim_samples_scen5[[i]],
                                                                cfr = cfr,
                                                                seq_prob = seq_prob,
                                                                times = times
                                               ) %>%
                                                 mutate(iter=i)
                                             }))
  prevented_by_sero_quant = prevented_by_sero_samples %>%
    group_by(scen, type, sero) %>%
    summarise(q025 = quantile(diff, 0.025),
              q975 = quantile(diff, 0.975))
  
  prevented_cases_by_sero = prevented_by_sero %>% left_join(prevented_by_sero_quant) %>%
    mutate(sero = factor(sero, levels = c("c", "awy", "o"),
                         labels = c("C", "AWY", "Other/B"))) %>%
    arrange(scen, sero, type) %>%
    select(scen, sero, type, prevented=diff, q025, q975)
  # NNV
  nnv = est_nnv(imd_sim_scen1,
                imd_sim_scen2,
                imd_sim_scen3,
                imd_sim_scen4,
                imd_sim_scen5,
                vacc_num_scen1,
                vacc_num_scen2,
                vacc_num_scen3,
                vacc_num_scen4,
                vacc_num_scen5,
                cfr,
                seq_prob,
                times = times)
  nnv = nnv %>% left_join(prevented_cases) %>%
    mutate(q975 = pmax(q975, 0),
           q025 = pmax(q025, 0),
           nnv_q025 = add_immu/q975,
           nnv_q975 = add_immu/q025,
           nnv_wo_o_q025 = add_immu/wo_o_q975,
           nnv_wo_o_q975 = add_immu/wo_o_q025) %>%
    select(scen, type, diff, diff_wo_o, add_immu, nnv, nnv_q025, nnv_q975, nnv_wo_o, nnv_wo_o_q025, nnv_wo_o_q975)
  # Results for sensitivity table
  # Expected IMD BL scenario (by sero and overall)
  scen_1_imd = imd_sim_scen1 %>%
    filter(time %in% times) %>%
    group_by(sero) %>%
    summarise(imd = sum(exp_inv)) %>%
    rbind(imd_sim_scen1 %>%
            filter(time %in% times) %>%
            ungroup %>%
            summarise(imd = sum(exp_inv)) %>%
            mutate(sero = "overall"))
  scen_1_imd_samples = do.call(
    rbind,
    lapply(1:length(imd_sim_samples_scen1),
           function(i) {
             imd_sim_samples_scen1[[i]] %>%
               filter(time %in% times) %>%
               group_by(sero) %>%
               summarise(imd = sum(exp_inv)) %>%
               rbind(imd_sim_samples_scen1[[i]] %>%
                       filter(time %in% times) %>%
                       ungroup %>%
                       summarise(imd = sum(exp_inv)) %>%
                       mutate(sero = "overall")) %>%
               ungroup() %>%
               mutate(iter = i)
           }))
  scen_1_imd_smry = scen_1_imd %>%
    left_join(
      scen_1_imd_samples %>%
        group_by(sero) %>%
        summarise(q025 = quantile(imd, 0.025),
                  q975 = quantile(imd, 0.975))
    ) %>%
    mutate(scen = "scen_1",
           name = "Exp. IMD") %>%
    rename(value = imd) %>%
    select(scen, name, sero, value, q025, q975)
  # Scen 2-5
  scen_2_5_imd_diff = imd_sim_scen2 %>% mutate(scen = "scen_2") %>%
    rbind(imd_sim_scen3 %>% mutate(scen = "scen_3")) %>%
    rbind(imd_sim_scen4 %>% mutate(scen = "scen_4")) %>%
    rbind(imd_sim_scen5 %>% mutate(scen = "scen_5")) %>%
    filter(time %in% times) %>%
    left_join(imd_sim_scen1 %>% rename(imd_1 = exp_inv)) %>%
    group_by(sero, scen) %>%
    summarise(imd = sum(exp_inv),
              imd_1 = sum(imd_1)) %>%
    mutate(imd_prev = imd_1 - imd)
    
  scen_2_5_imd_diff_samples = do.call(
    rbind,
    lapply(1:length(imd_sim_samples_scen1),
           function(i) {
             imd_sim_samples_scen2[[i]] %>% mutate(scen = "scen_2") %>%
               rbind(imd_sim_samples_scen3[[i]] %>% mutate(scen = "scen_3")) %>%
               rbind(imd_sim_samples_scen4[[i]] %>% mutate(scen = "scen_4")) %>%
               rbind(imd_sim_samples_scen5[[i]] %>% mutate(scen = "scen_5")) %>%
               filter(time %in% times) %>%
               left_join(imd_sim_samples_scen1[[i]] %>% rename(imd_1 = exp_inv)) %>%
               group_by(sero, scen) %>%
               summarise(imd = sum(exp_inv),
                         imd_1 = sum(imd_1)) %>%
               mutate(imd_prev = imd_1 - imd) %>%
               mutate(iter = i)
           }))
  
  scen_2_5_imd_smry = scen_2_5_imd_diff %>%
    select(sero, scen, imd_prev) %>%
    left_join(
      scen_2_5_imd_diff_samples %>%
        group_by(scen, sero) %>%
        summarise(q025 = quantile(imd_prev, 0.025),
                  q975 = quantile(imd_prev, 0.975))
    ) %>%
    mutate(name = "Exp. prevented IMD") %>%
    rename(value = imd_prev) %>%
    select(scen, name, sero, value, q025, q975) %>%
    rbind(prevented_cases %>% filter(type == "imd") %>%
            mutate(name = "Exp. prevented IMD",
                   sero = "overall") %>%
            select(scen, name, sero, value = diff, q025, q975)) %>%
    rbind(nnv %>% filter(type == "imd") %>%
            mutate(name = "NNV",
                   sero = "overall") %>%
            select(scen, name, sero, value = nnv, q025=nnv_q025, q975=nnv_q975)) %>%
    arrange(scen, sero)
  
  list(prevented_cases = prevented_cases,
       nnv = nnv,
       smry_sens_tab = scen_1_imd_smry %>%
         rbind(scen_2_5_imd_smry),
       prevented_cases_by_sero = prevented_cases_by_sero)
}



#' Integrate DTM (model 2B) during simulation period given specific parameters 
#' under different vaccination scenarios
#' @param par Vector of parameters governing transmission dynamics (e.g., obtained from calibration)
#' @param int Integrated model from the calibration period used for initialization of the model
#' at the beginning of the simulation period
#' @param vacc_scen integer in [0, ..., 5] specifying the simulated vaccination strategy
#' @param grid_len grid-length used during numerical integration of the ODEs
#' @params parms list of parameters of the ODE (specified from external data)
#' @params rootfun Function defining time points of step changes
forward_sim_lambda2_8agecat= function(par,
                                      int,
                                      vacc_scen = 0,
                                      grid_len=1/12,
                                      parms,
                                      rootfun) {
  # Get initial conditions
  init = unlist(int[nrow(int),-1])
  # set helper compartments to zero
  init[unlist(map(0:85, function(x) 37*x + 21:37))] = 0
  integrate = deSolve::lsoda(
    y = init,
    times = seq(0,30, by = grid_len),
    func = function(t, y, parms) {
      list(calc_comp_dt_2(state_full = y,
                          parms = parms,
                          prob_inf_a_m = matrix(
                            rep(par, times = rep(c(5, 5, 5, 5, 10, 20, 20, 16), 3)), 
                            nrow = 86),
                          age_groups = 86)
      )
    },
    parms,
    rootfun = rootfun,
    events = list(func = function(t, y, parms) {
      stepChanges_sim(t, y, parms, scen = vacc_scen)
    }, root = TRUE))
  integrate
}

#' Plot carriage during calibration and simulation period for selected birth-cohorts
#' @param comp_nb_lambda2_age8 Compartments during calibration period
#' @param comp_sim_nb_lambda2_age8_scen*1-5* compartments during simulation for different vaccination scenarios
#' @param scens vector selecting which vaccination scenarios to plot
#' @param years vector of length 2 with plotting window (start and end year)
#' @param birthyears vector with birthyears included in plot
plot_comp_scen_by = function(comp_nb_lambda2_age8,
                             comp_sim_nb_lambda2_age8_scen1,
                             comp_sim_nb_lambda2_age8_scen2,
                             comp_sim_nb_lambda2_age8_scen3,
                             comp_sim_nb_lambda2_age8_scen4,
                             comp_sim_nb_lambda2_age8_scen5,
                             scens = paste0("scen_", 1:3),
                             years = c(2005, 2030),
                             birthyears = c(2020, 2010, 2000)) {
  plot_dat = comp_nb_lambda2_age8 %>% mutate(scen = "scen_1") %>%
    rbind(comp_nb_lambda2_age8 %>% mutate(scen = "scen_2")) %>%
    rbind(comp_nb_lambda2_age8 %>% mutate(scen = "scen_3")) %>%
    rbind(comp_nb_lambda2_age8 %>% mutate(scen = "scen_4")) %>%
    rbind(comp_nb_lambda2_age8 %>% mutate(scen = "scen_5")) %>%
    rbind(comp_sim_nb_lambda2_age8_scen1 %>% mutate(scen="scen_1") %>%
            rbind(comp_sim_nb_lambda2_age8_scen2 %>% mutate(scen="scen_2")) %>%
            rbind(comp_sim_nb_lambda2_age8_scen3 %>% mutate(scen="scen_3")) %>%
            rbind(comp_sim_nb_lambda2_age8_scen4 %>% mutate(scen="scen_4")) %>%
            rbind(comp_sim_nb_lambda2_age8_scen5 %>% mutate(scen="scen_5")) %>%
            filter(time>0) %>%
            mutate(time = time + 15)) %>%
    filter(year >= min(years),
           year<=max(years),
           birthyear %in% birthyears,
           scen %in% scens) %>% 
    group_by(time, comp, vacc, birthyear, scen) %>%
    summarise(value = sum(value)) %>% 
    ungroup() %>%
    group_by(time, birthyear, scen) %>%
    mutate(n=sum(value)) %>%
    ungroup() %>%
    mutate(share=value/n) %>% mutate(comp = factor(comp, levels = c("s", "car_c", "car_awy", "car_o"))) %>%
    mutate(time = time + 2005,
           birthyear = factor(birthyear, levels = sort(birthyears, decreasing = T)),
           vacc = factor(vacc, levels = c("uv",
                                          "vw_pmACWYc",
                                          "vi_pmACWYc",
                                          "vw_mmCc",
                                          "vi_mmCc"),
                         labels = c("Unvacc.\n",
                                    "MenACWY:\nwaned",
                                    "MenACWY\n",
                                    "MenC:\nwaned",
                                    "MenC\n")),
           comp = factor(comp, levels = c("s", "car_c", "car_awy", "car_o"),
                         labels = c("Susc.", "Carr. C", "Carr. AWY", "Carr. O")),
           scen = factor(scen, levels = paste0("scen_", 1:5),
                         labels = paste0("Scen. ", 1:5)))
  
  plot_list = lapply(birthyears, function(by) {
    plot_dat %>% filter(birthyear==by) %>%
      ggplot() + 
      geom_area(aes(time, y=share, fill = vacc), col = "black", linewidth=0.5) + 
      # geom_blank(aes(y=share, group = comp), data = limits_data) +
      facet_grid(cols= vars(scen), rows = vars(comp), scales = "free_y") + 
      theme_bw() +
      scale_y_continuous(labels = scales::percent) +
      ylab("Fraction") +
      xlab("Year") +
      theme(legend.position = "top", legend.title = element_blank()) +
      # guides(fill=guide_legend(title="Vacc. status:")) +
      scale_fill_manual(values = c("MenC\n" = "#1B9E77", 
                                   "Unvacc.\n" = "#D95F02", 
                                   "MenC:\nwaned" = "#7570B3", 
                                   "MenACWY:\nwaned" = "#E7298A", 
                                   "MenACWY\n" = "#66A61E"),
                        guide = guide_legend(reverse = TRUE)) +
      geom_vline(aes(xintercept=x), data = tibble(x=2020), lty = 2) +
      expand_limits(x = c(2005))
  })
  names(plot_list) = paste0("birthyear_", birthyears)
  plot_list
}

#' plots_sim_scen Plots of expected IMD cases under different simulation scenarios
#' @param new_imd_main_scen_*0-5* expected IMD cases under different simulation (vaccination) scenarios
#' @param new_imd_main_scen_*0-5_samples* list of expected IMD cases under different simulation (vaccination) scenarios
#' for different simulation parameters
#' @params times vector of simulation years used for plotting (1:10 corresponds to ten years simulation, i.e., 2020-2029)

plots_sim_scen = function(new_imd_main_scen_0,
                          new_imd_main_scen_1,
                          new_imd_main_scen_2,
                          new_imd_main_scen_3,
                          new_imd_main_scen_4,
                          new_imd_main_scen_5,
                          new_imd_main_scen_0_samples,
                          new_imd_main_scen_1_samples,
                          new_imd_main_scen_2_samples,
                          new_imd_main_scen_3_samples,
                          new_imd_main_scen_4_samples,
                          new_imd_main_scen_5_samples,
                          times = 1:10) {

  scen_col = RColorBrewer::brewer.pal(3, "Dark2")[c(1:3, 2, 3)] 
  names(scen_col) = paste0("Scen. ", 1:5)
  scen_col[1] = "darkgrey"
  dat = new_imd_main_scen_0 %>% mutate(scen = "Scen. 0") %>%
    rbind(new_imd_main_scen_1 %>% mutate(scen = "Scen. 1")) %>%
    rbind(new_imd_main_scen_2 %>% mutate(scen = "Scen. 2")) %>%
    rbind(new_imd_main_scen_3 %>% mutate(scen = "Scen. 3")) %>%
    rbind(new_imd_main_scen_4 %>% mutate(scen = "Scen. 4")) %>%
    rbind(new_imd_main_scen_5 %>% mutate(scen = "Scen. 5")) %>%
    mutate(year=2020-1+time) %>%
    mutate(sero = factor(sero, levels = c("c", "awy", "o"))) %>%
    filter(time %in% times)
  
  
  plot_abs = dat %>% 
    mutate(age_group = cut(age, breaks=c(0,5,15,20,86),
                                 right=FALSE, include.lowest=T,
                                 labels = c("0-4", "5-14", "15-19", "20+")),
                 sero = factor(sero, levels = c("c", "awy", "o"),
                               labels = c("C", "AWY", "Other/B"))) %>%
    select(-time, -age) %>%
    # group_by(age_group, sero, year, scen) %>%
    group_by(sero, year, scen) %>%
    summarise(exp_imd=sum(exp_inv)) %>%
    ungroup() %>%
    filter(scen != c("Scen. 0")) %>%
    ggplot() + 
    geom_line(aes(year, exp_imd, col = scen, lty = scen)) +
    # facet_grid(rows = vars(sero), cols = vars(age_group), scales = "free") +
    facet_grid(rows = vars(sero), scales = "free") +
    scale_color_manual(values = scen_col) +
    scale_linetype_manual(values = c("Scen. 1" = 1,
                                     "Scen. 2" = 2,
                                     "Scen. 3" = 2,
                                     "Scen. 4" = 3,
                                     "Scen. 5" = 3)) + 
    expand_limits(y=0) +
    ylab("Expected yearly IMD") +
    theme_bw() +
    scale_x_continuous(breaks=c(2020, 2025, 2030)) +
    xlab("Year") +
    theme(legend.position = "top", legend.title = element_blank()) +
    guides(lty=guide_legend(nrow=1,byrow=TRUE),
           col = guide_legend(nrow=1,byrow=TRUE))
  
  plot_abs_1 = dat %>% 
    mutate(age_group = cut(age, breaks=c(0,5,15,20,86),
                           right=FALSE, include.lowest=T,
                           labels = c("0-4", "5-14", "15-19", "20+")),
           sero = factor(sero, levels = c("c", "awy", "o"),
                         labels = c("C", "AWY", "Other/B"))) %>%
    select(-time, -age) %>%
    # group_by(age_group, sero, year, scen) %>%
    group_by(sero, year, scen) %>%
    summarise(exp_imd=sum(exp_inv)) %>%
    ungroup() %>%
    filter(scen != c("Scen. 0")) %>%
    filter(scen %in% c("Scen. 1")) %>%
    ggplot() + 
    geom_line(aes(year, exp_imd, col = scen, lty = scen)) +
    # facet_grid(rows = vars(sero), cols = vars(age_group), scales = "free") +
    facet_grid(rows = vars(sero), scales = "free") +
    scale_color_manual(values = scen_col) +
    scale_linetype_manual(values = c("Scen. 1" = 1,
                                     "Scen. 2" = 2,
                                     "Scen. 3" = 2,
                                     "Scen. 4" = 3,
                                     "Scen. 5" = 3)) + 
    expand_limits(y=0) +
    ylab("Expected yearly IMD") +
    theme_bw() +
    scale_x_continuous(breaks=c(2020, 2025, 2030)) +
    xlab("Year") +
    theme(legend.position = "top", legend.title = element_blank()) +
    guides(lty=guide_legend(nrow=1,byrow=TRUE),
           col = guide_legend(nrow=1,byrow=TRUE))
  
  plot_abs_124 = dat %>% 
    mutate(age_group = cut(age, breaks=c(0,5,15,20,86),
                           right=FALSE, include.lowest=T,
                           labels = c("0-4", "5-14", "15-19", "20+")),
           sero = factor(sero, levels = c("c", "awy", "o"),
                         labels = c("C", "AWY", "Other/B"))) %>%
    select(-time, -age) %>%
    # group_by(age_group, sero, year, scen) %>%
    group_by(sero, year, scen) %>%
    summarise(exp_imd=sum(exp_inv)) %>%
    ungroup() %>%
    filter(scen != c("Scen. 0")) %>%
    filter(scen %in% c("Scen. 1", "Scen. 2", "Scen. 4")) %>%
    ggplot() + 
    geom_line(aes(year, exp_imd, col = scen, lty = scen)) +
    # facet_grid(rows = vars(sero), cols = vars(age_group), scales = "free") +
    facet_grid(rows = vars(sero), scales = "free") +
    scale_color_manual(values = scen_col) +
    scale_linetype_manual(values = c("Scen. 1" = 1,
                                     "Scen. 2" = 2,
                                     "Scen. 3" = 2,
                                     "Scen. 4" = 3,
                                     "Scen. 5" = 3)) + 
    expand_limits(y=0) +
    ylab("Expected yearly IMD") +
    theme_bw() +
    scale_x_continuous(breaks=c(2020, 2025, 2030)) +
    xlab("Year") +
    theme(legend.position = "top", legend.title = element_blank()) +
    guides(lty=guide_legend(nrow=1,byrow=TRUE),
           col = guide_legend(nrow=1,byrow=TRUE))
  
  
  calc_rel_change = function(dat) {
    dat %>% ungroup() %>%
      mutate(age_group = cut(age, breaks=c(0,5,15,20,86), 
                             right=FALSE, include.lowest=T,
                             labels = c("0-4", "5-14", "15-19", "20+"))) %>%
      filter(scen == "scen_1") %>%
      group_by(time, sero, age_group, scen) %>%
      summarise(exp_inv=sum(exp_inv)) %>%
      group_by(sero, age_group, scen) %>%
      arrange(scen, age_group, sero, time) %>%
      mutate(cum_exp_inv = cumsum(exp_inv)) %>%
      ungroup() %>%
      select(time, sero, age_group, exp_imd_scen_1 = cum_exp_inv) %>% 
      right_join(dat %>% ungroup() %>%
                   mutate(age_group = cut(age, breaks=c(0,5,15,20,86), 
                                          right=FALSE, include.lowest=T,
                                          labels = c("0-4", "5-14", "15-19", "20+"))) %>%
                   group_by(time, sero, age_group, scen) %>%
                   summarise(exp_inv=sum(exp_inv)) %>%
                   group_by(sero, age_group, scen) %>%
                   arrange(scen, age_group, sero, time) %>%
                   mutate(cum_exp_inv = cumsum(exp_inv)) %>%
                   ungroup(), 
                 by = c("time", "sero", "age_group")) %>%
      mutate(perc_diff = (cum_exp_inv-exp_imd_scen_1)/exp_imd_scen_1*100)
  }
  
  rel_change_samples = do.call(rbind, 
                               lapply(1:length(new_imd_main_scen_0_samples), 
                                      function(iter) {
                                        dat = new_imd_main_scen_0_samples[[iter]] %>% 
                                          mutate(scen = "scen_0") %>%
                                          rbind(new_imd_main_scen_1_samples[[iter]] %>% 
                                                  mutate(scen = "scen_1")) %>%
                                          rbind(new_imd_main_scen_2_samples[[iter]] %>% 
                                                  mutate(scen = "scen_2")) %>%
                                          rbind(new_imd_main_scen_3_samples[[iter]] %>% 
                                                  mutate(scen = "scen_3")) %>%
                                          rbind(new_imd_main_scen_4_samples[[iter]] %>% 
                                                  mutate(scen = "scen_4")) %>%
                                          rbind(new_imd_main_scen_5_samples[[iter]] %>% 
                                                  mutate(scen = "scen_5")) %>%
                                          mutate(year=2020-1+time) %>%
                                          mutate(sero = factor(sero, levels = c("c", "awy", "o"))) %>%
                                          filter(time %in% times) %>%
                                          ungroup() 
                                        rel_change = calc_rel_change(dat) %>% mutate(iter)
                                        rel_change
                                      }))
  
  quant = rel_change_samples %>% group_by(sero, time, age_group, scen) %>% 
    summarise(mean = round(mean(perc_diff),3), 
              q025 = round(quantile(perc_diff, 0.025),1),
              q975 = round(quantile(perc_diff, 0.975),1)) %>%
    ungroup() %>%
    mutate(year=2020-1+time) %>%
    mutate(scen=factor(scen, levels = paste0("scen_", 0:5),
                       labels = paste0("Scen. ", 0:5)),
           sero=factor(sero, levels = c("c", "awy", "o"),
                       labels = c("C", "AWY", "Other/B")))
    
  
  plot_rel_age_cat = quant %>% filter(scen!="Scen. 0") %>%
    droplevels()%>%
    ggplot() + 
    geom_line(aes(year, mean/100, col = scen, lty = scen)) +
    # geom_ribbon(aes(year, ymin=q025, ymax=q975, fill=scen), alpha=.9) +
    facet_grid(rows = vars(sero), cols = vars(age_group), scales = "free_y") +
    scale_color_manual(values = scen_col) +
    scale_linetype_manual(values = c("Scen. 1" = 1,
                                     "Scen. 2" = 2,
                                     "Scen. 3" = 2,
                                     "Scen. 4" = 3,
                                     "Scen. 5" = 3)) + 
    scale_x_continuous(breaks=c(2020, 2025, 2030)) +
    scale_y_continuous(n.breaks=3, labels = scales::percent) +
    xlab("Year") +
    ylab("Relative expected IMD") +
    theme_bw() +
    theme(legend.position = "top", legend.title = element_blank()) +
    expand_limits(y=c(-0.025, 0.05))
  
  
  
  ##### barplot
  dat = new_imd_main_scen_0 %>% mutate(scen = "scen_0") %>%
    rbind(new_imd_main_scen_1 %>% mutate(scen = "scen_1")) %>%
    rbind(new_imd_main_scen_2 %>% mutate(scen = "scen_2")) %>%
    rbind(new_imd_main_scen_3 %>% mutate(scen = "scen_3")) %>%
    rbind(new_imd_main_scen_4 %>% mutate(scen = "scen_4")) %>%
    rbind(new_imd_main_scen_5 %>% mutate(scen = "scen_5")) %>%
    mutate(year=2020-1+time) %>%
    mutate(sero = factor(sero, levels = c("c", "awy", "o"))) %>%
    filter(time %in% times)
  
  
  calc_rel_change = function(dat) {
    dat %>% group_by(scen, sero) %>%
      summarise(exp_imd = sum(exp_inv)) %>%
      filter(scen == "scen_1") %>%
      ungroup() %>%
      select(sero, exp_imd_scen_1 = exp_imd) %>% 
      right_join(dat %>% group_by(scen, sero) %>%
                   summarise(exp_imd = sum(exp_inv))) %>%
      mutate(perc_diff = (exp_imd-exp_imd_scen_1)/exp_imd_scen_1*100)
  }
  
  rel_change_samples = do.call(rbind, lapply(1:length(new_imd_main_scen_0_samples), 
                                             function(iter) {
                                               dat = new_imd_main_scen_0_samples[[iter]] %>% 
                                                 mutate(scen = "scen_0") %>%
                                                 rbind(new_imd_main_scen_1_samples[[iter]] %>% 
                                                         mutate(scen = "scen_1")) %>%
                                                 rbind(new_imd_main_scen_2_samples[[iter]] %>% 
                                                         mutate(scen = "scen_2")) %>%
                                                 rbind(new_imd_main_scen_3_samples[[iter]] %>% 
                                                         mutate(scen = "scen_3")) %>%
                                                 rbind(new_imd_main_scen_4_samples[[iter]] %>% 
                                                         mutate(scen = "scen_4")) %>%
                                                 rbind(new_imd_main_scen_5_samples[[iter]] %>% 
                                                         mutate(scen = "scen_5")) %>%
                                                 mutate(year=2020-1+time) %>%
                                                 mutate(sero = factor(sero, levels = c("c", "awy", "o"))) %>%
                                                 filter(time %in% times)
                                               rel_change = calc_rel_change(dat) %>% mutate(iter)
                                               rel_change
                                             }))
  
  
  quant = rel_change_samples %>% group_by(sero, scen) %>% 
    summarise(mean_change = round(mean(perc_diff),3), 
              q025_change = round(quantile(perc_diff, 0.025),1),
              q975_change = round(quantile(perc_diff, 0.975),1),
              mean = round(mean(exp_imd),3), 
              q025 = round(quantile(exp_imd, 0.025),1),
              q975 = round(quantile(exp_imd, 0.975),1)) %>% ungroup()
  
  plot_dat = dat %>% group_by(scen, sero) %>%
    summarise(exp_imd = sum(exp_inv)) %>%
    filter(scen == "scen_1") %>%
    ungroup() %>%
    select(sero, exp_imd_scen_1 = exp_imd) %>% 
    right_join(dat %>% group_by(scen, sero) %>%
                 summarise(exp_imd = sum(exp_inv)) %>%
                 filter(scen != "scen_0")) %>%
    left_join(quant) %>%
    droplevels()
  
  fill_col = scen_col[1:3]
  names(fill_col) = c("No booster", "MenC", "MenACWY")
  barplot = plot_dat %>% 
    mutate(scen = gsub(pattern = "scen_", "", scen),
           rel_change = ifelse(scen!=1,
                               paste0(ifelse(mean_change >0, "+", ""), round(mean_change ,1), 
                                      "% [",
                                      ifelse(abs(mean_change)>10, 
                                             round(q025_change,0),
                                             round(q025_change,1)),
                                      "; ", 
                                      ifelse(abs(mean_change)>10, 
                                             round(q975_change,0),
                                             round(q975_change,1)),"]"),
                               ""),
           sero = factor(sero, 
                         levels = c("c", "awy", "o"),
                         labels = c("C", "AWY", "O")),
           fill_col = factor(scen, levels = 0:5,
                             labels = c("No vac.",
                                        "No booster", 
                                        "MenC",
                                        "MenACWY",
                                        "MenC",
                                        "MenACWY"))) %>% 
    ggplot() +
    geom_col(aes(scen, exp_imd, fill = fill_col), alpha = .8, col = "darkgrey") +
    geom_errorbar(aes(scen, ymin = q025, ymax = q975), width=0.1) +
    geom_text(aes(scen, .4*exp_imd, label = rel_change),
              family = 'Times', size = 3.5, lineheight=0.9) +
    facet_wrap(~sero, ncol = 1, scales = "free_y") +
    theme_bw() +
    xlab("Simulation scenario") +
    ylab("Expected IMD cases (2020-2049)") +
    scale_fill_manual(values = fill_col)  +
    theme(legend.position = "top", legend.title = element_blank())
  ###
  list(plot_abs=plot_abs, 
       plot_rel_age_cat=plot_rel_age_cat,
       barplot = barplot,
       plot_abs_1 = plot_abs_1,
       plot_abs_124 = plot_abs_124)  
    
}


#' plot_carr_age_year_sim Plot carriage by age for a sequence of years during the simulation period
#' for simulation scenarios 1-3
plot_carr_age_year_sim = function(comp_sim_lambda2_age8_scen1, 
                                  comp_sim_lambda2_age8_scen2,
                                  comp_sim_lambda2_age8_scen3) {
  col_vec = RColorBrewer::brewer.pal(3, name = "Dark2")
  names(col_vec) = c("Scen. 1/\nno Booster",
                     "Scen. 2/\nMenC",
                     "Scen. 3/\nMenACWY")
  col_vec[1] = "darkgrey"
  years = seq(2020, 2030, by = 2)
  comp_sim_lambda2_age8_scen1 %>% mutate(scen="Scen. 1/\nno Booster") %>%
    rbind(comp_sim_lambda2_age8_scen2 %>% mutate(scen="Scen. 2/\nMenC")) %>%
    rbind(comp_sim_lambda2_age8_scen3 %>% mutate(scen="Scen. 3/\nMenACWY")) %>%
    filter((time %% 1)==0) %>%
    filter(year %in% years) %>%
    mutate(comp = factor(comp, levels = c("s", "car_c", "car_awy", "car_o"),
                         labels = c("s", "C", "AWY", "Other/B"))) %>%
    group_by(scen, year, age, comp) %>% 
    summarise(n=sum(value)) %>% 
    group_by(scen, year, age) %>% 
    mutate(N=sum(n), frac = n/N) %>% 
    filter(comp!="s") %>%
    ungroup() %>%
    ggplot() + 
    geom_line(aes(age, frac, col = scen, lty=scen)) +
    facet_grid(cols = vars(year),
               rows = vars(comp), scales = "free_y") +
    scale_y_continuous(labels = scales::percent) +
    theme_bw() +
    ylab("Fraction") +
    scale_color_manual(values = col_vec) +
    theme(legend.position = "bottom",
          legend.title = element_blank())
}


#' Formatting tables and figures
format_averted = function(x) {
  x_form = case_when(abs(x)>1000 ~ round(x, -2), 
            abs(x)>100 ~ round(x, -1),
            abs(x)>=10 ~ round(x, 0),
            abs(x)<10 ~ round(x, 1))
  case_when(abs(x_form)>=10 ~ format(x_form, big.mark = ",", decimal.mark = ".",  digits = 1, trim = T),
            abs(x_form)<10 ~ format(x_form, big.mark = ",", decimal.mark = ".",  digits = 1, trim = T))
}

format_nnv = function(x) {
  x_form = case_when(
    abs(x)>1e8 ~ round(x, -7),
    abs(x)>1e7 ~ round(x, -6),
    abs(x)>1e6 ~ round(x, -5),
    abs(x)>1e5 ~ round(x, -4),
    abs(x)>1e4 ~ round(x, -3),
    abs(x)>1e3 ~ round(x, -2), 
    abs(x)>100 ~ round(x, -1),
    abs(x)>=10 ~ round(x, 0),
    abs(x)<10 ~ round(x, 1))
  case_when(is.infinite(x_form) ~ "-",
            x_form>=0 ~ format(x_form, big.mark = ",", decimal.mark = ".",  
                               trim = T, 
                               scientific = FALSE, digits = 1))
}


format_plot_sim_res_30y = function(plot_sim_res_list) {
  ggpubr::ggarrange(
    plot_sim_res_list[[3]] + theme(axis.text=element_text(size=10)),
    ggpubr::ggarrange(plot_sim_res_list[[1]] + theme(axis.text=element_text(size=10)) +
                        xlab("Year") + ylab("Expected IMD cases\nper year") +
                        scale_x_continuous(breaks = c(2020, 2035, 2050)) +
                        geom_vline(aes(xintercept = x), data = tibble(x=2030), lty = 2),
                      plot_sim_res_list[[2]] + theme(axis.text.x = element_text(size=10),
                                                     axis.text.y = element_text(size=9)) +
                        scale_x_continuous(breaks = c(2020, 2035)) +
                        xlab("Year") +
                        ylab("Rel. change \ncumulative IMD") +
                        geom_vline(aes(xintercept = x), data = tibble(x=2030), lty = 2),
                      ncol = 2, labels = c("B", "C"), common.legend = TRUE),
    nrow = 2, labels = c("A", ""), heights = c(1.25, 1))
}

format_plot_sim_res_fig_3 = function(plot_sim_res_list,
                                     plot_sim_res_30y_list) {
  ggpubr::ggarrange(
    plot_sim_res_list[[3]] + theme(axis.text=element_text(size=10)),
    ggpubr::ggarrange(
      plot_sim_res_list[[1]] + theme(axis.text=element_text(size=10)) +
        xlab("Year") + ylab("Expected IMD cases\nper year") +
        scale_x_continuous(breaks = c(2020, 2025, 2030)),
      plot_sim_res_list[[2]] + theme(axis.text.x = element_text(size=10),
                                     axis.text.y = element_text(size=9)) +
        scale_x_continuous(breaks = c(2020, 2025)) +
        xlab("Year") +
        ylab("Rel. change \ncumulative IMD"),
      plot_sim_res_30y_list[[1]] + theme(axis.text=element_text(size=10)) +
        xlab("Year") + ylab("Expected IMD cases\nper year") +
        scale_x_continuous(breaks = c(2020, 2035, 2050)) +
        geom_vline(aes(xintercept = x), data = tibble(x=2030), lty = 2),
      plot_sim_res_30y_list[[2]] + theme(axis.text.x = element_text(size=10),
                                         axis.text.y = element_text(size=9)) +
        scale_x_continuous(breaks = c(2020, 2035)) +
        xlab("Year") +
        ylab("Rel. change \ncumulative IMD") +
        geom_vline(aes(xintercept = x), data = tibble(x=2030), lty = 2),
      ncol = 2, nrow = 2, labels = c("B", "C", "D", "E"), common.legend = TRUE, widths = c(1, 1.5)),
    nrow = 2, labels = c("A", ""), heights = c(1, 1.25))
}
