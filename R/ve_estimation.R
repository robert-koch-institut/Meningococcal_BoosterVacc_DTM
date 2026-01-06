library(tidyverse)

#' mean VE for exponential waning with rate r (in years)
#' during the time interval a, b in years
#' ve is VE at t=0

int = function(r, ve, a, b) {
  ((ve*(exp(-a*r) - exp(-b*r)))/r)/(b-a)
} 

#' target function for optimisation
#' @param r waning rate
#' @param ve0_imd VE against IMD at time t=0
#' @param ve0_arr VE against carriage at time t=0
#' @param obs_dat observed VE estimates with information on VE estimate, study 
#' period (duration of observation) and study target (VE against imd or carriage) 
#' in specific format (see below)

target = function(r,
                  ve0_imd,
                  ve0_carr,
                  obs_dat) {
  
  mod_ve_imd = int(r=r,
                   ve=ve0_imd,
                   a = obs_dat %>% filter(type=="imd") %>% pull(a),
                   b = obs_dat %>% filter(type=="imd") %>% pull(b))
  mod_ve_carr = int(r=r,
                    ve=ve0_carr,
                    a = obs_dat %>% filter(type=="carr") %>% pull(a),
                    b = obs_dat %>% filter(type=="carr") %>% pull(b))
  
  sum((obs_dat %>% filter(type=="imd") %>% pull(obs_ve) - mod_ve_imd)^2) +
    sum((obs_dat %>% filter(type=="carr") %>% pull(obs_ve) - mod_ve_carr)^2)
}


# MenC infant
obs_dat = tibble(obs_ve = c(0.96, .68, 0.60, .31,  # VE estimates against imd Update post-license
                            .75), # Maiden (2008) VE estimate against carriage
                 a = c(0, 1, 2, 3,
                       0),
                 b = c(1, 2, 3, 5,
                       1),
                 type = c(rep("imd", 4), 
                          "carr"))


obs_dat %>% ggplot() + 
  geom_segment(aes(x=a, xend=b, y=obs_ve, yend=obs_ve)) +
  expand_limits(y=c(0,1)) +
  theme_bw() +
  ylab("VE") +
  xlab("Time in years") +
  facet_wrap(~type)


opt_res = optim(c(log_exp_dur = 0,
                  par_car=0,
                  par_imd=0.1), 
                fn = function(par) {
                  r=1/exp(par[1])
                  ve0_carr = plogis(par[2])
                  ve0_imd = plogis(par[2] + par[3])
                  target(r=r,
                         ve0_carr = ve0_carr,
                         ve0_imd = ve0_imd,
                         obs_dat = obs_dat)
                },
                hessian = T)

tibble(par = c("mean_dur_prot", "ve0_carr", "ve0_imd"),
       est = c(exp(opt_res$par[1]),
               plogis(opt_res$par[2]),
               plogis(opt_res$par[2] + opt_res$par[3])))

n_pred = 200
ve_pred = tibble(t=rep(seq(0, 15, length.out=n_pred), 2),
                 prot = c(plogis(opt_res$par[2])*exp(-1/exp(opt_res$par[1])*seq(0,15, length.out=n_pred)),
                          plogis(opt_res$par[2] + opt_res$par[3])*exp(-1/exp(opt_res$par[1])*seq(0,15, length.out=n_pred))),
                 type = c(rep("carr", n_pred),
                          rep("imd", n_pred)))

obs_dat %>% ggplot() + 
  geom_segment(aes(x=a, xend=b, y=obs_ve, yend=obs_ve, col = type), lty = 2) +
  expand_limits(y=c(0,1)) +
  theme_bw() +
  ylab("VE") +
  xlab("Time in years") +
  geom_line(aes(t, prot, col = type), data = ve_pred)



# MenACWY infant
obs_dat = tibble(obs_ve = c(0.96, .68, 0.60, .31,  # VE estimates against imd Update post-license
                            .14), # Maiden (2008) VE estimate against carriage
                 a = c(0, 1, 2, 3,
                       1/12),
                 b = c(1, 2, 3, 5,
                       1),
                 type = c(rep("imd", 4), 
                          "carr"))

obs_dat %>% ggplot() + 
  geom_segment(aes(x=a, xend=b, y=obs_ve, yend=obs_ve)) +
  expand_limits(y=c(0,1)) +
  theme_bw() +
  ylab("VE") +
  xlab("Time in years") +
  facet_wrap(~type)


opt_res = optim(c(log_exp_dur = 0,
                  par_car=0,
                  par_imd=0.1), 
                fn = function(par) {
                  r=1/exp(par[1])
                  ve0_carr = plogis(par[2])
                  ve0_imd = plogis(par[2] + par[3])
                  target(r=r,
                         ve0_carr = ve0_carr,
                         ve0_imd = ve0_imd,
                         obs_dat = obs_dat)
                },
                hessian = T)

tibble(par = c("mean_dur_prot", "ve0_carr", "ve0_imd"),
       est = c(exp(opt_res$par[1]),
               plogis(opt_res$par[2]),
               plogis(opt_res$par[2] + opt_res$par[3])))

n_pred = 200
ve_pred = tibble(t=rep(seq(0, 15, length.out=n_pred), 2),
                 prot = c(plogis(opt_res$par[2])*exp(-1/exp(opt_res$par[1])*seq(0,15, length.out=n_pred)),
                          plogis(opt_res$par[2] + opt_res$par[3])*exp(-1/exp(opt_res$par[1])*seq(0,15, length.out=n_pred))),
                 type = c(rep("carr", n_pred),
                          rep("imd", n_pred)))

obs_dat %>% ggplot() + 
  geom_segment(aes(x=a, xend=b, y=obs_ve, yend=obs_ve, col = type), lty = 2) +
  expand_limits(y=c(0,1)) +
  theme_bw() +
  ylab("VE") +
  xlab("Time in years") +
  geom_line(aes(t, prot, col = type), data = ve_pred)

