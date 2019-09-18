# Set up ------------------------------------------------------------------

if(!require(here)) install.packages("here")
if(require(here)) {
  setwd(here::here())
}else{
  warning('make sur you have package "here" installed')
}

source("GlobalParameters.R")

list.files(glob_params$RFunc,
           pattern="*.R$",
           full.names=TRUE, ignore.case=TRUE)

pckToLoad = c('tidyverse', "rstan")
reloadpck()
# source(glob_params$f.RFunc("get_os.R"))
source(glob_params$f.RFunc("get_timestamp.R"))
tmp.timestamp = get_timestamp() ; tmp.timestamp

# Modelling ---------------------------------------------------------------
## See the corresponding Rmarkdown file.

# Compile stan code -------------------------------------------------------

tmp.stanfile.name = "step_by_step_uptoKE1_CIa_loglogistic" ; paste0("model_",tmp.stanfile.name) %>% cat

model_step_by_step_uptoKE1_CIa_loglogistic = 
  stan_model(stanc_ret = stanc(glob_params$Stan$f.CodePath(tmp.stanfile.name,".stan")))

# Prepare stan input ------------------------------------------------------
source(glob_params$f.RScript('load-and-reshape-Konstanz-data.R'))

# Fit for Rotenone --------------------------------------------------------

tmp.chemi = "Rotenone"

tmp.data.KE1 = KE1.both %>% filter(., chemi == tmp.chemi, concentration_MuMol > 0)

stan_input.Neuro = list(
  "is_rot_or_deg" = as.integer(tmp.chemi == "Rotenone"),
  "ndoses_CIa" = length(tmp.data.KE1$concentration_MuMol),
  "dose_CIa"   = tmp.data.KE1$concentration_MuMol,
  "CIa_CIa"    = tmp.data.KE1$c1_activity
)

tmp.parsName = c("CIa_min" ,
                 "CIa_max" ,
                 "Par_CIa_50" ,
                 "log_Par_CIa_50",
                 "k_CIa" ,
                 "sigma_CIa",
                 
                 "lp__")

nb.chains = 3
nb.iter = 10000
my.seed = tmp.timestamp

source(glob_params$f.RFunc("stanFit_Name.R"))
tmp.stanFit.name = f.stanFit.newName(chemi         = tmp.chemi,
                                     stanfile.name = tmp.stanfile.name,
                                     timestamp     = tmp.timestamp,
                                     nb.chains     = nb.chains,
                                     nb.iter       = nb.iter,
                                     seed          = my.seed) 

fit_step_by_step_uptoKE1_CIa_loglogistic = sampling(
  model_step_by_step_uptoKE1_CIa_loglogistic,
  data = stan_input.Neuro,
  pars = tmp.parsName,
  chains = nb.chains,
  refresh = min(nb.iter/10,50),
  # sample_file = paste0(tmp.stanFit.name,".csv"),
  iter = nb.iter,
  seed = my.seed
)

# source(glob_params$f.RScript("sound_simulation_finished.R"))

# shinystan::launch_shinystan(fit_step_by_step_uptoKE1_CIa_loglogistic)

fit_step_by_step_uptoKE1_CIa_loglogistic

tmp.summary = summary(fit_step_by_step_uptoKE1_CIa_loglogistic,
                      probs = c(0.05, 0.95),
                      pars = tmp.parsName)$summary %>% 
  as.data.frame() 

tmp.n_eff.min     = tmp.summary$n_eff[!is.na(tmp.summary$n_eff)] %>% min() %>% ceiling()
tmp.rhat.max      = tmp.summary$Rhat[is.finite(tmp.summary$Rhat)] %>% max() %>% round(digits = 2)

tmp.stanFit.name = f.stanFit.newName(chemi         = tmp.chemi,
                                     stanfile.name = tmp.stanfile.name,
                                     timestamp     = tmp.timestamp,
                                     nb.chains     = nb.chains,
                                     nb.iter       = nb.iter,
                                     seed          = my.seed,
                                     n_eff.min     = tmp.n_eff.min, 
                                     rhat.max      = tmp.rhat.max) 

write_rds(fit_step_by_step_uptoKE1_CIa_loglogistic,
          path = glob_params$Stan$f.FitPath(tmp.stanFit.name, ".stanFit"))
write_lines(get_stancode(fit_step_by_step_uptoKE1_CIa_loglogistic), 
            glob_params$Stan$f.FitPath(tmp.stanFit.name, ".stan"), append = FALSE)
file.copy(here::here("Rnotebook", "step_by_step_uptoKE1_CIa_loglogistic.R"),
          glob_params$Stan$f.FitPath(tmp.stanFit.name,".R"), overwrite = T)


## Posterior summary -------------------------------------------------------

if(F){
  tmp.fit = fit_step_by_step_uptoKE1_CIa_loglogistic
  
  la_all = as.data.frame(rstan::extract(tmp.fit, pars = tmp.parsName))
  
  la.map_all = la_all[which.max(la_all$lp__),]
  
  
  la.summary_all = summary(tmp.fit,
                           probs = c(0.05, 0.95),
                           pars = tmp.parsName)$summary %>% as.data.frame()
  
  la.summary_all$MAP = la.map_all %>% as.matrix() %>% c()
  la.mean_all = la.summary_all["mean"] %>% t() %>% as.data.frame()
  
  la.summary_all = la.summary_all %>% select(MAP, mean,se_mean, sd, `5%`, `95%`,n_eff,Rhat)
  
  # options("scipen"= 2,digits=5)
  # options("scipen"= 0,digits=7)
  la.summary_all %>% round(.,digits = 3) %>%  DT::datatable()
}

# Fit for Deguelin --------------------------------------------------------

if(exists("fit_step_by_step_uptoKE1_CIa_loglogistic")) rm(fit_step_by_step_uptoKE1_CIa_loglogistic)

tmp.chemi = "Deguelin"

tmp.data.KE1 = KE1.both %>% filter(., chemi == tmp.chemi, concentration_MuMol > 0)

stan_input.Neuro = list(
  "is_rot_or_deg" = as.integer(tmp.chemi == "Rotenone"),
  "ndoses_CIa" = length(tmp.data.KE1$concentration_MuMol),
  "dose_CIa"   = tmp.data.KE1$concentration_MuMol,
  "CIa_CIa"    = tmp.data.KE1$c1_activity
)

tmp.stanFit.name = f.stanFit.newName(chemi         = tmp.chemi,
                                     stanfile.name = tmp.stanfile.name,
                                     timestamp     = tmp.timestamp,
                                     nb.chains     = nb.chains,
                                     nb.iter       = nb.iter,
                                     seed          = my.seed) 

fit_step_by_step_uptoKE1_CIa_loglogistic = sampling(
  model_step_by_step_uptoKE1_CIa_loglogistic,
  data = stan_input.Neuro,
  pars = tmp.parsName,
  chains = nb.chains,
  refresh = min(nb.iter/10,50),
  # sample_file = paste0(tmp.stanFit.name,".csv"),
  iter = nb.iter,
  seed = my.seed
)

# shinystan::launch_shinystan(fit_step_by_step_uptoKE1_CIa_loglogistic)

fit_step_by_step_uptoKE1_CIa_loglogistic

tmp.summary = summary(fit_step_by_step_uptoKE1_CIa_loglogistic,
                      probs = c(0.05, 0.95),
                      pars = tmp.parsName)$summary %>% 
  as.data.frame() 

tmp.n_eff.min     = tmp.summary$n_eff[!is.na(tmp.summary$n_eff)] %>% min() %>% ceiling()
tmp.rhat.max      = tmp.summary$Rhat[is.finite(tmp.summary$Rhat)] %>% max() %>% round(digits = 2)

tmp.stanFit.name = f.stanFit.newName(chemi         = tmp.chemi,
                                     stanfile.name = tmp.stanfile.name,
                                     timestamp     = tmp.timestamp,
                                     nb.chains     = nb.chains,
                                     nb.iter       = nb.iter,
                                     seed          = my.seed,
                                     n_eff.min     = tmp.n_eff.min, 
                                     rhat.max      = tmp.rhat.max) 

write_rds(fit_step_by_step_uptoKE1_CIa_loglogistic,
          path = glob_params$Stan$f.FitPath(tmp.stanFit.name, ".stanFit"))
write_lines(get_stancode(fit_step_by_step_uptoKE1_CIa_loglogistic), 
            glob_params$Stan$f.FitPath(tmp.stanFit.name, ".stan"), append = FALSE)
file.copy(here::here("Rnotebook", "step_by_step_uptoKE1_CIa_loglogistic.R"),
          glob_params$Stan$f.FitPath(tmp.stanFit.name,".R"), overwrite = T)

source(glob_params$f.RScript("sound_simulation_finished.R"))


# Posterior summary -------------------------------------------------------

if(F){
  tmp.fit = fit_step_by_step_uptoKE1_CIa_loglogistic
  
  la_all = as.data.frame(rstan::extract(tmp.fit, pars = tmp.parsName))
  
  la.map_all = la_all[which.max(la_all$lp__),]
  
  
  la.summary_all = summary(tmp.fit,
                           probs = c(0.05, 0.95),
                           pars = tmp.parsName)$summary %>% as.data.frame()
  
  la.summary_all$MAP = la.map_all %>% as.matrix() %>% c()
  la.mean_all = la.summary_all["mean"] %>% t() %>% as.data.frame()
  
  la.summary_all = la.summary_all %>% select(MAP, mean,se_mean, sd, `5%`, `95%`,n_eff,Rhat)
  
  # options("scipen"= 2,digits=5)
  # options("scipen"= 0,digits=7)
  la.summary_all %>% round(.,digits = 3) %>%  DT::datatable()
}
