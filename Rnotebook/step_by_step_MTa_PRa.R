# Set up ------------------------------------------------------------------

if(!require(here)) install.packages("here")
if(require(here)) {
  setwd(here::here())
}else{
  warning('make sur you have package "here" installed')
}

source("GlobalParameters.R") 
pckToLoad = c('tidyverse', "rstan")
reloadpck()
source(glob_params$RFunc %>% paste0("get_os.R"))

tmp.timestamp =  gsub("[: -]", "" , Sys.time(), perl=TRUE) %>% as.numeric() %% 1E12 %/% 100

# Modelling ---------------------------------------------------------------
## See the corresponding Rmarkdown file.

# Compile stan code -------------------------------------------------------

tmp.stanfile.name = "step_by_step_MTa_PRa" ; paste0("model_",tmp.stanfile.name) %>% cat

tmp.stanCodeBU.name = paste0(tmp.stanfile.name,'_',tmp.timestamp,"_BU.stan")

model_step_by_step_MTa_PRa = 
  stan_model(stanc_ret = stanc(glob_params$Stan$f.CodePath(tmp.stanfile.name,".stan")))

file.copy(glob_params$Stan$f.CodePath(tmp.stanfile.name,".stan"), 
          glob_params$Stan$f.FitPath(tmp.stanCodeBU.name), overwrite = T)

# Prepare stan input ------------------------------------------------------
source(glob_params$f.RScript('load-and-reshape-Konstanz-data.R'))

# Fit for Rotenone --------------------------------------------------------

tmp.chemi = "Rotenone"

tmp.data.KE1 = KE1.both %>% filter(., chemi == tmp.chemi, concentration_MuMol > 0)
tmp.data.KE2 = KE2.both %>% filter(., chemi == tmp.chemi, concentration_MuMol > 0)
tmp.data.KE3 = KE3.both %>% filter(., chemi == tmp.chemi, concentration_MuMol > 0)
tmp.data.KE4 = KE4.both %>% filter(., chemi == tmp.chemi, concentration_MuMol > 0)

stan_input.Neuro = list(
  "ndoses_CIa" = length(tmp.data.KE1$concentration_MuMol),
  "dose_CIa"   = tmp.data.KE1$concentration_MuMol,
  "CIa_CIa"    = tmp.data.KE1$c1_activity
  ,
  "ndoses_MTa" = length(tmp.data.KE2$concentration_MuMol),
  "dose_MTa"   = tmp.data.KE2$concentration_MuMol,
  "MTa_MTa"    = tmp.data.KE2$mito_resp
  ,
  "ndoses_PRa" = length(tmp.data.KE3$concentration_MuMol),
  "dose_PRa"   = tmp.data.KE3$concentration_MuMol,
  "PRa_PRa"    = tmp.data.KE3$Prot_acti
  ,
  # "ndoses_NRa" = length(tmp.data.KE4$concentration_MuMol),
  # "dose_NRa"   = tmp.data.KE4$concentration_MuMol,
  # "NRa_NRa"    = tmp.data.KE4$Neurite_Area
  # ,
  "ndoses_CIa_new" = length(tmp.data.KE1$concentration_MuMol),
  "dose_CIa_new"   = tmp.data.KE1$concentration_MuMol
  ,
  "ndoses_MTa_new" = length(tmp.data.KE2$concentration_MuMol),
  "dose_MTa_new"   = tmp.data.KE2$concentration_MuMol
  ,
  "ndoses_PRa_new" = length(tmp.data.KE3$concentration_MuMol),
  "dose_PRa_new"   = tmp.data.KE3$concentration_MuMol
  # ,
  # "ndoses_NRa_new" = length(tmp.data.KE4$concentration_MuMol),
  # "dose_NRa_new"   = tmp.data.KE4$concentration_MuMol
)

tmp.parsName = c("PRa_min",
                 "PRa_max",
                 "Par_PRa_50",
                 "q_PRa",
                 "sigma_PRa",
                 
                 "lp__")

nb.chains = 3
nb.iter = 10000
my.seed = tmp.timestamp

tmp.stanFit.name = gsub("\\+", "", 
                        glob_params$Stan$f.FitPath(tmp.stanfile.name,
                                                   "_Chemi_",tmp.chemi,
                                                   "_Chain_",nb.chains,
                                                   "_iter_",nb.iter,
                                                   "_seed_",my.seed,
                                                   "_",tmp.timestamp))


fit_step_by_step_MTa_PRa = sampling(
  model_step_by_step_MTa_PRa,
  data = stan_input.Neuro,
  pars = tmp.parsName,
  chains = nb.chains,
  # refresh = min(nb.iter/10,50),
  # sample_file = paste0(tmp.stanFit.name,".csv"),
  iter = nb.iter,
  seed = my.seed
)

write_rds(fit_step_by_step_MTa_PRa,
          path = paste0(tmp.stanFit.name, ".stanFit"))

# shinystan::launch_shinystan(fit_step_by_step_MTa_PRa)



# Analysis ----------------------------------------------------------------


{
  tmp.stanfile.name = "step_by_step_X_CIa_MTa"
  list.files(glob_params$Stan$FitPath, pattern = tmp.stanfile.name)
  fit_step_by_step_X_CIa_MTa =
    readr::read_rds(glob_params$Stan$f.FitPath("step_by_step_X_CIa_MTa_Chemi_Rotenone_Chain_3_iter_4000_seed_1907121347_1907121347.stanFit"))
  
  tmp.fit = fit_step_by_step_X_CIa_MTa
  
  tmp.parsName = c("CIa_min" ,
                   "CIa_max" ,
                   "Par_CIa_50" ,
                   "k_CIa" ,
                   "sigma_CIa",
                   
                   "beta_MTa",
                   "beta_0_MTa",
                   "sigma_MTa",
                   
                   "log_Par_CIa_50" ,
                   "CIa_delta",
                   
                   # "CIa_pred_CIa",
                   # "CIa_pred_MTa",
                   # "MTa_pred_MTa",
                   
                   "lp__")
  
  ## MAP
  la_all = as.data.frame(rstan::extract(tmp.fit, pars = tmp.parsName))
  
  la.map_all = la_all[which.max(la_all$lp__),]
  la.summary_all = summary(tmp.fit,
                           probs = c(0.05, 0.95),
                           pars = tmp.parsName)$summary %>% as.data.frame()
  la.mean_all = la.summary_all["mean"] %>% t() %>% as.data.frame()
  
  tmp.est.ponct = la.map_all
  # tmp.est.ponct = la.mean_all
  
  # Parameters of stanFit for CIa
  tmp.pars_MTa = list(beta     = tmp.est.ponct$beta_MTa,
                      beta_0   = tmp.est.ponct$beta_0_MTa)
  # Parameters of stanFit for MTa
  tmp.pars_CIa = list(y_min     = tmp.est.ponct$CIa_min,
                      y_max     = tmp.est.ponct$CIa_max,
                      x_50      = tmp.est.ponct$log_Par_CIa_50,
                      k         = tmp.est.ponct$k_CIa)
}

# cleanup()

{
  tmp.stanfile.name = "step_by_step_MTa_PRa"
  list.files(glob_params$Stan$FitPath, pattern = tmp.stanfile.name)
  tmp.fit = readr::read_rds(glob_params$Stan$f.FitPath("step_by_step_MTa_PRa_Chemi_Rotenone_Chain_3_iter_10000_seed_1907131050_1907131050.stanFit"))
  
  tmp.parsName = c("PRa_min",
                   "PRa_max",
                   "Par_PRa_50",
                   "q_PRa",
                   "sigma_PRa",
                   
                   "lp__")
  
  la_all = as.data.frame(rstan::extract(tmp.fit, pars = tmp.parsName))
  
  la.map_all = la_all[which.max(la_all$lp__),]
  tmp.est.ponct = la.map_all
  
  tmp.pars_PRa  = list(y_min     = tmp.est.ponct$PRa_min,
                       y_max     = tmp.est.ponct$PRa_max,
                       x_50      = tmp.est.ponct$Par_PRa_50,
                       k         = exp(tmp.est.ponct$q_PRa)
                       )
}
source(glob_params$f.RScript("CPDs.R")) ;

pckToLoad = c('tidyverse', 'DT',  "threejs",
              "plyr","citr","readxl",'scales',"rstan","reshape2")
reloadpck()

{
  p.data = tmp.data.KE3
  
  # dim(p.data)
  # length(p.data$CIa_pred_MAP)
  # length(p.data$MTa_pred_MAP)
  
  p.data$CIa_pred_MAP = log(p.data$concentration_MuMol) %>% 
    f.link.lgstc(.,
                 params = tmp.pars_CIa)
  
  p.data$MTa_pred_MAP = 
    f.link.linear(x = p.data$CIa_pred_MAP,
                  params = tmp.pars_MTa) %>% 
    f.link.truncation(., lower = 0, upper = 120) 
  
  ggplot(data = p.data,
         aes(x = MTa_pred_MAP,
             y = Prot_acti)) +
    geom_point(aes(colour = replication)) +
    scale_x_continuous(limit=c(0,150)) +
    scale_y_continuous(limit=c(0,150),oob=squish) +
    stat_function(fun = function(x) 
      f.link.lgstc(x,
                   params = tmp.pars_PRa))
}

if(F){
  library(bayesplot)
  mcmc_pairs(tmp.fit)
}
