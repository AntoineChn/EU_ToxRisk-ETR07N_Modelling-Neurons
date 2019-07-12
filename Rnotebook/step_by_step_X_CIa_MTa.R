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

tmp.stanfile.name = "step_by_step_X_CIa_MTa" ; paste0("model_",tmp.stanfile.name) %>% cat

tmp.stanCodeBU.name = paste0(tmp.stanfile.name,'_',tmp.timestamp,"_BU.stan")

model_step_by_step_X_CIa = 
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
  # "ndoses_PRa" = length(tmp.data.KE3$concentration_MuMol),
  # "dose_PRa"   = tmp.data.KE3$concentration_MuMol,
  # "PRa_PRa"    = tmp.data.KE3$Prot_acti
  # ,
  # "ndoses_NRa" = length(tmp.data.KE4$concentration_MuMol),
  # "dose_NRa"   = tmp.data.KE4$concentration_MuMol,
  # "NRa_NRa"    = tmp.data.KE4$Neurite_Area
  # ,
  "ndoses_CIa_new" = length(tmp.data.KE1$concentration_MuMol),
  "dose_CIa_new"   = tmp.data.KE1$concentration_MuMol,
  "ndoses_MTa_new" = length(tmp.data.KE2$concentration_MuMol),
  "dose_MTa_new"   = tmp.data.KE2$concentration_MuMol
)

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
                 
                 "CIa_pred_CIa",
                 "CIa_pred_MTa",
                 "MTa_pred_MTa",
                 
                 "lp__")

nb.chains = 3
nb.iter = 4000
my.seed = tmp.timestamp

tmp.stanFit.name = gsub("\\+", "", 
                        glob_params$Stan$f.FitPath(tmp.stanfile.name,
                                                   "_Chemi_",tmp.chemi,
                                                   "_Chain_",nb.chains,
                                                   "_iter_",nb.iter,
                                                   "_seed_",my.seed,
                                                   "_",tmp.timestamp))


fit_step_by_step_X_CIa_MTa = sampling(
  model_step_by_step_X_CIa,
  data = stan_input.Neuro,
  pars = tmp.parsName,
  chains = nb.chains,
  refresh = min(nb.iter/10,50),
  # sample_file = paste0(tmp.stanFit.name,".csv"),
  iter = nb.iter,
  seed = my.seed
)

write_rds(fit_step_by_step_X_CIa_MTa,
          path = paste0(tmp.stanFit.name, ".stanFit"))



# Fit for Deguelin -------------------------------------------------------

# tmp.chemi = "Deguelin"
# 
# tmp.data.KE1 = KE1.both %>% filter(., chemi == tmp.chemi, concentration_MuMol > 0)
# tmp.data.KE2 = KE2.both %>% filter(., chemi == tmp.chemi, concentration_MuMol > 0)
# tmp.data.KE3 = KE3.both %>% filter(., chemi == tmp.chemi, concentration_MuMol > 0)
# tmp.data.KE4 = KE4.both %>% filter(., chemi == tmp.chemi, concentration_MuMol > 0)
# 
# stan_input.Neuro = list(
#   "ndoses_CIa" = length(tmp.data.KE1$concentration_MuMol),
#   "dose_CIa"   = tmp.data.KE1$concentration_MuMol,
#   "CIa_CIa"    = tmp.data.KE1$c1_activity
#   # ,
#   # "ndoses_MTa" = length(tmp.data.KE2$concentration_MuMol),
#   # "dose_MTa"   = tmp.data.KE2$concentration_MuMol,
#   # "MTa_MTa"    = tmp.data.KE2$mito_resp
#   # ,
#   # "ndoses_PRa" = length(tmp.data.KE3$concentration_MuMol),
#   # "dose_PRa"   = tmp.data.KE3$concentration_MuMol,
#   # "PRa_PRa"    = tmp.data.KE3$Prot_acti
#   # ,
#   # "ndoses_NRa" = length(tmp.data.KE4$concentration_MuMol),
#   # "dose_NRa"   = tmp.data.KE4$concentration_MuMol,
#   # "NRa_NRa"    = tmp.data.KE4$Neurite_Area
# )
# 
# tmp.stanFit.name = gsub("\\+", "", 
#                         glob_params$Stan$f.FitPath(tmp.stanfile.name,
#                                                    "_Chemi_",tmp.chemi,
#                                                    "_Chain_",nb.chains,
#                                                    "_iter_",nb.iter,
#                                                    "_seed_",my.seed,
#                                                    "_",tmp.timestamp))
# 
# 
# fit_step_by_step_X_CIa = sampling(
#   model_step_by_step_X_CIa,
#   data = stan_input.Neuro,
#   pars = tmp.parsName,
#   chains = nb.chains,
#   refresh = min(nb.iter/10,50),
#   # sample_file = paste0(tmp.stanFit.name,".csv"),
#   iter = nb.iter,
#   seed = my.seed
# )
# 
# write_rds(fit_step_by_step_X_CIa,
#           path = paste0(tmp.stanFit.name, ".stanFit"))
# 
