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

tmp.stanfile.name = "step_by_step_X_CIa_MTa_PRa_NRa" ; paste0("model_",tmp.stanfile.name) %>% cat
tmp.stanCodeBU.name = paste0(tmp.stanfile.name,'_',tmp.timestamp,"_BU.stan")

model_step_by_step_X_CIa_MTa_PRa_NRa = 
  stan_model(stanc_ret = stanc(glob_params$Stan$f.CodePath(tmp.stanfile.name,".stan")))

# Prepare stan input ------------------------------------------------------
source(glob_params$f.RScript('load-and-reshape-Konstanz-data.R'))

# Fit for Rotenone --------------------------------------------------------

tmp.chemi = "Deguelin"

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
  "ndoses_NRa" = length(tmp.data.KE4$concentration_MuMol),
  "dose_NRa"   = tmp.data.KE4$concentration_MuMol,
  "NRa_NRa"    = tmp.data.KE4$Neurite_Area
)

tmp.parsName = c("CIa_min" ,
                 "CIa_max" ,
                 "Par_CIa_50" ,
                 "log_Par_CIa_50",
                 "k_CIa" ,
                 "sigma_CIa",
                 
                 "beta_MTa",
                 "beta_0_MTa",
                 "sigma_MTa",
                 
                 "PRa_min",
                 "PRa_max",
                 "Par_PRa_50",
                 "q_PRa",
                 "sigma_PRa",
                 
                 "k_NRa",
                 "NRa_max",
                 "MTa_tox_min_NRa",
                 "sigma_NRa",

                 "lp__")

nb.chains = 3
nb.iter = 8000
my.seed = tmp.timestamp

source(glob_params$f.RFunc("stanFit_Name.R"))
tmp.stanFit.name = f.stanFit.newName(chemi         = tmp.chemi,
                                     stanfile.name = tmp.stanfile.name,
                                     timestamp     = tmp.timestamp,
                                     nb.chains     = nb.chains,
                                     nb.iter       = nb.iter,
                                     seed          = my.seed)

fit_step_by_step_X_CIa_MTa_PRa_NRa = sampling(
  model_step_by_step_X_CIa_MTa_PRa_NRa,
  data = stan_input.Neuro,
  pars = tmp.parsName,
  chains = nb.chains,
  refresh = min(nb.iter/10,50),
  # sample_file = paste0(tmp.stanFit.name,".csv"),
  iter = nb.iter,
  seed = my.seed
)

source(glob_params$f.RScript("sound_simulation_finished.R"))

write_rds(fit_step_by_step_X_CIa_MTa_PRa_NRa,
          path = glob_params$Stan$f.FitPath(tmp.stanFit.name, ".stanFit"))
file.copy(glob_params$Stan$f.CodePath(tmp.stanfile.name,".stan"), 
          glob_params$Stan$f.FitPath(tmp.stanCodeBU.name), overwrite = T)
