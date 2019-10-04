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

tmp.stanfile.name = "morechemi_uptoKE1_CIa_loglogistic" ; paste0("model_",tmp.stanfile.name) %>% cat

model_morechemi_uptoKE1_CIa_loglogistic = 
  stan_model(stanc_ret = stanc(glob_params$Stan$f.CodePath(tmp.stanfile.name,".stan")))

# Prepare stan input ------------------------------------------------------
source(glob_params$f.RScript('load-and-reshape-Konstanz-data.R'))
source(glob_params$f.RScript('load-and-reshape-morechemi.R'))

# Fit for all --------------------------------------------------------

tmp.data.KE1.ds1 = KE1.both %>% 
  inner_join(chemi.common.all, by = "chemi") %>% 
  filter(., concentration_MuMol > 0, c1_activity > 0) 
tmp.data.KE1.ds2 = KE1.other.chemi %>% 
  inner_join(chemi.common.all, by = "chemi") %>% 
  filter(., concentration_MuMol > 0, c1_activity > 0) 

tmp.data.KE1 = bind_rows(tmp.data.KE1.ds1,
                         tmp.data.KE1.ds2) %>% 
  arrange(chemiID, concentration_MuMol)
tmp.data.KE1

stan_input.Neuro = list(
  "nchemi"             = tmp.data.KE1$chemi %>% n_distinct(),
  "ndata_CIa"          = nrow(tmp.data.KE1),
  "dose_CIa"           = tmp.data.KE1$concentration_MuMol,
  "CIa_CIa"            = tmp.data.KE1$c1_activity,
  "chemiID_CIa"        = tmp.data.KE1$chemiID 
)

p.pars.CIa   = c(
  "y_min_CIa_Rot",
  "y_min_CIa_Deg",
  "y_max_CIa_Rot",
  "y_max_CIa_Deg",
  "k_CIa_Rot",
  "k_CIa_Deg",
  "y_min_CIa_All",
  "y_max_CIa_All",
  "x_50_CIa_All",
  "k_CIa_All",
  "log_x_50_CIa_All",
  "sigma_CIa")

tmp.parsName = c(p.pars.CIa,
                 "lp__")

nb.chains = 3
nb.iter = 5000
my.seed = tmp.timestamp

source(glob_params$f.RFunc("stanFit_Name.R"))
tmp.stanFit.name = f.stanFit.newName(chemi         = "All",
                                     stanfile.name = tmp.stanfile.name,
                                     timestamp     = tmp.timestamp,
                                     nb.chains     = nb.chains,
                                     nb.iter       = nb.iter,
                                     seed          = my.seed) 

fit_morechemi_uptoKE1_CIa_loglogistic = sampling(
  model_morechemi_uptoKE1_CIa_loglogistic,
  data = stan_input.Neuro,
  pars = tmp.parsName,
  chains = nb.chains,
  refresh = min(nb.iter/10,50),
  # sample_file = paste0(tmp.stanFit.name,".csv"),
  iter = nb.iter,
  seed = my.seed
)

# source(glob_params$f.RScript("sound_simulation_finished.R"))

# shinystan::launch_shinystan(fit_morechemi_uptoKE1_CIa_loglogistic)

fit_morechemi_uptoKE1_CIa_loglogistic

tmp.summary = summary(fit_morechemi_uptoKE1_CIa_loglogistic,
                      probs = c(0.05, 0.95),
                      pars = tmp.parsName)$summary %>% 
  as.data.frame() 

tmp.n_eff.min     = tmp.summary$n_eff[!is.na(tmp.summary$n_eff)] %>% min() %>% ceiling()
tmp.rhat.max      = tmp.summary$Rhat[is.finite(tmp.summary$Rhat)] %>% max() %>% round(digits = 2)

tmp.stanFit.name = f.stanFit.newName(chemi         = "all",
                                     stanfile.name = tmp.stanfile.name,
                                     timestamp     = tmp.timestamp,
                                     nb.chains     = nb.chains,
                                     nb.iter       = nb.iter,
                                     seed          = my.seed,
                                     n_eff.min     = tmp.n_eff.min, 
                                     rhat.max      = tmp.rhat.max) 

write_rds(fit_morechemi_uptoKE1_CIa_loglogistic,
          path = glob_params$Stan$f.FitPath(tmp.stanFit.name, ".stanFit"))
write_lines(get_stancode(fit_morechemi_uptoKE1_CIa_loglogistic), 
            glob_params$Stan$f.FitPath(tmp.stanFit.name, ".stan"), append = FALSE)
file.copy(here::here("Rnotebook", "morechemi_uptoKE1_CIa_loglogistic.R"),
          glob_params$Stan$f.FitPath(tmp.stanFit.name,".R"), overwrite = T)
