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

tmp.stanfile.name = "morechemi_uptoKE2_MTa_linear" ; paste0("model_",tmp.stanfile.name) %>% cat

model_morechemi_uptoKE2_MTa_linear = 
  stan_model(stanc_ret = stanc(glob_params$Stan$f.CodePath(tmp.stanfile.name,".stan")))

# Prepare stan input ------------------------------------------------------
source(glob_params$f.RScript('load-and-reshape-Konstanz-data.R'))
source(glob_params$f.RScript('load-and-reshape-morechemi.R'))

# Fit for all --------------------------------------------------------

{# join two data sets
  tmp.data.KE1.ds1 = KE1.both %>% 
    inner_join(chemi.common.all, by = "chemi") %>% 
    filter(., concentration_MuMol > 0, c1_activity > 0) 
  tmp.data.KE1.ds2 = KE1.other.chemi %>% 
    inner_join(chemi.common.all, by = "chemi") %>% 
    filter(., concentration_MuMol > 0, c1_activity > 0)
  
  tmp.data.KE1 = bind_rows(tmp.data.KE1.ds1,
                           tmp.data.KE1.ds2) %>% 
    arrange(chemiID, concentration_MuMol)
  
  tmp.data.KE2.ds1 = KE2.both %>% 
    inner_join(chemi.common.all, by = "chemi") %>% 
    filter(., concentration_MuMol > 0) 
  tmp.data.KE2.ds2 = KE2.other.chemi %>% 
    inner_join(chemi.common.all, by = "chemi") %>% 
    filter(., concentration_MuMol > 0)
  
  tmp.data.KE2 = bind_rows(tmp.data.KE2.ds1,
                           tmp.data.KE2.ds2) %>% 
    arrange(chemiID, concentration_MuMol)
}

# Choose chemicals 
(chemiID_chosen = (1:10)[-c(5,7,10)])
chemiID_chosen = chemi.common.all %>% filter(chemiID %in% chemiID_chosen) %>% 
  mutate(chemiID_New = row_number())

# train data set
tmp.train.data.KE1 = tmp.data.KE1 %>% inner_join(chemiID_chosen, by = c("chemiID","chemi"))
tmp.train.data.KE2 = tmp.data.KE2 %>% inner_join(chemiID_chosen, by = c("chemiID","chemi"))

stan_input.Neuro = list(
  "nchemi"             = tmp.train.data.KE1$chemiID_New %>% n_distinct(),
  "ndata_CIa"          = nrow(tmp.train.data.KE1),
  "dose_CIa"           = tmp.train.data.KE1$concentration_MuMol,
  "CIa_CIa"            = tmp.train.data.KE1$c1_activity,
  "chemiID_CIa"        = tmp.train.data.KE1$chemiID_New
  ,
  "ndata_MTa"          = nrow(tmp.train.data.KE2),
  "dose_MTa"           = tmp.train.data.KE2$concentration_MuMol,
  "MTa_MTa"            = tmp.train.data.KE2$mito_resp,
  "chemiID_MTa"        = tmp.train.data.KE2$chemiID_New 
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
  "sigma_CIa"
  )

p.pars.MTa = c(
  "beta_MTa",
  "beta_0_MTa",
  "sigma_MTa"
)

tmp.parsName = c(p.pars.CIa,
                 p.pars.MTa,
                 "lp__")

nb.chains = 3
nb.iter = 5000
my.seed = tmp.timestamp

source(glob_params$f.RFunc("stanFit_Name.R"))
tmp.stanFit.name = f.stanFit.newName(chemi         = paste(chemiID_chosen$chemiID,collapse = '.'),
                                     stanfile.name = tmp.stanfile.name,
                                     timestamp     = tmp.timestamp,
                                     nb.chains     = nb.chains,
                                     nb.iter       = nb.iter,
                                     seed          = my.seed) 

fit_morechemi_uptoKE2_MTa_linear = sampling(
  model_morechemi_uptoKE2_MTa_linear,
  data = stan_input.Neuro,
  pars = tmp.parsName,
  chains = nb.chains,
  refresh = min(nb.iter/10,50),
  # sample_file = paste0(tmp.stanFit.name,".csv"),
  iter = nb.iter,
  seed = my.seed
)

# source(glob_params$f.RScript("sound_simulation_finished.R"))

# shinystan::launch_shinystan(fit_morechemi_uptoKE2_MTa_linear)

fit_morechemi_uptoKE2_MTa_linear

tmp.summary = summary(fit_morechemi_uptoKE2_MTa_linear,
                      probs = c(0.05, 0.95),
                      pars = tmp.parsName)$summary %>% 
  as.data.frame() 

tmp.n_eff.min     = tmp.summary$n_eff[!is.na(tmp.summary$n_eff)] %>% min() %>% ceiling()
tmp.rhat.max      = tmp.summary$Rhat[is.finite(tmp.summary$Rhat)] %>% max() %>% round(digits = 2)

tmp.stanFit.name = f.stanFit.newName(chemi         = paste(chemiID_chosen$chemiID,collapse = '.'),
                                     stanfile.name = tmp.stanfile.name,
                                     timestamp     = tmp.timestamp,
                                     nb.chains     = nb.chains,
                                     nb.iter       = nb.iter,
                                     seed          = my.seed,
                                     n_eff.min     = tmp.n_eff.min, 
                                     rhat.max      = tmp.rhat.max) 

write_rds(fit_morechemi_uptoKE2_MTa_linear,
          path = glob_params$Stan$f.FitPath(tmp.stanFit.name, ".stanFit"))
write_lines(get_stancode(fit_morechemi_uptoKE2_MTa_linear), 
            glob_params$Stan$f.FitPath(tmp.stanFit.name, ".stan"), append = FALSE)
file.copy(here::here("Rnotebook", "morechemi_uptoKE2_MTa_linear.R"),
          glob_params$Stan$f.FitPath(tmp.stanFit.name,".R"), overwrite = T)
 