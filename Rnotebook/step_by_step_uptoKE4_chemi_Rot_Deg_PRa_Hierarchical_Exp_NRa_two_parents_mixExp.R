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

tmp.stanfile.name = "step_by_step_uptoKE4_chemi_Rot_Deg_PRa_Hierarchical_Exp_NRa_two_parents_mixExp" ; paste0("model_",tmp.stanfile.name) %>% cat

model_step_by_step_uptoKE4_chemi_Rot_Deg_PRa_Hierarchical_Exp_NRa_two_parents_mixExp = 
  stan_model(stanc_ret = stanc(glob_params$Stan$f.CodePath(tmp.stanfile.name,".stan")))

# Prepare stan input ------------------------------------------------------
source(glob_params$f.RScript('load-and-reshape-Konstanz-data.R'))
source(glob_params$f.RScript('load-and-reshape-morechemi.R'))

{# joint All data concerning KE1
  tmp.data.KE1.ds1 = KE1.both %>%
    inner_join(chemi.common.all, by = "chemi") %>%
    filter(., concentration_MuMol > 0, c1_activity > 0)
  tmp.data.KE1.ds2 = KE1.other.chemi %>%
    inner_join(chemi.common.all, by = "chemi") %>%
    filter(., concentration_MuMol > 0, c1_activity > 0)
  
  tmp.data.KE1 = bind_rows(tmp.data.KE1.ds1,
                           tmp.data.KE1.ds2) %>%
    arrange(chemiID, concentration_MuMol)
}

{# joint All data concerning KE2
  tmp.data.KE2.ds1 = KE2.both %>%
    inner_join(chemi.common.all, by = "chemi") %>%
    filter(., concentration_MuMol > 0, mito_resp > 0)
  tmp.data.KE2.ds2 = KE2.other.chemi %>%
    inner_join(chemi.common.all, by = "chemi") %>%
    filter(., concentration_MuMol > 0, mito_resp > 0)
  
  tmp.data.KE2 = bind_rows(tmp.data.KE2.ds1,
                           tmp.data.KE2.ds2) %>%
    arrange(chemiID, concentration_MuMol)
}

{# joint All data concerning KE3
  tmp.data.KE3.ds1 = KE3.both %>%
    inner_join(chemi.common.all, by = "chemi") %>%
    filter(., concentration_MuMol > 0, Prot_acti > 0)
  tmp.data.KE3.ds2 = KE3.other.chemi %>%
    inner_join(chemi.common.all, by = "chemi") %>%
    filter(., concentration_MuMol > 0, Prot_acti > 0)
  
  tmp.data.KE3 = bind_rows(tmp.data.KE3.ds1,
                           tmp.data.KE3.ds2) %>%
    arrange(chemiID, concentration_MuMol)
}

{# joint All data concerning KE4
  tmp.data.KE4.ds1 = KE4.both %>%
    inner_join(chemi.common.all, by = "chemi") %>%
    filter(., concentration_MuMol > 0, Neurite_Area > 0)
  tmp.data.KE4.ds2 = KE4.other.chemi %>%
    inner_join(chemi.common.all, by = "chemi") %>%
    filter(., concentration_MuMol > 0, Neurite_Area > 0)
  
  tmp.data.KE4 = bind_rows(tmp.data.KE4.ds1,
                           tmp.data.KE4.ds2) %>%
    arrange(chemiID, concentration_MuMol)
}

# Choose chemicals 
(chemiID_chosen = (1:2))
chemiID_chosen = chemi.common.all %>% filter(chemiID %in% chemiID_chosen) %>% 
  mutate(chemiID_Param = row_number())

# train data set
{
  tmp.train.data.KE1 = tmp.data.KE1 %>% inner_join(chemiID_chosen, by = c("chemiID","chemi"))
  tmp.train.data.KE2 = tmp.data.KE2 %>% inner_join(chemiID_chosen, by = c("chemiID","chemi"))
  tmp.train.data.KE3 = tmp.data.KE3 %>% inner_join(chemiID_chosen, by = c("chemiID","chemi"))
  tmp.train.data.KE4 = tmp.data.KE4 %>% inner_join(chemiID_chosen, by = c("chemiID","chemi"))
  
  tmp.train.data.KE3 = # add replication ID
    tmp.train.data.KE3 %>% 
    mutate(repID = str_sub(replication, -1) %>% as.numeric() )
  tmp.repID.KE3 = # creat new replication ID
    tmp.train.data.KE3 %>% 
    distinct(chemiID, repID) %>% 
    mutate(repID_new = (1:dim(.)[1])[row_number()] )
  
  tmp.train.data.KE3 = # joint new replication ID back to original dataset
    tmp.train.data.KE3 %>% right_join(tmp.repID.KE3)
  rm(tmp.repID.KE3)
  
  tmp.train.data.KE4 = # add replication ID
    tmp.train.data.KE4 %>% 
    mutate(repID = str_sub(replication, -1) %>% as.numeric() )
  tmp.repID.KE4 = # creat new replication ID
    tmp.train.data.KE4 %>% 
    distinct(chemiID, repID) %>% 
    mutate(repID_new = (1:dim(.)[1])[row_number()] )
  
  tmp.train.data.KE4 = # joint new replication ID back to original dataset
    tmp.train.data.KE4 %>% right_join(tmp.repID.KE4)
  rm(tmp.repID.KE4)
}



# Fit for Rotenone --------------------------------------------------------

stan_input.Neuro = list(
  "nchemi"             = tmp.train.data.KE1$chemiID_Param %>% n_distinct(),
  "ndata_CIa"          = nrow(tmp.train.data.KE1),
  "dose_CIa"           = tmp.train.data.KE1$concentration_MuMol,
  "CIa_CIa"            = tmp.train.data.KE1$c1_activity,
  "chemi_id_CIa"        = tmp.train.data.KE1$chemiID_Param
  ,
  "ndata_MTa"          = nrow(tmp.train.data.KE2),
  "dose_MTa"           = tmp.train.data.KE2$concentration_MuMol,
  "MTa_MTa"            = tmp.train.data.KE2$mito_resp,
  "chemi_id_MTa"        = tmp.train.data.KE2$chemiID_Param 
  ,
  "ndata_PRa"          = nrow(tmp.train.data.KE3),
  "nrep_PRa"           = tmp.train.data.KE3$repID_new %>% length(),
  "dose_PRa"           = tmp.train.data.KE3$concentration_MuMol,
  "PRa_PRa"            = tmp.train.data.KE3$Prot_acti,
  "rep_id_PRa"         = tmp.train.data.KE3$repID_new,
  "chemi_id_PRa"       = tmp.train.data.KE3$chemiID_Param
  ,
  "ndata_NRa"          = nrow(tmp.train.data.KE4),
  "nrep_NRa"           = tmp.train.data.KE4$repID_new %>% length(),
  "dose_NRa"           = tmp.train.data.KE4$concentration_MuMol,
  "NRa_NRa"            = tmp.train.data.KE4$Neurite_Area,
  # "rep_id_NRa"         = tmp.train.data.KE4$repID_new,
  "chemi_id_NRa"       = tmp.train.data.KE4$chemiID_Param
)

p.pars.CIa = c("y_min_CIa_Rot",
               "y_max_CIa_Rot",
               "k_CIa_Rot",
               "y_min_CIa_Deg",
               "y_max_CIa_Deg",
               "k_CIa_Deg",
               "x_50_CIa_All",
               "sigma_CIa")

p.pars.MTa = c("beta_MTa",
               "beta_0_MTa",
               "sigma_MTa")

p.pars.PRa = c("mean_y_max_PRa",
               "sd_y_max_PRa",
               "y_max_PRa",
               "x_toxMin_PRa",
               "k_PRa",
               "sigma_PRa")

p.pars.NRa = c("x_toxMin_MTa_NRa",
               "x_toxMin_PRa_NRa",
               "y_max_frac_NRa",
               "k_MTa_NRa",     
               "k_PRa_NRa",     
               "sigma_max_NRa")

tmp.parsName = c(p.pars.CIa,
                 p.pars.MTa,
                 p.pars.PRa,
                 p.pars.NRa,
                 
                 "lp__")

nb.chains = 1
nb.iter = 5000
my.seed = tmp.timestamp

source(glob_params$f.RFunc("stanFit_Name.R"))
tmp.stanFit.name = f.stanFit.newName(chemi         = paste(chemiID_chosen$chemiID,collapse = '.'),
                                     stanfile.name = tmp.stanfile.name,
                                     timestamp     = tmp.timestamp,
                                     nb.chains     = nb.chains,
                                     nb.iter       = nb.iter,
                                     seed          = my.seed) 

fit_step_by_step_uptoKE4_chemi_Rot_Deg_PRa_Hierarchical_Exp_NRa_two_parents_mixExp = sampling(
  model_step_by_step_uptoKE4_chemi_Rot_Deg_PRa_Hierarchical_Exp_NRa_two_parents_mixExp,
  data = stan_input.Neuro,
  pars = tmp.parsName,
  chains = nb.chains,
  refresh = min(nb.iter/10,50),
  # sample_file = paste0(tmp.stanFit.name,".csv"),
  iter = nb.iter,
  seed = my.seed
)

# shinystan::launch_shinystan(fit_step_by_step_uptoKE4_chemi_Rot_Deg_PRa_Hierarchical_Exp_NRa_two_parents_mixExp)

# source(glob_params$f.RScript("sound_simulation_finished.R"))

fit_step_by_step_uptoKE4_chemi_Rot_Deg_PRa_Hierarchical_Exp_NRa_two_parents_mixExp

tmp.summary = summary(fit_step_by_step_uptoKE4_chemi_Rot_Deg_PRa_Hierarchical_Exp_NRa_two_parents_mixExp,
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

write_rds(fit_step_by_step_uptoKE4_chemi_Rot_Deg_PRa_Hierarchical_Exp_NRa_two_parents_mixExp,
          path = glob_params$Stan$f.FitPath(tmp.stanFit.name, ".stanFit"))
write_lines(get_stancode(fit_step_by_step_uptoKE4_chemi_Rot_Deg_PRa_Hierarchical_Exp_NRa_two_parents_mixExp), 
            glob_params$Stan$f.FitPath(tmp.stanFit.name, ".stan"), append = FALSE)
file.copy(here::here("Rnotebook", "step_by_step_uptoKE4_chemi_Rot_Deg_PRa_Hierarchical_Exp_NRa_two_parents_mixExp.R"),
          glob_params$Stan$f.FitPath(tmp.stanFit.name,".R"), overwrite = T)

## Posterior summary -------------------------------------------------------

if(F){
  tmp.fit = fit_step_by_step_uptoKE4_chemi_Rot_Deg_PRa_Hierarchical_Exp_NRa_two_parents_mixExp
  
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

if(exists("fit_step_by_step_uptoKE4_chemi_Rot_Deg_PRa_Hierarchical_Exp_NRa_two_parents_mixExp")) rm(fit_step_by_step_uptoKE4_chemi_Rot_Deg_PRa_Hierarchical_Exp_NRa_two_parents_mixExp)

tmp.chemi = "Deguelin"

tmp.data.KE1 = KE1.both %>% filter(., chemi == tmp.chemi, concentration_MuMol > 0)
tmp.data.KE2 = KE2.both %>% filter(., chemi == tmp.chemi, concentration_MuMol > 0)
tmp.data.KE3 = KE3.both %>% filter(., chemi == tmp.chemi, concentration_MuMol > 0)
tmp.data.KE4 = KE4.both %>% filter(., chemi == tmp.chemi, concentration_MuMol > 0)

stan_input.Neuro = list(
  "is_rot_or_deg" = as.integer(tmp.chemi == "Rotenone"),
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
  "PRa_PRa"    = tmp.data.KE3$Prot_acti,
  "nrep_PRa"   = tmp.data.KE3$replication %>% str_sub(., -1) %>% unique() %>% length(),
  "rep_ID_PRa" = tmp.data.KE3$replication %>% str_sub(., -1) %>% as.numeric()
  ,
  "ndoses_NRa" = length(tmp.data.KE4$concentration_MuMol),
  "dose_NRa"   = tmp.data.KE4$concentration_MuMol,
  "NRa_NRa"    = tmp.data.KE4$Neurite_Area,
  "nrep_NRa"   = tmp.data.KE4$replication %>% str_sub(., -1) %>% unique() %>% length(),
  "rep_ID_NRa" = tmp.data.KE4$replication %>% str_sub(., -1) %>% as.numeric()
)


tmp.stanFit.name = f.stanFit.newName(chemi         = tmp.chemi,
                                     stanfile.name = tmp.stanfile.name,
                                     timestamp     = tmp.timestamp,
                                     nb.chains     = nb.chains,
                                     nb.iter       = nb.iter,
                                     seed          = my.seed) 

fit_step_by_step_uptoKE4_chemi_Rot_Deg_PRa_Hierarchical_Exp_NRa_two_parents_mixExp = sampling(
  model_step_by_step_uptoKE4_chemi_Rot_Deg_PRa_Hierarchical_Exp_NRa_two_parents_mixExp,
  data = stan_input.Neuro,
  pars = tmp.parsName,
  chains = nb.chains,
  refresh = min(nb.iter/10,50),
  # sample_file = paste0(tmp.stanFit.name,".csv"),
  iter = nb.iter,
  seed = my.seed
)

# shinystan::launch_shinystan(fit_step_by_step_uptoKE4_chemi_Rot_Deg_PRa_Hierarchical_Exp_NRa_two_parents_mixExp)

source(glob_params$f.RScript("sound_simulation_finished.R"))

fit_step_by_step_uptoKE4_chemi_Rot_Deg_PRa_Hierarchical_Exp_NRa_two_parents_mixExp

tmp.summary = summary(fit_step_by_step_uptoKE4_chemi_Rot_Deg_PRa_Hierarchical_Exp_NRa_two_parents_mixExp,
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

write_rds(fit_step_by_step_uptoKE4_chemi_Rot_Deg_PRa_Hierarchical_Exp_NRa_two_parents_mixExp,
          path = glob_params$Stan$f.FitPath(tmp.stanFit.name, ".stanFit"))
write_lines(get_stancode(fit_step_by_step_uptoKE4_chemi_Rot_Deg_PRa_Hierarchical_Exp_NRa_two_parents_mixExp), 
            glob_params$Stan$f.FitPath(tmp.stanFit.name, ".stan"), append = FALSE)
file.copy(here("Rnotebook", "step_by_step_uptoKE4_chemi_Rot_Deg_PRa_Hierarchical_Exp_NRa_two_parents_mixExp.R"),
          glob_params$Stan$f.FitPath(tmp.stanFit.name,".R"), overwrite = T)

## Posterior summary -------------------------------------------------------

if(F){
  tmp.fit = fit_step_by_step_uptoKE4_chemi_Rot_Deg_PRa_Hierarchical_Exp_NRa_two_parents_mixExp
  
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

