# Load-and-reshape-Konstanz-data.R


# KE1 ---------------------------------------------------------------------

if(exists("pckToLoad")){
  pckToLoad = unique(c(pckToLoad, 'tidyverse','rstan','magrittr','readxl'))
}else{
  pckToLoad = c('tidyverse','rstan','magrittr','readxl')
}
  
reloadpck()

{# KE1
  data.file.name = "KE1 rotenone deguelin c1 activity glucose.xlsx"
  KE1.Rotenone = 
    read_xlsx(glob_params$DataPath %>% paste0(data.file.name),
              sheet = "Rotenone")
  
  KE1.Deguelin =   read_xlsx(glob_params$DataPath %>% paste0(data.file.name),
                             sheet = "Deguelin")
  {
    KE1.Rotenone.rs = fill(KE1.Rotenone, `concentration [µM]`) %>% 
      gather(.,
             starts_with("c1 activity"),
             key = "replication", value = "c1_activity") %>% 
      drop_na(., "c1_activity")
    
    names(KE1.Rotenone.rs)[1] = "concentration_MuMol"
    KE1.Rotenone.rs$concentration_MuMol [KE1.Rotenone.rs$concentration_MuMol == "DMSO control"] = 0
    
    KE1.Rotenone.rs$concentration_MuMol = as.numeric(KE1.Rotenone.rs$concentration_MuMol)
    
    KE1.Rotenone.rs = arrange(KE1.Rotenone.rs,concentration_MuMol)
  }
  {
    KE1.Deguelin
    KE1.Deguelin.rs = fill(KE1.Deguelin, `concentration [µM]`) %>% 
      gather(.,
             starts_with("c1 activity"),
             key = "replication", value = "c1_activity") %>% 
      drop_na(., "c1_activity")
    
    names(KE1.Deguelin.rs)[1] = "concentration_MuMol"
    KE1.Deguelin.rs$concentration_MuMol [KE1.Deguelin.rs$concentration_MuMol == "DMSO control"] = 0
    
    KE1.Deguelin.rs$concentration_MuMol = as.numeric(KE1.Deguelin.rs$concentration_MuMol)
    
    KE1.Deguelin.rs = arrange(KE1.Deguelin.rs,concentration_MuMol)
  }
  
  KE1.Rotenone.rs$chemi = "Rotenone"
  KE1.Deguelin.rs$chemi = "Deguelin"
  
  KE1.both = 
    bind_rows(KE1.Rotenone.rs,
              KE1.Deguelin.rs)
  rm(data.file.name)
}


# KE2 ---------------------------------------------------------------------


{#KE2
  data.file.name = "KE2 rotenone deguelin mitochondrial respiration glucose.xlsx"
  KE2.Rotenone = 
    read_xlsx(glob_params$DataPath %>% paste0(data.file.name),
              sheet = "Rotenone")
  
  KE2.Deguelin =   read_xlsx(glob_params$DataPath %>% paste0(data.file.name),
                             sheet = "Deguelin")
  
  {
    KE2.Rotenone.rs = fill(KE2.Rotenone, `concentration [µM]`) %>% 
      gather(.,
             starts_with("mito. resp."),
             key = "replication", value = "mito_resp") %>% 
      drop_na(., "mito_resp")
    
    names(KE2.Rotenone.rs)[1] = "concentration_MuMol"
    KE2.Rotenone.rs$concentration_MuMol [KE2.Rotenone.rs$concentration_MuMol == "DMSO control"] = 0
    
    KE2.Rotenone.rs$concentration_MuMol = as.numeric(KE2.Rotenone.rs$concentration_MuMol)
    
    KE2.Rotenone.rs = arrange(KE2.Rotenone.rs,concentration_MuMol)
  }
  
  {
    KE2.Deguelin.rs = fill(KE2.Deguelin, `concentration [µM]`) %>% 
      gather(.,
             starts_with("mito. resp."),
             key = "replication", value = "mito_resp") %>% 
      drop_na(., "mito_resp")
    
    names(KE2.Deguelin.rs)[1] = "concentration_MuMol"
    KE2.Deguelin.rs$concentration_MuMol [KE2.Deguelin.rs$concentration_MuMol == "DMSO control"] = 0
    
    parse_double(KE2.Deguelin.rs$concentration_MuMol)
    
    KE2.Deguelin.rs$concentration_MuMol = as.numeric(KE2.Deguelin.rs$concentration_MuMol)
    
    KE2.Deguelin.rs = arrange(KE2.Deguelin.rs,concentration_MuMol)
  }
  
  KE2.Rotenone.rs$chemi = "Rotenone"
  KE2.Deguelin.rs$chemi = "Deguelin"
  
  KE2.both = 
    bind_rows(KE2.Rotenone.rs,
              KE2.Deguelin.rs)
  rm(data.file.name)
}


# KE3 ---------------------------------------------------------------------


{#KE3
  data.file.name = "KE3 rotenone deguelin Proteasomal activity glucose.xlsx"
  KE3.Rotenone = 
    read_xlsx(glob_params$DataPath %>% paste0(data.file.name),
              sheet = "Rotenone")
  
  KE3.Deguelin =   read_xlsx(glob_params$DataPath %>% paste0(data.file.name),
                             sheet = "Deguelin")
  
  {
    KE3.Rotenone
    KE3.Rotenone.rs = fill(KE3.Rotenone, `concentration [µM]`) %>% 
      gather(.,
             starts_with("Prot. acti."),
             key = "replication", value = "Prot_acti") %>% 
      drop_na(., "Prot_acti")
    
    names(KE3.Rotenone.rs)[1] = "concentration_MuMol"
    KE3.Rotenone.rs$concentration_MuMol [KE3.Rotenone.rs$concentration_MuMol == "DMSO control"] = 0
    
    KE3.Rotenone.rs$concentration_MuMol = as.numeric(KE3.Rotenone.rs$concentration_MuMol)
    
    KE3.Rotenone.rs = arrange(KE3.Rotenone.rs,concentration_MuMol)
  }
  {
    KE3.Deguelin.rs = fill(KE3.Deguelin, `concentration [µM]`) %>% 
      gather(.,
             starts_with("Prot. acti."),
             key = "replication", value = "Prot_acti") %>% 
      drop_na(., "Prot_acti")
    
    names(KE3.Deguelin.rs)[1] = "concentration_MuMol"
    KE3.Deguelin.rs$concentration_MuMol [KE3.Deguelin.rs$concentration_MuMol == "DMSO control"] = 0
    
    KE3.Deguelin.rs$concentration_MuMol = as.numeric(KE3.Deguelin.rs$concentration_MuMol)
    
    KE3.Deguelin.rs = arrange(KE3.Deguelin.rs,concentration_MuMol)
  }
  
  KE3.Rotenone.rs$chemi = "Rotenone"
  KE3.Deguelin.rs$chemi = "Deguelin"
  
  KE3.both = 
    bind_rows(KE3.Rotenone.rs,
              KE3.Deguelin.rs)
  rm(data.file.name)
}



# KE4 ---------------------------------------------------------------------


{#KE4
  data.file.name = "Example data set rotenone deguelin UKN4 NeuriTox glucose.xlsx"
  KE4.Rotenone = 
    read_xlsx(glob_params$DataPath %>% paste0(data.file.name),
              sheet = "Rotenone")
  
  KE4.Deguelin =   read_xlsx(glob_params$DataPath %>% paste0(data.file.name),
                             sheet = "Deguelin")
  {
    KE4.Rotenone.rs = fill(KE4.Rotenone, `concentration [µM]`) %>% 
      select(.,`concentration [µM]`,starts_with("Neurite Area")) %>% 
      gather(.,
             starts_with("Neurite Area"),
             key = "replication", value = "Neurite_Area") %>% 
      drop_na(., "Neurite_Area")
    
    names(KE4.Rotenone.rs)[1] = "concentration_MuMol"
    KE4.Rotenone.rs$concentration_MuMol [KE4.Rotenone.rs$concentration_MuMol == "DMSO control"] = 0
    
    KE4.Rotenone.rs$concentration_MuMol = as.numeric(KE4.Rotenone.rs$concentration_MuMol)
    
    KE4.Rotenone.rs = arrange(KE4.Rotenone.rs,concentration_MuMol)
  }
  {
    KE4.Deguelin
    KE4.Deguelin.rs = fill(KE4.Deguelin, `concentration [µM]`) %>% 
      select(.,`concentration [µM]`,starts_with("Neurite Area")) %>% 
      gather(.,
             starts_with("Neurite Area"),
             key = "replication", value = "Neurite_Area") %>% 
      drop_na(., "Neurite_Area")
    
    names(KE4.Deguelin.rs)[1] = "concentration_MuMol"
    KE4.Deguelin.rs$concentration_MuMol [KE4.Deguelin.rs$concentration_MuMol == "DMSO control"] = 0
    
    KE4.Deguelin.rs$concentration_MuMol = as.numeric(KE4.Deguelin.rs$concentration_MuMol)
    
    KE4.Deguelin.rs = arrange(KE4.Deguelin.rs,concentration_MuMol)
  }
  
  KE4.Rotenone.rs$chemi = "Rotenone"
  KE4.Deguelin.rs$chemi = "Deguelin"
  
  KE4.both = 
    bind_rows(KE4.Rotenone.rs,
              KE4.Deguelin.rs)
  rm(data.file.name)
}

rm(list = ls(pattern = 'KE\\d.Rotenone') )
rm(list = ls(pattern = 'KE\\d.Deguelin') )
