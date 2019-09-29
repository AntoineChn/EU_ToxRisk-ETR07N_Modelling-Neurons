# load-and-reshape-morechemi.R


# set up ------------------------------------------------------------------

if(!exists("reloadpck")) source('GlobalParameters.R')
if(exists("pckToLoad")){
  pckToLoad = unique(c(pckToLoad, 'tidyverse','rstan','magrittr','readxl'))
}else{
  pckToLoad = c('tidyverse','readxl')
}

reloadpck()

# KE1 ---------------------------------------------------------------------

{# KE1
  data.file.name = "2 KE1 other CS4 substances c1 activity glucose.xlsx"
  
  sheet.names = excel_sheets(glob_params$DataPath %>% paste0(data.file.name))
  
  KE1.other.chemi = 
    map(sheet.names, function(sn){
      read_xlsx(glob_params$DataPath %>% paste0(data.file.name),
                skip = 1,
                sheet = sn) %>%
        dplyr::rename(concentration_MuMol = `concentration [µM]`) %>% 
        filter(concentration_MuMol != "DMSO control") %>% 
        select(c(starts_with("concentration"),starts_with("c1 inhibition"))) %>% 
        gather(.,
               starts_with("c1 inhibition"),
               key = "replication", value = "c1_inhibition") %>% 
        drop_na(., "c1_inhibition") %>% 
        mutate(chemi = sn,
               c1_activity = 100. - c1_inhibition) %>% 
        select(chemi, concentration_MuMol, replication, c1_activity) %>% 
        return()
    }) %>% bind_rows() %>% 
    mutate(concentration_MuMol = as.numeric(concentration_MuMol))

  rm(data.file.name, sheet.names)
}


# KE2 ---------------------------------------------------------------------

{#KE2
  data.file.name = "2 KE2 other CS4 substances mitochondrial respiration glucose.xlsx"
  
  sheet.names = excel_sheets(glob_params$DataPath %>% paste0(data.file.name))

  KE2.other.chemi = 
    map(sheet.names, function(sn){
      read_xlsx(glob_params$DataPath %>% paste0(data.file.name),
                range = cell_rows(c(2,5)),
                sheet = sn) %>%
        dplyr::rename(concentration_MuMol = `concentration [µM]`) %>% 
        filter(concentration_MuMol != "DMSO control") %>% 
        select(c(starts_with("concentration"),starts_with("total"), starts_with("non-mito")) ) %>%
        gather(.,
               c(starts_with("total"), starts_with("non-mito")),
               key = "replication", value = "resp") %>% 
        drop_na(., resp) %>% 
        mutate(., repID = str_sub(replication, -1),
               replication = str_sub(replication, start = 1, 
                                     end = nchar(replication) - 5)) %>% 
        mutate(replication = gsub("\\s+", " ", str_trim(replication))) %>% 
        spread(., key = replication, value = resp) %>% 
        mutate(replication = paste0("mito. resp., N", repID),
               mito_resp = `total resp` - `non-mito. resp`) %>% 
        mutate(chemi = sn) %>% 
        select(chemi, concentration_MuMol, replication, mito_resp) %>% 
        return()
    }) %>% bind_rows() %>% 
    mutate(concentration_MuMol = as.numeric(concentration_MuMol))
  
  rm(data.file.name, sheet.names)
}

# KE3 ---------------------------------------------------------------------

{#KE3
  data.file.name = "2 KE3 other CS4 substances proteasomal activity glucose.xlsx"
  sheet.names = excel_sheets(glob_params$DataPath %>% paste0(data.file.name))
  
  KE3.other.chemi = 
    map(sheet.names, function(sn){
      read_xlsx(glob_params$DataPath %>% paste0(data.file.name),
                skip = 1,
                sheet = sn) %>%
        dplyr::rename(concentration_MuMol = `concentration [µM]`) %>% 
        filter(concentration_MuMol != "DMSO control") %>% 
        select(c(starts_with("concentration"), starts_with("Prot")) ) %>% 
        gather(.,
               starts_with("Prot"),
               key = "replication", value = "Prot_acti") %>% 
        drop_na(., "Prot_acti") %>% 
        mutate(chemi = sn) %>% 
        select(chemi, concentration_MuMol, replication, Prot_acti) %>% 
        return()
    }) %>% bind_rows() %>% 
    mutate(concentration_MuMol = as.numeric(concentration_MuMol))
  rm(data.file.name, sheet.names)
}

# KE4 ---------------------------------------------------------------------


{#KE4
  data.file.name = "2 KE4 other CS4 substances NeuriTox glucose.xlsx"
  
  sheet.names = excel_sheets(glob_params$DataPath %>% paste0(data.file.name))
  
  sn = sheet.names[1]
  
  KE4.other.chemi = 
    map(sheet.names, function(sn){
      read_xlsx(glob_params$DataPath %>% paste0(data.file.name),
                skip = 1,
                sheet = sn) %>%
        dplyr::rename( concentration_MuMol = `concentration [µM]`) %>% 
        filter(concentration_MuMol != "DMSO control") %>% 
        select(c(starts_with("concentration"), starts_with("Neurite")) ) %>% 
        gather(.,
               starts_with("Neurite"),
               key = "replication", value = "Neurite_Area") %>% 
        drop_na(., "Neurite_Area") %>% 
        mutate(chemi = sn) %>% 
        select(chemi, concentration_MuMol, replication, Neurite_Area) %>% 
        return()
    }) %>% bind_rows() %>% 
    mutate(concentration_MuMol = as.numeric(concentration_MuMol))
  
  rm(data.file.name, sheet.names)
}


# common KE ---------------------------------------------------------------
# KE1.other.chemi %>% DT::datatable()
# KE2.other.chemi %>% DT::datatable()
# KE3.other.chemi %>% DT::datatable()
# KE4.other.chemi %>% DT::datatable()


chemi.common =
  reduce(list(KE1.other.chemi$chemi,
              KE2.other.chemi$chemi,
              KE3.other.chemi$chemi,
              KE4.other.chemi$chemi),
         intersect) %>%
  tibble(chemi = .) %>%
  rowid_to_column(var = "chemiID")

chemi.common.all =
  c("Rotenone",
    "Deguelin",
    reduce(
      list(
        KE1.other.chemi$chemi,
        KE2.other.chemi$chemi,
        KE3.other.chemi$chemi,
        KE4.other.chemi$chemi
      ),
      intersect
    )) %>%
  tibble(chemi = .) %>%
  rowid_to_column(var = "chemiID")

common.chemi.uptoKE1 =
  c("Rotenone",
    "Deguelin",
    KE1.other.chemi$chemi %>% unique())

common.chemi.uptoKE2 =
  c("Rotenone",
    "Deguelin",
    reduce(list(common.chemi.uptoKE1,
                KE2.other.chemi$chemi),
           intersect)
    ) 

common.chemi.uptoKE3 =
  c("Rotenone",
    "Deguelin",
    reduce(list(common.chemi.uptoKE2,
                KE3.other.chemi$chemi),
           intersect)
  ) 

common.chemi.uptoKE4 =
  c("Rotenone",
    "Deguelin",
    reduce(list(common.chemi.uptoKE3,
                KE4.other.chemi$chemi),
           intersect)
  )

# list(KE1.other.chemi,
#      KE2.other.chemi,
#      KE3.other.chemi,
#      KE4.other.chemi) %>% 
#   map(., function(ds){
#     ds %>% 
#       filter(chemi %in% chemi.common$chemi) %>% 
#       return()
#   })


