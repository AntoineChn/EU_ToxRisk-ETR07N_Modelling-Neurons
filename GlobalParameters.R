if(!require(here)) install.packages("here")

# Global parameters
glob_params = NULL
glob_params$rootPath = here::here()
glob_params$RPath = paste0(glob_params$rootPath,"/R/")
glob_params$RFunc = paste0(glob_params$RPath, "function/")
glob_params$PlotPath = paste0(glob_params$rootPath,"plot/")
glob_params$RScript = paste0(glob_params$RPath, "script/")
glob_params$f.RScript = function(..., list = F) {
  file = .Internal(paste(list(...),sep = "",collapse = ""))
  res = paste0(glob_params$RScript, file)
  if(file != ""){
    return(res)
  }else{
    cat("The file name is not specified\n")
    if(list) print(list.files(res))
    return(glob_params$RScript)
  }
}

glob_params$f.RFunc = function(..., list = F) {
  file = .Internal(paste(list(...),sep = "",collapse = ""))
  res = paste0(glob_params$RFunc, file)
  if(file != ""){
    return(res)
  }else{
    cat("The file name is not specified\n")
    if(list) print(list.files(res))
    return(glob_params$RFunc)
  }
}

glob_params$MCSim_DataImage = paste0(glob_params$DataPath, "/MCSim_DataImage/")

glob_params$Stan = NULL
glob_params$Stan$CodePath = paste0(glob_params$rootPath,"/stan/code/")
glob_params$Stan$f.CodePath = function(..., list = F) {
  file = .Internal(paste(list(...),sep = "",collapse = ""))
  res = paste0(glob_params$Stan$CodePath, file)
  if(file != ""){
    return(res)
  }else{
    cat("The file name is not specified\n")
    if(list) print(list.files(res))
    return(here::here("/stan/code/"))
  }
}


glob_params$Stan$ProgPath = paste0(glob_params$rootPath,"/stan/prog/") # stan program after compilation
glob_params$Stan$f.ProgPath = function(..., list = F) {
  file = .Internal(paste(list(...),sep = "",collapse = ""))
  res = paste0(glob_params$Stan$ProgPath, file)
  if(file != ""){
    return(res)
  }else{
    cat("The file name is not specified\n")
    if(list) print(list.files(res))
    return(here::here("stan/prog/"))
  }
}

glob_params$Stan$FitPath  = paste0(glob_params$rootPath,"/stan/fit/") # fit result of stan program
glob_params$Stan$f.FitPath  = function(..., list = F) {
  file = .Internal(paste(list(...),sep = "",collapse = ""))
  res = paste0(glob_params$Stan$FitPath, file)
  if(file != ""){
    return(res)
  }else{
    cat("The file/folder name or is not specified")
    if(list) return(list.files(res))
    print(here::here("stan/fit/"))
  }
}

glob_params$DataPath      = paste0(glob_params$rootPath, "/data/") 
# glob_params$MCsimRegeneration = FALSE
# glob_params$DataRegeneration = FALSE

dir.create("./Rnotebook/img")
file.copy("./img","./Rnotebook", recursive = T)

{
  ipak <- function(pkg){
    new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
    if (length(new.pkg)) 
      install.packages(new.pkg, dependencies = TRUE)
    sapply(pkg, require, character.only = TRUE)
  }
  reloadpck = function(reload = TRUE){
    if(reload){
      print(ipak(pckToLoad))
      if("rstan" %in% pckToLoad) {
        rstan_options(auto_write = TRUE)
        options(mc.cores = parallel::detectCores() - 1)
      }
    }
  }
  writeLines(c("use 'pckToLoad = c(...)' to specify r packages in use,",
               "then 'reloadpck()' to install / load them"))
}

{
  sapply(list.files(glob_params$RFunc,
                    pattern="*.R$", 
                    full.names=TRUE, ignore.case=TRUE), 
         source, .GlobalEnv)
}

cleanup = function(){
  if(exists("fit_step_by_step_X_CIa")) {
    rm(fit_step_by_step_X_CIa)
    cat("fit_step_by_step_X_CIa removed\n")
  }
  rm(list = ls(pattern="tmp.",pos=.GlobalEnv), pos=.GlobalEnv) ; cat("tmp.* variables removed\n")
  rm(list = ls(pattern="la"  ,pos=.GlobalEnv), pos=.GlobalEnv) ; cat("la* variables removed\n")
}