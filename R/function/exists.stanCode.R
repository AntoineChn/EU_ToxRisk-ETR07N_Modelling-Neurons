f.exists.stanCode = function(stanfile.name,
                             path = "stan/code"){
  if(!file.exists(paste0(path,'/',stanfile.name,".stan"))){
    message(stanfile.name, ".stan", 
            " does NOT exist in ./",path,"/ folder\n",
            "  Please check the name of your file !!!")
    source(glob_params$f.RFunc("get_os.R"))
    os = get_os()
    rm(get_os,pos=.GlobalEnv)
    if(os == "osx"){
      cat("Do you want to check the folder ./",path, " ? [y/N]",
          sep = "")
      res = readline()
      if( tolower(res) == "y" ) {
        system(paste0("open -a Finder ./", path))
      }
    }
    return(FALSE)
  }else{
    return(TRUE)
  }
}

cat("One function is imported\n",
    "f.exists.stanCode() : verify if a *.stan exist in ./stan/code/ folder \n")
