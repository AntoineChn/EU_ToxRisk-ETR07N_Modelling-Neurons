f.open.wd = function(){
  os = get_os()
  if(os == "osx"){
    system("open -a Finder ./")
  }else{
    message('This only works on MacOS')
  }
}

f.open.stanFit = function(path = "stan/fit"){
  os = get_os()
  if(os == "osx"){
    system(paste0("open -a Finder ./", path))
  }else{
    message('This only works on MacOS')
  }
}

f.open.stanCode = function(path = "stan/code"){
  os = get_os()
  if(os == "osx"){
    system(paste0("open -a Finder ./", path))
  }else{
    message('This only works on MacOS')
  }
}

cat("\nThree functions are imported\n",
    "f.open.wd() : open working directory on MacOS \n",
    "f.open.stanFit() : open ./stan_fit/ folder on MacOS \n",
    "f.open.stanCode() : open ./stan_code/ folder on MacOS \n")
