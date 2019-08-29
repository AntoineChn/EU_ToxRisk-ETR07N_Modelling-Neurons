
source(glob_params$f.RFunc("get_os.R"))
if(get_os() == "osx"){
  system("say Simulation finished!")
}else{
  if(require(beepr)) install.packages('beepr')
  beepr::beep()
}
rm(get_os, pos=.GlobalEnv)