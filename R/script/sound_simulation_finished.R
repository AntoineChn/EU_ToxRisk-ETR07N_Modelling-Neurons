if(get_os() == "osx"){
  system("say Simulation finished!")
}else{
  if(require(beepr)) install.packages('beepr')
  beepr::beep()
}