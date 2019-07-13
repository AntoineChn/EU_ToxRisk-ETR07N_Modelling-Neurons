f.stanFit.newName = function(){
  newName = paste0(tmp.stanfile.name,
                   "_Chemi_",tmp.chemi,
                   tmp.timestamp,
                   "_Chain_",nb.chains,
                   "_iter_",nb.iter,
                   "_seed_",my.seed)
  cat(newName)
  newName =
    gsub("\\+", "", 
         glob_params$Stan$f.FitPath(newName))
  return(newName)
}
