f.stanFit.newName = function(chemi,
                             stanfile.name,
                             timestamp,
                             nb.chains,
                             nb.iter,
                             seed,
                             n_eff.min = 0,
                             rhat.max = 0){
  newName = paste0(stanfile.name,
                   "-ts_",timestamp,
                   "-chemi_",chemi,
                   "-chain_",nb.chains,
                   "-iter_",nb.iter,
                   "-seed_",seed)
  
  if(n_eff.min != 0){
    newName = paste0(newName,
           "-nEffMin_",n_eff.min,
           "-rHatMax100_", rhat.max * 100)
  }
  cat("The following is the name of the stan fit that will be exported to ./stan/fit/ folder\n", 
      newName,"\n")
  newName =
    gsub("\\+", "", 
         paste0(newName))
  return(newName)
}

# Read different timestamp ------------------------------------------------

f.stanFit.filter = function(pattern,chemi = "Rotenone", extension = "stanFit"){
  tmp.fitnames1 = list.files("stan/fit", 
                            pattern = glob2rx(paste0("*",pattern,
                                                     "*.stanFit"),)
                            )
  tmp.fitnames2 = list.files("stan/fit", 
                            pattern = glob2rx(paste0("*",chemi,
                                                     "*.stanFit"),)
  )
  tmp.fitnames = intersect(tmp.fitnames1, tmp.fitnames2)
  return(tmp.fitnames)
}

f.stanFit.filter.timestamp = function(pattern, chemi = "Rotenone",
                                      newest.first = T,
                                      as.num = F){
  # tmp.fitnames = list.files("stan/fit", pattern = glob2rx(paste0("*",pattern,"*.stanFit*")))
  tmp.fitnames = f.stanFit.filter(pattern = pattern, chemi = chemi)
  pos.begin_old = as.numeric ( regexpr('-timestamp_', tmp.fitnames))
  nb.char = nchar("-timestamp_")
  pos.begin_old2 = pos.begin_old + nb.char
  tmp.ts = NA
  tmp.ts = substr(tmp.fitnames[pos.begin_old >= 0], 
                  pos.begin_old2[pos.begin_old >= 0], 
                  pos.begin_old2[pos.begin_old >= 0] + 9)
  
  pos.begin_new = as.numeric ( regexpr('-ts_', tmp.fitnames))
  nb.char = nchar("-ts_")
  pos.begin_new2 = pos.begin_new + nb.char
  tmp.ts = c(tmp.ts,
             substr(tmp.fitnames[pos.begin_new >= 0], 
                    pos.begin_new2[pos.begin_new >= 0], 
                    pos.begin_new2[pos.begin_new >= 0] + 9)
  )
  
  tmp.ts = unique(tmp.ts)
  if(length(tmp.ts) > 1) {
  	message("multiple timestamps are found in ./stan/fit/ folder")
  	if(newest.first == T) message("Newest first")
  	if(newest.first == F) message("Oldest first")
  }
  	if(as.num) tmp.ts = as.numeric(tmp.ts)
  	return(sort(tmp.ts ,decreasing = newest.first))
}

cat("\nThree functions are imported\n",
        "f.stanFit.newName() : useful for simulation\n",
        "f.stanFit.filter() and f.stanFit.filter.timestamp : useful for posterior analysis\n")


