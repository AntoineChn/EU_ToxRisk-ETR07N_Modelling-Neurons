f.stanFit.newName = function(chemi,
                             stanfile.name,
                             timestamp,
                             nb.chains,
                             nb.iter,
                             seed){
  newName = paste0(stanfile.name,
                   "-chemi_",chemi,
                   "-timestamp_",timestamp,
                   "-chain_",nb.chains,
                   "-iter_",nb.iter,
                   "-seed_",seed)
  cat("The following is the name of the stan fit that will be exported to ./stan/fit/ folder\n", 
      newName)
  newName =
    gsub("\\+", "", 
         paste0(newName))
  return(newName)
}

# Read different timestamp ------------------------------------------------

f.stanFit.filter.timestamp = function(pattern,
                                      newest.first = T,
                                      as.num = F){
  tmp.fitnames = list.files("stan/fit", pattern = pattern)
  pos.begin = as.numeric ( regexpr('-timestamp_', tmp.fitnames))
  nb.char = nchar("-timestamp_")
  pos.begin = pos.begin + nb.char
  tmp.ts = substr(tmp.fitnames, pos.begin, pos.begin + 9)
  tmp.ts = unique(tmp.ts)
  if(length(tmp.ts) > 1) {
  	message("multiple timestamps are found in ./stan/fit/ folder")
  	if(newest.first == T) message("Newest first")
  	if(newest.first == F) message("Oldest first")
  }
  	if(as.num) tmp.ts = as.numeric(tmp.ts)
  	return(sort(tmp.ts ,decreasing = newest.first))
}

cat("\nTwo functions are imported\n",
        "f.stanFit.newName() : useful for simulation\n",
        "f.stanFit.filter.timestamp() : useful for posterior analysis\n")


