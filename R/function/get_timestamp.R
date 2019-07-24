get_timestamp = function(){
  options(scipen=0, digits=7)
  as.numeric(gsub("[: -]", "" , 
                  Sys.time(), 
                  perl=TRUE)
  ) %% 1E12 %/% 100
}

cat("One function is imported\n",
    " get_timestamp() : get a timestamp as `numeric` of format \n  yy-mm-dd-hh-mm :  ",
    sep = "")
