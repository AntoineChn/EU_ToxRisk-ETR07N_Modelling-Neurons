# CPDs.R


# Logistic function -------------------------------------------------------

f.link.lgstc = function(x,
                   params # = list(k = 1, y_max = 1, y_min = 0, x_50 = 0)
) {
  res = (params$y_max - params$y_min) / (1 + exp(-params$k * (x - params$x_50))) + params$y_min
  return(res)
}

# Truncation function -----------------------------------------------------

f.link.truncation = function(y, lower = 0, upper = 120) {
  require(tidyverse)
  case_when(y > upper ~ upper,
            y < lower ~ lower,
            TRUE      ~ y) %>%
    return()
}


# Trancated logistic function ---------------------------------------------

f.link.trancated.lgstc = function(x,
                             params,
                             lower = 0,
                             upper = 1) {
  require(tidyverse)
  f.link.lgstc(x, params) %>%
    f.link.truncation(., lower=lower, upper=upper) %>%
    return()
}