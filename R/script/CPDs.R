# CPDs.R


# Logistic function -------------------------------------------------------

f.link.lgstc = function(x,
                   params # = list(k = 1, y_max = 1, y_min = 0, x_50 = 0)
) {
  res = (params$y_max - params$y_min) / (1 + exp(-params$k * (x - params$x_50))) + params$y_min
  return(res)
}

# Logistic function -------------------------------------------------------

f.link.linear = function(x,
                         params # = list(k = 1, y_max = 1, y_min = 0, x_50 = 0)
) {
  res = params$beta * x + params$beta_0
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

# truncated tox_min
f.link.exp_tox_min = function(x,
                              params) {
  require(tidyverse)
  lapply(x, function(x.loc){
    if(x.loc > params$x_tox_min){
      res = params$y_max * (1 - exp(-params$k * (x.loc - params$x_tox_min)))
      return(res)
    } else {
      res = 0
      return(res)
    }
  }) %>% unlist() %>% return()
}
