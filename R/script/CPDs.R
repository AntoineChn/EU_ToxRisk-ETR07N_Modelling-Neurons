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
  if_else(x < params$x_tox_min,
          true = 0,
          false = params$y_max * (1 - exp(-params$k * (x - params$x_tox_min)))
          ) %>% 
  return()
}




