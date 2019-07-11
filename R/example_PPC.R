# Load the relevant libraries
library(reshape2); library(arm); library(ggplot2); library(dplyr)

# Set some known parameter values

# We'll just assume a single covariate x

# parameters for the logit
beta11 <- -1.5
beta12 <- 0.5

# parameters for the beta
beta21 <- 5
beta22 <- 0.5
beta31 <- 2
beta32 <- -0.3

# Number of observations
N <- 1000

# Generate our predictor x
x <- rnorm(N)

# Generate some draws from the beta distribution for the continuous outcome
continuous <- rbeta(N, beta21 + beta22*x, beta31 + beta32*x)

# Generate the probabilities for the binary outcome
prob_binary <- invlogit(beta11 + beta12*x)

# Generate the draws for the binary outcome
binary <- rbinom(N, size = 1, prob = prob_binary)

# The outcome as defined
outcome <- ifelse(binary==1, binary, continuous)

data_frame(outcome) %>% 
  ggplot(aes(x = outcome)) +
  geom_histogram(alpha = 0.3) +
  ggthemes::theme_economist() +
  ggtitle("One draw from our generative model")


expected_value <- (1 - prob_binary)*((beta21 + beta22*x)/(beta21 + beta22*x+ beta31 + beta32*x)) + prob_binary

data_frame(expected_value) %>% 
  ggplot(aes(x = expected_value)) +
  geom_histogram(alpha = 0.3) +
  ggthemes::theme_economist() +
  ggtitle("Predicted values given x") +
  xlim(0, 1)



# Stan model --------------------------------------------------------------

beta_logit_model <- "
data {
  int N; // number of observations
  int P; // number of explanatory variables
  int N2; // number of observations to predict
  vector[N] Y; // the proportion paid of expected at some period ahead
  matrix[N, P] X; // The explanatory variables (make sure no variables known in the future!)
  matrix[N2, P] X_new; // explanatory variables for out-of-sample
}
transformed data {
// we an integer Y and a real Y (Y_adj is only here for flexibility of extension)
  vector[N] Y_adj;
  int Y_int[N];
  
  for(i in 1:N) {
    if(Y[i] == 1.0) {
      Y_adj[i] = 1.0;
      Y_int[i] = 1;
    } else {
      Y_adj[i] = Y[i];
      Y_int[i] = 0;
    }
  }
}
parameters {
  matrix[3, P] beta;
  vector[3] mu;
}
model {
  // priors 
  to_vector(beta) ~ normal(0, .2);
  mu ~ normal(0, 1);

  // likelihood
  Y_int ~ bernoulli_logit(mu[1] + X*beta[1]');
  
  for(i in 1:N) {
    if(Y_int[i]==0) {
      // We constrain alpha and beta to be positive by taking exp()

      Y_adj[i] ~ beta(exp(mu[2] + X[i]*beta[2]'),exp(mu[3] + X[i]*beta[3]'));
    }
  }
}
generated quantities {
  vector[N2] Y_predict;
  int tmp;
  for(i in 1:N2) {
    tmp = bernoulli_rng(inv_logit(mu[1] + X_new[i]*beta[1]'));
    if(tmp ==1) {
      Y_predict[i] = 1.0;
    } else {
      Y_predict[i] = beta_rng(exp(mu[2] + X_new[i]*beta[2]'), exp(mu[3] + X_new[i]*beta[3]'));
    }
    if(is_nan(Y_predict[i])) {
      Y_predict[i] = 1.1;
    }
  }
}

"



# fit ---------------------------------------------------------------------

library(rstan)
options(mc.cores = parallel::detectCores() - 1)

beta_logit_estimation <- stan(model_code = beta_logit_model, 
                              data = list(N = N,
                                          P = 1, 
                                          N2 = N,
                                          Y = outcome,
                                          X = matrix(x, N, 1),
                                          X_new = matrix(x, N, 1)),
                              iter = 600)



# Post analysis -----------------------------------------------------------

predictions <- extract(beta_logit_estimation,"Y_predict" ) %>% 
  melt()

predictions %>% 
  ggplot() +
  geom_line(aes(x = value, group = Var2), colour = "orange", stat = "density", alpha = 0.1, adjust = .8) +
  geom_density(data = data_frame(outcome), aes(x = outcome), colour = "black", adjust = 0.8) +
  ggthemes::theme_economist() +
  ggtitle("Actual outcomes and posterior predictive replications") +
  annotate("text", x = 0.2, y = 5, label = "Density of actual outcomes", hjust = 0) +
  annotate("text", x = 0.2, y = 3.5, label = "Posterior replications", colour = "orange", hjust = 0) 


shinystan::launch_shinystan_demo(
)