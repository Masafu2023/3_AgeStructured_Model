
########
# SIR Model (Susceptible-Infected-Recovered) with 3 Age Groups
library(deSolve)
library(tidyverse)
library(coda)
library(MCMCvis)

# SIR model function with age compartments
SIR_model_age <- function(time, state, parameters) {
  # Children (age group 1)
  S1 <- state[1]
  I1 <- state[2]
  R1 <- state[3]
  
  # Adults (age group 2)
  S2 <- state[4]
  I2 <- state[5]
  R2 <- state[6]
  
  # Elderly (age group 3)
  S3 <- state[7]
  I3 <- state[8]
  R3 <- state[9]
  
  # Total population for each age group
  N1 <- S1 + I1 + R1
  N2 <- S2 + I2 + R2
  N3 <- S3 + I3 + R3
  
  # Transmission rates (force of infection) for each age group
  beta1 <- parameters["beta1"]
  beta2 <- parameters["beta2"]
  beta3 <- parameters["beta3"]
  gamma <- parameters["gamma"]
  
  # Equations for age group 1 (Children)
  dS1 <- -beta1 * S1 * I1 / N1 
  dI1 <- beta1 * S1 * I1 / N1 - gamma * I1
  dR1 <- gamma * I1
  
  # Equations for age group 2 (Adults)
  dS2 <- -beta2 * S2 * I2 / N2
  dI2 <- beta2 * S2 * I2 / N2 - gamma * I2
  dR2 <- gamma * I2
  
  # Equations for age group 3 (Elderly)
  dS3 <- -beta3 * S3 * I3 / N3
  dI3 <- beta3 * S3 * I3 / N3 - gamma * I3
  dR3 <- gamma * I3
  
  # Returning the rates of change for all compartments
  list(c(dS1, dI1, dR1, dS2, dI2, dR2, dS3, dI3, dR3))
}

# Parameters
beta1 <- 0.4   # Force of infection for children
beta2 <- 0.3   # Force of infection for adults
beta3 <- 0.25  # Force of infection for elderly
gamma <- 0.1   # Recovery rate

# Initial conditions for each age group
initial_state <- c(S1 = 50, I1 = 1, R1 = 0,  # Children
                   S2 = 40, I2 = 1, R2 = 0,  # Adults
                   S3 = 30, I3 = 1, R3 = 0)  # Elderly

# Time sequence (one year with daily time steps)
times <- seq(0, 365, by = 1)

# Simulation of the SIR model with age structure
parameters <- c(beta1 = beta1, beta2 = beta2, beta3 = beta3, gamma = gamma)
simulated_data <- ode(y = initial_state, times = times, func = SIR_model_age, parms = parameters) %>% data.frame()

# Ensuring that the observed data columns are numeric
simulated_data[, 2:10] <- lapply(simulated_data[, 2:10], as.numeric)

# Adding a Time column for plotting
simulated_data$Time <- times

# Reshaping data to long format for plotting
simulated_data_long <- simulated_data %>%
  pivot_longer(cols = -Time, names_to = "Compartment", values_to = "Count")

# Plotting all states on one graph with correct color mapping and y-axis limits
ggplot(simulated_data_long, aes(x = Time, y = Count, color = Compartment)) +
  geom_line() +
  labs(title = "SIR Model Simulation with 3 Age Groups", x = "Time (days)", y = "Count") +
  scale_color_manual(values = c("blue", "green", "red",      # Children
                                "blue4", "green4", "red4",  # Adults
                                "blue2", "green2", "red2",  # Elderly
                                "black")) +                  # Recovered (R1, R2, R3)
  coord_cartesian(ylim = c(0, 100)) +  #  y-axis limits from 0 to 100
  theme_minimal()


# Log-likelihood function for the model with age compartments
log_likelihood_age <- function(params, observed_data, initial_state, times) {
  
  # Parameters with the current estimate of beta for each age group
  parameters <- c(params["beta1"], params["beta2"], params["beta3"], gamma = 1/10)
  
  # Simulating the model
  out <- ode(y = initial_state, times = times, func = SIR_model_age, parms = parameters)
  model_data <- as.data.frame(out)
  
  #checking observed data(numeric form)
  observed_data[, 2:10] <- lapply(observed_data[, 2:10], as.numeric)
  
  # Log-likelihood for each compartment in each age group
  log_likelihood <- 0
  for (col in 2:10) {
    obs <- ceiling(observed_data[[col]])
    model <- model_data[[col]]
    log_likelihood <- log_likelihood + sum(dpois(obs, lambda = model, log = TRUE))
  }
  
  return(log_likelihood)
}

# MCMC settings
n_iter <- 5000
beta_init <- c(0.1, 0.15, 0.2)  # Initial beta values for MCMC
n_chains <- 3
sd_prop <- rep(0.001, n_chains)  # Initial SD for each chain
target_accept_rate <- 0.234  # Target acceptance rate for adaptive MCMC
adapt_rate <- 0.01  # Rate of adaptation for the proposal SD

# Incorporating RHAT logic: storage for multiple chains
beta_chains <- matrix(NA, ncol = n_chains, nrow = n_iter)

# Initializing chains with different initial values
beta_chains[1, ] <- beta_init

# Prior distribution: Beta(2, 2)
prior <- function(beta) {
  return(dbeta(beta, 2, 2, log = TRUE))
}

# MCMC loop for each chain
for (chain in 1:n_chains) {
  
  acceptance_counter <- 0 # Reset acceptance counter for each chain
  
  # Initializing log likelihood for current chain
  loglik_curr <- log_likelihood_age(params = c(beta1 = beta_chains[1, chain], beta2 = beta_chains[1, chain], beta3 = beta_chains[1, chain]),
                                    observed_data = simulated_data,
                                    initial_state = initial_state,
                                    times = times) +
    prior(beta_chains[1, chain])
  
  for (i in 2:n_iter) {
    # Proposing new beta from a normal distribution centered around the current value
    beta_proposed <- rnorm(1, mean = beta_chains[i - 1, chain], sd = sd_prop[chain])
    
    if (beta_proposed > 0) {  # Ensure beta is positive
      loglik_prop <- log_likelihood_age(params = c(beta1 = beta_proposed, beta2 = beta_proposed, beta3 = beta_proposed),
                                        observed_data = simulated_data,
                                        initial_state = initial_state,
                                        times = times) + 
        prior(beta_proposed)
    } else {
      loglik_prop <- -1E6  # Penalize invalid proposals
    }
    
    # Calculating acceptance probability
    acceptance_prob <- loglik_prop - loglik_curr
    
    # Metropolis-Hastings acceptance step
    if (log(runif(1)) < acceptance_prob) {
      beta_chains[i, chain] <- beta_proposed
      loglik_curr <- loglik_prop
      acceptance_counter <- acceptance_counter + 1
    } else {
      beta_chains[i, chain] <- beta_chains[i - 1, chain]
    }
    
    # Adaptive adjustment of proposal standard deviation (sd_prop)
    if (i > 100) {  # Adapt every 100 iterations
      acceptance_rate <- acceptance_counter / i
      sd_prop[chain] <- sd_prop[chain] * exp(adapt_rate * (acceptance_rate - target_accept_rate))
    }
  }
}

# Assigning parameter names to the columns in beta_chains
colnames(beta_chains) <- rep("Beta", n_chains)

#  mcmc.list with consistent variable names for each chain
burnin <- 1000
mcmc_out <- mcmc.list(
  as.mcmc(beta_chains[burnin:n_iter, 1, drop = FALSE]),
  as.mcmc(beta_chains[burnin:n_iter, 2, drop = FALSE]),
  as.mcmc(beta_chains[burnin:n_iter, 3, drop = FALSE])
)

# MCMC Trace for Beta
mcmc_df <- data.frame(Iteration = burnin:n_iter, Beta = c(beta_chains[burnin:n_iter, ]))
ggplot(mcmc_df, aes(x = Iteration, y = Beta)) +
  geom_line(color = "blue", linewidth = 1) +
  labs(title = "MCMC Trace for Beta (Adaptive MCMC)", x = "Iteration", y = "Beta") +
  theme_minimal()

# Histogram of posterior distribution
ggplot(data.frame(Beta = c(beta_chains[burnin:n_iter, ])), aes(x = Beta)) +
  geom_histogram(bins = 50, fill = "skyblue", color = "black") +
  labs(title = "Posterior Distribution of Beta", x = "Beta", y = "Frequency") +
  theme_minimal()

# R_hat calculations
summary <- MCMCsummary(mcmc_out, Rhat = TRUE)                     

