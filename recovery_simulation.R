######################## 100 Parameter Recovery Simulations (Parallelization) ################
rm(list=ls())
library(rethinking)
library(rstan)
library(lubridate)
library(msm)
library(RColorBrewer)
source("helpers/game_simulation.R")

extract_estimates <- function(model, data, seed, chains, iter, warmup){
  
  samples <- sampling(object = model, data = data, seed = seed,
                      chains = chains, iter = iter, warmup = warmup,
                      control = list(adapt_delta = 0.99, max_treedepth = 11))
  samples <- extract.samples(samples)
  
  # estimated
  estimates <- list(
    mu_k = exp(samples$mu_k),
    mu_m = exp(samples$mu_m),
    b = samples$b,
    sigma_k = samples$sigma_k,
    sigma_m = samples$sigma_m,
    beta_k = samples$beta_k,
    beta_m = samples$beta_m,
    m = samples$skill_parameters[, , 1],
    k = samples$skill_parameters[, , 2],
    age_peak = samples$age_peak
  )
  
  return(estimates)
  
}


set.seed(123)
mu_k <- exp(rnorm(150, -2.5, 0.5))
mu_m <- exp(rnorm(150, -4.5, 0.5))
b <- rtnorm(150, 2, 0.3, lower = 0)
sigma_k <- rtnorm(150, 0, 0.5, lower = 0)
sigma_m <- rtnorm(150, 0, 0.5, lower = 0)
beta_k <- rnorm(150, 0.15, 0.2)
beta_m <- rnorm(150, -0.15, 0.2)
true_values <- cbind(mu_k, mu_m, b, sigma_k, sigma_m, beta_k, beta_m)
true_values <- data.frame(true_values, stringsAsFactors = FALSE)
true_values$m <- NA
true_values$k <- NA
true_values$age_peak <- NA

estimated_values <- data.frame(matrix(NA, ncol = ncol(true_values),
                                      nrow = nrow(true_values)),
                               stringsAsFactors = FALSE)
colnames(estimated_values) <- colnames(true_values)


model <- stan_model(file = "stan/skill_model_birthyear_regression.stan")

########## Simulate Data
set.seed(123)
birth_years <- sample(seq(as.Date('1950/01/01'),
                          as.Date('2018/01/01'), by="day"), 150)
all_data_list <- list()
all_data <- list()
for (i in 1:100) {
  sim <- sim_games_regression_years(n_players = 150, n_games = 0, n_years = 15,
                                    scale_year = 10, location_year = 1950,
                                    beta_k = true_values[i, "beta_k"],
                                    beta_m = true_values[i, "beta_m"], b = true_values[i, "b"],
                                    Mu =  c(true_values[i, "mu_k"], true_values[i, "mu_m"]),
                                    game_draw = FALSE,
                                    sigma = c(true_values[i, "sigma_k"], true_values[i, "sigma_m"]),
                                    birth_years = birth_years)
  
  data_list <- list(
    N = nrow(sim$games),
    P = length(unique(sim$characteristics$players)),
    K = length(unique(sim$games$result)),
    years = (year(sim$characteristics[!duplicated(sim$characteristics$players), "birth_years"])-1950)/10,
    player1 = match(sim$games$player1, sim$characteristics$players),
    player2 = match(sim$games$player2, sim$characteristics$players),
    p1_age = sim$games$player1_age,
    p2_age = sim$games$player2_age,
    outcome = sim$games$result
  )
  
  all_data_list[[i]] <- data_list
  all_data[[i]] <- sim
}

library(parallel)

# dispatch to cores
num_cores <- parallel::detectCores()


Niter = 2000
chains = 2
warmup = 1000

result <- mclapply(
  X = all_data_list,
  FUN = function(i) extract_estimates(model = model, data = i, seed = 123,
                                      chains = chains, iter = Niter, warmup = warmup),
  mc.cores = 50, mc.silent = FALSE
)

save(result, all_data, all_data_list, true_values,
     file = "output/parameter_recovery_100.RData")

################################# Plot ####################################
load("output/parameter_recovery_100.RData")

#png("output/recovery_mcmc100.png", width = 5, height = 10, units = "in", res = 800)
postscript("output/recovery_mcmc100.eps", onefile = TRUE, horizontal = FALSE, width = 5, height = 10)
par(mfrow = c(4, 2))

plot(true_values$mu_k, lapply(result, function(x) mean(x$mu_k)),
     xlim = c(0, 0.3), ylim = c(0, 0.3), 
     ylab = "estimated", xlab = "true",
     main = "mu k",
     pch = 16)
abline(coef = c(0,1), lwd = 2)

plot(true_values$mu_m, lapply(result, function(x) mean(x$mu_m)),
     xlim = c(0, 0.1), ylim = c(0, 0.1), 
     ylab = "estimated", xlab = "true",
     main = "mu m",
     pch = 16)
abline(coef = c(0,1), lwd = 2)

plot(true_values$sigma_k, lapply(result, function(x) mean(x$sigma_k)),
     xlim = c(0, 1), ylim = c(0, 1), 
     ylab = "estimated", xlab = "true",
     main = "sigma k",
     pch = 16)
abline(coef = c(0,1), lwd = 2)

plot(true_values$sigma_m, lapply(result, function(x) mean(x$sigma_m)),
     xlim = c(0, 1), ylim = c(0, 1), 
     ylab = "estimated", xlab = "true",
     main = "sigma m",
     pch = 16)
abline(coef = c(0,1), lwd = 2)


plot(true_values$beta_k, lapply(result, function(x) mean(x$beta_k)),
     xlim = c(-0.5, 0.5), ylim = c(-0.5, 0.5), 
     ylab = "estimated", xlab = "true",
     main = "beta k",
     pch = 16)
abline(coef = c(0,1), lwd = 2)

plot(true_values$beta_m, lapply(result, function(x) mean(x$beta_m)),
     xlim = c(-0.5, 0.5), ylim = c(-0.5, 0.5), 
     ylab = "estimated", xlab = "true",
     main = "beta m",
     pch = 16)
abline(coef = c(0,1), lwd = 2)

plot(true_values$b, lapply(result, function(x) mean(x$b)),
     xlim = c(0, 3), ylim = c(0, 3),
     ylab = "estimated", xlab = "true",
     main = "b",
     pch = 16)
abline(coef = c(0,1), lwd = 2)
dev.off()
