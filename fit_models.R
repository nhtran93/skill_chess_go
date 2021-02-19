rm(list=ls())
library(rethinking)
library(rstan)
library(lubridate)
library(msm)
library(RColorBrewer)
source("helpers/game_simulation.R")
############ Prior Predictives ##############

# Knowledge
a  <- 1 # upper bound of the knowledge a > 1
k <- 0.15 # rate of growth, k > 0
curve(a-exp(-k*x), from=0, to=100, xlab="age", ylab="knowledge", bty = "l")

# Senescence
m <- 0.01 # decline rate due to senescence, m > 0
curve(exp(-m*x), from=0, to=100, xlab="age", ylab="knowledge", bty = "l")

# Skill
b = 2.5 # elasticity: b controls the relative importance of Knowledge
curve((exp(-m*x))*((a-exp(-k*x))^b), from=0, to=100,
      xlab="age", ylab="Skill", bty = "l", ylim = c(0, 1))

######## Prior Predictives Simple
mu_k <- exp(rnorm(500, -1.9, 0.25))
mu_m <- exp(rnorm(500, -4.6, 0.15))
b <- rtnorm(500, 2, 0.1, lower = 0)
curve((exp(-mean(mu_m)*x))*((1-exp(-mean(mu_k)*x))^mean(b)), lwd = 2,
      col = "black", from = 0, to = 100, ylim = c(0,1),
      xlab = "age", ylab = "Skill")
sigma_k <- rexp(500, 50)
sigma_m <- rexp(500, 50)
for (i in 1:100) {
  var <- rmvnorm2(100, Mu = log(c(mu_m[i], mu_k[i])),
                  sigma = c(sigma_m[i], sigma_k[i]),
                  matrix(c(1, 0.3, 0.3, 1), nrow = 2))
  var <- exp(var)
  m <- var[, 1]
  k <- var[, 2]
  curve((exp(-m[i]*x))*((1-exp(-k[i]*x))^b[i]), add = TRUE, col = col.alpha("#800026", 0.2))
}
######## Prior Predictives Complex
set.seed(1)
mu_k <- exp(rnorm(500, -2.5, 0.5))
mu_m <- exp(rnorm(500, -4.5, 0.5))
b <- rtnorm(500, 2, 0.3, lower = 0)
sigma_k <- rnorm(500, 0.5, 0.5)
sigma_m <- rnorm(500, 0.5, 0.5)
beta_k <- rnorm(500, 0.15, 0.2)
beta_m <- rnorm(500, -0.15, 0.2)
birth_years <- sample(seq(as.Date('1970/01/01'),
                          as.Date('2018/01/01'), by="day"), 500)
i = 1
var <- rmvnorm2(500, Mu = log(c(mu_m[i], mu_k[i])),
                sigma = c(sigma_m[i], sigma_k[i]),
                matrix(c(1, 0.3, 0.3, 1), nrow = 2))
k <- var[, 2] + beta_k*((year(birth_years) - 10)/2020)
m <- var[, 1] + beta_m*((year(birth_years) - 10)/2020)
k <- exp(k)
m <- exp(m)

png("output/prior_predictives.png", res = 800, height = 5, width = 5, units = "in")
cairo_ps(file = "output/prior_predictives.eps",
         onefile = TRUE, fallback_resolution = 800, width = 7, height = 7)
curve(skill_fun(age= x, k = mean(k), m = mean(m), b = mean(b)), lwd = 3,
      col = "black", from = 0, to = 100, ylim = c(0,1),
      xlab = "Age", ylab = "Skill")
peaks <- c()
for (i in 1:300) {
  curve(skill_fun(age= x, k = k[i], m = m[i], b = b[i]),
        add = TRUE, col = col.alpha("#800026", 0.2), from = 0, to = 100)
  peaks <- c(peaks, find_peaks(skill_fun(age= 1:100, k = 0.15, m = 0.01, b = 2), m = 0.1))
}
dev.off()

############### Fit Stats Model and Simulations #############################
# Go Game Simulation. Only Ind offsets without birth cohort effects. 
sim <- sim_games(n_players = 50, n_years = 10, n_games = 100,
                 Mu =  c(0.15, 0.01), sigma = c(0.35, 0.05), b = 2,
                 rho = 0.3, seed = 123)
range(year(sim$characteristics$birth_years))
summary(sim$games$p1_win)
summary(sim$characteristics$age)
nrow(sim$games)
table(sim$characteristics$age)


#### Birth year offset
# Go games
sim <- sim_games_regression_years(n_players = 100, n_years = 5,
                                  scale_year = 10, location_year = 2020,
                                  beta_k = 0.15, beta_m = -0.15,
                                  sigma = c(0.35, 0.05), b = 2,
                                  Mu =  c(0.15, 0.01), game_draw = FALSE,
                                  start_birth_year = c('1950/01/01', '2015/01/01'),
                                  seed = 123)
range(year(sim$characteristics$birth_years))
summary(sim$characteristics$age)
summary(sim$characteristics$k)
nrow(sim$games)
table(sim$characteristics$age)

data_list <- list(
  N = nrow(sim$games),
  P = length(unique(sim$characteristics$players)),
  K = length(unique(sim$games$result)),
  years = (year(sim$characteristics[!duplicated(sim$characteristics$players), "birth_years"])-2020)/10,
  player1 = match(sim$games$player1, sim$characteristics$players),
  player2 = match(sim$games$player2, sim$characteristics$players),
  p1_age = sim$games$player1_age,
  p2_age = sim$games$player2_age,
  outcome = sim$games$result
)

Niter = 1500
Nchains = 2
Ncores = 4
samples <- stan("stan/skill_model_birthyear_regression.stan",
                data = data_list, chains = Nchains, iter = Niter,
                warmup = 500,
                control = list(adapt_delta = 0.99),
                cores = Ncores, seed = 123)

precis(samples)
pairs(samples, pars = c("k", "b", "m"))
posterior <- extract(samples)

# Chess Games
sim <- sim_games_regression_years(n_players = 10, n_years = 5,
                                  scale_year = 10, location_year = 2020, b = 2,
                                  beta_k = 0.1, beta_m = -0.1,
                                  Mu =  c(0.15, 0.01), game_draw = TRUE,
                                  ordered_logistic_cutpoints = c(-0.5, 0.5),
                                  start_birth_year = c('1950/01/01', '2015/01/01'),
                                  seed = 123)
range(year(sim$characteristics$birth_years))
summary(sim$characteristics$age)
summary(sim$characteristics$k)
nrow(sim$games)
table(sim$characteristics$age)

data_list <- list(
  N = nrow(sim$games),
  P = length(unique(sim$characteristics$players)),
  K = length(unique(sim$games$result)),
  years = (year(sim$characteristics[!duplicated(sim$characteristics$players), "birth_years"])-2020)/10,
  player1 = match(sim$games$player1, sim$characteristics$players),
  player2 = match(sim$games$player2, sim$characteristics$players),
  p1_age = sim$games$player1_age,
  p2_age = sim$games$player2_age,
  outcome = sim$games$result
)

Niter = 1500
Nchains = 2
Ncores = 4
samples <- stan("stan/skill_chess_model_birthyear_regression.stan",
                data = data_list, chains = Nchains, iter = Niter,
                warmup = 500,
                control = list(adapt_delta = 0.99),
                cores = Ncores, seed = 123)

precis(samples)
pairs(samples, pars = c("k", "b", "m"))
posterior <- extract(samples)
