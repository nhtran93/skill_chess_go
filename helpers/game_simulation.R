library(rethinking)
library(rstan)
library(lubridate)
library(msm)
library(RColorBrewer)
source("helpers/helper_functions.R")

##### Simulation of games with individual offsets (multivariate normal) for k and m #########
# n_players = 100; n_games = 10; n_years = 0;
# Mu =  c(0.15, 0.01); rho = 0.3; sigma = c(0.05, 0.01);
# b = NA; a = 1; k = NA; m = NA;
# start_year = '2020/09/08';
# start_birth_year = c('1970/01/01', '2018/01/01');
# seed = 123

sim_games <- function(n_players = 100, n_games = 10, n_years = 0,
                      Mu =  c(0.15, 0.01), rho = 0.3, sigma = c(0.05, 0.01),
                      b = NA, a = 1, k = NA, m = NA,
                      start_year = '2020/09/08',
                      start_birth_year = c('1970/01/01', '2018/01/01'),
                      seed = 123){
  
  if (is.na(b)) {
    # b parameter should not have individual variation.
    # The importance of the knowlege should be a property of the game of Chess or Go.
    b <- 2
  }
  # building in individual variation --> individual offsets
  if (is.na(k) & is.na(m)) {
    # sample individual deviations from means.
    # Need to be log-transformed and exponentiated to make k & m positive
    var <- rmvnorm2(n_players, Mu = log(Mu),
                    sigma = sigma,
                    matrix(c(1, rho, rho, 1), nrow = 2))
    var <- exp(var)
    k <- var[, 1]
    m <- var[, 2]
  }
  
  all_games <- data.frame(stringsAsFactors = FALSE)
  all_players <- data.frame(stringsAsFactors = FALSE)
  
  # This loop, loops over the number of years. 
  # n_year = 0 implies cross-sectional data of the same year.
  for (t in 0:n_years) {
    current_year <- as.Date(start_year) %m+% years(t)
    
    # Player Characteristics
    if (t > 0 ) {
      # Update the Current Year, Age and Skill
      age <- day(as.period(interval(birth_years, current_year), "days"))/365
      characteristics$age <- age
      skill <- logit(skill_fun(a = a, k = characteristics$k,
                               m = characteristics$m, b = characteristics$b,
                               age = age))
      characteristics$skill <- skill
      characteristics$current_year <- current_year
    }
    else{
      players <- id_maker(n = n_players, seed = seed, nchars = 3)
      birth_years <- sample(seq(as.Date(start_birth_year[1]),
                                as.Date(start_birth_year[2]), by="day"), n_players)
      age <- day(as.period(interval(birth_years, current_year), "days"))/365
      characteristics <- cbind.data.frame(players, birth_years, age)
      characteristics$players <- as.character(characteristics$players)
      
      # Check whether k & m are fixed or have individual variation
      if (length(k) == n_players & length(m) == n_players) {
        characteristics$k <- k
        characteristics$m <- m
      }
      else{
        # In this case, there is no individual variationa and k & m is fixed across ind.
        characteristics$k <- rep(k, n_players)
        characteristics$m <- rep(m, n_players)
      }
      characteristics$b <- rep(b, n_players)
      skill <- logit(skill_fun(a = a, k = characteristics$k,
                               m = characteristics$m, b = characteristics$b,
                               age = age))
      characteristics$skill <- skill
      characteristics$current_year <- current_year
    }
    
    # Game Dataframe
    games <- t(combn(x = players, m = 2)) # choose(10, 2) = 45 games
    games <- cbind.data.frame(current_year, 1:nrow(games),
                              games, rep(NA, nrow(games)))
    colnames(games) <- c("current_year", "match idx", "player1", "player2", "result")
    games$player1 <- as.character(games$player1)
    games$player2 <- as.character(games$player2)
    
    # samples evenly extra among all uniqe pairs of games, additional games
    idx <- sample(x = 1:nrow(games[!duplicated(games$player1), ][!duplicated(games$player2), ]),
                  size = n_games, replace = TRUE)
    add_games <- games[!duplicated(games$player1), ][!duplicated(games$player2), ][idx, ]
    games <- rbind(games, add_games)
    
    idx <- characteristics[match(games$player1, characteristics$players), "birth_years"]
    games$player1_age <- day(as.period(interval(idx, current_year), "days"))/365
    idx <- characteristics[match(games$player2, characteristics$players), "birth_years"]
    games$player2_age <- day(as.period(interval(idx, current_year), "days"))/365
    games$player1_skill <- logit(skill_fun(a = a,
                                           m = characteristics[match(games$player1,
                                                                     characteristics$players), "m"],
                                           k = characteristics[match(games$player1,
                                                                     characteristics$players), "k"],
                                           b = b, age = games$player1_age))
    games$player2_skill <- logit(skill_fun(a = a,
                                           m = characteristics[match(games$player2,
                                                                     characteristics$players), "m"],
                                           k = characteristics[match(games$player2,
                                                                     characteristics$players), "k"],
                                           b = b, age = games$player2_age))
    
    games$p1_win <-  exp(games$player1_skill)/(exp(games$player2_skill) + exp(games$player1_skill))
    games$result <- rbern(nrow(games), games$p1_win)
    
    all_games <- rbind(all_games, games)
    all_players <- rbind(all_players, characteristics)
  }
  
  
  return(list(games = all_games, characteristics = all_players))
}

#### Simulate the birth year effect as a regression and draw individual variation from multivariate normal #####
# n_players = 100; n_games = 10; n_years = 0;
# location_year = year(Sys.Date()); scale_year = 5;
# beta_k = 0.01; beta_m = -0.015;
# Mu =  c(0.15, 0.01); rho = 0.3; sigma = c(0.05, 0.01);
# b = NA; a = 1; k = NA; m = NA;
# start_year = '2020/07/31';
# start_birth_year = c('1950/01/01', '2018/01/01');
# seed = 123; game_draw = FALSE;
# ordered_logistic_cutpoints = c(-0.5, 0.5)

sim_games_regression_years <- function(n_players = 100, n_games = 10, n_years = 0,
                                       location_year = year(Sys.Date()), scale_year = 5,
                                       beta_k = 0.01, beta_m = -0.015,
                                       Mu =  c(0.15, 0.01), rho = 0.3, sigma = c(0.05, 0.01),
                                       b = NA, a = 1, k = NA, m = NA,
                                       start_year = '2020/07/31',
                                       start_birth_year = c('1950/01/01', '2018/01/01'),
                                       seed = 123, game_draw = FALSE, 
                                       ordered_logistic_cutpoints = c(-0.5, 0.5),
                                       birth_years = NA){
  set.seed(seed)
  all_games <- data.frame(stringsAsFactors = FALSE)
  all_players <- data.frame(stringsAsFactors = FALSE)
  
  # This loop, loops over the number of years. 
  # n_year = 0 implies cross-sectional data of the same year.
  for (t in 0:n_years) {
    current_year <- as.Date(start_year) %m+% years(t)
    
    # Player Characteristics
    if (t > 0) {
      # Update the Current Year, Age and Skill
      age <- day(as.period(interval(characteristics$birth_years, current_year), "days"))/365
      characteristics$age <- age
      skill <- logit(skill_fun(a = a, k = characteristics$k,
                               m = characteristics$m, b = characteristics$b,
                               age = age))
      characteristics$skill <- skill
      characteristics$current_year <- current_year
    }
    else{
      players <- id_maker(n = n_players, seed = seed, nchars = 3)
      
      if (all(is.na(birth_years))) {
        birth_years <- sample(seq(as.Date(start_birth_year[1]),
                                  as.Date(start_birth_year[2]), by="day"), n_players)
      }
      
      if (is.na(b)) {
        # b parameter should not have individual variation.
        # The importance of the knowlege should be a property of the game of Chess or Go.
        b <- 2
      }
      # building in individual variation --> individual offsets
      # sample individual deviations from means.
      # Need to be log-transformed and exponentiated to make k & m positive
      var <- rmvnorm2(n_players, Mu = log(Mu),
                      sigma = sigma,
                      matrix(c(1, rho, rho, 1), nrow = 2))
      # Compute k, m and b with the effect of birth years
      # scale year will influence the interpretation of the beta.
      # if the scale is 1 then, beta will give the effect of birth years by year.
      # if the scale is 10 then, beta will give the effect of birth years by decade.
      # location of the year will influence the meaning of the individual intercepts.
      # if location is set 2020, then the individual intercepts are the mean individual intercepts for the birth year 2020
      # Setting the location to a specific year, makes this comparable across Go and Chess
      # because pure normalizing will make this regression mean different things (different scale).
      k <- var[, 1] + beta_k*((year(birth_years) - location_year)/scale_year)
      m <- var[, 2] + beta_m*((year(birth_years) - location_year)/scale_year)
      k <- exp(k)
      m <- exp(m)
      
      
      age <- day(as.period(interval(birth_years, current_year), "days"))/365
      characteristics <- cbind.data.frame(players, birth_years, age)
      characteristics$players <- as.character(characteristics$players)
      characteristics$age_peak <- maximum_skill_age(b, k, m)
      
      # Check whether k & m are fixed or have individual variation
      if (length(k) == n_players & length(m) == n_players) {
        characteristics$k <- k
        characteristics$m <- m
      }
      else{
        # In this case, there is no individual variation and k & m is fixed across ind.
        characteristics$k <- rep(k, n_players)
        characteristics$m <- rep(m, n_players)
      }
      characteristics$b <- rep(b, n_players)
      skill <- logit(skill_fun(a = a, k = characteristics$k,
                               m = characteristics$m, b = characteristics$b,
                               age = age))
      characteristics$skill <- skill
      characteristics$current_year <- current_year
    }
    
    # Game Dataframe
    games <- t(combn(x = players, m = 2)) # choose(10, 2) = 45 games
    games <- cbind.data.frame(current_year, 1:nrow(games),
                              games, rep(NA, nrow(games)))
    colnames(games) <- c("current_year", "match idx", "player1", "player2", "result")
    games$player1 <- as.character(games$player1)
    games$player2 <- as.character(games$player2)
    
    # samples evenly extra among all uniqe pairs of games, additional games
    idx <- sample(x = 1:nrow(games[!duplicated(games$player1), ][!duplicated(games$player2), ]),
                  size = n_games, replace = TRUE)
    add_games <- games[!duplicated(games$player1), ][!duplicated(games$player2), ][idx, ]
    games <- rbind(games, add_games)
    
    idx <- characteristics[match(games$player1, characteristics$players), "birth_years"]
    games$player1_age <- day(as.period(interval(idx, current_year), "days"))/365
    idx <- characteristics[match(games$player2, characteristics$players), "birth_years"]
    games$player2_age <- day(as.period(interval(idx, current_year), "days"))/365
    games$player1_skill <- logit(skill_fun(a = a,
                                           m = characteristics[match(games$player1,
                                                                     characteristics$players), "m"],
                                           k = characteristics[match(games$player1,
                                                                     characteristics$players), "k"],
                                           b = b, age = games$player1_age))
    games$player2_skill <- logit(skill_fun(a = a,
                                           m = characteristics[match(games$player2,
                                                                     characteristics$players), "m"],
                                           k = characteristics[match(games$player2,
                                                                     characteristics$players), "k"],
                                           b = b, age = games$player2_age))
    
    if (game_draw) {
      games$p1_win <- exp(games$player1_skill)/(exp(games$player2_skill) + exp(games$player1_skill))
      
      games$result <- ologis_rng(logit(games$p1_win), ordered_logistic_cutpoints)
    }
    else{
      games$p1_win <-  exp(games$player1_skill)/(exp(games$player2_skill) + exp(games$player1_skill))
      games$result <- rbern(nrow(games), games$p1_win)
    }
    
    all_games <- rbind(all_games, games)
    all_players <- rbind(all_players, characteristics)
  }
  
  
  return(list(games = all_games, characteristics = all_players))
}
