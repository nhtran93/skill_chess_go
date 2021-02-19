maximum_skill_age <- function(b, k, m){
  out <- log(1 + (b * k)/m) / k
  return(out)
}

id_maker <- function(n, reserved = "", seed = NA, nchars = NA){
  my_let <- letters 
  my_num <- 0:9 
  if(is.na(seed) | !is.numeric(seed)) set.seed(as.numeric(as.POSIXlt(Sys.time())))
  if(!is.na(seed) & is.numeric(seed)) set.seed(seed)
  output <- replicate(n, paste(sample(c(my_let, my_num), nchars, replace=TRUE), 
                               collapse=''))
  rejected <- duplicated(output) | output %in% reserved |
    substr(output, 1, 1) %in% my_num
  while (any(rejected)) {
    output <- output[-which(rejected)]
    remaining <- n - length(output)
    output <- c(output, replicate(remaining, paste(sample(c(my_let, my_num), nchars, 
                                                          replace=TRUE), collapse="")))
    rejected <- duplicated(output) | output %in% reserved |
      substr(output, 1, 1) %in% my_num
  }
  output
}

inv_logit <- function (x) {
  p <- 1/(1 + exp(-x))
  p <- ifelse(x == Inf, 1, p)
  p
}

logit <- function (x){log(x) - log(1 - x)}

skill_fun <- function(a = 1, k, m, b, age){
  (exp(-m*age))*((a-exp(-k*age))^b)
}

rmvnorm2 <- function (n, Mu = rep(0, length(sigma)), sigma = rep(1, length(Mu)), 
                      Rho = diag(length(Mu)), method = "chol") 
{
  DS <- diag(sigma)
  SIGMA <- DS %*% Rho %*% DS
  rmvnorm(n = n, mean = Mu, sigma = SIGMA, method = method)
}

ordered_logistic_stan_code <- 
  '
functions {
  int[] ologis_rng(vector eta, vector c) {
    int l = num_elements(eta);
    int out[l];
    
    for(i in 1:l) out[i] = ordered_logistic_rng(eta[i], c);
    
    return out;
  }
}
'

rstan::expose_stan_functions(rstan::stanc(model_code = ordered_logistic_stan_code))

find_peaks <- function (x, m = 3){
  shape <- diff(sign(diff(x, na.pad = FALSE)))
  pks <- sapply(which(shape < 0), FUN = function(i){
    z <- i - m + 1
    z <- ifelse(z > 0, z, 1)
    w <- i + m + 1
    w <- ifelse(w < length(x), w, length(x))
    if(all(x[c(z : i, (i + 2) : w)] <= x[i + 1])) return(i + 1) else return(numeric(0))
  })
  pks <- unlist(pks)
  pks
}

jeffreys <- function(p, q){
  sum((p - q)*(log(p)*log(q)))
}

kl <- function(p, q){
  sum(p*log(p/q))
}