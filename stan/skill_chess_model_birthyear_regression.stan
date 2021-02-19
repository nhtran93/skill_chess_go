functions{
    // Skill function
    real skill_fun(real age, real m, real k, real b){
        real score = exp(-m * age)*((1-exp(-k * age))^b);
        return(score);
    } 
}
data {
    int<lower=0> N; // N games
    int<lower=0> P; // N Player
    int<lower=0> K; // Numer of outcomes = 3: Win, Draw, Lose
    vector[P] years; // Indicator for year of the game
    
     
    // Each player is referred to by an integer that acts as an index for the skill vector. 
    int<lower=1, upper = P> player1[N]; // Indicator for player 1
    int<lower=1, upper = P> player2[N]; // Indicator for player 2
    vector[N] p1_age;
    vector[N] p2_age;
    int<lower = 0, upper = 3> outcome[N]; // Results. 1 if player 1 won, 0 if player 2 won.
}

parameters {
    matrix[2, P] z_skill;
    real<lower=0> sigma_m;
    real<lower=0> sigma_k;
    cholesky_factor_corr[2] L_Rho_skill;
    real<lower=0> b;
    real mu_k;
    real mu_m;
    real beta_k;
    real beta_m;
    ordered[K-1] cut_points;
}

transformed parameters{
    matrix[P, 2] skill_parameters;
    vector[N] p1_skill;
    vector[N] p2_skill; 
    vector[2] sigma = [sigma_m, sigma_k]';
    matrix[P, 2] mu = append_col(rep_vector(mu_m, P), rep_vector(mu_k, P));
    
    // the Mu estimates for k and m are estimated on the log scale
    skill_parameters = mu + (diag_pre_multiply(sigma, L_Rho_skill) * z_skill)' + append_col(beta_m*years, beta_k*years);
    skill_parameters = exp(skill_parameters);
    
    for (i in 1: N){
        p1_skill[i] = logit(skill_fun(p1_age[i], skill_parameters[player1[i], 1],
        skill_parameters[player1[i], 2], b));
        p2_skill[i] = logit(skill_fun(p2_age[i], skill_parameters[player2[i], 1],
        skill_parameters[player2[i], 2], b));
    }
}
 
model {
    b ~ normal(2, 0.3);
    mu_m ~ normal(-4.5, 0.5);
    mu_k ~ normal(-2.5, 0.5);
    beta_k ~ normal(0.15, 0.2);
    beta_m ~ normal(-0.15, 0.2);
    
    L_Rho_skill ~ lkj_corr_cholesky(3);
    sigma_m ~ normal(0, 0.5);
    sigma_k ~ normal(0, 0.5);
    to_vector(z_skill) ~ normal(0, 1);
    
    for (i in 1:N){
        // Using the Categorical_logit allows us to directly compute the outcome without the additional step of using the inverse_logit.
        outcome[i] ~ ordered_logistic(p1_skill[i] - p2_skill[i], cut_points);
    }
}
generated quantities{
    vector[P] age_peak;
    for(i in 1:P){
      age_peak[i] = log(1 + (b * skill_parameters[i, 2])/skill_parameters[i, 1]) / skill_parameters[i, 2];
    }
}
