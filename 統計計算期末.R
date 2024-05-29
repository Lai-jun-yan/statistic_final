library(LearnBayes)
data("birthweight")
library(lattice)
mygender <- with(birthweight, ifelse(birthweight$gender==0,"Male","Female"))
dotplot(mygender ~ weight, data = birthweight,
        xlab = "Birthweight", ylab = "Gender")

library(rjags)

model_string <- "
model {
  for (i in 1:N) {
    weight[i] ~ dnorm(mu[gender[i]], tau)
  }
  
  for (j in 1:2) {
    mu[j] ~ dnorm(eta, tau_mu)
  }
  
  tau <- 1 / sigma2
  sigma2 ~ dgamma(a, b)
  
  eta ~ dnorm(0, 1e-6)
  tau_mu <- 1 / tau2
  tau2 ~ dgamma(a_tau, b_tau)
}
"

data_jags <- list(
  weight = birthweight$weight,
  gender = birthweight$gender + 1,
  N = nrow(birthweight),
  a = 0.3,
  b = 0.5,
  a_tau = 0.1,
  b_tau = 0.01
)

parameters <- c("mu", "eta", "sigma2", "tau2")

model <- jags.model(textConnection(model_string), data = data_jags, n.chains = 3)
update(model, 1000)  # Burn-in

samples <- coda.samples(model, variable.names = parameters, n.iter = 10000)

library(coda)
gelman.diag(samples)
autocorr.diag(samples)
traceplot(samples)
                  
plot(samples)                       
                        
# Acceptance-Reject Sampling
acceptance_reject <- function(a, b, b_sigma) {
  repeat {
    sigma2_star <- 1 / rgamma(1, a, b)
    Cmax <- 2 / (pi * b_sigma * (1 + (0 / b_sigma)^2))
    acceptance_prob <- (2 / (pi * b_sigma * (1 + (sigma2_star / b_sigma)^2))) / Cmax
    if (runif(1) < acceptance_prob) return(sigma2_star)
  }
}

# Independence Metropolis-Hastings
independence_mh <- function(a, b, b_sigma, sigma2_current) {
  sigma2_star <- 1 / rgamma(1, a, b)
  acceptance_prob <- min((2 / (pi * b_sigma * (1 + (sigma2_star / b_sigma)^2))) / 
                           (2 / (pi * b_sigma * (1 + (sigma2_current / b_sigma)^2))), 1)
  if (runif(1) < acceptance_prob) return(sigma2_star)
  else return(sigma2_current)
}

# Metropolis-Hastings with Chi-square proposal
mh_with_chisq <- function(a, b, b_sigma, sigma2_current) {
  sigma2_star <- sigma2_current * rchisq(1, df = 2)
  acceptance_prob <- min((2 / (pi * b_sigma * (1 + (sigma2_star / b_sigma)^2))) * 
                           dchisq(sigma2_current, df = 2, ncp = sigma2_star) / 
                           ((2 / (pi * b_sigma * (1 + (sigma2_current / b_sigma)^2))) * 
                              dchisq(sigma2_star, df = 2, ncp = sigma2_current)), 1)
  if (runif(1) < acceptance_prob) return(sigma2_star)
  else return(sigma2_current)
}

# 初始化參數
a <- 0.3
b <- 0.5
b_sigma <- 10
sigma2_current <- 1

# Acceptance-Reject Sampling
sigma2_ar <- acceptance_reject(a, b, b_sigma)

# Independence Metropolis-Hastings
sigma2_imh <- independence_mh(a, b, b_sigma, sigma2_current)

# Metropolis-Hastings with Chi-square proposal
sigma2_mh <- mh_with_chisq(a, b, b_sigma, sigma2_current)

# 輸出結果
cat("Acceptance-Reject:", sigma2_ar, "\n")
cat("Independence Metropolis-Hastings:", sigma2_imh, "\n")
cat("Metropolis-Hastings with Chi-square proposal:", sigma2_mh, "\n")

                        
                        
                        