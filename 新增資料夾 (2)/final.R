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
    weight[i] ~ dnorm(mu[gender[i]], sigma2)
  }
  
  for (j in 1:2) {
    mu[j] ~ dnorm(eta, tau2)
  }
  
  eta ~ dnorm(0, 1e-6)
  
  sigma2 <- abs(1 / (sigma+0.00001))
  sigma ~ dt(0,b_s,1)
  
  tau2 <- abs(1 / (tau+0.00001))
  tau ~ dt(0,b_r,1)
}
"

data_jags <- list(
  weight = birthweight$weight,
  gender = birthweight$gender + 1,
  N = nrow(birthweight),
  b_r = 0.00001,
  b_s = 0.00001
)

parameters <- c("mu", "eta", "sigma2", "tau2")

model <- jags.model(textConnection(model_string), data = data_jags, n.chains = 3)
update(model, 2000)  # Burn-in

samples <- coda.samples(model, variable.names = parameters, n.iter = 10000)

library(coda)
gelman.diag(samples)
autocorr.diag(samples)
traceplot(samples)
plot(samples)          
