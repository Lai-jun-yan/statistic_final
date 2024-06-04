rigamma = function(a,b){
  x = 1/rgamma(1,shape = a,rate = b)
  return(x)
}

truncated_Cauchy = function(x,b){
  value = 2/(pi*b*(1+(x/b)^2))
  return(value)
}

AR = function(mu1,mu2,y){
  
  a = 11
  b = (sum((y[y$gender==0,3]-mu1)^2)+sum((y[y$gender==1,3]-mu2)^2))/2
  
  b_sigma = 10000000
  Cmax = truncated_Cauchy(0,b_sigma)
  resample =TRUE
  
  while(resample){
    sigma_star = rigamma(a,b)
    acceptance_probability = truncated_Cauchy(sigma_star,b_sigma)/Cmax
    y = runif(1)
    if(y<acceptance_probability){
      value = sigma_star
      resample = FALSE
    }
  }
  return(value)
  
}

mu_male = c()
mu_female = c()
mu_male[1] = mean(birthweight[birthweight$gender==0,3])
mu_female[1] = mean(birthweight[birthweight$gender==1,3])

sigma = c()
sigma[1] = (var(birthweight[birthweight$gender==0,3])+var(birthweight[birthweight$gender==1,3]))/2

eta = c()
eta[1] = (mu_male[1]+mu_female[1])/2

tau = c()
tau[1] = ((mu_male[1]-eta[1])^2+(mu_female[1]-eta[1])^2)/2

n = 10000
for (i in 2:n) {
  mu_male_normal_mean = (1/(1/(sigma[i-1]/12)+1/tau[i-1]))*(mean(birthweight[birthweight$gender==0,3])/(sigma[i-1]/12)+eta[i-1]/tau[i-1])
  mu_female_normal_mean = (1/(1/(sigma[i-1]/12)+1/tau[i-1]))*(mean(birthweight[birthweight$gender==1,3])/(sigma[i-1]/12)+eta[i-1]/tau[i-1])
  mu_normal_sigma = sqrt((1/(1/(sigma[i-1]/12)+1/tau[i-1])))
  mu_male[i] = rnorm(1,mu_male_normal_mean,mu_normal_sigma)
  mu_female[i] = rnorm(1,mu_female_normal_mean,mu_normal_sigma)
  
  mu_mean = ((mu_male[i-1]+mu_female[i-1])/2)
  eta[i] = rnorm(1,mu_mean,sqrt(tau[i-1]/2))
  
  sigma[i] = AR(mu_male[i-1],mu_female[i-1],birthweight)
  
  a_r = 0.00001
  a = a_r+1
  b_r = 1
  b = ((mu_male[i-1]-eta[i-1])^2+(mu_female[i-1]-eta[i-1])^2)/2 + b_r
  tau[i] = rigamma(a,b)
}

burn_in = 2000

plot(mu_male, type = "l", col = "blue", lwd = 2, main = "Trace of µ1", xlab = "Iterations", ylab = "µ1", xlim = c(burn_in+1, length(mu_male)), ylim = c(0, max(mu_male) * 1.1))

hist(mu_male, breaks = seq(floor(min(mu_male)/1000)*1000, ceiling(max(mu_male)/1000)*1000, by = 100), xlim = c(2000, 4000), ylim = c(0, 0.005), main = "Density of µ2", xlab = "µ1", ylab = "Density",freq = FALSE)

plot(mu_female, type = "l", col = "blue", lwd = 2, main = "Trace of µ2", xlab = "Iterations", ylab = "µ2", xlim = c(burn_in+1, length(mu_female)), ylim = c(0, max(mu_female) * 1.1))

hist(mu_female, breaks = seq(floor(min(mu_female)/1000)*1000, ceiling(max(mu_female)/1000)*1000, by = 100), xlim = c(2000, 4000), ylim = c(0, 0.005), main = "Density of µ2", xlab = "µ2", ylab = "Density",freq = FALSE)

plot(eta, type = "l", col = "blue", lwd = 2, main = "Trace of η", xlab = "Iterations", ylab = "η", xlim = c(burn_in+1, length(eta)), ylim = c(0, max(eta) * 1.1))

hist(eta, breaks = seq(floor(min(eta)/1000)*1000, ceiling(max(eta)/1000)*1000, by = 100), xlim = c(2000, 4000), ylim = c(0, 0.005), main = "Density of η", xlab = "η", ylab = "Density",freq = FALSE)

sigma2 = 1/sigma

plot(sigma2, type = "l", col = "blue", lwd = 2, main = "Trace of 1/σ^2", xlab = "Iterations", ylab = "1/σ^2", xlim = c(burn_in+1, length(sigma2)), ylim = c(0, max(sigma2) * 1.1))

hist(sigma2, breaks = seq(min(sigma2)-1.0e-06, max(sigma2)+1.0e-06, by = 5.0e-06), xlim = c(min(sigma2), max(sigma2)), ylim = c(0, 1.5e+05), main = "Density of 1/σ^2", xlab = "sigma2", ylab = "Density",freq = FALSE)

tau2 = 1/tau

plot(tau2, type = "l", col = "blue", lwd = 2, main = "Trace of 1/τ^2", xlab = "Iterations", ylab = "1/τ^2", xlim = c(burn_in+1, length(tau2)), ylim = c(0, max(tau2) * 1.1))

hist(tau2, breaks = seq(min(tau2)-sd(tau2), max(tau2)+sd(tau2), by = sd(tau2)), xlim = c(mean(tau2)-3*sd(tau2), max(tau2)), ylim = c(0, 200), main = "Density of 1/τ^2", xlab = "1/τ^2", ylab = "Density",freq = FALSE)
