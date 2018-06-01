library(numDeriv)
library(statmod)

getirt_prob <- function(theta,a,b,guess){
  p <- guess+(1-guess)*(1 / (1 + exp(a*outer(b,theta,"-"))))
  return (p)
}

irt_group <- function(mu,sigma,a,b,guess){
  p <- guess + (1-guess)*pnorm(outer(mu,b,"-")/sqrt((1/a)^2 + (1.7*sigma)^2))
  return (as.vector(p))
}

likperson <- function(theta,a,b,guess,y){
  P <- getirt_prob(theta,a,b,guess)
  temp <- y * log(P) + (1-y) * log(1-P)
  loglik_vec <- apply(temp,2,sum)
  return (exp(loglik_vec))
}

integrand <- function(theta,mu,sigma,a,b,guess,y){
  g <- likperson(theta,a,b,guess,y)
  out <- g * dnorm(theta,mu,sigma)
  return (out)
}

marginal_likelihood_person <- function(mu,sigma,a,b,guess,y){
  out <- rep(NA,dim(y)[2])
  for (i in 1:dim(y)[2]){
    out[i] <- integrate(integrand,
                        mu=mu,
                        sigma=sigma,
                        a=a,b=b,
                        guess=guess,
                        y=y[,i],
                        lower=-10,upper=10)$value
  }
  return (out)
}

# marginal_likelihood_person <- function(mu,sigma,a,b,guess,y){
#   Ns <- gauss.quad.prob(17, dist = "normal", mu = mu, sigma = sigma)$nodes
#   Ws <- gauss.quad.prob(17, dist = "normal", mu = mu, sigma = sigma)$weights
#   P <- getirt_prob(theta = Ns, a = a, b = b,guess = guess)
#   out <- apply(P, 2,
#                function(x){exp(colSums(y*log(x)+(1-y)*log(1-x)))}
#                )
#   return (colSums(t(out)*Ws))
# }

marginal_loglik_person <- function(mu,sigma,a,b,guess,y){
  out <- log(integrate(integrand,mu=mu,sigma=sigma,a=a,b=b,guess=guess,y=y,lower=-10,upper=10)$value)
  return (out)
}


marginal_likelihood <- function(mu,sigma,a,b,guess,y){
  temp <- marginal_likelihood_person(mu,sigma,a,b,guess,y)
  out <- sum(log(temp))
  return (out)
}

marginal_likelihood2 <- function(par,a,b,guess,y){
  temp <- marginal_likelihood_person(par[1],par[2],a,b,guess,y)
  out <- sum(log(temp))
  return (out)
}

getmuhat <- function(a, b, guess, y, sigmahat){
  out <- optimize(f = marginal_likelihood, interval = c(-3,3),sigma=sigmahat,a=a,b=b,guess=guess,y=y,maximum = T)$maximum
  return(out)
}

getsigmahat <- function(a, b, guess, y){
  out <- optim(par = c(0,1), fn = marginal_likelihood2, 
               lower=c(-4, 0.01), upper=c(4, 5), method="L-BFGS-B", 
               a=a, b=b, guess=guess, y=y, 
               control=list(fnscale=-1))$par[2]
  return(out)
}

getWald_bypp <- function(muhat, sigma, a1, a2, b1, b2, guess1, guess2, y1, y2){
  l1p <- grad(marginal_likelihood, x = muhat, sigma=sigma,a=a1,b=b1,guess=guess1,y=y1)
  l1pp <- hessian(marginal_likelihood, x = muhat, sigma=sigma,a=a1,b=b1,guess=guess1,y=y1)
  l2p <- grad(marginal_likelihood, x = muhat, sigma=sigma,a=a2,b=b2,guess=guess2,y=y2)
  l2pp <- hessian(marginal_likelihood, x = muhat, sigma=sigma,a=a2,b=b2,guess=guess2,y=y2)
  I <- dim(y1)[2]
  l1ps <- rep(NA, I)
  l2ps <- rep(NA, I)
  # l1pps <- rep(NA, I)
  # l2pps <- rep(NA, I)
  for(i in 1:I){
    l1ps[i] <- grad(marginal_loglik_person, x=muhat,sigma=sigma,a=a1,b=b1,guess=guess1,y=y1[,i])
    l2ps[i] <- grad(marginal_loglik_person, x=muhat,sigma=sigma,a=a2,b=b2,guess=guess2,y=y2[,i])
    # l1pps[i] <- hessian(marginal_loglik_person, x=muhat,sigma=sigma,a=a1,b=b1,guess=guess1,y=y1[,i])
    # l2pps[i] <- hessian(marginal_loglik_person, x=muhat,sigma=sigma,a=a2,b=b2,guess=guess2,y=y2[,i])
  }
  Num <- -(l2p/l2pp) + (l1p/l1pp)
  # Den <- sqrt(I*var(l2ps)/(l2pp)^2 + I*var(l1ps)/(l1pp)^2 - I*cov(l1ps,l2ps)/(l1pp*l2pp))
  Den <- sqrt((-l2pp)^(-1) + (-l1pp)^(-1) - 2*I*cov(l1ps,l2ps)/(l1pp*l2pp))
  Wald <- Num/Den
  return(as.numeric(Wald))
}

getWald_byemp <- function(sigma, a1, a2, b1, b2, guess1, guess2, y1, y2){
  muhat <- getmuhat(a1, b1, guess1, y1, sigma)
  l1p <- grad(marginal_likelihood, x = muhat, sigma=sigma,a=a1,b=b1,guess=guess1,y=y1)
  l2p <- grad(marginal_likelihood, x = muhat, sigma=sigma,a=a2,b=b2,guess=guess2,y=y2)
  I <- dim(y1)[2]
  l1ps <- rep(NA, I)
  l2ps <- rep(NA, I)
  # l1pps <- rep(NA, I)
  # l2pps <- rep(NA, I)
  for(i in 1:I){
    l1ps[i] <- grad(marginal_loglik_person, x=muhat,sigma=sigma,a=a1,b=b1,guess=guess1,y=y1[,i])
    l2ps[i] <- grad(marginal_loglik_person, x=muhat,sigma=sigma,a=a2,b=b2,guess=guess2,y=y2[,i])
    # l1pps[i] <- hessian(marginal_loglik_person, x=muhat,sigma=sigma,a=a1,b=b1,guess=guess1,y=y1[,i])
    # l2pps[i] <- hessian(marginal_loglik_person, x=muhat,sigma=sigma,a=a2,b=b2,guess=guess2,y=y2[,i])
  }
  Num <- -(l2p/(-I*var(l2ps))) + (l1p/(-I*var(l1ps)))
  Den <- sqrt((I*var(l2ps))^(-1) + (I*var(l1ps))^(-1) - 2*cov(l1ps,l2ps)/(I*var(l2ps)*var(l1ps)))
  Wald <- Num/Den
  return(as.numeric(Wald))
}

getWald_byemp_withsigma <- function(sigma, a1, a2, b1, b2, guess1, guess2, y1, y2){
  muhat1 <- getmuhat(a1, b1, guess1, y1, sigma)
  muhat2 <- getmuhat(a2, b2, guess2, y2, sigma)
  l1p <- grad(marginal_likelihood, x = muhat1, sigma=sigma,a=a1,b=b1,guess=guess1,y=y1)
  l2p <- grad(marginal_likelihood, x = muhat1, sigma=sigma,a=a2,b=b2,guess=guess2,y=y2)
  I <- dim(y1)[2]
  l1ps <- rep(NA, I)
  l2ps <- rep(NA, I)
  # l1pps <- rep(NA, I)
  # l2pps <- rep(NA, I)
  for(i in 1:I){
    l1ps[i] <- grad(marginal_loglik_person, x=muhat1,sigma=sigma,a=a1,b=b1,guess=guess1,y=y1[,i])
    l2ps[i] <- grad(marginal_loglik_person, x=muhat1,sigma=sigma,a=a2,b=b2,guess=guess2,y=y2[,i])
    # l1pps[i] <- hessian(marginal_loglik_person, x=muhat,sigma=sigma,a=a1,b=b1,guess=guess1,y=y1[,i])
    # l2pps[i] <- hessian(marginal_loglik_person, x=muhat,sigma=sigma,a=a2,b=b2,guess=guess2,y=y2[,i])
  }
  Num <- muhat2 -muhat1
  Den <- sqrt((I*var(l2ps))^(-1) + (I*var(l1ps))^(-1) - 2*cov(l1ps,l2ps)/(I*var(l2ps)*var(l1ps)))
  Wald <- Num/Den
  return(as.numeric(Wald))
}

getWald_bypp_withsigma <- function(sigma, a1, a2, b1, b2, guess1, guess2, y1, y2){
  muhat1 <- getmuhat(a1, b1, guess1, y1, sigma)
  muhat2 <- getmuhat(a2, b2, guess2, y2, sigma)
  l1p <- grad(marginal_likelihood, x = muhat1, sigma=sigma,a=a1,b=b1,guess=guess1,y=y1)
  l2p <- grad(marginal_likelihood, x = muhat1, sigma=sigma,a=a2,b=b2,guess=guess2,y=y2)
  l1pp <- hessian(marginal_likelihood, x = muhat1, sigma=sigma,a=a1,b=b1,guess=guess1,y=y1)
  l2pp <- hessian(marginal_likelihood, x = muhat1, sigma=sigma,a=a2,b=b2,guess=guess2,y=y2)
  I <- dim(y1)[2]
  l1ps <- rep(NA, I)
  l2ps <- rep(NA, I)
  # l1pps <- rep(NA, I)
  # l2pps <- rep(NA, I)
  for(i in 1:I){
    l1ps[i] <- grad(marginal_loglik_person, x=muhat1,sigma=sigma,a=a1,b=b1,guess=guess1,y=y1[,i])
    l2ps[i] <- grad(marginal_loglik_person, x=muhat1,sigma=sigma,a=a2,b=b2,guess=guess2,y=y2[,i])
    # l1pps[i] <- hessian(marginal_loglik_person, x=muhat1,sigma=sigma,a=a1,b=b1,guess=guess1,y=y1[,i])
    # l2pps[i] <- hessian(marginal_loglik_person, x=muhat1,sigma=sigma,a=a2,b=b2,guess=guess2,y=y2[,i])
  }
  Num <- muhat2 -muhat1
  Den <- sqrt(I*var(l2ps)*(l2pp)^(-2) + I*var(l1ps)*(l1pp)^(-2) - 2*I*cov(l1ps,l2ps)/(l1pp*l2pp))
  Wald <- Num/Den
  return(as.numeric(Wald))
}

sim_t1e_emp <- function(nrep = 100, nstu = 100, nitem = 40, nct = 20, mu, sigma1, sigma2){
  out <- rep(NA, nrep)
  J <- nitem
  I <- nstu
  guess <- rep(0, J)
  for (i in 1:nrep){
    # mu <- runif(1, min = -2, max = 2)
    # sigma <- runif(1, min = .8, max = 1.2)
    theta1 <- rnorm(I,mu,sigma1)
    theta2 <- ((theta1 - mu)/sigma1)*sigma2 + mu
    a <- runif(J,0.5,2)
    b <- runif(J,-2,2)
    a2 <- a[1:nct]
    a1 <- a[-(1:nct)]
    b2 <- b[1:nct]
    b1 <- b[-(1:nct)]
    guess2 <- guess[1:nct]
    guess1 <- guess[-(1:nct)]
    P1 <- getirt_prob(theta1,a1,b1,guess1)
    y1 <- (P1 > runif(I*(J-nct),0,1)) * 1  
    P2 <- getirt_prob(theta2,a2,b2,guess2)
    y2 <- (P2 > runif(I*nct,0,1)) * 1 
    
    # dimindex <- expand.grid(sample(nct,nct*cip), sample(nct,nct*cjp))
    # apply(dimindex, 1, function(x){y2[x[2], x[1]] <- 1})
    sigmahat <- getsigmahat(a, b, guess, y=rbind(y2, y1))
    temp <- getWald_byemp_withsigma(sigma = sigmahat, a1, a2, b1, b2, guess1, guess2, y1, y2) 
    out[i] <- temp
    # out[i] <- ifelse(temp > 1.96 || temp < -1.96, 1, 0)
  }
  return(out)
}

sim_pow_thetashift_emp <- function(nrep = 100, nstu = 100, nitem = 40, nct = 20, shift = .1, mu, sigma1, sigma2){
  out <- rep(0, nrep)
  J <- nitem
  I <- nstu
  guess <- rep(0, J)
  for (i in 1:nrep){
    # mu <- runif(1, min = -2, max = 2)
    # sigma <- runif(1, min = .8, max = 1.2)
    theta1 <- rnorm(I,mu,sigma1)
    theta2 <- ((theta1 - mu)/sigma1)*sigma2 + mu
    a <- runif(J,0.5,2)
    b <- runif(J,-2,2)
    a2 <- a[1:nct]
    a1 <- a[-(1:nct)]
    b2 <- b[1:nct]
    b1 <- b[-(1:nct)]
    guess2 <- guess[1:nct]
    guess1 <- guess[-(1:nct)]
    P1 <- getirt_prob(theta1,a1,b1,guess1)
    y1 <- (P1 > runif(I*(J-nct),0,1)) * 1  
    P2 <- getirt_prob(theta2,a2,b2,guess2)
    y2 <- (P2 > runif(I*nct,0,1)) * 1
    sigmahat <- getsigmahat(a, b, guess, y = rbind(y2, y1))
    # sigmahat <- sigma
    # sigmahat <- getsigmahat(a1, b1, guess1, y=y1)
    temp <- getWald_byemp_withsigma(sigma = sigmahat, a1, a2, b1, b2, guess1, guess2, y1, y2) 
    out[i] <- temp
    # out[i] <- ifelse(temp > 1.96 || temp < -1.96, 1, 0)
  }
  # return(mean(out))
  return(out)
}

##################################################################################################

# library(parallel)
# 
# simCon_pow1 <- expand.grid(effsize = c(0.05, 0.15, 0.25), mu = c(-1,0,1), sigma1 = 1, sigma2 = .5, nct = c(5, 10, 20))
# simCon_pow2 <- expand.grid(effsize = c(0.05, 0.15, 0.25), mu = c(-1,0,1), sigma1 = .5, sigma2 = .5, nct = c(5, 10))
# simCon_pow3 <- expand.grid(effsize = c(0.05, 0.15, 0.25), mu = c(-1,0,1), sigma1 = 1, sigma2 = 1, nct = c(5, 10))
# 
# cl <- makeCluster(2)
# clusterExport(cl,ls())
# clusterExport(cl,"grad")
# clusterExport(cl,"hessian")
# 
# sim_pow_thetashift_emp(nrep=2 ,shift=0.05, nstu = 100, 
#                        nitem = 40, nct = 5,
#                        mu=0, sigma1 = 1, sigma2 = 2)
# 
# pow_emp_shift_1_5 <- parRapply(cl, simCon[c(1, 10),],
#                                FUN = function(x)
#                                  sim_pow_thetashift_emp(nrep=2 ,shift=x[1], nstu = 100, 
#                                                         nitem = 40, nct = x[5],
#                                                         mu=x[2], sigma1 = x[3], sigma2 = x[4]))
# 
# stopCluster(cl)
# 
# cl <- makeCluster(5)
# clusterExport(cl,ls())
# clusterExport(cl,"grad")
# clusterExport(cl,"hessian")
# 
# pow_emp_shift_6_10 <- parRapply(cl, simCon[21:25,],
#                                 FUN = function(x)
#                                   sim_pow_thetashift_emp(nrep=300 ,shift=x[1], nstu = 100, 
#                                                          nitem = 40, nct = 20,
#                                                          mu=x[2], sigma = x[3]))
# 
# 
# stopCluster(cl)


##################################################################################################

simCon_t1e <- expand.grid(nct = c(5, 10, 20), mu = c(-1,0,1), sigma2 = c(.75, .5, .25))

cl <- makeCluster(10)
clusterExport(cl,ls())
clusterExport(cl,"grad")
clusterExport(cl,"hessian")

t1e_sim1_10 <- parRapply(cl, simCon_t1e[1:10,],
                               FUN = function(x)
                                 sim_t1e_emp(nrep=3000, nstu = 100, 
                                             nitem = 40, nct = x[1],
                                             mu=x[2], sigma1 = 1, sigma2 = x[3]))

stopCluster(cl)

cl <- makeCluster(10)
clusterExport(cl,ls())
clusterExport(cl,"grad")
clusterExport(cl,"hessian")

t1e_sim11_20 <- parRapply(cl, simCon_t1e[11:20,],
                     FUN = function(x)
                       sim_t1e_emp(nrep=3000, nstu = 100, 
                                   nitem = 40, nct = x[1],
                                   mu=x[2], sigma1 = 1, sigma2 = x[3]))

stopCluster(cl)
    
# cl <- makeCluster(7)
# clusterExport(cl,ls())
# clusterExport(cl,"grad")
# clusterExport(cl,"hessian")
# 
# t1e_sim21_27 <- parRapply(cl, simCon_t1e[21:27,],
#                           FUN = function(x)
#                             sim_t1e_emp(nrep=3000, nstu = 100, 
#                                         nitem = 40, nct = x[1],
#                                         mu=x[2], sigma1 = 1, sigma2 = x[3]))
# 
# stopCluster(cl)



######################################################################################################


