sim_pow_parc_emp <- function(nrep = 100, nstu = 100, nitem = 40, nct = 20, cip =.9, cjp =.9){
  out <- rep(NA, nrep)
  J <- nitem
  I <- nstu
  guess <- rep(0, J)
  for (i in 1:nrep){
    mu <- runif(1, min = -2, max = 2)
    sigma <- runif(1, min = .8, max = 1.2)
    theta <- rnorm(I,mu,sigma)
    a <- runif(J,0.5,2)
    b <- runif(J,-2,2)
    P <- getirt_prob(theta,a,b,guess)
    response <- (P > runif(I*J,0,1)) * 1 
    a2 <- a[1:nct]
    a1 <- a[-(1:nct)]
    b2 <- b[1:nct]
    b1 <- b[-(1:nct)]
    guess2 <- guess[1:nct]
    guess1 <- guess[-(1:nct)]
    y2 <- response[1:nct,]
    y1 <- response[-(1:nct),]
    dimindex <- expand.grid(sample(I,I*cip), sample(nct,nct*cjp))
    for (j in 1:NROW(dimindex)){
      index <- dimindex[j,]
      y2[as.numeric(index[2]),as.numeric(index[1])] <- 1
    }
    sigmahat <- getsigmahat(a, b, guess, y=response)
    temp <- getWald_byemp_withsigma(sigma = sigmahat, a1, a2, b1, b2, guess1, guess2, y1, y2) 
    out[i] <- ifelse(temp >= 1.96 || temp <= -1.96, 1, 0)
  }
  return(mean(out))
}

sim_pow_thetashift_emp <- function(nrep = 100, nstu = 100, nitem = 40, nct = 20, shift = .1, mu, sigma){
  out <- rep(NA, nrep)
  J <- nitem
  I <- nstu
  guess <- rep(0, J)
  for (i in 1:nrep){
    # mu <- runif(1, min = -2, max = 2)
    # sigma <- runif(1, min = .8, max = 1.2)
    theta <- rnorm(I,mu,sigma)
    a <- runif(J,0.5,2)
    b <- runif(J,-2,2)
    a2 <- a[1:nct]
    a1 <- a[-(1:nct)]
    b2 <- b[1:nct]
    b1 <- b[-(1:nct)]
    guess2 <- guess[1:nct]
    guess1 <- guess[-(1:nct)]
    theta1 <- theta 
    theta2 <- theta + shift
    P1 <- getirt_prob(theta1,a1,b1,guess1)
    P2 <- getirt_prob(theta2,a2,b2,guess2)
    y1 <- (P1 > runif(I*(J-nct),0,1)) * 1 
    y2 <- (P2 > runif(I*nct,0,1)) * 1
    sigmahat <- getsigmahat(a, b, guess, y=cbind(y1,y2))
    temp <- getWald_byemp_withsigma(sigma = sigmahat, a1, a2, b1, b2, guess1, guess2, y1, y2) 
    out[i] <- ifelse(temp >= 1.96 || temp <= -1.96, 1, 0)
  }
  return(mean(out))
}


sim_t1e_emp <- function(nrep = 100, nstu = 100, nitem = 40, nct = 20){
  out <- rep(NA, nrep)
  J <- nitem
  I <- nstu
  guess <- rep(0, J)
  for (i in 1:nrep){
    mu <- runif(1, min = -2, max = 2)
    sigma <- runif(1, min = .8, max = 1.2)
    theta <- rnorm(I,mu,sigma)
    a <- runif(J,0.5,2)
    b <- runif(J,-2,2)
    P <- getirt_prob(theta,a,b,guess)
    response <- (P > runif(I*J,0,1)) * 1 
    a2 <- a[1:nct]
    a1 <- a[-(1:nct)]
    b2 <- b[1:nct]
    b1 <- b[-(1:nct)]
    guess2 <- guess[1:nct]
    guess1 <- guess[-(1:nct)]
    y2 <- response[1:nct,]
    y1 <- response[-(1:nct),]
    # dimindex <- expand.grid(sample(nct,nct*cip), sample(nct,nct*cjp))
    # apply(dimindex, 1, function(x){y2[x[2], x[1]] <- 1})
    sigmahat <- getsigmahat(a, b, guess, y=response)
    temp <- getWald_byemp_withsigma(sigma = sigmahat, a1, a2, b1, b2, guess1, guess2, y1, y2) 
    out[i] <- ifelse(temp > 1.96 || temp < -1.96, 1, 0)
  }
  return(mean(out))
}

sim_t1e_pp <- function(nrep = 100, nstu = 100, nitem = 40, nct = 20){
  out <- rep(NA, nrep)
  J <- nitem
  I <- nstu
  guess <- rep(0, J)
  for (i in 1:nrep){
    mu <- runif(1, min = -2, max = 2)
    sigma <- runif(1, min = .8, max = 1.2)
    theta <- rnorm(I,mu,sigma)
    a <- runif(J,0.5,2)
    b <- runif(J,-2,2)
    P <- getirt_prob(theta,a,b,guess)
    response <- (P > runif(I*J,0,1)) * 1 
    a2 <- a[1:nct]
    a1 <- a[-(1:nct)]
    b2 <- b[1:nct]
    b1 <- b[-(1:nct)]
    guess2 <- guess[1:nct]
    guess1 <- guess[-(1:nct)]
    y2 <- response[1:nct,]
    y1 <- response[-(1:nct),]
    # dimindex <- expand.grid(sample(nct,nct*cip), sample(nct,nct*cjp))
    # apply(dimindex, 1, function(x){y2[x[2], x[1]] <- 1})
    sigmahat <- getsigmahat(a, b, guess, y=response)
    temp <- getWald_bypp_withsigma(sigma = sigmahat, a1, a2, b1, b2, guess1, guess2, y1, y2) 
    out[i] <- ifelse(temp > 1.96 || temp < -1.96, 1, 0)
  }
  return(mean(out))
}

shifts <- 0.1*c(1:10)
pow_shift_100 <- sapply(shifts, function(x){sim_pow_thetashift_emp(shift = x)})

t1e1000 <- sim_t1e_pp(nrep = 1000)

