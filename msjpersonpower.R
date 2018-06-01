getirt_prob <- function(theta,a,b,guess){
  p <- guess+(1-guess)*(1 / (1 + exp(a*outer(b,theta,"-"))))
  return (p)
}

getloglik <-  function(theta, y, a, b, guess)
{ #:: x is the theta
  prob <- getirt_prob(theta, a, b, guess)
  return(colSums(y*log(prob) + (1 - y)*log(1-prob)))
}

getIRTprobprime <- function(a, b, guess, theta) {
  probprime <- (1-guess)*(a*exp(a*(outer(b, theta, "-"))))/((1 + exp(a*(outer(b, theta, "-"))))^2)
  return(probprime)
}

getThetaHat <- function(y, a, b, guess){
  #:: pattern could be a matrix
  theta_hat <- NULL
  for (n in 1 : NCOL(y))
  {
    temp <- optim(par = 0, fn = getloglik,
                  y = y[, n],
                  a = a, b = b, guess = guess,
                  lower = -3.5, upper = 3.5,
                  method = "L-BFGS-B",
                  control = list(fnscale = -1))
    theta_hat <- c(theta_hat,temp$par)
  }
  return(theta_hat)
}

library(PerFit)

getNaivelz_PerFit <- function(res, itempar){
  lzstat <- PerFit::lz(matrix = res, IP = itempar)
  lzstat <- lzstat$PFscores
  return(lzstat)
}

# getNaivelz <- function(p_hat, y){
#   prob <- p_hat
#   responsePat <- y
#   w <- (-1)*log(prob/(1-prob))
#   l <- colSums((prob - responsePat) * w)
#   var_l <- colSums((w^2)*prob*(1-prob))
#   lz <- l / (var_l ^ (1/2))
#   return(lz)
# }
# 
# getlz <- function(p_hat, y, prime, nitem){
#   l_pat <- nitem
#   responsePat <- y
#   prob <- p_hat
#   probprime <- prime
#   w <- (-1)*log(prob/(1-prob))
#   # /
#   # t(outer(colSums((prob^2)*((1-prob)^2)*((prob - t(outer(colMeans(prob), rep(1, l_pat))))^2))^(1/2), rep(1, l_pat)))
#   l <- colSums((prob - responsePat) * w)
#   r <- probprime / (prob*(1-prob))
#   c_n <- t(outer(colSums(probprime * w) / colSums(probprime * r), rep(1, l_pat)))  
#   w_t <- w - c_n*r
#   var_l <- colSums((w_t^2)*prob*(1-prob))
#   lz <- l / (var_l ^ (1/2))
#   return(lz)
# }

sim_pow_parc_emp_bonefor <- function(nrep = 100, nstu = 100, nitem = 40, nct = 20, cip =.9, cjp =.9,
                                     mu = 0, sigma = 1){
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
    y <- rbind(y2, y1)
    #thetahat <- getThetaHat(y=y1, a1, b1, guess1)
    # p_hat <- getirt_prob(thetahat, a, b, guess)
    # p_hat <- P
    cfinv <- qnorm(1-0.05/I)
    # prime <- getIRTprobprime(a, b, guess, theta)
    # lzstat <- getlz(p_hat, y, prime, nitem = J)
    # lzstat <- getNaivelz(p_hat, y)
    lzstat <- getNaivelz_PerFit(res = t(y), itempar = cbind(a, b, guess))
    out[i] <- ifelse(any(lzstat < -cfinv), 1, 0)
  }
  outlist <- list(out, mean(out))
  return(outlist)
}

library(parallel)

simCon <- expand.grid(cip = c(.5,.7,.9), cjp = (1:6)*0.05, mu = c(-1, 0), 
                      sigma = 1)



cl <- makeCluster(4)
clusterExport(cl, ls())
junk <- clusterEvalQ(cl, library(PerFit))


pow_emp_prac_bonefor <- parRapply(cl, simCon,
                          FUN = function(x)
                            sim_pow_parc_emp_bonefor(nrep=3000, nstu = 100, 
                                             nitem = 100, nct = 40,
                                             mu=x[3], sigma = x[4], 
                                             cip = x[1], cjp = x[2]))
stopCluster(cl)
saveRDS(pow_emp_prac_bonefor, "/Users/zhuangzhuanghan/Rwd/Rscripts/GroupDetect/pow_emp_prac_bonefor2.rds")

powpractable <- NULL
for(i in 1: 36){
  temppow <- pow_emp_prac_bonefor[[i]][[2]]
  powpractable <- c(powpractable, temppow)
}
powpractable <- cbind(simCon, powpractable)

powprac <-  cbind(powpractable[1:9, 1:2], powpractable[1:9, 5], powpractable[10:18, 5])


cl <- makeCluster(5)
clusterExport(cl,ls())


pow_emp_prac_bonefor2 <- parRapply(cl, simCon,
                                  FUN = function(x)
                                    sim_pow_parc_emp_bonefor(nrep=300, nstu = 100, 
                                                             nitem = 40, nct = 20,
                                                             mu=x[3], sigma = x[4], 
                                                             cip = x[1], cjp = x[2]))
stopCluster(cl)

sim_lz <- function(nstu = 100, nitem = 40, nct = 20, mu = 0, sigma = 1){
  J <- nitem
  I <- nstu
  guess <- rep(0, J)
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
  # dimindex <- expand.grid(sample(I,I*cip), sample(nct,nct*cjp))
  # for (j in 1:NROW(dimindex)){
  #   index <- dimindex[j,]
  #   y2[as.numeric(index[2]),as.numeric(index[1])] <- 1
  # }
  thetahat <- getThetaHat(y=y1, a1, b1, guess1)
  # p_hat <- getirt_prob(thetahat, a, b, guess)
  p_hat <- P
  cfinv <- qnorm(1-0.025/I)
  # prime <- getIRTprobprime(a, b, guess, thetahat)
  # lz <- getlz(p_hat, response, prime, nitem = J)
  lz <- getNaivelz(p_hat, response)
  return(lz)
}