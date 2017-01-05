prep.models <- function(sub, cond){
  
  I <- length(unique(sub))
  R <- length(sub)
  X.full <- matrix(nrow = R, ncol = 2 * I + 2, 0)
  for (r in 1:R){
    X.full[r, 1] <- 1
    X.full[r, sub[r] + 1] <- 1
    if (cond[r] == 2) {
      X.full[r, I + 2] <- 1
      X.full[r, I + 2 + sub[r]] <- 1}}
  
  gMap.full <- c(rep(0, I), 1, rep(2, I))
  
  X.one <- matrix(nrow=R,ncol=I+2,0)
  for (r in 1:R){
    X.one[r, 1] <- 1
    X.one[r, sub[r] + 1] <- 1
    if (cond[r] == 2) {
      X.one[r, I + 2] <- 1
    }}
  
  gMap.one <- c(rep(0, I), 1)
  
  X.null <- matrix(nrow = R, ncol = I + 1, 0)
  for(r in 1:R){
    X.null[r, 1] <- 1
    X.null[r, sub[r] + 1] <- 1
  }
  
  gMap.null <- rep(0, I)
  
  if(length(cond) != length(sub)) return(print("oops, your condition vector does not match your subject vector"))
  return(list(X.full = X.full
              , gMap.full = gMap.full
              , X.one = X.one
              , gMap.one = gMap.one
              , X.null = X.null
              , gMap.null = gMap.null
              , R = R
              , I= I))
}

ratio.greater <- function(M, I, a = alpha, b = beta, sd_mu = sigma_mu){
  
  x <- matrix(nrow = M, ncol = I)
  s2 <- rinvgamma(M, a, b)
  mu <- rcauchy(M, 0, sd_mu)
  for (i in 1:I)
    x[, i] <- rnorm(M, mu, sqrt(s2))
  
  all.greater <- function(x) as.integer(mean(x > 0) == 1)
  return(1/mean(apply(x, 1, all.greater)))
}

makeBF <- function(y, meanScale, effectScale, prep = prep.1, keep = 1001:10000)
{
  mcmc.full <- nWayAOV(y
                       , prep$X.full
                       , prep$gMap.full
                       , rscale = c(1, meanScale, effectScale)
                       , posterior = T
                       , method = "auto"
                       , iterations = max(keep))
  bf.full <- nWayAOV(y
                     , prep$X.full
                     , prep$gMap.full
                     , rscale = c(1, meanScale, effectScale)
                     , posterior = F
                     , method = "auto"
                     , iterations = max(keep))
  bf.one <- nWayAOV(y
                    , prep$X.one
                    , prep$gMap.one
                    , rscale = c(1, meanScale)
                    , posterior = F
                    , method = "auto"
                    , iterations = max(keep))
  mcmc.one <- nWayAOV(y
                      , prep$X.one
                      , prep$gMap.one
                      , rscale = c(1, meanScale)
                      , posterior = T
                      , method = "auto"
                      , iterations = max(keep))
  bf.null <- nWayAOV(y
                       , prep$X.null
                       , prep$gMap.null
                       , rscale = 1
                       , posterior = F
                       , method = "auto"
                     , iterations = max(keep))
  
  i.theta0 <- prep$I + 2
  i.theta <- (prep$I + 3):(2 * prep$I + 2)
  
  myTheta <- mcmc.full[keep, i.theta] + mcmc.full[keep, i.theta0]
  good <- myTheta > 0
  all.good <- apply(good, 1, mean)
  PostCount <- mean(all.good == 1)
  
  #prior settings
  R <- max(keep)
  beta <- .5 * effectScale^2
  alpha <- .5
  mu.theta.sd <- .5 * meanScale
  
  x <- matrix(nrow = R, ncol = prep$I)
  s2 <- rinvgamma(R, alpha, beta)
  mu <- rcauchy(R, 0, mu.theta.sd)
  for (i in 1:prep$I)
    x[,i] <- rnorm(R, mu, sqrt(s2))
  
  all.greater <- function(x) as.integer(mean(x > 0) == 1)
  PriorCount <- mean(apply(x, 1, all.greater))
  bf.FP <- PriorCount/PostCount
  bf.F0 <- exp(bf.full$bf - bf.null$bf)
  bf.F1 <- exp(bf.full$bf - bf.one$bf)
  
  m <- apply(myTheta, 2, mean)
  new.sd <- sd(m)
  new.mean <- mean(m)
  return(list(mean = new.mean, sd = new.sd, bf.1f = 1/ bf.F1, bf.pf =  1/ bf.FP, bf.0f = 1 / bf.F0, bf.full <- bf.full, bf.one = bf.one, bf.null = bf.null, est.full = mcmc.full, est.one = mcmc.one, prior.c = PriorCount, post.c = PostCount))
}