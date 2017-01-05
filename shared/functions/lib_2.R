MeanNormalTruncated<- function(mu = 0, sigma = 1, a = -Inf, b = Inf){
  mu + sigma*(dnorm((a-mu)/sigma) - dnorm((b-mu)/sigma)) / (pnorm((b-mu)/sigma) - pnorm((a-mu)/sigma))
}

VarianceNormalTruncated<- function(mu = 0, sigma = 1, a = -Inf, b = Inf){
  sigma^2*(
    1
    + (((a-mu)/sigma)*dnorm((a-mu)/sigma) - ((b-mu)/sigma)*dnorm((b-mu)/sigma)) / (pnorm((b-mu)/sigma) - pnorm((a-mu)/sigma))
    - ((dnorm((a-mu)/sigma) - dnorm((b-mu)/sigma)) / (pnorm((b-mu)/sigma) - pnorm((a-mu)/sigma)))^2
    )
}

var.tnorm<-function(mu,sd,lower,upper){
  ##return the variance of a truncated normal distribution
  lower.std=(lower-mu)/sd
  upper.std=(upper-mu)/sd
  variance=sd^2*(1+(lower.std*dnorm(lower.std)-upper.std*dnorm(upper.std))/
                   (pnorm(upper.std)-pnorm(lower.std))-((dnorm(lower.std)-dnorm(upper.std))/
                                                          (pnorm(upper.std)-pnorm(lower.std)))^2)
  return(variance)
}

mean.tnorm<-function(mu,sd,lower,upper){
  ##return the expectation of a truncated normal distribution
  lower.std=(lower-mu)/sd
  upper.std=(upper-mu)/sd
  mean=mu+sd*(dnorm(lower.std)-dnorm(upper.std))/
    (pnorm(upper.std)-pnorm(lower.std))
  return(mean)
}

search.var <- function(boundaries, optimum, cutoff, mu = 0, a, b = Inf){
  n <- 1
  var.opt <- optimum^2
  repeat{
    val <- mean( boundaries )
    x <- VarianceNormalTruncated(mu, val, a, b)
    if(abs(x - var.opt) < cutoff){break}
    if(x < var.opt){ boundaries <- c(val, boundaries[2]) }
    else { boundaries <- c(boundaries[1], val) }
    n <- n+1
    if(n > 100000){
      print("blob")
      break}
  }
  return(
    val)
}


search.means <- function(boundaries, optimum, cutoff, sigma = 1, a, b = Inf){
  n <- 1
  repeat{
    val <- mean(boundaries)
    x <- MeanNormalTruncated(val, sigma, a, b)
    if(abs(x - optimum) < cutoff){break}
    if(x < optimum){boundaries <- c(val, boundaries[2])}
    else{boundaries <- c(boundaries[1], val)}
    n <- n+1
    if(n > 100000){
      print("blob")
      break}
  }
  return(val)
}


produce.gif <- function(lower){
  for(i in 1:length(lower)){
    png(filename = paste0("trunc", sprintf("%04d", i), ".png"))
    means <- search.means(c(0,1), .3, cutoff = .00001, a = lower[i], b = Inf, sigma = .2)
    sd.t <- search.var(c(0, 1), .2, .00001, a = 0, b = 100000000000000, mu = .3)
    y.t.norm <- dtnorm(x, means, sd.t, lower = lower[i])
    y.t.m.norm <- dnorm(x, .3, .02)/5
    plot(x
         , y.t.m.norm
         , type = 'l'
         , lty = 2
         , frame.plot = FALSE
         , ylab = "Density"
         , xlab = "Effect")
    lines(x
          , y.t.norm)
    abline(v = 0
           , col = "gray80")
    lines(x = c(means, means), y = c(0, max(y.t.norm)), lwd = 1.5)
    lines(x = c(.3, .3), y = c(0, max(y.t.m.norm)), lwd = 1.5)
    legend("topleft", legend = paste("Truncation at", round(lower[i], 2)), bty = "n")
    legend(x = .25
           , y = 3.8
           , legend = expression(x5)
           , bty = "n")
    dev.off()
  }
}

#qinvgamma=function(x,a,b,lower.tail=T) qgamma(1/x,a,rate=b,lower.tail=!lower.tail)
r.t.invgamma <- function(n, shape, scale, lower, upper){
  low <- pinvgamma(lower, shape, scale = scale, lower.tail = FALSE)
  up <- pinvgamma(upper, shape, scale = scale, lower.tail = FALSE)
  qs <- runif(n, low, up)
  qinvgamma(qs, shape, scale = scale)
}
