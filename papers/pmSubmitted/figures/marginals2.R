alpha <- seq(-1.1,1.1,.01)
theta <- seq(-1.1,1.1,.01)

# dcauchy(1)
# temp <- seq(-1000, 1000, .001)
# test <- dmvc(cbind(temp, rep(0, length(temp))), rep(0, 2) ,diag(2))
# summtest <- sum(test * .001)
# dmvc(c(1, 0), c(0, 0), diag(2))/summtest
# dcauchy(1)

layout(mat = matrix(c(1, 2, 3, 3), nrow = 2, byrow = TRUE), height = c(.9, .1), width = c(1/2, 1/2))

percentiles <- qcauchy(c(.25, .35, .45))

dist <- function(alpha1, alpha2, sigma2, cov = 0){dmvc(cbind(alpha1, alpha2), rep(0, 2), diag(2)*sigma2 + cov)}
alpha.mar <- outer(alpha, alpha, dist, sigma2 = 1)
levels.p.a <- dist(percentiles, percentiles, 1)

par(mar=c(3,3,3,1), mgp = c(2,1,0))
image(alpha
      , alpha
      , alpha.mar > 0.03456156
      , col = grey((256:200)/256)
      # , frame.plot = F
      , xlab = expression(alpha[1]*"*")
      , ylab = expression(alpha[2]*"*"))
image(alpha
      , alpha
      , alpha.mar > 0.09590448
      , col = adjustcolor(grey((256:200)/256), alpha.f = .5)
      , add = TRUE)
image(alpha
      , alpha
      , alpha.mar > 0.16687259
      , col = adjustcolor(grey((256:200)/256), alpha.f = .5)
      , add = TRUE)
contour(alpha
        , alpha
        , alpha.mar
        , levels = levels.p.a
        , drawlabels = F
        , col = 'gray40'
        #, lwd = 2
        , add = TRUE
)
title("A.", line = 2)

axis(side = 3, at = seq(-1, 1, 1), labels = seq(-300, 300, 300))

percentiles.theta <- qcauchy(c(.05, .15, .25, .35, .45), 0, .1)
theta.mar <- outer(theta, theta, dist, sigma2 = (1/5)^2, cov = .1)
levels.p.t <- dist(percentiles.theta, percentiles.theta, sigma2 = .1)

par(mar=c(3,2,3,3), mgp = c(2,1,0))
image(theta
      , theta
      , theta.mar > levels.p.t[3]
      , col = grey((256:200)/256)
      , yaxt = 'n'
      # , frame.plot = F
      , xlab = expression(theta[1]*"*")
      , ylab = expression(theta[2]*"*"))
image(theta
      , theta
      , theta.mar > levels.p.t[4]
      , col = adjustcolor(grey((256:200)/256), alpha.f = .5)
      , add = TRUE)
image(theta
      , theta
      , theta.mar > levels.p.t[5]
      , col = adjustcolor(grey((256:200)/256), alpha.f = .5)
      , add = TRUE)
contour(theta
        , theta
        , theta.mar
        , levels = levels.p.t[3:5]
        , col = 'gray40'
        #, lwd = 2
        , add = TRUE
        , drawlabels = FALSE
)
contour(theta
        , theta
        , theta.mar
        , levels = levels.p.t[1:2]
        , labels = c("10 %", "30 %")
        , col = 'gray40'
        #, lwd = 2
        , add = TRUE
        , drawlabels = TRUE
        , labcex = 0.8
)
title("B.", line = 2)

axis(side = 4, at = seq(-1, 1, 1), labels = seq(-300, 300, 300))
axis(side = 3, at = seq(-1, 1, 1), labels = seq(-300, 300, 300))
mtext(expression(theta[2]*"*"), side = 2, line = .5, cex = .85, las = 0)

par(mar=c(0.1,0.1,0.1,0.1), mgp = c(1,1,0))
plot(1, type = "n", axes=FALSE, xlab="", ylab="")
par(xpd=TRUE)
legend("center", legend = paste("upper", seq(50, 90, 20), "%"), fill = c('gray90', 'gray80', 'gray70'), cex = .9, horiz = TRUE)
par(xpd=FALSE)


# di <- function(a1, a2) a1^2 + a2^2
# out <- outer(alpha, alpha, di)
# selector <- out < 2
# 
# int <- function(a1, a2, s2){
#   if(sqrt(a1^2 + a2^2) < sqrt(2)){
#     dist(a1, a2, sigma2 = s2)
#   }
# }
# 
# p <- 0
# for(i in 1:length(alpha)){
#   for(j in 1:length(alpha)){
#     if(selector[i, j] == TRUE) p <- p + dist(alpha[i], alpha[j], sigma2 = 1)
#   }
# }
# p*.01^2
