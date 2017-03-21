
########################
#####DELTA PLOTS########
########################

#plot Stroop
layout(mat = matrix(1:2, nrow = 1, byrow = TRUE))
par(cex=1, las=1, mar=c(4, 3, 2, 2), mgp=c(2,.6,0)) 
prob <- seq(.1,.9,.1)
quants.1 <- tapply(stroop$rt, list(stroop$sub, stroop$cond), probs=prob, quantile)
quants.2 <- tapply(dat.stroop.p1$rt, list(dat.stroop.p1$sub, dat.stroop.p1$cond), probs=prob, quantile)
quants.3 <- tapply(dat.stroop.p2$rt, list(dat.stroop.p2$sub, dat.stroop.p2$cond), probs=prob, quantile)
quants.1 <- array(unlist(quants.1), dim=c(length(prob), 121, 2))
quants.2 <- array(unlist(quants.2), dim=c(length(prob), 38, 2))
quants.3 <- array(unlist(quants.3), dim=c(length(prob), 38, 2))

d.effect.1 <- apply(quants.1[,,2] - quants.1[,,1], 1, mean)
d.effect.2 <- apply(quants.2[,,1] - quants.2[,,2], 1, mean)
d.effect.3 <- apply(quants.3[,,1] - quants.3[,,2], 1, mean)
d.effect <- cbind(d.effect.1, d.effect.2, d.effect.3)
average.1 <- apply((quants.1[,,1] + quants.1[,,2])/2, 1, mean)
average.2 <- apply((quants.2[,,1] + quants.2[,,2])/2, 1, mean)
average.3 <- apply((quants.3[,,1] + quants.3[,,2])/2, 1, mean)
average <- cbind(average.1, average.2, average.3)

matplot(average
        , d.effect
        , t='p'
        , col=c("black")
        , pch = 20
        , axes=F
        , ylab="Effect (ms)"
        , xlab="Average RT (ms)"
        , main = "A.")
axis(1, at=seq(.4,1,.2), labels=seq(400,1000,200))
axis(2, at=seq(0, .15,.05), labels= seq(0,150,50))
lines(average[,1], d.effect[,1], lty= 1, lwd = 1.5)
lines(average[,2], d.effect[,2], lty= 2, lwd = 1.5)
lines(average[,3], d.effect[,3], lty= 3, lwd = 1.5)
abline(h = 0, col = "gray70")
legend('topleft'
       , legend = c("Data Set 1", "Data Set 2", "Data Set 3")
       , lty = c(1, 2, 3)
       , bty = "n"
       , lwd = 1.5
       , cex = .9)




#plot Simon
quants.1 <- tapply(simon$rt, list(simon$sub, simon$cond), probs=prob, quantile)
quants.2 <- tapply(dat.simon.p1$rt, list(dat.simon.p1$sub, dat.simon.p1$cond), probs=prob, quantile)
quants.3 <- tapply(dat.simon.p2$rt, list(dat.simon.p2$sub, dat.simon.p2$cond), probs=prob, quantile)
quants.1 <- array(unlist(quants.1), dim=c(length(prob), 121, 2))
quants.2 <- array(unlist(quants.2), dim=c(length(prob), 38, 2))
quants.3 <- array(unlist(quants.3), dim=c(length(prob), 38, 2))

d.effect.1 <- apply(quants.1[,,2] - quants.1[,,1], 1, mean)
d.effect.2 <- apply(quants.2[,,1] - quants.2[,,2], 1, mean)
d.effect.3 <- apply(quants.3[,,1] - quants.3[,,2], 1, mean)
d.effect <- cbind(d.effect.1, d.effect.2, d.effect.3)
average.1 <- apply((quants.1[,,1] + quants.1[,,2])/2, 1, mean)
average.2 <- apply((quants.2[,,1] + quants.2[,,2])/2, 1, mean)
average.3 <- apply((quants.3[,,1] + quants.3[,,2])/2, 1, mean)
average <- cbind(average.1, average.2, average.3)

matplot(average
        , d.effect
        , t='p'
        , col=c("black")
        , pch = 20
        , axes=F
        , ylab="Effect (ms)"
        , xlab="Average RT (ms)"
        , main = "B."
        , ylim = c(-.02, .1))
axis(1, at=seq(.4,1,.2), labels=seq(400,1000,200))
axis(2, at=seq(0, .1,.02), labels= seq(0,100,20))
lines(average[,1],d.effect[,1], lty= 1, lwd = 1.5)
lines(average[,2], d.effect[,2], lty= 2, lwd = 1.5)
lines(average[,3], d.effect[,3], lty= 3, lwd = 1.5)
abline(h = 0, col = "gray70")
legend('topright'
       , legend = c("Data Set 4", "Data Set 5", "Data Set 6")
       , lty = c(1, 2, 3)
       , bty = "n"
       , lwd = 1.5
       , cex = .9)


