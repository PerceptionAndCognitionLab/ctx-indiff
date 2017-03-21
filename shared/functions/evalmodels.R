evalchain <- function(alpha
                      , theta
                      , s2
                      , M = M
                      , keep = 51:M #keep <- 51:M
                      , means
                      , comp){
  
  h <- apply(alpha[keep,], 2, mean)
  plot(means[, 1]
       , h
       , main = "Estimate of mean of congruent trials"
       , xlab = "Observed mean values per person"
       , ylab = "Estimate of alpha"
       , pch = 20
       , col = "azure4") #shrinkage?
  
  if(is.matrix(theta)){
    th <- apply(theta[keep,], 2, mean)
    plot(means[,2] - means[,1]
         , th
         , main = "Estimate of difference scores"
         , xlab = "Observed mean differences per person"
         , ylab = "Estimate of theta"
         , pch = 20
         , col = "azure4") #shrinkage?
  }else{
    th <- mean(theta)
    plot(theta
         , main = "Estimate of difference scores"
         , xlab = "Iteration"
         , ylab = "Estimate of theta"
         , pch = 20
         , col = "azure4" 
         )
    abline(h = th
           , col = "red")
  }
  
  apa_table(matrix(c("Estimates", round(c(mean(h), mean(h) + mean(th), mean(th)), 3), "Observed", round(comp, 3)), ncol = 4, byrow = TRUE)
            , caption = "Observed and estimated condition means and difference scores"
            , align = c("l", "r", "r", "r"))
}

##################################################################################
#######################ESTIMATE PLOT############################################
#################################################################################
plot_freq_bayes <- function(data
                            , sub_var = "sub"
                            , cond_var = "cond"
                            , congruent = 0
                            , rt_var = "rt"
                            , trial_var = "trial"
                            , thet.n = theta.n
                            , thet.1 = theta.one
                            , main = ""
){
  
  #####################FREQUENTIST PLOT#############################################
  sub_num <- unique(data[, sub_var])
  delta <- data.frame(sub = sub_num, delta = NA, lower = NA, upper = NA)
  if(!is.factor(data[, cond_var])){
    data[, cond_var] <- ifelse(data[, cond_var] == congruent, 0, 1)
    data[, cond_var] <- factor(data[, cond_var], labels = c("c","i"))
  }
  for(i in sub_num){
    data_i <- subset(data, data[,sub_var] == i)
    rts <- split(data_i, data_i[,cond_var])
    t <- t.test(rts$i[, rt_var], rts$c[, rt_var])
    delta[delta$sub == i, 2:4] <- c((t$estimat[1]-t$estimat[2]), t$conf.int[1:2])
  }
  delta$th1 <- mean(thet.1) * 1000
  delta$thn <- colMeans(thet.n) * 1000
  delta[, ncol(delta)+(1:2)] <- t(apply(thet.n, 2, function(x) quantile(x, c(.025, .975))))
  delta_i <- delta[order(delta$delta),]
  delta_i$delta <- delta_i$delta *1000
  
  plot(x = 1: nrow(delta_i)
       , y = delta_i$delta
       , ylim = c(min(delta_i$delta), max(delta_i$delta))
       , type = "l"
       , col = "gray50"
       , ylab = expression(Effect ~ parameter ~ theta[i])
       , xlab = "Participants"
       , frame.plot = FALSE
       , xaxt = 'n'
       , yaxt = 'n'
  )
  title(main, line = -1)
  axis(side = 1
       , at = c(1, nrow(delta_i))
  )
  axis(side = 2
  )
  polyCI(upper = delta_i[,ncol(delta)] * 1000, lower = delta_i[,ncol(delta)-1] * 1000, col = "gray85")
#   polygon(c(rev(1: nrow(delta_i)), 1: nrow(delta_i))
#           , c(rev(delta_i[,ncol(delta)]), delta_i[,ncol(delta)-1]), col = 'gray60', border = NA)
  
  abline(h = 0, col = "gray40")
  
  lines(delta_i$delta
         , col = "gray65"
        , lwd = 2)
#   points(delta_i$delta
#         , col = "gray30"
#         , pch = 20)
  abline(h = delta_i$th1
         , lty = 2
         , col = "indianred"
         , lwd = 1.2
  )
  
  points(delta_i$thn
         , col = "darkslateblue"
         , pch = 19
         )
 

}
  #####################BAYESIAN PLOT#############################################
  


##################################################################################
#######################PLOT BAYES FACTORS############################################
#################################################################################
plot.bf <- function(SD.trunc, SD.norm, SD.one){
  log.BFtrunc <- log10(SD.trunc * 1/SD.one)
  log.BFnorm <- log10(SD.norm * 1/SD.one)
  log.Null <- log10(1/SD.one)
  
  bfs <- c(log.Null, log.BFnorm, log.BFtrunc, 0)
  bfs <- ifelse(bfs == Inf, 1000, bfs)
  bfs <- ifelse(bfs == -Inf, -1000, bfs)
  plot(x = 1
       , col = "white"
       , ylim = c(min(bfs)+min(bfs)/10, max(bfs))
       , xlim = c(0, 5)
       , ylab = expression(log[10] ~ BF)
       , axes = FALSE
       , xlab = ""
       , frame.plot = FALSE
  )
  axis(side = 2
       , at = seq(round(max(bfs), -1)
                  , min(bfs)+min(bfs)/10
                  , -(round(max(bfs), -1) - round(min(bfs), -1))/4)
       , labels = TRUE)
  Arrows(x0 = 1:3
         , y0 = rep(0, 3)
         , x1 = 1:3
         , y1 = bfs[1:3]
         , arr.type = "triangle"
         , lwd = 2
         , col = "gray40"
         , arr.adj = 1
         , arr.length = .15
         , arr.width = .15)
  points(4, 0, pch = 15, col = "gray40")
  axis(side = 1
       , at = 1:4
       , labels = FALSE)
  text(x = 1:4
       , par("usr")[3]+min(bfs)/15
       , labels = c("Null", "Normal", "Truncated", "One effect")
       , srt = 45
       , xpd = TRUE
       , adj = 1)
}


##########BAYES FACTORS#####################
return.bf <- function(SD.trunc, SD.norm, SD.one){
  BFtrunc <- SD.trunc * 1/SD.one
  BFnorm <- SD.norm * 1/SD.one
  Null <- 1/SD.one
  
  bfs <- c(Null, BFnorm, BFtrunc, 1)
  return(bfs)
}


##########################plotestimates for internal use###############################
plot_freq_bayes_intern <- function(data
                            , sub_var = "sub"
                            , cond_var = "cond"
                            , congruent = 0
                            , rt_var = "rt"
                            , trial_var = "trial"
                            , thet.n = theta.n
                            , thet.t = theta.t
                            , thet.1 = theta.one
                            , keep = keep.min:M
                            , experiment
                            , bfs = simon.com$bfs[[2, 1:3]]
                            , f = f.stat
                            , trunccount
                            , poscount
                            , eval
){
  
  #####################FREQUENTIST PLOT#############################################
  sub_num <- unique(data[, sub_var])
  delta <- data.frame(sub = sub_num, delta = NA, lower = NA, upper = NA)
  if(!is.factor(data[, cond_var])){
    data[, cond_var] <- ifelse(data[, cond_var] == congruent, 0, 1)
    data[, cond_var] <- factor(data[, cond_var], labels = c("c","i"))
  }
  for(i in sub_num){
    data_i <- subset(data, data[,sub_var] == i)
    rts <- split(data_i, data_i$cond)
    t <- t.test(rts$i[, rt_var], rts$c[, rt_var])
    delta[i, 2:4] <- c((t$estimat[1]-t$estimat[2]), t$conf.int[1:2])
  }
  delta$th1 <- mean(thet.1)
  delta$tht <- apply(thet.t[keep,], 2, mean)
  delta$thn <- apply(thet.n[keep,], 2, mean)
  delta[, c(8:9)] <- t(apply(thet.n[keep, ], 2, function(x) quantile(x, c(.025, .975))))
  delta_i <- delta[order(delta$delta),]
  
  plot(x = 1: nrow(delta_i)
       , y = delta_i$delta
       , ylim = c(-.05, .18)
       , pch = 20
       , col = "gray60"
       , ylab = expression(Effect ~ parameter ~ theta[i])
       , xlab = "Participants"
       , frame.plot = FALSE
       , xaxt = 'n'
       , main = experiment
  )
  axis(side = 1
       , at = c(1, nrow(delta_i))
  )
  axis(side = 1
       , at = c(seq(1, nrow(delta_i), 2), nrow(delta_i))
       , tck = -.025
       , labels = FALSE)
  
  polygon(c(rev(1: nrow(delta_i)), 1: nrow(delta_i))
          , c(rev(delta_i[,9]), delta_i[,8]), col = 'grey90', border = NA)
  
  abline(h = 0, col = "gray40")

  points(delta_i$delta
         , col = "gray60"
         , pch = 20)
  points(delta_i$thn
         # , col = adjustcolor("#DD4814", alpha.f = .5)
         , col = "#DD4814"
         , pch = 15)
  #   lines(delta_i$tht
  #         , col = "steelblue")
  points(delta_i$tht
         # , col = adjustcolor("royalblue4", alpha.f = .5)
         , col = "royalblue4"
         , pch = 19)
  abline(h = delta_i$th1[1]
         , col = "darkgreen"
         , lty=2)
  
  lbf=round(log10(bfs/bfs[1]),3)
  lab=paste(c("One","Trunc","Gen"),lbf,sep=": ")
  legend("topleft"
         , legend = lab
         , col = c('darkgreen', 'royalblue4', '#DD4814')
         , pch = c(NA, 15, 19)
         , lty = c(2, NA, NA)
         , title = "log10 BF")
  legend("bottomright",legend=c(
    paste("F: ",round(f,3))
#     paste("Naive ES:",round(mean(effect)/sd(effect),3)),
#     paste("Model ES: ",round(mean(g.theta)/sd(g.theta),3)),
    , paste("Trunc Count:", round(trunccount[1]/M, 2), round(trunccount[2]/M, 2))
    , paste("Pos count:", round(poscount/M, 2))
    , paste("evaluated at", round(eval, 3))
))
  
  
}

