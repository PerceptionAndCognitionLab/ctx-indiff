plot_delta_i <- function(data
                         , sub_var
                         , cond_var
                         , rt_var
                         , trial_var
                         , congruent
                         , main = ""){
  sub_num <- unique(data[, sub_var])
  delta <- data.frame(sub = sub_num
                      , delta = rep(NA, length(sub_num))
                      , lower = rep(NA, length(sub_num))
                      , upper = rep(NA, length(sub_num))
                      , lower.within = rep(NA, length(sub_num))
                      , upper.within = rep(NA, length(sub_num))
                      )
  if(!is.factor(data[, cond_var])){
    data[, cond_var] <- ifelse(data[, cond_var] == congruent, 0, 1)
    data[, cond_var] <- factor(data[, cond_var], labels = c("c","i"))
  }
  for(i in sub_num){
    data_i <- subset(data, data[,sub_var] == i)
    rts <- split(data_i, data_i[,cond_var])
    t <- t.test(rts$i[, rt_var], rts$c[, rt_var])
    K <- tapply(data_i[, rt_var], data_i[, cond_var], length)
    ci <- (sqrt(sum(tapply(data_i[, rt_var], data_i[, cond_var], var)/K))*qt(.975,min(K)-1))
    sd_within <- sqrt(sum(tapply(data_i[, rt_var], data_i[, cond_var], var)/tapply(data_i[, rt_var], data_i[, cond_var], length)))
    diff <- t$estimat[1]-t$estimat[2]
    delta[delta$sub == i, 2:6] <- c(diff
                       , t$conf.int[1:2]
                       , diff - sd_within
                       , diff + sd_within)
  }
  delta_i <- delta[order(delta$delta),]
  
  #PLOT DELTA WITH CONFIDENCE INTERVALS
  plot(x = 1: nrow(delta_i)
       , y = delta_i$delta
       , ylim = c(-.15, .25)
       , pch = 19
       , col = "gray40"
       , ylab = expression(Observed ~ effect ~ d[i] ~ (ms))
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
       , at = seq(-.1, .2, .1)
       , labels = seq(-100, 200, 100)
  )
#   axis(side = 1
#        , at = c(seq(1, nrow(delta_i), 2), nrow(delta_i))
#        , tck = -.025
#        , labels = FALSE)
  
#   plotCI(x = 1: nrow(delta_i)
#          , y = delta_i$delta
#          , ui= delta_i$upper
#          , li=delta_i$lower
#          , add = TRUE
#          , sfrac = 1/150
#          , pch = 19
#          , col = "gray30")
  
  polygon(c(rev(1: nrow(delta_i)), 1: nrow(delta_i))
          , c(rev(delta_i$upper), delta_i$lower), col = 'gray85', border = NA)
  
  abline(h = mean(delta_i$delta)
         , col = "gray50"
         , lty = 2)
  
  points(x = 1: nrow(delta_i)
         , y = delta_i$delta
         , col = "gray40"
         , pch = 19)
  neg <- delta_i$delta < 0
  points(x = 1: sum(neg)
         , y = delta_i$delta[neg]
         , pch = 19
         , col = "darkorchid4")
  
  abline(h = 0, col = "gray40")
  
  ######PLOT DELTA IN HISTOGRAM###########
#   hist(delta_i$delta
#        , xlab = "Delta"
#        , main = "Histogram of delta")
}

plot_delta_q <- function(data
                         , sub_var
                         , cond_var
                         , rt_var
                         , trial_var
                         , congruent){
  sub_num <- unique(data[, sub_var])
  q <- qnorm(sub_num/(length(sub_num)+1))
  delta <- data.frame(sub = sub_num
                      , delta = rep(NA, length(sub_num))
                      , lower = rep(NA, length(sub_num))
                      , upper = rep(NA, length(sub_num))
                      , lower.within = rep(NA, length(sub_num))
                      , upper.within = rep(NA, length(sub_num))
                      , var.within = rep(NA, length(sub_num))
  )
  if(!is.factor(data[, cond_var])){
    data[, cond_var] <- ifelse(data[, cond_var] == congruent, 0, 1)
    data[, cond_var] <- factor(data[, cond_var], labels = c("c","i"))
  }
  for(i in sub_num){
    data_i <- subset(data, data[,sub_var] == i)
    rts <- split(data_i, data_i$cond)
    t <- t.test(rts$i[, rt_var], rts$c[, rt_var])
    var_within <- sum(tapply(data_i[, rt_var], data_i[, cond_var], var)/tapply(data_i[, rt_var], data_i[, cond_var], length))
    sd_within <- sqrt(var_within)
    diff <- t$estimat[1]-t$estimat[2]
    delta[i, 2:7] <- c(diff
                       , t$conf.int[1:2]
                       , diff - sd_within
                       , diff + sd_within
                       , var_within)
  }
  delta_i <- delta[order(delta$delta),]
  
  #PLOT DELTA WITH CONFIDENCE INTERVALS
  plot(x = q
       , y = delta_i$delta
       , ylim = c(-.2, .3)
       , pch = 19
       , col = "gray40"
       , ylab = "Mean difference"
       , xlab = "Participants"
       , frame.plot = FALSE
       , xaxt = 'n'
  )
  
  axis(side = 1
       , at = c(q[1], q[length(q)])
       , labels = c("1", paste(nrow(delta_i))))
  axis(side = 1
       , at = c(q[seq(1, nrow(delta_i), 2)], nrow(delta_i))
       , tck = -.025
       , labels = FALSE)
  
#   plotCI(x = q
#          , y = delta_i$delta
#          , ui= delta_i$upper
#          , li=delta_i$lower
#          , add = TRUE
#          , sfrac = 1/150
#          , pch = 19
#          , col = "gray30")
  
  polygon(c(rev(q), q)
          , c(rev(delta_i$upper), delta_i$lower), col = 'gray90', border = NA)
  
  points(x = q
         , y = delta_i$delta
         , col = "gray40"
         , pch = 19)
  
  abline(h = 0, col = "gray40")
  
  sd.bar <- sqrt(mean(delta$var.within))
  abline(a = mean(delta$delta), b = sd.bar, lwd = 2, col = "red")
  abline(a = 0, b = sd.bar, lwd = 2, col = "gray")
  
  ######PLOT DELTA IN HISTOGRAM###########
  #   hist(delta_i$delta
  #        , xlab = "Delta"
  #        , main = "Histogram of delta")
}
