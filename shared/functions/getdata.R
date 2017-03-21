
########################Simulate######################################
mkdata <- function(I = 100, J = 2, K = 20, t.d, t.mean, t.sd){
  N <- I*J*K #Total number of observations
  sub <- rep(1:I, each = J*K)
  condition <- rep(0:(J-1), times = I, each = K)
  replicate <- rep(1:K, times = I*J)
  
  #########################
  #Set true values
  cond <- rep(0:(J-1), I)
  get.truth <- function(x){
    rnorm(n = 1
          , mean = t.mean
          , sd = t.sd
    )
  }
  t.alpha <- sapply(1:(I), get.truth)
  
  #################################
  #make data
  #Lets use a vector rather than matrix format.
  sd.d <- t.sd - .3*t.sd
  
  get.data <- function(x){
    rtruncnorm(n = K
               , mean = (rep(t.alpha, each = 2)[x] + cond[x] * t.d)
               , sd = sd.d
               , a = 0
    )
  }
  y <- sapply(1:(I*J), get.data)
  y <- as.vector(y)
  
  dat <- cbind(sub = sub, cond = condition, trial = replicate, y = y)
  dat <- as.data.frame(dat)
  
  means <- tapply(dat$y, list(dat$sub, dat$cond), mean)
  comp <- apply(means, 2, mean)
  comp[3] <- comp[2] - comp[1]
  
  return(list(dat = dat, means = means, compares = comp))
}



########################readdata######################################
readdata <- function(filename
                     , clnames #colnames of dataset
                     , v.y = "y" #Y/rt
                     , v.cond = "cond" #condition variable
                     , cond.congruent = 0 #label of congruent cond
                     , v.sub = "sub" #Subject variable
                     , unit = "sec" #Which unit is  it in? Standard sec
                     , trial = FALSE #Remove fist trials/block -> 1:x?
                     , exclude = FALSE #other excluding conditions like accuracy, in matrix
                     , exp.pratte = TRUE
                     ){
  
  if(exp.pratte == TRUE){
    dat <- read.table(filename)
    colnames(dat) <- clnames
    
    I <- length(levels(as.factor(dat[, v.sub])))
  }else{
    dat <- filename
    I <- length(unique(dat[, v.sub]))
  }

  if(trial[1] != FALSE){
    dat <- dat[!(dat$trial %in% trial), ]}  #REMOVE FIRST 5 OF EVERY BLOCK

  if(unit == "sec"){
    a <- paste("dat$", v.y, " > .2 & dat$", v.y, " < 2", sep = "")
    dat <- dat[eval(parse(text = a)), ]
  }else{
    a <- paste0("dat$", v.y, " > 200 & dat$", v.y, " < 2000")
    dat <- dat[eval(parse(text = a)), ]
  }

  N = nrow(exclude)
  for(n in 1:N){
    a <- paste("dat$"
               , exclude[n, 1]
               , " "
               , exclude[n,2]
               , " "
               , exclude[n, 3]
               , sep = "")
    dat <- dat[eval(parse(text = a)), ]
  } 
  
  J <- length(unique(dat[, v.cond]))# Number of conditions
  N <- nrow(dat) #Total number of observations
  
  sub.o <- dat[, v.sub]
  sub <- c(1)
  for(i in 2:N){
    if(sub.o[i] == sub.o[i-1]) sub[i] <- sub[i-1]
    else sub[i] <- sub[i-1] + 1
  }
  
  
  condition <- ifelse(dat[, v.cond] == cond.congruent, 0, 1)
  cond <- rep(0:(J-1), I)
  y <- dat[, v.y]

  means <- tapply(y, list(sub, condition), mean)
  comp <- tapply(y, condition, mean)
  comp[3] <- comp[2] - comp[1]
  
  return(list(I = I, J = J, N = N, comp = comp, means = means, y = y, cond = cond, condition = condition, sub = sub, dat = dat))
}


  