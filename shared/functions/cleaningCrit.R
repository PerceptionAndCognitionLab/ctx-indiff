cleaning <- function(filename, clean.dat, cond.exclude = NA, exp.include = NA, blktype.include = NA, rts = c(.2, 2), acc.exclude = 0, colnms = NA, trials = 1){
  if(!is.na(colnms[1])){
    dat <- read.table(filename)
    colnames(dat) <- colnms} else{
      dat <- read.csv2(filename, header=TRUE, dec=".")
    }
  if(is.factor(dat$congruency)) dat$cond <- as.numeric(dat$congruency)
  if(is.vector(dat$RT)) dat$rt <- dat$RT / 1000
  if(is.numeric(cond.exclude)) dat <- subset(dat, cond != cond.exclude)
  if(is.numeric(exp.include)) dat <- subset(dat, exp == exp.include)
  if(is.numeric(blktype.include)) dat <- subset(dat, blktype == blktype.include)
  
  dat.time <- dat[dat$rt > rts[1] & dat$rt < rts[2], ]
  percent.time <- round((1 - nrow(dat.time) / nrow(dat)) * 100, 2)
  
  if(is.vector(dat$RT)){
    dat.acc <- subset(dat, accuracy != acc.exclude )} else{
      dat.acc <- subset(dat, acc != acc.exclude )
    }
  percent.acc <- round((1 - nrow(dat.acc) / nrow(dat)) * 100, 2)
  
  if(is.vector(dat$trial)){
    dat.trial <- dat[dat$trial >= trials, ]
    percent.trial <- round((1 - nrow(dat.trial) / nrow(dat))*100, 2)}else{
      percent.trial = 0
    }
  
  percent.total <- round((1 - nrow(clean.dat) / nrow(dat)) * 100, 2)
return(list(N = nrow(dat), N.clean = nrow(clean.dat), p.time = percent.time, p.acc = percent.acc, p.trial = percent.trial, p.total = percent.total))
}