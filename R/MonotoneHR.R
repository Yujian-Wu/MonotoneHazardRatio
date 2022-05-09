library(survival)
library(fdrtool)
library(plyr)
library(tidyverse)
library(kedd)
library(KernSmooth)
library(twostageTE)

data("chernoff_realizations")

#########
NA.est = function(surv.data){
  surv.model <- survfit(Surv(time, status) ~ 1, data=surv.data, ctype = 1)
  est = data.frame(surv.model$time, surv.model$cumhaz)
  if (all(est[1, ] != c(0, 0))) est <- rbind(c(0, 0), est)
  colnames(est) = c('time', 'cumhaz')
  return(est)
}

# gcmlcm function with duplicates removal

gcm.unique <- function(x, y){
  # x <- Lambda.T$cumhaz
  x.unique <- unique(x)
  y.unique <- rep(NA, length(x.unique))
  for (i in 1:length(x.unique)){
    y.unique[i] <- y[max(which(x == x.unique[i]))]
  }
  #### GCM
  logcm <- gcmlcm(x.unique, y.unique, type = 'gcm')
  logcm
}

# derivative and scaling parameter estimator
tau.est <- function(time.grid, f.data, g.data, ci.lvl=0.05){
  N <- dim(f.data)[1] + dim(g.data)[1]
  prop <- dim(f.data)[1] / N
  n.size <- N %/% 2
  
  Lambda.S <- NA.est(f.data)
  Lambda.T <- NA.est(g.data)
  
  Lambda.S.fn <- stepfun(Lambda.S$time, c(0, Lambda.S$cumhaz))
  Lambda.T.fn <- stepfun(Lambda.T$time, c(0, Lambda.T$cumhaz))
  
  mapped.t <- Lambda.T.fn(time.grid)
  
  #### inverse function of Lambda.T
  Lambda.T.inv <- sapply(Lambda.T$cumhaz, function(x) Lambda.T$time[min(which(Lambda.T$cumhaz >= x))])
  
  #### Lambda.S compose Lambda.T inverse
  LambdaS.LambdaT.inv <- Lambda.S.fn(Lambda.T.inv)
  logcm <- gcm.unique(Lambda.T$cumhaz, LambdaS.LambdaT.inv)
  
  #### left derivative of GCM
  left.deriv <- logcm$slope.knots
  left.deriv.time <- logcm$x.knots[2:length(logcm$y.knots)]
  
  smooth.x <- seq(0, max(Lambda.T$cumhaz), length.out = n.size^(2/3))
  #### monotone hazard ratio estimates
  smooth.y <- sapply(smooth.x, function(x) ifelse(x < min(left.deriv.time), 0, left.deriv[min(which(left.deriv.time >= x))]))
  
  # direct derivative estimation from primitive functions (with kernel smoothing)
  prim.drv.ker <- dkde(smooth.x, smooth.y, deriv.order = 1, kernel = 'gaussian')
  prim.drv.fit <- locpoly(smooth.x, smooth.y, bandwidth = prim.drv.ker$h, drv = 1)
  
  prim.drv.ker.indices <- unlist(sapply(mapped.t, function(x) which.min(abs(prim.drv.fit$x - x))))
  mapped.drv <- prim.drv.fit$y[prim.drv.ker.indices]
  
  pro.s <- sapply(time.grid, function(x) mean(f.data$time >= x))
  pro.t <- sapply(time.grid, function(x) mean(g.data$time >= x))
  
  theta.t <- sapply(mapped.t, function(x) ifelse(x < min(left.deriv.time), 0, left.deriv[min(which(left.deriv.time >= x))]))
  
  tau.hat <- (4 * mapped.drv * (theta.t / (prop*pro.s) + theta.t^2 / ((1- prop) * pro.t)))^(1/3)
  
  # confidence interval
  chernoff.cdf <- unique(c(rev(1 - chernoff_realizations$DF), chernoff_realizations$DF))
  chernoff.x <- unique(c(rev(-1 * chernoff_realizations$xcoor), chernoff_realizations$xcoor))
  
  inter.cher <- approx(chernoff.x, chernoff.cdf, xout = seq(-2, 2, 0.0001))
  
  ub <- 1 - ci.lvl / 2
  lb <- ci.lvl / 2
  x.ub <- inter.cher$x[which.min(abs(inter.cher$y - ub))]
  x.lb <- inter.cher$x[which.min(abs(inter.cher$y - lb))]
  
  ci.lower <- x.lb * tau.hat / N^(1/3) + theta.t
  ci.upper <- x.ub * tau.hat / N^(1/3) + theta.t
  
  return(list(
    tau = tau.hat,
    theta = theta.t,
    ci.upper = ci.upper,
    ci.lower = ci.lower
  ))
}

