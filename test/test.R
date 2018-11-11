#test standardization
source("test/prep_data.R")
out <- .Fortran("standardize", PACKAGE = "MTCox",
             ntasks = as.integer(ntasks),
             ii = as.integer(ii),
             io = as.integer(io),
             p = as.integer(p),
             w = as.double(data[,wid]),
             x = as.double(X),
             n = as.integer(n),
             iex = as.integer(rep(1,n)),
             wsum = double(ntasks),
             xmean = double(ntasks*p),
             xsd = double(ntasks*p),
             flnullsd = integer(2))
t(sapply(seq(ntasks), function(k){
apply(X[ii[k]:io[k],], 2, function(col) sqrt(sum((col-mean(col))^2)/length(col)))
}))-
matrix(out$xsd, nrow = ntasks, ncol = p, byrow = T)

#test grad_null
source("test/prep_data.R")
out <- .Fortran("gradnull", PACKAGE = "MTCox",
                ntasks = as.integer(ntasks),
                p = as.integer(p),
                n = as.integer(n),
                ns = as.integer(ns),
                iski = as.integer(iski),
                isko = as.integer(isko),
                nski = as.integer(nski),
                nsko = as.integer(nsko),
                w = as.double(data[,wid]),
                x = as.double(X),
                delta = as.double(data[, did]),
                d = as.double(ds),
                grad = double(p*ntasks)
                )
matrix(out$grad,nrow = ntasks,ncol = p, byrow=T)
#test saturated llk
source("test/prep_data.R")
out <- .Fortran("llk_saturated", PACKAGE = "MTCox",
                ntasks = as.integer(ntasks),
                ns = as.integer(ns),
                nski = as.integer(nski),
                nsko = as.integer(nsko),
                d = as.double(ds),
                llk_sat = double(ntasks)
)
out$llk_sat
sapply(seq(K), function(k){
   -sum(ds[nski[k]:nsko[k]] * log(ds[nski[k]:nsko[k]]))
})
#test lambda_max
out <- .Fortran("lambda_max", PACKAGE = "MTCox",
                p = as.integer(p),
                ntasks = as.integer(ntasks),
                reg = as.integer(0),
                grad = as.double(grad),
                pf = as.double(rep(1,p)),
                al = double(1),
                iex = as.integer(rep(1,p))
)
out$al
out$grad
#test loglik
source("test/prep_data.R")
eta <- rep(0,n)

t0 <- proc.time()
out <- .Fortran("loglik", PACKAGE = "MTCox",
                ntasks = as.integer(ntasks),
                n = as.integer(n),
                p = as.integer(p),
                ns = as.integer(ns),
                iski = as.integer(iski),
                isko = as.integer(isko),
                nski = as.integer(nski),
                nsko = as.integer(nsko),
                w = as.double(data[,wid]),
                x = as.double(X),
                delta = as.integer(data[,did]),
                d = as.double(ds),
                eta = as.double(eta),
                llk = double(ntasks)
)
proc.time()-t0
out$llk

llk <- rep(0,ntasks)
t0 <- proc.time()
for(k in seq(ntasks)){
   SE =0
   for(isk in seq(nsko[k], nski[k])){
      SI = 0
      for(i in seq(iski[isk],isko[isk])){
         #eta = sum(X[i,] * b[,k])
         SE = SE + data[i,wid]*exp(eta[i])
         if(data[i,did] == 0)SI = SI +  data[i,wid] * eta[i]
      }
      llk[k] = llk[k] + SI - ds[isk] * SE
   }
}
proc.time()-t0
llk

#test proximal
source("test/prep_data.R")
b <- seq(ntasks)^2*0.1
g <- -seq(ntasks)
out <- .Fortran("prox", PACKAGE = "MTCox",
                lam = as.double(0.1),
                pf = as.double(1),
                sig = as.double(0.1),
                alpha = as.double(0.0),
                reg = as.integer(2),
                ntasks = as.integer(ntasks),
                beta = as.double(b),
                grad = as.double(g)
)
b
out$beta


#test grad_hess
# source("test/prep_data.R")
# b <- matrix(seq(p*ntasks)*0, nrow = p, ncol = ntasks)
# out <- .Fortran("grad_hess", PACKAGE = "MTCox",
#                 ntasks = as.integer(ntasks),
#                 p = as.integer(p),
#                 n = as.integer(n),
#                 ns = as.integer(ns),
#                 iski = as.integer(iski),
#                 isko = as.integer(isko),
#                 nski = as.integer(nski),
#                 nsko = as.integer(nsko),
#                 w = as.double(data[,wid]),
#                 x = as.double(X),
#                 delta = as.double(data[, did]),
#                 d = as.double(ds),
#                 beta = as.double(b),
#                 grad = double(p*ntasks),
#                 hess = double(p*ntasks)
# )
# out$grad
# out$hess


#test grad_hessj
source("test/prep_data.R")
eta <- rep(0,n)
j <- 1
out <- .Fortran("grad_hessj", PACKAGE = "MTCox",
                ntasks = as.integer(ntasks),
                n = as.integer(n),
                ns = as.integer(ns),
                iski = as.integer(iski),
                isko = as.integer(isko),
                nski = as.integer(nski),
                nsko = as.integer(nsko),
                w = as.double(data[,wid]),
                x = as.double(X[,j]),
                delta = as.double(data[, did]),
                d = as.double(ds),
                eta = as.double(eta),
                grad = double(ntasks),
                hess = double(ntasks)
)
out$grad
out$hess
#
# mu <- exp(eta)
# grad <- rep(0,ntasks)
# hess <- rep(0,ntasks)
# for(k in seq(ntasks)){
#    SE=0;SR=0;SR2=0;SI=0;hhat=0
#    for(isk in seq(nsko[k],nski[k])){
#       ind <- iski[isk]:isko[isk]
#       SE <- SE + sum(data[ind,wid]*mu[ind])
#       hhat <- 1/SE
#       SR <- SR + sum(data[ind,wid]*mu[ind]*X[ind,j])
#       SR2 <- SR2 + sum(data[ind,wid]*mu[ind]*(X[ind,j])^2)
#       SI <- sum(X[ind,j]*data[ind,wid]*(1-data[ind,did]))
#       grad[k] <- grad[k] - SI + ds[isk]*hhat*SR
#       hess[k] <- hess[k] - (ds[isk] * hhat * SR)^2 + ds[isk] * hhat * SR2
#    }
# }
# grad
# hess

source("test/prep_data.R")
plot(NA, xlim = c(0,1)*1, ylim = c(0,1)*1)
abline(h=0, v=0, lty=3)

for(bx in seq(0,1,0.1)){
   for(by in seq(0,1,0.1)){
      beta <- matrix(rep(c(bx,by,rep(0,p-2)),ntasks),p,ntasks)
      eta <- predict.MTcox(X,data$task,beta,T)
      gd <- matrix(0,p,ntasks)
      for(j in seq(p)){
         out <- .Fortran("grad_hessj", PACKAGE = "MTCox",
                         ntasks = as.integer(ntasks),
                         n = as.integer(n),
                         ns = as.integer(ns),
                         iski = as.integer(iski),
                         isko = as.integer(isko),
                         nski = as.integer(nski),
                         nsko = as.integer(nsko),
                         w = as.double(data[,wid]),
                         x = as.double(X[,j]),
                         delta = as.double(data[, did]),
                         d = as.double(ds),
                         #beta = as.double(beta[j,]),
                         eta=as.double(eta),
                         grad = double(ntasks),
                         hess = double(ntasks)
         )
         gd[j,] <- out$grad
      }
      segments(beta[1,], beta[2,],beta[1,] + gd[1,], beta[2,] + gd[2,], col=1:ntasks)
   }
}

#test gpg_majorized
#source("test/prep_data.R")
#plot(NA, xlim = c(0,1)*1, ylim = c(0,1)*1)
beta <- matrix(rep(0,ntasks*p), nrow = p, ncol = ntasks)
eta <- rep(0,n)
points(beta[1,], beta[2,], col=1:ntasks, pch=16)
#abline(h=0, v=0, lty=3)

for(i in seq(20)){
   beta0<-beta
   out <- .Fortran("gpg_cycle", PACKAGE = "MTCox",
                   ntasks = as.integer(ntasks),
                   p = as.integer(p),
                   n = as.integer(n),
                   ii = as.integer(ii),
                   io = as.integer(io),
                   ns = as.integer(ns),
                   iski = as.integer(iski),
                   isko = as.integer(isko),
                   nski = as.integer(nski),
                   nsko = as.integer(nsko),
                   w = as.double(data[,wid]),
                   x = as.double(X),
                   delta = as.double(data[, did]),
                   d = as.double(ds),
                   beta = as.double(beta),
                   eta = as.double(eta),
                   lam = as.double(0.1),
                   pf = as.double(rep(1,p)),
                   reg = as.integer(0),
                   alpha = as.double(0),
                   grad = double(p*ntasks),
                   sig = double(p),
                   ierr= double(5)
   )
   beta <- matrix(out$beta, nrow = p, byrow=F)
   eta <- out$eta
   grad <- matrix(out$grad, nrow = p, byrow=F)
   points(beta[1,], beta[2,], col=1:ntasks, pch=16)
   #segments(beta0[1,], beta0[2,],beta0[1,] + grad[1,], beta0[2,] + grad[2,], col=1:ntasks)
}
beta
out$grad
#test gpg__cycle_backtracking
#source("test/prep_data.R")

#plot(NA, xlim = c(0,1)*1, ylim = c(0,1)*1)
beta <- matrix(rep(0,ntasks*p), nrow = p, ncol = ntasks)
eta <- rep(0,n)
points(beta[1,], beta[2,], col=1:ntasks, pch=16)
#abline(h=0, v=0, lty=3)
llk<-llk0
for(i in seq(10)){
   sig = rep(0,p)
   if(i>1)sig=out$sig
   beta0<-beta
   out <- .Fortran("gpg_cycle_backtracking", PACKAGE = "MTCox",
                   ntasks = as.integer(ntasks),
                   p = as.integer(p),
                   n = as.integer(n),
                   ii = as.integer(ii),
                   io = as.integer(io),
                   ns = as.integer(ns),
                   iski = as.integer(iski),
                   isko = as.integer(isko),
                   nski = as.integer(nski),
                   nsko = as.integer(nsko),
                   w = as.double(data[,wid]),
                   x = as.double(X),
                   delta = as.double(data[, did]),
                   d = as.double(ds),
                   beta = as.double(beta),
                   eta = as.double(eta),
                   lam = as.double(0.1),
                   pf = as.double(rep(1,p)),
                   reg = as.integer(0),
                   alpha = as.double(0),
                   llka = as.double(llk),
                   frac = as.double(0.5),
                   grad = double(p*ntasks),
                   hess = double(p*ntasks),
                   sig = as.double(sig),
                   ierr = integer(5),
                   ntries = integer(p)
   )
   beta <- matrix(out$beta, nrow = p, byrow=F)
   eta <- out$eta
   grad <- matrix(out$grad, nrow = p, byrow=F)
   points(beta[1,], beta[2,], col=1:ntasks, pch=17)
   text(beta[1,], beta[2,], labels = out$ntries, pos = 3, col=  out$ntries>1)
   llk = out$llka
   #segments(beta0[1,], beta0[2,],beta0[1,] + grad[1,], beta0[2,] + grad[2,], col=1:ntasks)
}
beta
out$ntries
out$sig
out$grad


#test gpg_descent
source("test/prep_data.R")
beta <- matrix(rep(0,ntasks*p), nrow = p, ncol = ntasks)
llk0
t0 <- proc.time()
out <- .Fortran("gpg_descent", PACKAGE = "MTCox",
                     ntasks = as.integer(ntasks),
                     p = as.integer(p),
                     n = as.integer(n),
                     ii = as.integer(ii),
                     io = as.integer(io),
                     iex = as.integer(rep(1,p)),
                     ns = as.integer(ns),
                     iski = as.integer(iski),
                     isko = as.integer(isko),
                     nski = as.integer(nski),
                     nsko = as.integer(nsko),
                     w = as.double(data[,wid]),
                     x = as.double(X),
                     delta = as.double(data[, did]),
                     d = as.double(ds),
                     beta = as.double(beta),
                     eta = as.double(eta),
                     lam = as.double(0.05),
                     pf = as.double(rep(1,p)),
                     reg = as.integer(0),
                     alpha = as.double(0),
                     llkb = as.double(llk0),
                     ierr = integer(5),
                     alg = as.integer(1),
                     eps = as.double(1e-5),
                     frac = as.double(0.1),
                     ncycles = integer(1),
                     nupdates = integer(1),
                     iactive = as.integer(1)
   )
beta <- matrix(out$beta, nrow = p, byrow=F)
beta
out$ncycles
out$nupdates
out$llkb
proc.time()-t0


######################################
## SOLUTION PATH TEST
######################################
source("test/prep_data.R")
nlam <- 100
t0 <- proc.time()
out <- .Fortran("mtcox_solutionpath", PACKAGE = "MTCox",
                ntasks = as.integer(ntasks), ii = as.integer(ii), io = as.integer(io),
                p = as.integer(p), ns = as.integer(ns), delta = as.double(data[, did]),
                w = as.double(data[,wid]), x = as.double(X), d = as.double(ds),
                iski = as.integer(iski), isko = as.integer(isko),
                nski = as.integer(nski), nsko = as.integer(nsko),
                n = as.integer(n),
                iex = as.integer(rep(1,p)),
                str = as.integer(1),
                dfmax = as.integer(p),
                pmax = as.integer(p),
                nlam = as.integer(nlam),
                earlystop = as.integer(1),
                alg = as.integer(1),
                lamfrac = as.double(1e-3),
                lam = double(nlam),
                maxit = as.integer(1e4),
                eps = as.double(1e-4),
                frac = as.double(0.1),
                alpha = as.double(0),
                reg = as.integer(0),
                pf = as.double(rep(1,p)),
                hhat = double(sum(ns)*nlam),
                llk_path = double(ntasks*nlam),
                dev = double(nlam),
                pdev_path = double(nlam),
                eta_path = double(n*nlam),
                null_dev = double(1),
                ierr = integer(5),
                beta = double(p* ntasks* nlam),
                ncycles = integer(1),
                nupdates = integer(1),
                llk_null = double(ntasks),
                llk_sat = double(ntasks),
                grad = double(p*ntasks),
                alf=double(1)
                )
proc.time()-t0
out$ncycles
out$nupdates
out$ierr
out$lam
matrix(out$llk_path,ntasks,nlam)
out$dev
round(out$pdev_path*100,1)
array(out$beta, c(p, ntasks, nlam))

###############################################################
# test for MTcox.R

########################
# non-simultaneous tests
#parameters
p <- 5
K <- 3
n <- 100
#raw data
X = round(matrix(rnorm(n*p,0,1),n,p),2)
beta <- c(1,-2,3,rep(0,p-3))
#eventually do weibull distribution or something else
y <- round(rexp(n,exp(X %*% beta)),1)
data <- data.frame(
   y = y,
   d = rbinom(n,1,0.3),
   task = sample.int(K,n, replace=TRUE),
   X = X
)
yid <- 1
did <- 2
taskid <- 3

fit <- MTcox(data, yid, taskid, did)

#status: need to test arguments
########################
# simultaneous tests

#parameters
p <- 3
K <- 5
n <- 100/K
#raw data
X = round(matrix(rnorm(n*p,0,1),n,p),2)
beta <- c(1,-2,3,rep(0,p-3))
#eventually do weibull distribution or something else
y <- matrix(round(rexp(n*K,exp(X %*% beta)),1),n,K)
data <- data.frame(
   y = y,
   d = rbinom(n,1,0.3),
   X = X
)
yid <- 1:K
did <- K+1
taskid <- K+2
