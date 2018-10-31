library(MTcox)
# mtcox_solutionpath(ntasks,ii,io,p,ns,delta,w,x,d,iski,isko, &
#                       nski,nsko,n,iex,str,dfmax,pmax,nlam,earlystop,alg,lamfrac, &
#                       lam,maxit,eps,frac,alpha,reg,pf,hhat,llk_path,dev,pdev_path, &
#                       eta_path,null_dev,ierr,beta,betanorm,kkt,kkt0,&
#                       ncycles,nupdates,nbeta,idvars,llk_null,llk_sat,grad,alf)

source("test/prep_data.R")
t0 <- proc.time()
nlam <- 30
out <- .Fortran("mtcox_solutionpath", PACKAGE = "MTCox",
                ntasks = as.integer(ntasks), ii = as.integer(ii), io = as.integer(io),
                p = as.integer(p), ns = as.integer(ns), delta = as.double(data[, did]),
                w = as.double(data[,wid]), x = as.double(X), d = as.double(ds),
                iski = as.integer(iski), isko = as.integer(isko),
                nski = as.integer(nski), nsko = as.integer(nsko),
                n = as.integer(n), iex = as.integer(rep(1,p)), str = as.integer(1),
                dfmax = as.integer(p), pmax = as.integer(p), nlam = as.integer(nlam),
                earlystop = as.integer(1), alg = as.integer(1), lamfrac = as.double(1e-3),
                lam = double(nlam), maxit = as.integer(1e4), eps = as.double(1e-4),
                frac = as.double(0.1), alpha = as.double(0), reg = as.integer(0),
                pf = as.double(rep(1,p)),
                hhat = double(sum(ns)*nlam),
                llk_path = double(ntasks*nlam),
                dev = double(nlam),
                pdev_path = double(nlam),
                eta_path = double(n*nlam),
                null_dev = double(1),
                ierr = integer(5),
                beta = double(p* ntasks* nlam),
                betanorm = double(p* nlam),
                kkt = double(p*nlam),
                kkt0 = double(p*nlam),
                ncycles = integer(1),
                nupdates = integer(1),
                nbeta = integer(nlam),
                entered = integer(p),
                llk_null = double(ntasks),
                llk_sat = double(ntasks),
                grad = double(p*ntasks),
                alf=double(1)
)
proc.time()-t0

# output checks
out$ncycles
out$nupdates
out$ierr
cbind(lambda = out$lam,
      nbeta = out$nbeta,
      deviance = out$dev,
      percent.deviance = round(out$pdev_path*100,12))
out$entered
matplot(t(matrix(out$llk_path,ntasks,nlam)))
matplot(t(matrix(out$betanorm,p,nlam)))
matplot(t(matrix(out$kkt,p,nlam)))
lines(out$lam)
matplot(t(matrix(out$kkt0,p,nlam)))
lines(out$lam)
