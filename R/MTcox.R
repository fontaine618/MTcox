MTcox <- function(data, yid, taskid, did, ...){
   call <- match.call()
   #catch all optional arguments
   opt <- list(...)
   for (name in names(opt) ) assign(name, opt[[name]])
   optionalParamNames <- c('wid','xid','pf','lambda','lambda.frac',
                           'n.lambda','str','maxit','earlystop','alg','dfmax',
                           'pmax','reg','eps','frac','alpha')
   unusedParams <- setdiff(names(opt),optionalParamNames)
   if(length(unusedParams))
      stop('Unused parameters: ',paste(unusedParams,collapse = ', '))
   #-----------------------------------
   # Add weight if not supplied
   if(is.null(opt$wid)){
      n <- nrow(data)
      data$w <- rep(1/(n*length(yid)),n)
      wid <- which(names(data) == 'w')
   }else{
      wid <- opt$wid
   }
   #-----------------------------------
   # Simultaneous or not
   if(length(yid)>1){
      nc <- ncol(data)
      # From tidyr package
      data <- tidyr::gather(data, key = 'task', value = 'y', -seq(nc)[!seq(nc) %in% yid])
      yid <- which(names(data) == 'y')
      taskid <- which(names(data) == 'task')
      wid <- which(names(data) == 'w')
      did <- which(names(data) == 'd')
   }
   #-----------------------------------
   # prepare dataset additional informations (ordering, risk sets, etc.)
   if(is.null(opt$xid))  D <- prepare.dataset(data, yid, taskid, wid, did)
   else  D <- prepare.dataset(data, yid, taskid, wid, did, opt$xid)
   #-----------------------------------
   # penalty factors (in order of apparition in original data)
   if(is.null(opt$pf)){
      pf <- as.double(rep(1,p))
   }else{
      if(length(opt$pf)!=p) stop("The size of penalty factor must be same as the number of input variables")
      if(any(opt$pf<0)) stop('pf must be non-negative')
      pf <- as.double(opt$pf)
   }
   #-----------------------------------
   # regularization parameters
      if(is.null(opt$lambda))
      {
         if(is.null(opt$lambda.frac)) opt$lambda.frac <- 1e-3
         if(is.null(opt$n.lambda)) opt$n.lambda <- 100
         if(opt$lambda.frac>=1) stop("lambda.frac should be less than 1 when no lambda values are supplied")
         if(opt$lambda.frac<1.0E-6) stop("lambda.frac is too small")
         lamfrac <- as.double(opt$lambda.frac)
         if(opt$n.lambda<5) stop("n.lambda should be at least 5")
         if(opt$n.lambda>1000) stop("n.lambda should be at most 1000")
         nlam <- as.integer(opt$n.lambda)
         userlam <- double(opt$n.lambda) #userlam=0 if lambda is is.null
      }
      else
      {
         #user defined lambda
         lamfrac <- as.double(1)#will trigger use of user lambda in fortran core
         if(any(opt$lambda<0)) stop("lambdas should be non-negative")
         userlam <- as.double(rev(sort(opt$lambda))) #lambda is declining
         nlam <- as.integer(length(userlam))
      }
   #-----------------------------------
   # options (integers)
      #strong rule (default 1= do strong rule)
         str <- ifelse(is.null(opt$str), 1, opt$str)
         if(! (str %in% 0:1)) stop("str must be either 0 (no strong rule) or 1 (strong rule) .")
      #maximum number of variables in model
         dfmax <- ifelse(is.null(opt$dfmax), p, opt$dfmax)
         dfmax <- min(max(dfmax,0),p)
         pmax <- ifelse(is.null(opt$pmax), dfmax*1.2, opt$dfmax)
         pmax <- min(max(pmax,0),p)
      #early stop trigger
         earlystop <- ifelse(is.null(opt$earlystop), 1, opt$earlystop)
         if(! (earlystop %in% 0:1)) stop("earlystop must be either 0 (no early stop) or 1 (early stop) .")
      #which algorithm to use
         alg <- ifelse(is.null(opt$alg), 1, opt$alg)
         if(! (alg %in% 0:2)) stop("alg must be either 0 (IRLS), 1 (majorized) or 2 (backtracking).")
      #maximum number if iterations(number of updates)
         maxit <- ifelse(is.null(opt$maxit), 1e5, opt$maxit)
         maxit <- min(max(maxit,1e3),1e7)
      #type of regularization
         reg <- ifelse(is.null(opt$reg), 2, opt$reg)
         if(! (reg %in% c(0,2))) stop("reg must be either 0 (L_infinity) or 2 (L_2).")
      #group all options in a vector
         optInt <- as.integer(c(str=str,maxit=maxit,earlystop=earlystop,alg=alg,
                                dfmax=dfmax,pmax=pmax,reg=reg))
   #-----------------------------------
   # options (doubles)
      #convergence threshold
         eps <- ifelse(is.null(opt$eps), 1e-3, opt$eps)
         eps <- min(max(eps,1e-9),1e-1)
      #fraction for backtracking
         frac <- ifelse(is.null(opt$frac), 0.1, opt$frac)
         frac <- min(max(frac,1e-3),0.99)
      #sparsity parameter (0 if no L1, 1 is only L1 i.e. Lasso)
         alpha <- ifelse(is.null(opt$alpha), 0, opt$alpha)
         alpha <- min(max(alpha,0),1)
      #fraction for regularization parameter
         lamfrac <- ifelse(is.null(opt$lamfrac), 1e-3, opt$lamfrac)
         lamfrac <- min(max(lamfrac,1e-6),0.99)
      #group all options in a vector
         optDbl <- as.double(c(eps=eps,frac=frac,alpha=alpha,lamfrac=lamfrac))
   #-----------------------------------
   # Call Fortran core
   D
   # Prepare output
}


