#' Title
#'
#' @param x
#' @param taskids
#' @param beta
#' @param log
#'
#' @return
#' @export
#'
#' @examples
predict.MTcox <- function(x, taskids, beta, log=F){
   out <- 0*taskids
   for(k in unique(taskids)){
      out[taskids == k] <- x[taskids == k,] %*% t(t(beta[,k]))
   }
   if(!log) out <- exp(out)
   out
}

#' Title
#'
#' @param data
#' @param yid
#' @param taskid
#' @param wid
#' @param did
#' @param xid
#'
#' @return
#' @export
#'
#' @examples
prepareDataset <- function(data, yid, taskid, wid, did, xid){
   # Extract dimensions
   if(missing(xid)) xid <- seq(ncol(data))[-c(yid,taskid,did,wid)]
   vars <- names(data)[-c(yid,taskid,did,wid)]
   p <-length(vars)

   # Prepare dataset
   ord <- order(data[,taskid],data[,yid],data[,did])
   data[,taskid] <- as.factor(data[,taskid])
   tasks <- levels(data[,taskid])
   K <- length(tasks)
   data <- data[ord,]
   ntasks <- length(tasks)
   n <- nrow(data)
   ii <- sapply(tasks, function(k) min(which(data[,taskid] == k)))
   io <- sapply(tasks, function(k) max(which(data[,taskid] == k)))
   ns <- rep(0,K)
   iskid <- rep(0,n)
   for(k in seq(K)){
      index <- seq(ii[k],io[k])[data[seq(ii[k],io[k]),did] == 0]
      yevent <- unique(data[index,yid])
      ids <- sapply(yevent, function(y) which.max(y == data[seq(ii[k],io[k]),yid] & data[seq(ii[k],io[k]),did] == 0))
      ns[k] <- length(ids)
      iskid[ids+ii[k]-1] <- 1
   }
   sns <- sum(ns)
   ds <- rep(0, sns)
   nski <- cumsum(c(1,ns[-K]))
   nsko <- cumsum(ns)
   iskid <- iskid * cumsum(iskid)
   iski <- rep(0, sns)
   isko <- rep(0, sns)
   for(isk in seq(sns)){
      iski[isk] <- which(iskid == isk)
      k <- which(data[iski[isk], taskid] == tasks)
      isko[isk] <- min(which(iskid == isk+1)-1, io[k])
   }
   # sum of weights for each distinct failure time
   ds <- sapply(seq(sns), function(isk){
      sum( data[seq(iski[isk], isko[isk]), wid] * (1-data[seq(iski[isk], isko[isk]), did]) )
   })
   Xmat <- apply(as.matrix.noquote(data[, -c(yid,did,taskid,wid)]),2,as.numeric)
   colnames(Xmat) <- vars
   iex <- vars %in% names(data)[xid]
   data$isk <- iskid
   out <- list(data=data, yid=yid, taskid=taskid, wid=wid, did=did,
               tasks=tasks, K=K, vars=vars, p=p, ntasks=ntasks, n=n,
               ii=ii, io=io, ns=ns, sns=sns, iski=iski, isko=isko,
               nski=nski, nsko=nsko, ds=ds, X=Xmat, iex=iex)
   class(out) <- "MTcox.dataset"
   out
}


#' dataCox
#'
#' @param lambda
#' @param rho
#' @param x
#' @param beta
#' @param censRate
#'
#' @return
#' @export
#'
#' @examples
dataCox <- function(lambda, rho, x, beta, censRate){
   N <- nrow(x)
   # real Weibull times
   u <- stats::runif(N)
   Treal <- (- log(u) / (lambda * exp(x %*% beta)))^(1 / rho)
   # censoring times
   Censoring <- stats::rexp(N, censRate)
   # follow-up times and event indicators
   time <- pmin(Treal, Censoring)
   status <- as.numeric(Treal <= Censoring)
   # data set
   data.frame(y=time, d=1-status, w=1/N, task= as.factor('t1'), X=x)
}
