predict.MTcox <- function(x, taskids, beta, log=F){
   out <- 0*taskids
   for(k in unique(taskids)){
      out[taskids == k] <- x[taskids == k,] %*% t(t(beta[,k]))
   }
   if(!log) out <- exp(out)
   out
}

prepare.dataset <- function(data, yid, taskid, wid, did, xid){
   # Extract dimensions
   tasks <- unique(data[,taskid])
   K <- length(tasks)
   if(missing(wid)) wid <- NA
   if(missing(xid)) xid <- seq(ncol(data))[-c(yid,taskid,did,wid)]
   vars <- names(data)[-c(yid,taskid,did,wid)]
   p <-length(vars)

   # Prepare dataset
   ord <- order(data[,taskid],data[,yid],data[,did])
   data <- data[ord,]
   ntasks <- length(tasks)
   n <- nrow(data)
   ii <- sapply(tasks, function(k) min(which(data[,taskid] == k)))
   io <- sapply(tasks, function(k) max(which(data[,taskid] == k)))
   ns <- sapply(tasks, function(k) length(unique(data[ii[k]:io[k],yid][data[ii[k]:io[k],did]==0])))
   sns <- sum(ns)
   nsk <- rep(0,sns)
   iski <- rep(0, sns)
   isko <- rep(0, sns)
   isk<-1
   for(k in seq(K)){
      i <- ii[k]
      #skip first censored which are in no risk set
      while(data[i,did]==1)i<-i+1
      while(i<=io[k]){
         #first will always be a failure
         ytmp <- data[i,yid]
         iski[isk] <- i
         nsk[isk] <- 1
         #check next values
         mx <- i<n
         #add all failures with same time
         while(data[min(i+1,n), yid]==ytmp &
               data[min(i+1,n), did]==0 &
               data[min(i+1,n), taskid]==k & mx){
            nsk[isk]<-nsk[isk]+1
            if(i+1==n) mx=F
            i<-i+1
         }
         #add all censored before next failure time
         while(data[min(i+1,n), did]==1 &
               data[min(i+1,n), taskid]==k & mx){
            nsk[isk]<-nsk[isk]+1
            if(i+1==n) mx=F
            i<-i+1
         }
         isko[isk] <- min(i,n)
         isk <- isk+1
         i<-i+1
      }
   }
   nski <- cumsum(c(1,ns[-K]))
   nsko <- cumsum(ns)
   ds <- sapply(seq(sns), function(isk){
      sum(data[seq(iski[isk], isko[isk]), wid] *
             (1-data[seq(iski[isk], isko[isk]), did]))
   })
   X <- apply(as.matrix.noquote(data[, -c(yid,did,taskid,wid)]),2,as.numeric)
   colnames(X) <- vars
   iex <- vars %in% names(data)[xid]
   out <- list(data=data, yid=yid, taskid=taskid, wid=wid, did=did,
               tasks=tasks, K=K, vars=vars, p=p, ntasks=ntasks, n=n,
               ii=ii, io=io, ns=ns, sns=sns, iski=iski, isko=isko,
               nski=nski, nsko=nsko, ds=ds, X=X, iex=iex)
   class(out) <- "MTcox.dataset"
   out
}
