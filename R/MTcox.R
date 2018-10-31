MTcox <- function(data, yid, taskid, wid, did){

   # Add weight if not supplied
   if(missing(wid)){
      n <- nrow(data)
      data$w <- rep(1/n,n)
      wid <- which(names(data) == 'w')
   }

   # Simultaneous or not
   if(length(yid)>1){
      nc <- ncol(data)
      # From tidyr package
      data <- gather(data, key = 'task', value = 'y', -seq(nc)[!seq(nc) %in% yid])
      yid <- which(names(data) == 'y')
      taskid <- which(names(data) == 'task')
      wid <- which(names(data) == 'w')
      did <- which(names(data) == 'd')
   }

   # Extract dimensions
   tasks <- unique(data[,taskid])
   K <- length(tasks)
   vars <- names(data)[-c(yid,taskid,did,wid)]
   p<-length(vars)

   # Prepare dataset
   ord <- order(data[,taskid],data[,yid],data[,did])
   data <- data[ord,]
   tasks <- unique(data[,taskid])
   ntasks <- length(tasks)
   n <- nrow(data)
   ii <- sapply(tasks, function(k) min(which(data[,taskid] == k)))
   io <- sapply(tasks, function(k) max(which(data[,taskid] == k)))
   ns <- sapply(tasks, function(k) length(unique(data[ii[k]:io[k],yid][data[ii[k]:io[k],did]==0])))
   nsi <- cumsum(c(1,ns[-K]))
   nso <- cumsum(ns)
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

   # Check input
   X <- apply(as.matrix.noquote(data[, -c(yid,did,taskid,wid)]),2,as.numeric)

   # Call Fortran core

   # Prepare output
}
data <- data.frame(y1=1:10, y2=111:120, X1=rnorm(10),X2=rnorm(10),X3=rnorm(10),d=rep(0,10))
yid <- c(1,2)
did <- 6


