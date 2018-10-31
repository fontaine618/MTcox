source("utils.R")
#parameters
p <- 10
K <- 5
n <- 1000
#raw data
X = matrix(rnorm(n*p,0,1),n,p)
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
wid <- NA
if(is.na(wid)){
   data$w <- 1/n
   wid <- ncol(data)
}
#pre processing
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

X <- apply(as.matrix.noquote(data[, -c(yid,did,taskid,wid)]),2,as.numeric)
#X <- apply(X, 2, function(col) (col-mean(col))/sd(col))

eta <- rep(0,n)
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
                #b = as.double(beta),
                d = as.double(ds),
                eta = as.double(eta),
                llk = double(ntasks)
)
llk0 <- out$llk
