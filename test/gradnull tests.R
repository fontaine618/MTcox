
library(MTcox)

p <- 4
K <- 1
n <- 1000
#raw data
X = round(matrix(rnorm(n*p,0,1),n,p),2)
beta <- c(2,3,-1,rep(0,p-3))
#eventually do weibull distribution or something else
y <- round(rexp(n,exp(X %*% beta)),3)+0.001
data <- data.frame(
   y = y,
   d = rbinom(n,1,0.3),
   w=rep(1,n)/n,
   task = paste("t",seq(K),sep="")[sample.int(K,n, replace=TRUE)],
   X = X
)
yid <- 1
did <- 2
wid <- 3
taskid <- 4
nlam=20

D <- prepare.dataset(data,yid=yid,taskid=taskid,did=did,wid=wid)

# gradient at 0 from grad null

out <- .Fortran("gradnull", PACKAGE = "MTCox",
                ntasks = as.integer(D$K),
                p = as.integer(D$p),
                n = as.integer(D$n),
                ns = as.integer(D$ns),
                iski = as.integer(D$iski),
                isko = as.integer(D$isko),
                nski = as.integer(D$nski),
                nsko = as.integer(D$nsko),
                w = as.double(D$data[,wid]),
                x = as.double(D$X),
                delta = as.double(D$data[, did]),
                d = as.double(D$ds),
                grad = double(D$p*D$K)
)
# gradient at zero from grad_hessj

eta <- rep(0,D$n)
grad <- rep(0,D$p)
for(j in seq(D$p)){
   tmp <- .Fortran("grad_hessj", PACKAGE = "MTCox",
                   ntasks = as.integer(D$K),
                   n = as.integer(D$n),
                   ns = as.integer(D$ns),
                   iski = as.integer(D$iski),
                   isko = as.integer(D$isko),
                   nski = as.integer(D$nski),
                   nsko = as.integer(D$nsko),
                   w = as.double(D$data[,wid]),
                   x = as.double(D$X[,j]),
                   delta = as.double(D$data[, did]),
                   d = as.double(D$ds),
                   eta = as.double(eta),
                   grad = double(D$K),
                   hess = double(D$K)
   )
   grad[j] <- tmp$grad
}
rbind(null=out$grad, hessj=grad, prop=out$grad/ grad)
