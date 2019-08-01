library(dplyr)

############################################################################
prepareData <- function(formula, data) {
   # Parameter identification as in  `survival::coxph()`.
   Call <- match.call()
   indx <- match(c("formula", "data"),
                 names(Call), nomatch = 0)
   if (indx[1] == 0)
      stop("A formula argument is required")
   temp <- Call[c(1, indx)]
   temp[[1]] <- as.name("model.frame")

   mf <- eval(temp, parent.frame())
   Y <- model.extract(mf, "response")

   if (!inherits(Y, "Surv"))
      stop("Response must be a survival object")
   type <- attr(Y, "type")

   if (type != "right" && type != "counting")
      stop(paste("Cox model doesn't support \"", type, "\" survival data",
                 sep = ""))

   # collect times, status, variables and reorder samples
   # to make the algorithm more clear to read and track
   cbind(event = unclass(Y)[,2], # 1 indicates event, 0 indicates cens
         times = unclass(Y)[,1],
         mf[, -1]) %>%
      arrange(times)
}


##################################################################
set.seed(456)
dataCox <- function(N, lambda, rho, x, beta, censRate){

   # real Weibull times
   u <- runif(N)
   Treal <- (- log(u) / (lambda * exp(x %*% beta)))^(1 / rho)

   # censoring times
   Censoring <- rexp(N, censRate)

   # follow-up times and event indicators
   time <- pmin(Treal, Censoring)
   status <- as.numeric(Treal <= Censoring)

   # data set
   data.frame(id=1:N, time=time, status=status, x=x)
}

#x <- matrix(sample(-1:4, size = 3000, replace = TRUE), ncol = 3)
n <- 10
p <- 3
x <- matrix(rnorm(n*p)^2, ncol=p)
dataCox(n, lambda = 5, rho = 1.5, x, beta = rep(2,p), censRate = 10) -> dCox

##################################################################
batchData <- prepareData(formula = survival::Surv(time, status)~x.1+x.2+x.3,
                         data = dCox)
batchData <- arrange(batchData, -times)
beta = rep(0,p)
scores <- apply(batchData[, -c(1, 2)], MARGIN = 1,
                function(element) {
                   exp(element %*% beta)
                })

nominator <- apply(batchData[, -c(1, 2)], 2, function(x) cumsum(scores*x))
denominator <- cumsum(scores)
partial_sum <- (batchData[, -c(1, 2)] - nominator/denominator) * batchData[, "event"]

U_batch <- colSums(partial_sum)


# batchData to MTcox format
batchData2 <- batchData
batchData2$event <- 1-batchData2$event
batchData2$task <- "t1"
batchData2$weight <- 1

D <- prepareDataset(batchData2,yid=2,taskid=p+3,did=1,wid=p+4)

# gradient at 0 from grad null
#D <- fit$D
#D$X <- matrix(fit$x, n, p)

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
out$grad


# eta <- rep(0,n)
eta <-  D$X %*% beta
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
grad


round(rbind(U_batch,grad,U_batch/grad),3)

# ########################################################
# library(numDeriv)
# llk <- function(beta) {
#    eta <- D$X %*% beta
#    .Fortran("loglik", PACKAGE = "MTCox",
#                    ntasks = as.integer(D$K),
#                    n = as.integer(D$n),
#                    p = as.integer(D$p),
#                    ns = as.integer(D$ns),
#                    iski = as.integer(D$iski),
#                    isko = as.integer(D$isko),
#                    nski = as.integer(D$nski),
#                    nsko = as.integer(D$nsko),
#                    w = as.double(D$data[,wid]),
#                    x = as.double(D$X),
#                    delta = as.integer(D$data[,did]),
#                    d = as.double(D$ds),
#                    eta = as.double(eta),
#                    llk = double(D$K)
#    )$llk
# }
#
# grad_num <- numDeriv::grad(llk, beta, method= "simple", method.args=list(eps=1e-10))


#D$data[,did]
