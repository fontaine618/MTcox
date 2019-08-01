dataCox2 <- function(N, lambda, rho, x, beta, censRate){

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
llk <- function(beta) {
   eta <- D$X %*% beta
   .Fortran("loglik", PACKAGE = "MTCox",
            ntasks = as.integer(D$K),
            n = as.integer(D$n),
            p = as.integer(D$p),
            ns = as.integer(D$ns),
            iski = as.integer(D$iski),
            isko = as.integer(D$isko),
            nski = as.integer(D$nski),
            nsko = as.integer(D$nsko),
            w = as.double(D$data[,wid]),
            x = as.double(D$X),
            delta = as.integer(D$data[,did]),
            d = as.double(D$ds),
            eta = as.double(eta),
            llk = double(D$K)
   )$llk
}


x <- matrix(sample(0:1, size = 3000, replace = TRUE), ncol = 3)
x <- apply(x, 2, function(col) (col-mean(col))/sd(col))
#x <- matrix(rnorm(3000)^2, ncol=3)
true.beta <- c(1,2,3)
dataCox2(10^3, lambda = 5, rho = 1.5, x, beta = true.beta, censRate = 0.0001) -> dCox

DCox <- data.frame(y=dCox$time,d=1-dCox$status,   w=1/1000,task =as.factor("t1"),
                   X.1=dCox$x.1, X.2=dCox$x.2, X.3=dCox$x.3)
wid=3
did=2
D<- prepare.dataset(DCox,yid=1,taskid=4,did=did,wid=wid)
Y <- survival::Surv(time=D$data$y, event= 1-D$data$d, type="right")

# FROM MTcox
# FROM mle by coxph
fit <- survival::coxph(Y ~ D$X, ties="breslow", weights = rep(1/1000,1000))
llk.MT <- c(null=llk(rep(0,3)),mle=llk(coef(fit)))
beta <- coef(fit)
beta
rbind(MTcox = c(llk.MT, diff = diff(llk.MT)),
      coxph = c(fit$loglik,diff(fit$loglik)))

# null likelihood
-sum(log(D$isko))/1000#w=1
-sum(D$ds*log(D$isko/1000))#w=1/1000
