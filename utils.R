predict.MTcox <- function(x, taskids, beta, log=F){
   out <- 0*taskids
   for(k in unique(taskids)){
      out[taskids == k] <- x[taskids == k,] %*% t(t(beta[,k]))
   }
   if(!log) out <- exp(out)
   out
}

