
Y <- survival::Surv(time=data$y, event= 1-data$d, type="right")

sorted <- order(-Y[, 1], Y[, 2])
newstrat <- rep.int(0L, nrow(Y))
newstrat[1] <- 1L

if (storage.mode(Y) != "double")
   storage.mode(Y) <- "double"
counts <- .Call(survival:::Ccoxcount1, Y[sorted, ], as.integer(newstrat))


counts$nrisk

D <- prepare.dataset(data,yid=yid,taskid=taskid,did=did,wid=wid)
cumsum(rev(D$isko - D$iski + 1))-counts$nrisk

D$data



