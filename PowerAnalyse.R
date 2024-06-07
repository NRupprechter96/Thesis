N <- 600
X <- runif(N, min = 1, max = 7)
M <- abs(a * X + rnorm(N))
Y <- abs(b * M + c * X + rnorm(N))
exoData <- data.frame(X, M, Y)

N <- nrow(exoData)
a <- .4 #effect of X on M
b <- .4 #effect of M on Y, controlling for X
c <- .1 #effect of X on Y, controlling for M
set.seed(777)

kappa.free <- kappa.pop <- matrix(NA, nrow = 2, ncol = 1, dimnames = list(c("M","Y"), "X"))
kappa.free["M","X"] <- "a"
kappa.free["Y","X"] <- "c"
kappa.pop["M","X"] <- a
kappa.pop["Y","X"] <- c

exoPaths <- bind(free = kappa.free, popParam = kappa.pop)
beta.free <- beta.pop <- matrix(0, nrow = 2, ncol = 2, dimnames = list(c("M","Y"), c("M","Y")))
beta.free["Y", "M"] <- "b"
beta.pop["Y", "M"] <- b
endoPaths <- bind(free = beta.free, popParam = beta.pop)
residCor <- binds(free = diag(as.numeric(NA), 2), popParam = diag(2))
userParams <- ' ind := a * b
total := ind + c '
simMod1 <- model.path(BE = endoPaths, RPS = residCor, KA = exoPaths, con = userParams,
                      indLab = rownames(kappa.free), covLab = colnames(kappa.free))

rejectMCCI <- function(object) {
  CIs <- semTools::monteCarloCI(object)
  apply(CIs, 1, function(CI) 0 < CI["ci.lower"] | 0 > CI["ci.upper"])
}
out1 <- sim(nRep = 5000, model = simMod1, covData = exoData,
            n = nrow(exoData), seed = 777, outfun = rejectMCCI)
summaryParam(out1, matchParam = TRUE, digits = 3)
testMCCI <- do.call(rbind, getExtraOutput(out1))
colMeans(testMCCI) #empirical estimate of power for Monte Carlo CIs

print(exoData)


