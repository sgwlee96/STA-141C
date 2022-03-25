rm(list = ls())

# Load ddirmult.fit.R function 
source("hw4/ddirmult.fit.R")
# Load ddirmult.R function
source("hw4/ddirmult.R")

#---------------------------------------#
# Please read in the training data here #
traindata <- read.table("hw4/data/optdigits.tra", sep = ",")
traindata <- as.matrix(traindata)
#---------------------------------------# 

# pre-allocate variables
alphahat <- matrix(0, nrow = 10, ncol = 64)
phat <- matrix(0, nrow = 10, ncol = 64)
loglmn <- rep(0, 10)
logldirmn <- rep(0, 10)
digitCount <- rep(0, 10)
iters <- rep(0, 10)
runtime <- rep(0, 10)
# loop over digit
for (dig in 0:9) {
  # MLE for Dirichlet-multinomial
  digitdata <- traindata[traindata[, 65] == dig, -65]
  digitCount[dig + 1] <- nrow(digitdata)
  ptm <- proc.time() # start timing
  digMLE <- dirmultfit(digitdata, tolfun = 1e-8, display = TRUE)
  runtime[dig + 1] <- (proc.time() - ptm)[3] # stop timing
  alphahat[dig + 1, ] <- digMLE$estimate
  iters[dig + 1] <- digMLE$iterations
  logldirmn[dig + 1] <- digMLE$maximum
  # MLE and log-likelihood for multinomial
  phat[dig + 1, ] <- colSums(digitdata) / sum(digitdata)
  cidx <- (colSums(digitdata) != 0)
  loglmn[dig + 1] <- sum(lgamma(rowSums(digitdata[, cidx]) + 1)) -
    sum(lgamma(digitdata[, cidx] + 1)) +
    sum(colSums(digitdata[, cidx]) * log(phat[dig + 1, cidx]))
}
