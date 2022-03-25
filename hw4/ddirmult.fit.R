# This is the main function - calculate MLE for the Newton's method

#' MLE of Dirichlet-multinomial distribution by Newton's method
#'
#' \code{dirmultfit} finds the MLE of Dirichlet-multinomial by Newton method
#'
#' @param X an n-by-d matrix of counts
#' @param weights an n vector of observation weights
#' @param alpha0 an n vector of starting point
#' @param tolfun convergence tolerance in log-likelihood values
#' @param maxiters maximum number of iterations
#' @param display verbose mode (TRUE) or not (FALSE)
#' @return a list that contains
#' alphaHat MLE
#' se standard error of MLE
#' maximum log-likelihood at MLE
#' iterations number of iterations performed
#' gradient gradient at MLE
#' obsinfo observed information at MLE
#' obsinfoInv inverse of observed information at MLE
#' @seealso
dirmultfit <- function(X, weights = NULL, alpha0 = NULL,
                       tolfun = 1e-6, maxiters = 100, display = FALSE) {
  # check observation weights
  if (!is.null(weights)) {
    weights <- weight[rowSums(X) != 0]
    if (any(weights < 0))
      stop("Error: observation weigths should be positive")
  }
  # remove data points with batch size 0
  rsum <- rowSums(X)
  if (any(rsum == 0)) {
    rmv <- sum(rsum == 0)
    message(paste( "Warning: ", rmv,
                   " rows are removed because the row sums are 0"))
  }
  # remove the bins with no observations
  csum <- colSums(X)
  if (any(csum == 0)) {
    rmv <- sum(csum == 0)
    message(paste("Warning: ", rmv,
                  " columns are removed because the column sums are 0"))
  }
  # cleaned up data
  data <- X[rsum != 0, csum != 0]
  N <- nrow(data) ## Sample size
  d <- ncol(data) ## Number of parameters
  m <- rowSums(data) ## batch sizes
  # set default obs weights to be all 1s
  if (is.null(weights))
    weights <- rep(1, N)
  weightsum = sum(weights)
  # set starting points
  if (is.null(alpha0)) {
    # method of moment estimate
    rho <- sum(colSums(weights * (data / m)^2) / (colSums(weights * data / m)))
    alpha0 <- as.vector(colSums(weights * data / m) * (d - rho) / (rho - 1) / N)
    alpha0[alpha0 <= 0] = 1e-6
  } else {
    alpha0 <- alpha0[csum != 0]
    if (!is.vector(alpha0) && length(alpha0) != d) {
      stop("Error: dimension of alpha0 does not match data")
    } else if (any(alpha0 <= 0)) {
      # user provided starting values
      stop("Error: starting values should be positive")
    }
  }
  # prepare some variables before looping
  alphaHat <- alpha0
  alphaSum <- sum(alphaHat)
  loglIter <- sum(weights * ddirmult(data, alpha0, log = TRUE))
  # Print initial log-like if asked.
  if (display)
    print(paste("iteration ", 1, ", logL = ", loglIter, sep=""))
  # backtrack max iterations
  backtrackMaxiters <- 10
  ##----------------------------------------##
  ## The Newton loop
  if (maxiters == 1) iter <- 1
  else {
    for (iter in 2:maxiters) {
      # score vector
      alphaMat <- matrix(alphaHat, nrow = N, ncol = d, byrow = TRUE)
      score <- colSums(weights * digamma(data + alphaMat)) -
        weightsum * (digamma(alphaHat) - digamma(alphaSum)) -
        sum(weights * digamma(alphaSum + m))
      # observed info. matrix = diag(obsinfoDvec) - obsinfoC
      obsinfoDvec <- weightsum * trigamma(alphaHat) -
        colSums(weights * trigamma(data + alphaMat))
      obsinfoDvecInv <- 1 / obsinfoDvec
      obsinfoC <- weightsum * trigamma(alphaSum) -
        sum(weights * trigamma(m + alphaSum))
      # shrink c if necessary to make obs. info. pos def
      if (obsinfoC * sum(obsinfoDvecInv) >= 1) {
        if (display) print("shrink c")
        obsinfoC <- 0.95 / sum(obsinfoDvecInv)
      }
      # compute Newton direction
      newtondir <- obsinfoDvecInv * score
      newtondir <- newtondir +
        (sum(newtondir) / (1 / obsinfoC - sum(obsinfoDvecInv))) * obsinfoDvecInv
      # line search by step halving
      if (any(newtondir < 0)) {
        # make sure Newton iterate always lands within boundary
        stepsize <- min(- alphaHat[newtondir < 0] / newtondir[newtondir < 0])
        stepsize <- min(0.95 * stepsize, 1)
      } else {
        stepsize <- 1
      }
      for (btiter in 1:backtrackMaxiters) {
        alphaNew <- alphaHat + stepsize * newtondir
        loglNew <- sum(weights * ddirmult(data, alphaNew, log = TRUE))
        # line search successful if improving log-L
        if (loglNew > loglIter) break
        else if (btiter == backtrackMaxiters) {
          warning("line search failed")
        } else {
          if (display) print("step halving")
          stepsize <- stepsize / 2
        }
      }
      alphaHat <- alphaNew
      alphaSum <- sum(alphaHat)
      loglOld <- loglIter
      loglIter <- loglNew
      # Print the iterate log-like if requested
      if (display)
        print(paste("iteration ", iter, ", logL = ", loglIter, sep=""))
      # check convergence criterion
      if (abs(loglIter - loglOld) < tolfun * (abs(loglOld) + 1)) break
    }
  }
  ##----------------------------------------##
  ## End of Newton loop
  # score, i.e., gradient
  alphaMat <- matrix(alphaHat, nrow = N, ncol = d, byrow = TRUE)
  score <-
    colSums(weights * digamma(data + alphaMat)) -
    weightsum * (digamma(alphaHat) - digamma(alphaSum)) -
    sum(weights * digamma(alphaSum + m))
  # diagonal part of the observed information matrix
  obsinfoDvec <- weightsum * trigamma(alphaHat) -
    colSums(weights * trigamma(data + alphaMat))
  obsinfoDvecInv <- 1 / obsinfoDvec
  # the constant c in the observed information matrix
  obsinfoC <- weightsum * trigamma(alphaSum) -
    sum(weights * trigamma(m + alphaSum))
  # compute standard errors
  obsinfo <- diag(obsinfoDvec) - obsinfoC
  obsinfoInv <- diag(obsinfoDvecInv) +
    outer(obsinfoDvecInv, obsinfoDvecInv) / (1 / obsinfoC - sum(obsinfoDvecInv))
  se <- sqrt(diag(obsinfoInv))
  # restore to original data size
  if (any(csum == 0)) {
    colidx <- (csum != 0)
    # parameter estimate
    tmp <- alphaHat
    alphaHat <- rep(0, ncol(X))
    alphaHat[colidx] <- tmp
    # gradient/score vector
    tmp <- score
    score <- rep(0, ncol(X))
    score[colidx] <- tmp
    score[!colidx] <- weightsum * digamma(alphaSum) -
      sum(weights * digamma(alphaSum + m))
    # obs info matrix
    tmp <- obsinfo
    obsinfo <- matrix(0, ncol(X), ncol(X))
    obsinfo[colidx, colidx] <- tmp
    obsinfo[!colidx, ] <- - obsinfoC
    obsinfo[, !colidx] <- - obsinfoC
    # inverse of obs info matrix
    tmp <- obsinfoInv
    obsinfoInv <- matrix(0, ncol(X), ncol(X))
    obsinfoInv[colidx, colidx] <- tmp
  }
  # output
  list(estimate = alphaHat, se = se, maximum = loglIter, iterations = iter,
       gradient = score, obsinfo = obsinfo, obsinfoInv = obsinfoInv)
}