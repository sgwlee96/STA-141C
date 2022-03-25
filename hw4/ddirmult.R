# This is the function to evaluate the Dirichlet-multinomial density. 
# Please do not change anything inside the function!


#' Density of Dirichlet-multinomial
#'
#' \code{ddirmult} evaluates the Dirichlet-multinomial density
#'
#' @param x an n-by-d matrix of counts
#' @param alpha an n-by-d or n-by-1 matrix of parameters
#' @param log log-density (TRUE) or density (FALSE)
#' @return an n vector of (log)-densities
#' @seealso
ddirmult <- function(x, alpha, log = FALSE) {
  # check x is count data
  if (any(x < 0))
    stop("Error: x should be nonnegative count data")
  # check positivity of alpha
  if (any(alpha < 0))
    stop("Error: Dirichlet-multinomial parameter alpha should be nonnegative")
  # check dimensions of input arguments
  if (is.vector(alpha) && length(alpha) > 1) {
    if (is.vector(x) && length(x) > 1) {
      if (length(x) != length(alpha)) {
        stop("Error: sizes of x and alpha do not match.")
      } else {
        # expand x to a matrix of matching size with alpha
        x <- matrix(x, 1, length(x))
      }
    } else if (is.vector(x) && length(x) <= 1) {
      stop("Error: x can not be a scalar")
    }
    # expand alpha to a matrix of matching size with x
    alpha <- matrix(alpha, nrow = nrow(x), ncol = length(alpha), byrow = TRUE)
  }
  if (any(dim(alpha) != dim(x)))
    stop("Error: dimensions of alpha and x do not match")
  # compute log-densities
  alphaRowSums <- rowSums(alpha)
  xRowSums <- rowSums(x)
  # lgamma(0 + 0) - lgamma(0) will produce NaNs.
  # This fixes x_{ij} = alpha_{ij} = 0 cases.
  # Is there better way to deal with this?
  alpha[(x == 0) & (alpha == 0)] <- 1
  # assemble log-likelihood
  # lgamma(0) throws a lot of warnings but they are valid
  logl <- suppressWarnings(
    lfactorial(xRowSums) + rowSums(lgamma(x + alpha)) +
      lgamma(alphaRowSums) - (rowSums(lfactorial(x)) +
                                rowSums(lgamma(alpha)) + lgamma(alphaRowSums + xRowSums))
  )
  # Deal with alphaRowSums == 0 cases
  # fix the lgamma(0 + 0) - lgamma(0) produces NaN here
  logl[(xRowSums == 0) & (alphaRowSums == 0)] <- 0
  logl[(xRowSums > 0) & (alphaRowSums == 0)] <- - Inf
  # output
  if (log)
    return(logl)
  else
    return(exp(logl))
}