#' Title
#'
#' description
#'
#' details
#'
#' @param x is a matrix containing the observations in columns
#' @param mean mean vector of the multivariate gaussian distribution
#' @param varcovM variance covariance matrix
#' @param Log logical parameter (default TRUE)
#'
#' @return list of the results
#' @export
#'
#' @examples
mvnpdf <- function(x, mean =  rep(0, nrow(x)),
                   varcovM = diag(nrow(x)), Log = TRUE) {
  n <- ncol(x)
  p <- nrow(x)
  x0 <- x - mean
  Rinv <- solve(varcovM)
  LogDetvarcovM <- log(det(varcovM))

  y <- NULL
  for (j in 1:n) {
    yj <- - p/2 * log(2*pi) - 0.5 * LogDetvarcovM -
      0.5 * t(x0[, j]) %*% Rinv %*% x0[, j]
    y <- c(y, yj)
  }

  if (!Log) {
    y <- exp(y)
  }

  res <- list(x = x, y = y)
  class(res) <- "mvnpdf"
  return(res)
}

#' Plot of the mvnpdf function
#'
#' @param x an object of class \code{mvnpdf} resulting from a call of
#' \code{mnvpdf()} function.
#' @param ... graphical parameters passed to \code{plot()} function.
#'
#' @return Nothing is returned, only a plot is given.
#' @export
#'
#' @examples
#' pdfvalues <- mvnpdf(x=matrix(seq(-3, 3, by = 0.1), nrow = 1), Log=FALSE)
#' plot(pdfvalues)
plot.mvnpdf <- function(x, ...) {
  plot(x$x, x$y, type = "l", ...)
}


# mvnpdf <- function(x, mean, varcovM, Log =TRUE) {
#  p = nrow(x)
#  (2* pi)^(-p/2) * det(varcovM)^(-1/2) * exp((-1/2)*t(x-mean) %*% solve(varcovM) %*% (x-mean))
#}

#p = nrow(varcovM)

#X <- matrix(c(-0.5,1.5, 0,1,-1,1), nrow=2)
#results<- c()

#for(j in 1:ncol(X)){
#  results <- c(results, mvnpdf(X[,j, drop=FALSE], mean, varcovM))
#}
#results



#M <- matrix(c(1,2,2,3), nrow=2,byrow=TRUE)
#det(M)

# Um in einer Matrix zu multiplizieren, braucht man %*%
#X <- matrix(c(-0.5, 1.5), nrow=2)

#solve(M)
#M %*% solve (M)
#t(x)
#t(x) %*% solve (M) %*% x

#f <- function(x) {
#  return(x^2)
#}

# runif gibt mir random numbers
#y <-runif(10)
#y

#z<- c()
#z

#for (i in 1:10) {
# z <- c(z,f(y[i]))
#}

# Um in einer Matrix zu multiplizieren, braucht man %*%

# transpose matrix:

#mean <- matrix(c(0, 0), nrow=2)
#varcovM <- matrix(c(1, 0, 0,1), nrow=2)

#mvnpdf(X, mean, varcovM)
#mvtnorm::dmvnorm(t(X), mean, varcovM)


