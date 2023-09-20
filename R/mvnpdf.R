################################################################################
##############################   HeavyR   ######################################
############################### Chapter 1 ######################################
################################################################################


# si on met un = dans une fonction, ça donne les variables par défaut
# donc ici, si Log = T, par défaut le paramètre Log sera associé à T, il faudra préciser F si on veut pas avoir Log = T


#' My package
#'
#'  description
#'
#' détails
#'
#' @param x x must be a matrix, containing the observations of interest
#' @param mean mean is the vector of the multivariate gaussian distribution
#' @param varcovM variance covariance matrix (squared matrix)
#' @param Log by default, Log = T
#'
#' @return A list with the matrix and the value of the density of a multivariate normal distribution
#' @export
#'
#' @examples x <- 1


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
