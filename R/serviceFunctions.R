#' Moving average for a vector
#'
#' This is a wrapper for the 'filter()' function that uses the convolution method
#' to calculate the moving averages of the terms in a vector, i.e. a series
#' of averages of different subsets of the full vector
#'
#' @param x a numeric vector representing a time series
#' @param n an integer, representing the number of values that compose each subset to be averaged
#' @param sides can be either 1 or 2. If sides = 1 the n values to be averaged are taken before
#' the preent values; If sides = 2 the n values are taken before and after the present value; If n
#' is odd, (n - 1)/2 values are taken before the present value and (n - 1)/2 are
#' taken after the present value, while, if n is even, more values are taken after
#' the present value
#'
#' @return This function returns a vector of moving averages
#' @author Andrea Onofri
#'
#' @examples
#' series <- c(319, 317, 332, 271, 301, 292, 351, 358, 259, 270)
#' ma(series, n = 4, sides = 2)
#'
ma <- function(x, n = 5, sides = 2){
  # moving average for a vector with the convolution method
  res <- stats::filter(x, rep(1 / n, n),
                       method = "convolution", sides = sides)
  as.numeric(res)
  }


