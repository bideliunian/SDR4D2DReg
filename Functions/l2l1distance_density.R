#######################################
## l2 distance between two density ##
######################################
l2density <- function (x, y, nSup = 201) {
  # Param
  #   x, y: two vectors, each contains iid observations from a density
  # Return
  #   l2 distance between the estimated density curves
  lower <- range(c(x, y))[1]
  upper <- range(c(x, y))[2]
  a <- density(x,from = lower, to = upper, n = nSup, kernel = 'gaussian')$y 
  b <- density(y,from = lower, to = upper, n = nSup, kernel = 'gaussian')$y
  sup <- density(x,from = lower, to = upper, n = nSup, kernel = 'gaussian')$x
  m <- length(a)
  n <- length(b)
  stopifnot(m > 0 && n > 0)
  return(sqrt(pracma::trapz(sup, (a - b)^2)))
}

#######################################
## l1 distance between two density ##
######################################
l1density <- function (x, y, nSup = 201) {
  # Param
  #   x, y: two vectors, each contains iid observations from a density
  # Return
  #   l1 distance between the estimated density curves
  lower <- range(c(x, y))[1]
  upper <- range(c(x, y))[2]
  a <- density(x,from = lower, to = upper, n = nSup, kernel = 'gaussian')$y 
  b <- density(y,from = lower, to = upper, n = nSup, kernel = 'gaussian')$y
  sup <- density(x,from = lower, to = upper, n = nSup, kernel = 'gaussian')$x
  m <- length(a)
  n <- length(b)
  stopifnot(m > 0 && n > 0)
  return(pracma::trapz(sup, abs(a - b)))
}
