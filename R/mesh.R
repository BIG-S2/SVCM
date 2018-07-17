meshgrid<- function(x, y) {
  
  
  x <- c(x); y <- c(y)
  n <- length(x)
  m <- length(y)
  
  X <- matrix(rep(x, each = m),  nrow = m, ncol = n)
  Y <- matrix(rep(y, times = n), nrow = m, ncol = n)
  
  return(list(x = X, y = Y))
}
