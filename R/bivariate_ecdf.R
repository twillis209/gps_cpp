bivariate_ecdf <- function(u, v) {
  n <- length(u)
  
  ecdf <- numeric(n)
  
  for(i in seq(1, n)) {
    for(j in seq(1, n)) {
      if(u[j] <= u[i] && v[j] <= v[i]) {
        ecdf[i] <- ecdf[i] + 1
      }
    }

  }

  ecdf/n
}

a <- c(.1, .2, .3, .4, .5, .5)
b <- c(.5, .4, .3, .2, .1, .1)
