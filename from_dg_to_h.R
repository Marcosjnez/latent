f <- function(x, y, w, X) {

  res <- y-X%*%x
  return(sum(w*res^2))

}
g <- function(x, y, w, X) {

  res <- y-X%*%x
  return(-2*t(X) %*% (w*res))

}


set.seed(2025)
w <- runif(10)
y <- rnorm(10)
X <- replicate(5, rnorm(10))
x <- rnorm(5)
numDeriv::grad(f, x, y=y, w=w, X=X)
g(x,y,w,X)
numDeriv::hessian(f, x, y=y, w=w, X=X)

dx <- rnorm(5)
-2*t(X) %*% (w*(-X%*%dx))
eps <- 1e-04
(g(x+eps*dx,y,w,X) - g(x-eps*dx,y,w,X)) / (2*eps)
repl <- c(0,1,0,0,0)
-2*t(X) %*% (w*(-X%*%repl))
