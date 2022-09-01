f <- function(x){
  exp(-(abs(x-3)/2)^0.4)
}

curve(f,- 20, 20)
x <- seq(-150, 150, 0.01)
f.x <- f(x)
w <- f.x / sum(f.x)

n <- 1000
v.sim <- sample(x, n, replace=T, prob=w)

K <- 1/integrate(f, -Inf, Inf)$value

f.f <- function(x){
  K*f(x)
}
hist(v.sim, freq=F)
curve(f.f, add=T, col=2, lwd=2)

