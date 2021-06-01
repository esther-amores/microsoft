## ----setup, include=FALSE------------------------------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, warning = FALSE)


## ------------------------------------------------------------------------------------------------------
# Paquets
library(poweRlaw)

# Dades
data("danish", package="evir")


## ------------------------------------------------------------------------------------------------------
ePL <- function(xdt){
  xm <- min(xdt)
  xi <- mean(log(xdt/xm))
  n <- length(xdt)
  al <- 1/xi
  lpl <- n*log(al)+n*al*log(xm)-(al+1)*sum(log(xdt))
  list(min=xm, alpha=1/xi, lPL=lpl)
}

ePL(danish)


## ------------------------------------------------------------------------------------------------------
rgpl <- function(mu, alpha, n, seed){
  set.seed(seed)
  y <- runif(n)
  x <- mu*((1-y)^(-1/alpha))
  return(x)
}


## ------------------------------------------------------------------------------------------------------
mu <- 1
al <- 1

# Simulació d'un vector de dades aleatòries
x <- rgpl(mu=mu, alpha=al, n=10000, seed=1)

# Estadístic de contrast
X <- sort(x)
n <- length(X)
F <- 1-(X/mu)^(-al)
E <- seq(1:n)/n
D <- ks.test(x=F, y=E)$statistic[[1]]
est.con <- sqrt(n)*D 
 
print(sprintf("est.con=%.3f, D=%.3f", est.con, D))


## ------------------------------------------------------------------------------------------------------
mu <- 1
al <- 3

# Simulació d'un vector de dades aleatòries
x <- rgpl(mu=mu, alpha=al, n=10000, seed=1)

# Estadístic de contrast
X <- sort(x)
n <- length(X)
F <- 1-(X/mu)^(-al) 
E <- seq(1:n)/n
D <- max(abs(E-F))
est.con <- sqrt(n)*D 
print(sprintf("est.con=%.8f, D=%.8f", est.con, D))

