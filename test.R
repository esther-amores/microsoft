install.packages("poweRlaw")
library(poweRlaw)

data("danish", package="evir")
danish


###1###
#Estimar parametres power-law de dades per mázima verzemblança
ePL <- function(xdt){
  xm <- min(xdt)
  xi <- mean(log(xdt/xm))
  n <- length(xdt)
  al <- 1/xi
  lpl <- n*log(al)+n*al*log(xm)-(al+1)*sum(log(xdt))
  list(min=xm, alpha=1/xi, lPL=lpl)
}

ePL(danish)


###2###
#Funció generadora de valors d'una Power-law amb parametres mu i alpha
rgpl <- function(mu, alpha, n, seed){
  set.seed(seed)
  y <- runif(n)
  x <- mu*((1-y)^(-1/alpha))
  return(x)
}


###3###
#Generació de valors d'una power-law amb parametres determinats
mu <- 1
alpha <- 1
x <- rgpl(mu=mu, alpha=alpha, n=1000, seed=1)

#Calcul de l'estadistic de contrast D
X <- sort(x)
n <- length(X)
F <- 1 - (X/mu)^(-alpha) 
E <- seq(1:n)/n
D <- max(abs(E-F))


sqrt(n)*D
#0.77<1.22 No tenim evidencies per rebutjar la hiposesis nula. Per tant, 
#les dades generades segueixen una distribució Power-law amb parametres mu=1, alpha=1


###4###


data = danish
mus = unique(data)
dat = numeric(length(mus))
z = sort(data)

for(i in 1:length(mus)){
  mu = mus[i] # escull la mu candidata
  z1 = z[z >= mu] # trunca les dades per sota d’aquest valor mu
  n = length(z1)
  a = ((1/n)* sum(log(z1/mu)))^(-1) # estima l’alpha mitjancant l’EMV
  cx = 1-(n:1)/n # construeix la CDF empirica
  cf = 1-((mu/z1)^(a)) # construeix la CDF teorica
  dat[i] = max(abs(cf-cx)) # calcula l’estadistic KS
}

D  = min(dat[dat >0], na.rm = TRUE) # troba el valor mes petit de D
mu = mus[which(dat==D)] # troba el corresponent valor mu
z  = data[data >= mu]
z  = sort(z)
n  = length(z)
alpha = ((1/n)*sum(log(z/mu)))^(-1) # obte la corresponent estimacio d’alpha
mu
alpha


plmag <- conpl$new(sort(danish))
xmin <- estimate_xmin(plmag)$xmin
plmag$setXmin(xmin)
estimate_pars(plmag)$pars - 1


