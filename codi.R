###################################
#  ESTIMACIÓ D'UN QUANTIL ELEVAT  #
###################################

# Paquets
library(poweRlaw)

# Carreguem les dades
data("danish", package="evir")
data("nidd.thresh", package="evir")

# Construcció de la distribució Power-Law
dPL <- function(xdt, mu, alpha){
  return((alpha/mu)*((mu/xdt)^(alpha+1)))
}

pPL <- function(xdt, mu, alpha){
  return(1 - (xdt/mu)^(-alpha))
}

qPL <- function(xdt, mu, alpha){
  return(mu*((1-xdt)^(-1/alpha)))
}

rPL <- function(mu, alpha, n, seed){
  set.seed(seed)
  y <- runif(n)
  x <- qPL(y, mu, alpha)
  return(x)
}

# Construcció de la funció ePL (estimació per MLE)
ePL <- function(xdt){
  xm <- min(xdt)
  xi <- mean(log(xdt/xm))
  n <- length(xdt)
  al <- 1/xi
  lpl <- n*log(al)+n*al*log(xm)-(al+1)*sum(log(xdt))
  list(min=xm, alpha=1/xi, lPL=lpl)
}

# Avaluació ePL a nidd.thresh i a danish (validació ePL)
pars_nidd <- ePL(nidd.thresh)
pars_danish <- ePL(danish)

# Construcció de la funció per calcular l'estadistic de contrast de Kolmogrov-Smirnov
ks_statistic <- function(xdt, mu = min(xdt)){
  al <- ePL(xdt)$alpha
  X <- sort(xdt)
  n <- length(xdt)
  F <- pPL(X, mu, al)
  E <- seq(1:n)/n
  D <- max(abs(E-F))
  est.con <- sqrt(n)*D
  list(est.con = est.con, D = D)
}

k_al <- c(1.22, 1.36, 1.63)

nidd_ks <- ks_statistic(nidd.thresh)
print(sprintf("est.con=%.8f, D=%.8f", nidd_ks$est.con, nidd_ks$D))
print(nidd_ks$est.con>k_al)

danish_ks <- ks_statistic(danish)
print(sprintf("est.con=%.8f, D=%.8f", danish_ks$est.con, danish_ks$D))
print(danish_ks$est.con>k_al)

# Validació de la funció rPL
mus <- c(1, 4, 5.5, 10, 48.4)
als <- c(1, 3.5, 5, 10, 12.4)

mu_est <- NULL
al_est <- NULL
val_mu_est <- NULL
val_al_est <- NULL
x <- NULL
i <- 1

for(mu in mus){
  for(al in als){
    x[[i]] <- rPL(mu=mu, alpha=al, n=10000, seed=1)
    mu_est[i] <- round(ePL(x[[i]])$min,3)
    al_est[i] <- round(ePL(x[[i]])$alpha,3)
    val_mu_est[i] <- round(ePL(x[[i]]/mu)$min,3)
    val_al_est[i] <- round(ePL(x[[i]]/mu)$alpha,3)
    i <- i + 1
  }
}

tabla <- data.frame(expand.grid(als, mus)[2:1], mu_est, al_est, val_mu_est, val_al_est)
colnames(tabla)[1:2] <- c("mu", "al")
tabla

for(i in 1:length(x)){
  x_ks <- ks_statistic(x[[i]])
  print(sprintf("est.con=%.8f, D=%.8f", x_ks$est.con, x_ks$D))
  print(x_ks$est.con>k_al)
}

# Ja que hem obtingut tots els estadístics de contrast inferiors als diferents punts crítics podem dir que hem validat la nostra funció rPL.

# Selecció del llindar
selmin <- function(xd){
  nn <- length(xd)
  xd <- sort(xd)
  dista <- numeric(nn)
  for(k in 1:nn){
    xdt <- xd[k:nn]
    dista[k] <- ks_statistic(xdt)$D
  }
  omin <- which.min(dista)
  xmin <- xd[omin]
  KS <- dista[omin]
  list(n=nn, min=xmin, KS=KS, ordre=omin)
}

parset <- function(xdat){
  nn <- length(xdat)
  xdat <- sort(xdat)
  st <- selmin(xdat)
  xm <- st$min
  nb <- st$ordre-1
  na <- nn-nb
  al <- ePL(xdat[(nb+1):nn])$alpha
  list(n=nn, alpha=al, min=xm, ntail=na)
}

parset(danish)

plmag <- conpl$new(sort(danish))
xmin_PL <- estimate_xmin(plmag)$xmin
plmag$setXmin(xmin_PL)
estimate_pars(plmag)$pars - 1
xmin_PL


#############################################################
#  ANÀLISI DE LES DADES D'UNA EMPRESA QUE COTITZA AL NASDAQ #
#############################################################

# Empresa: MICROSOFT

# Paquets
library(quantmod)
library(evir)
library(ercv)

# Carreguem les dades
getSymbols("MSFT", from="2016-01-01", to="2020-12-31")

names(MSFT)
plot(MSFT$MSFT.Close, col="blue")

# Preus
preus <- MSFT$MSFT.Close
preus <- as.numeric(preus)
preus <- preus[!is.na(preus)]
length(preus)
identical(as.numeric(na.omit(MSFT$MSFT.Close)), preus)

ts.plot(preus, main="MICROSOFT prices", col="blue")
length(preus)/5

prices.ts <- ts(preus, frequency=252, start=c(2016, 1), end=c(2020, 250))
plot(prices.ts, main="MICROSOFT (2016-2021)", col="blue")

# Rendibilitats
returns <- diff(log(preus))

mu <- mean(returns)
sg <- sd(returns)

# Gràfic rendibilitats
ts.plot(returns, main="Returns", col="red")

# Càlcul rendiment anual
ma <- 250*mu

# Càlcul volatilitat
sa <- sqrt(250)*sg
sprintf("Rendibilitat anual: %.4f (%.2f %%) \n Volatilitat anual: %.4f (%.2f %%)", ma, ma*100, sa, sa*100)

# Avaluem el nombre d’observacions que estan més enllà de la mitjana i tres desviacions estàndard
ts.plot(returns, main="Returns", col="red")
abline(h=mu-3*sg, col="blue")
abline(h=mu+3*sg, col="blue")

paste("Hi ha", length(returns[returns<mu-3*sg | returns>mu+3*sg]),
      "observacions que están més enllá de la mitjana i 3 desviacions estàndard")

# Comparem amb els valors que esperaríem sota el model normal
normal <- rnorm(length(returns), mu, sg)
paste("Hi ha", length(normal[normal<mu-3*sg | normal>mu+3*sg]),
      "observacions que están més enllá de la mitjana i 3 desviacions estàndard")

# Ràtio de sharpe
sharpe <- ma/sa
sharpe

# Cues
pos.ret <- returns[returns>0]
neg.ret <- -returns[returns<0]
par(mfrow = c(1,2))
hist(pos.ret, border="blue", col="white")
hist(neg.ret, border="blue", col="white")
dev.off()
hist(returns, breaks="Scott", main="Histogram of returns",
     freq=F, xlim=c(-0.07, 0.07), border="blue", lwd=2)

curve(dnorm(x, mu, sg), add=TRUE, col="red", lwd=2)

boxplot(returns, main="Boxplot of returns", horizontal=TRUE, col="cyan", border="blue")

# DESCRIPTIVA CUES
# Cua positiva
np <- length(pos.ret)
np
mp <- mean(pos.ret)
1/mp
cvp <- sd(pos.ret)/mp
paste("El coeficient de variació per la cua positiva és" , round(cvp,3))

# Cua negativa
nn <- length(neg.ret)
nn
mn <- mean(neg.ret)
1/mn
cvn <- sd(neg.ret)/mp
paste("El coeficient de variació per la cua negativa és" , round(cvn,3))

# Cua positiva
alp <- sd(pos.ret)/mean(pos.ret)
alp
xip <- 1/alp
paste("El paràmetre xi per la cua positiva és" , round(xip,3))

s1 <- gpd(pos.ret, threshold=0)
s1$par.ests
s1$par.ses
ci.pos <- c(s1$par.ests[1]-2*s1$par.ses[1], s1$par.ests[1]+2*s1$par.ses[1])
paste0("El paràmetre estimat per evir és ", round(s1$par.ests[1],3), " amb un interval de confiança de [", round(ci.pos[1],3), ",", round(ci.pos[2],3), "]")

# Cua negativa
aln <- sd(neg.ret)/mean(neg.ret)
aln
xin <- 1/aln
paste("El paràmetre xi per la cua negativa és" , round(xin,3))

s2 <- gpd(neg.ret, threshold=0)
s2$par.ests
s2$par.ses
ci.neg <- c(s2$par.ests[1]-2*s2$par.ses[1], s2$par.ests[1]+2*s2$par.ses[1])
paste0("El paràmetre estimat per evir és ", round(s2$par.ests[1],3), " amb un interval de confiança de [", round(ci.neg[1],3), ",", round(ci.neg[2],3), "]")

# Trobar el llindar que dona el millor ajust de PL

parset(pos.ret)
parset(neg.ret)

# Cua positiva

# Estimació PL
fit1 <- parset(pos.ret); fit1
a1 <- fit1$alpha; a1
xm <- fit1$min; xm
n <- fit1$ntail
llik1 <- n*log(a1)+n*a1*log(xm)-(a1+1)*sum(log(pos.ret[pos.ret>xm])); llik1

# Estimació LPD
fit2 <- gpd(pos.ret, threshold = xm); fit2
fit2$par.ests
fit2$nllh.final
par <- as.numeric(fit2$par.ests); par
ah <- 1/par[1]; ah # alpha
dh <- par[2]/par[1]-xm; dh # sigma
llik2 <- -fit2$nllh.final; llik2

# Valors crítics
qchisq(0.95, df=1) # Xi^2 amb 1 gl 95% confiança
qchisq(0.99, df=1) # Xi^2 amb 1 gl 99% confiança

# Likelihood Ratio Test (LRT)
LRT <- 2*(llik2-llik1); LRT
# El test de raó de versembança diu que LRT=11.88
# Ja que 17.22842>3.84, el test és significatiu al 95% de confiança.
# Ja que 17.22842>6.63, el test és significatiu al 99% de confiança.
# Per tant, el LRT rebutja que les dades segueixin una distribució Power-Law.
# En conclusió, ens quedem amb la distribució LPD.

# Cua negativa

# Estimació PL
fit1 <- parset(neg.ret); fit1
a1 <- fit1$alpha; a1
xm <- fit1$min; xm
n <- fit1$ntail
llik1 <- n*log(a1)+n*a1*log(xm)-(a1+1)*sum(log(neg.ret[neg.ret>xm])); llik1

# Estimació LPD
fit2 <- gpd(neg.ret, threshold=xm); fit2
fit2$par.ests
fit2$nllh.final
par <- as.numeric(fit2$par.ests); par
ah <- 1/par[1]; ah # alpha
dh <- par[2]/par[1]-xm; dh # sigma
llik2 <- -fit2$nllh.final; llik2 # canvi signe

# Likelihood Ratio Test (LRT)
LRT <- 2*(llik2-llik1); LRT
# El test de raó de versembança diu que LRT=11.88
# Ja que 18.84034>3.84, el test és significatiu al 95% de confiança.
# Ja que 18.84034>6.63, el test és significatiu al 99% de confiança.
# Per tant, el LRT rebutja que les dades segueixin una distribució Power-Law.
# En conclusió, ens quedem amb la distribució LPD.

# Càlcul del llindar per la cua negativa
cvplot(neg.ret, main="Cua negativa")
Tm(neg.ret, m=10)
Tm(neg.ret, evi = 0)
thr <- thrselect(neg.ret, m=30, nsim=1000, conf.level=0.95)
thr$solution # accepto que les dades són GPD perquè el p-valor>=0.05
cievi(nextrem=thr$solution$nextremes, evi=thr$solution$evi,  nsim=1000)

fitpot(neg.ret, evi=abs(thr$solution$evi), nextremes=thr$solution$nextremes)

llindar_final <- parset(neg.ret)$min

fitpot1 <- fitpot(neg.ret, evi=abs(thr$solution$evi), nextremes=thr$solution$nextremes); fitpot1 # ercv
fitpot2 <- fitpot(neg.ret, evi=as.numeric(fit2$par.ests[[1]]), nextremes=fit2$n.exceed); fitpot2  # LPD

ccdfplot(neg.ret, pars=c(fitpot1, fitpot2), log="xy", main="MICROSOFT")

# Gràfics segons els paràmetres obtinguts amb threshold (fitpot1)
xi <- thr$solution$evi; xi
mu <- thr$solution$threshold; mu
beta <- fitpot1$coeff[[2]]; beta
n <- thr$solution$nextremes

X <- sort(neg.ret)
E <- seq(1:n)/(n+1) # distribució empírica
F <- pgpd(X, xi=xi, mu=mu, beta=beta) # distribució teòrica
Q <- qgpd(E, xi=xi, mu=mu, beta=beta)

# Funcions de distribució
plot(X, F, main = "Funcions de distribució", col = "red", lwd=2, type="l")
lines(X, E, col = "blue", lwd = 2, type = "l")
legend("bottomright", legend=c("Empírica", "Teòrica"), col=c("red", "blue"), lty=1, lwd=2)

# Les funcions de distribució estan solapades. Això significa que la distribució teòrica s'ajusta a la distribució empírica.

# PP-plot
plot(E, F, main="PP-plot Dow Jones", col="red", type="l", lwd=2)
abline(0, 1, col="green")

# QQ-plot
plot(Q, X, main="QQ-plot Dow Jones", col="red", lwd=2, type="l")
abline(0, 1, col="green")

# Gràfics segons els paràmetres obtinguts amb GPD (fitpot2)
xi <- as.numeric(fit2$par.ests[[1]])
mu <- fit2$threshold
beta <- as.numeric(fit2$par.ests[[2]])
n <- fit2$n.exceed

X <- sort(neg.ret[neg.ret>mu])
E <- seq(1:n)/(n+1) # distribució empírica
F <- pgpd(X, xi=xi, mu=mu, beta=beta) # distribució teòrica
Q <- qgpd(E, xi=xi, mu=mu, beta=beta)

# Funcions de distribució
plot(X, F, main = "Funcions de distribució", col = "red", lwd=2, type="l")
lines(X, E, col = "blue", lwd = 2, type = "l")
legend("bottomright", legend=c("Empírica", "Teòrica"), col=c("red", "blue"), lty=1, lwd=2)

# Les funcions de distribució estan solapades. Això significa que la distribució teòrica s'ajusta a la distribució empírica.

# PP-plot
plot(E, F, main="PP-plot Dow Jones", col="red", type="l", lwd=2)
abline(0, 1, col="green")

# QQ-plot
plot(Q, X, main="QQ-plot Dow Jones", col="red", lwd=2, type="l")
abline(0, 1, col="green")

# Decisió del llindar escollit


# Càlcul dels VARs
llindar_final <- 0.02820915 # mu (codi propi PL)

p_u <- length(neg.ret[neg.ret>llindar_final])/length(neg.ret)
mod <- gpd(neg.ret, threshold=llindar_final)
chi_u <- mod$par.ests[1]
beta_u <- mod$par.ests[2]
epsilon <- c(0.01, 0.001)

vars <- llindar_final + (beta_u/chi_u)*((p_u/epsilon)^chi_u-1)
vars
