library(ggplot2)

# ----------------------------
#     Capítulo 1
# ----------------------------

# Bernoulli

# Previa no informativa Jeffreys Beta(1/2,1/2)

p <- ggplot(data = data.frame(x = 0), mapping = aes(x = x)) +
	ggtitle("Densidad Beta(1/2,1/2)") + ylab("f(x)")
Beta.Jeffreys <- function(x) dbeta(x,0.5,0.5)
p + stat_function(fun = Beta.Jeffreys) + xlim(0.01,0.99)


# Primera figura del ejemplo electoral

p <- ggplot(data = data.frame(x = 0), mapping = aes(x = x)) + ylab("f(x)")
previa <- function(x) dbeta(x,350, 650)
posterior <- function(x) dbeta(x,6710, 6290)
p + layer(stat="function", fun = previa, mapping = aes(linetype="solid"))+
	layer(stat="function", fun = posterior, mapping = aes(linetype="dashed")) +
	xlim(0.3,0.6) +	
	theme(legend.position="none")
	
# Segunda figura del ejemplo electoral

p <- ggplot(data = data.frame(x = 0), mapping = aes(x = x)) +
		ylab("f(x)")
posterior2 <- function(x) dbeta(x,6360.5, 5640.5)
p + layer(stat="function", fun = Beta.Jeffreys, mapping = aes(linetype="solid"))+
	layer(stat="function", fun = posterior2, mapping = aes(linetype="dashed")) +
	xlim(0.3,0.6)	 +	
	theme(legend.position="none")
	
# Binomial

# Figura 1.4
p <- ggplot(data = data.frame(x = 0), mapping = aes(x = x)) + ylab("") +xlab(expression(theta))
previa <- function(x) dbeta(x, 2,5)
veros <- function(x) {x^8*(1-x)^2/beta(8,2)}
posterior <- function(x) dbeta(x,10, 7)
p + layer(stat="function", fun = previa, mapping = aes(linetype="solid"))+
	layer(stat="function", fun = veros, mapping = aes(linetype="dashed")) +
	layer(stat="function", fun = posterior, mapping = aes(linetype="dotdash")) +
	xlim(0,1)	 +
	scale_linetype_manual(name="",values=c("solid","dashed","dotdash"),
		label=c("Verosimilitud","A posterior","A priori")) +
		theme(legend.text = element_text(size=20))

# Figura 1.5, estimación bayesiana para diferentes valores de n
alfa <- 5; beta <- 5
n <- 3:100
theta <- n/(n+alfa+beta)*0.33 + (alfa+beta)/(n+alfa+beta)*0.5
qplot(3:100,theta)+ ylab("Estimación a posteriori")+ xlab("n")+
	ylim(0.3,0.5) +
	geom_hline(yintercept=0.33)+
	geom_point(shape=3)


# Figura 1.6
p <- ggplot(data = data.frame(x = 0), mapping = aes(x = x)) + ylab("") +xlab(expression(theta))
previa <- function(x) dbeta(x, 215,595)
posterior <- function(x) dbeta(x,2033, 5426)
p + layer(stat="function", fun = previa, mapping = aes(linetype="solid"))+
	layer(stat="function", fun = posterior, mapping = aes(linetype="dashed")) +
	xlim(0.2,0.35)	 +
	scale_linetype_manual(name="",values=c("solid","dashed"),
		label=c("A posterior","A priori")) +
		theme(legend.text = element_text(size=20))

# Figura 1.7
n <- 70; s<- 70*0.2
alp <-7; bet <- 38; n.ast <- 400

predictiva <- rep(NA,n.ast+1)
for(k in 0:n.ast)
{
predictiva[(k+1)] <-
choose(n.ast,k)*beta(k+s+alp,bet-k-s+n.ast+n)/beta(s+alp,bet-s+n)
}

sum(predictiva)

predic <- data.frame(time=c(1:(n.ast+1)), y=predictiva)
ggplot(predic, aes(time, y)) + geom_line() + 
  xlab(expression(tilde(y))) + ylab("Probabilidad predictiva")
  
# Figura 1.8

predict<-function(y,alfa,beta,s,n,k){
 choose(y-1,k-1)*exp(lbeta(alfa+k+s,beta+y-k+n-s)-lbeta(alfa+s,beta+n-s))
 }

 alfa <- beta <- 0.5
 s <- 152;n <- 29620;k <- 5
 fun <- rep(NA)
 for(y in 5:5000){
 fun[y-4] <- predict(y,alfa,beta,s,n,k)
 }
 sum(fun)
 
 func <- data.frame(time=5:5000, y=fun)
 ggplot(func, aes(time, y)) + geom_line() + 
  xlab(expression(tilde(y))) + ylab("Probabilidad predictiva")


# Figura 1.9,  Poisson
Trans <- c(22, 9, 9, 20, 10, 14, 11, 14, 11, 11, 19, 12, 8, 9, 16, 8, 13, 8, 14, 12, 14, 11, 14, 13, 11, 14, 13, 11, 7, 12 )

qplot(Trans, geom="histogram", binwidth=2, xlab="Accidentes")

# Figura 1.10
library(gridExtra)
p <- ggplot(data = data.frame(x = 0), mapping = aes(x = x)) + ylab("") +xlab(expression(theta))
previa_J <- function(x){
	x^{-1/2}/20
}
posterior_J <- function(x) dgamma(x,370.5, 30)
previa_I <- function(x) dgamma(x,38,9)
posterior_I <- function(x) dgamma(x,408, 39)
p1 <- p + layer(stat="function", fun = posterior_J, mapping = aes(linetype="dashed"))+
	layer(stat="function", fun = previa_J, mapping = aes(linetype="solid")) +
	xlim(0, 20)	 +
	scale_linetype_manual(name="",values=c("solid","dashed"),
		label=c("Posterior","Previa")) +
		theme(legend.text = element_text(size=14))

p2 <- p + layer(stat="function", fun = posterior_I, mapping = aes(linetype="dashed"))+
	layer(stat="function", fun = previa_I, mapping = aes(linetype="solid")) +
	xlim(0, 20)	 +
	scale_linetype_manual(name="",values=c("solid","dashed"),
		label=c("Posterior","Previa")) +
		theme(legend.text = element_text(size=14))
		
grid.arrange(p1, p2, ncol=2)

# Figura 1.11. Distribución predictiva posterior para el ejemplo de tránsitos
Trans <- c(22, 9, 9, 20, 10, 14, 11, 14, 11, 11, 19, 12, 8, 9, 16, 8, 13, 8, 14, 12, 14, 11, 14, 13, 11, 14, 13, 11, 7, 12 )
alfa <- 38; beta <- 9; n <- length(Trans)
pre.Transito <- function(s){
if(s>0){
val <- gamma(s)*((n+beta)/(n+beta+1))^(sum(Trans)+alfa)/(beta(s,sum(Trans)+alfa)*(1+n+beta)^s*prod(1:s))}
if(s==0){
val <- ((n+beta)/(n+beta+1))^(sum(Trans)+alfa)
}
val
}
pre.Transito.NoInf <- function(s){
if(s>0){
	val <- gamma(s)*(n/(n+1))^(sum(Trans)+0.5)/(beta(s,sum(Trans)+0.5)*prod(1:s)*(n+1)^s)}
if(s==0){
	val <- (n/(n+1))^(sum(Trans)+0.5)
}
val
}
s.max <- 40; s.val <- 0:s.max
pre.NoInf.val <- pre.Inf.val<- c()
for(i in 1:length(s.val)){
	pre.NoInf.val[i] <- pre.Transito.NoInf(s.val[i])
	pre.Inf.val[i] <- pre.Transito(s.val[i])
}
sum(pre.NoInf.val)
sum(pre.Inf.val)
func <- data.frame(time=0:s.max, y=pre.NoInf.val)
func2 <- data.frame(time=0:s.max, y=pre.Inf.val)
 ggplot(func, aes(time, y)) + geom_line() + 
  xlab(expression(tilde(y))) + ylab("Probabilidad predictiva") + 
   geom_line(data=func2, aes(x=time, y=pre.Inf.val), linetype="dashed") 
   
 
 
 
 
 
# Figura 1.12. Tiempo de sobrevivencia para los pacientes 
 pred_exp <- function(x){
 ((s/(s+x*n.mono))^n)*((x*n.mono/(s+x*n.mono))^n.mono)/(x*beta(n,n.mono))
 }

 alfa <- beta <- 0
 s <- 25998.5
 n <- 69
 n.mono <- 5

p <- ggplot(data = data.frame(x = 0), mapping = aes(x = x)) +
	ggtitle("Densidad predictiva para el tiempo promedio de sobreviviencia") + 
	xlab("Tiempo promedio") + ylab("")
p + stat_function(fun = pred_exp) + xlim(0.01,1500)


# Figura 1.13 Gráfica comparativa del estimador Bayesiano
mu <- 5
tau2 <- 0.01
n <- c(5,10,50,200)
sigma2 <- 1
y.bar <- 2
mu.n <- ((n/sigma2)*y.bar + mu/tau2)/(n/sigma2+1/tau2)
tau2.n <- (n/sigma2+1/tau2)^-1
previa <- function(x){
	dnorm(x,mu,sqrt(tau2))
}
verosi1 <- function(x){
	(2*pi*sigma2)^(-n[1]/2)*exp(-(n[1]*x^2-2*n[1]*y.bar*x)/(2*sigma2))
}
veros1 <- function(x){
	verosi1(x)/integrate(verosi1,-20,20)$value
}
posterior1 <- function(x){
	dnorm(x,mu.n[1],sqrt(tau2.n[1]))
}
verosi2 <- function(x){
	(2*pi*sigma2)^(-n[2]/2)*exp(-(n[2]*x^2-2*n[2]*y.bar*x)/(2*sigma2))
}
veros2 <- function(x){
	verosi2(x)/integrate(verosi2,-20,20)$value
}
posterior2 <- function(x){
	dnorm(x,mu.n[2],sqrt(tau2.n[2]))
}
verosi3 <- function(x){
	(2*pi*sigma2)^(-n[3]/2)*exp(-(n[3]*x^2-2*n[3]*y.bar*x)/(2*sigma2))
}
veros3 <- function(x){
	verosi3(x)/integrate(verosi3,-20,20)$value
}
posterior3 <- function(x){
	dnorm(x,mu.n[3],sqrt(tau2.n[3]))
}
verosi4 <- function(x){
	(2*pi*sigma2)^(-n[4]/2)*exp(-(n[4]*x^2-2*n[4]*y.bar*x)/(2*sigma2))
}
veros4 <- function(x){
	verosi4(x)/integrate(verosi4,-20,20)$value
}
posterior4 <- function(x){
	dnorm(x,mu.n[4],sqrt(tau2.n[4]))
}


p <- ggplot(data = data.frame(x = 0), mapping = aes(x = x)) + ylab("") +xlab(expression(theta))
p1 <- p + layer(stat="function", fun = previa, mapping = aes(linetype="dashed"))+
	layer(stat="function", fun = posterior1, mapping = aes(linetype="solid")) +
		layer(stat="function", fun = veros1, mapping = aes(linetype="dotdash"))+
	xlim(1, 6)	  +ggtitle("n=5")+ theme(legend.position="none")

p2 <- p + layer(stat="function", fun = previa, mapping = aes(linetype="dashed"))+
	layer(stat="function", fun = posterior2, mapping = aes(linetype="solid")) +
		layer(stat="function", fun = veros2, mapping = aes(linetype="dotdash"))+
	xlim(1, 6)	 +
	scale_linetype_manual(name="",values=c("solid","dashed","dotdash"),
		label=c("Previa","Verosimilitud","Posterior")) + ggtitle("n=10")+
		theme(legend.text = element_text(size=10))

p3 <- p + layer(stat="function", fun = previa, mapping = aes(linetype="dashed"))+
	layer(stat="function", fun = posterior3, mapping = aes(linetype="solid")) +
		layer(stat="function", fun = veros3, mapping = aes(linetype="dotdash"))+ xlim(1, 6)	  + ggtitle("n=50") + theme(legend.position="none")

p4 <- p + layer(stat="function", fun = previa, mapping = aes(linetype="dashed"))+
	layer(stat="function", fun = posterior4, mapping = aes(linetype="solid")) +
		layer(stat="function", fun = veros4, mapping = aes(linetype="dotdash"))+
	xlim(1, 6)	 +
	scale_linetype_manual(name="",values=c("solid","dashed","dotdash"),
		label=c("Previa","Verosimilitud","Posterior")) +ggtitle("n=200")+
		theme(legend.text = element_text(size=10))
		
library(gridExtra)
grid.arrange(p1, p2, p3, p4, ncol=2)

# Figura 1.14. Distribución previas gamma inversa para sigma^2
sigma_02 <- c(1,10)
n0 <- c(50,100)
previa1 <- function(x){
	a <- n0[1]/2
	b <- n0[1]*sigma_02[1]/2
	densigamma(x,alpha=a,beta=b)
}

previa2 <- function(x){
	a <- n0[2]/2
	b <- n0[2]*sigma_02[1]/2
	densigamma(x,alpha=a,beta=b)
}
previa3 <- function(x){
	a <- n0[1]/2
	b <- n0[1]*sigma_02[2]/2
	densigamma(x,alpha=a,beta=b)
}

previa4 <- function(x){
	a <- n0[2]/2
	b <- n0[2]*sigma_02[2]/2
	densigamma(x,alpha=a,beta=b)
}


p1 <- ggplot(data=data.frame(x = c(0.3, 2)), aes(x)) + ylab("")+ xlab("")+
  stat_function(fun = previa1, geom = "line",aes(linetype="n = 50")) +
  stat_function(fun = previa2, geom = "line", aes(linetype="n = 100"))+
  scale_linetype_manual("", values = c("solid","dashed")) + ggtitle(expression(paste(sigma[0]^2,"=1")))

p2 <- ggplot(data=data.frame(x = c(0.5,25)), aes(x)) + xlab("")+ylab("")+
  stat_function(fun = previa3, geom = "line",aes(linetype="n = 50")) +ylab("")+
  stat_function(fun = previa4, geom = "line", aes(linetype="n = 100"))+
  scale_linetype_manual("", values = c("solid","dashed")) + ggtitle(expression(paste(sigma[0]^2,"=10")))

library(gridExtra)
grid.arrange(p1, p2, ncol=2)


# Figura 1.15 Gráfica comparativa de sigma2
library(pscl)
n0 <- 20
sigma2_0 <- 10
n <- c(5,20,50,100)
sigma2_C <- 50

previa <- function(x){
	a <- n0/2
	b <- n0*sigma2_0/2
	densigamma(x, a, b)
}
verosi1 <- function(x){
	x^(-n[1]/2)*exp(-n[1]*sigma2_C/(2*x))
}
veros1 <- function(x){
	verosi1(x)/integrate(verosi1,0,500)$value
}
posterior1 <- function(x){
	densigamma(x,(n[1]+n0)/2, (n0*sigma2_0+n[1]*sigma2_C)/2)
}

verosi2 <- function(x){
	x^(-n[2]/2)*exp(-n[2]*sigma2_C/(2*x))
}
veros2 <- function(x){
	verosi2(x)/integrate(verosi2,0,500)$value
}
posterior2 <- function(x){
	densigamma(x,(n[2]+n0)/2, (n0*sigma2_0+n[2]*sigma2_C)/2)
}

verosi3 <- function(x){
	x^(-n[3]/2)*exp(-n[3]*sigma2_C/(2*x))
}
veros3 <- function(x){
	verosi3(x)/integrate(verosi3,0,500)$value
}
posterior3 <- function(x){
	densigamma(x,(n[3]+n0)/2, (n0*sigma2_0+n[3]*sigma2_C)/2)
}

verosi4 <- function(x){
	x^(-n[4]/2)*exp(-n[4]*sigma2_C/(2*x))
}
veros4 <- function(x){
	verosi4(x)/integrate(verosi4,0,500)$value
}
posterior4 <- function(x){
	densigamma(x,(n[4]+n0)/2, (n0*sigma2_0+n[4]*sigma2_C)/2)
}


p1 <- ggplot(data=data.frame(x = c(0.01, 120)), aes(x)) + ylab("")+ xlab("")+
  stat_function(fun = previa, geom = "line",aes(linetype="Previa")) +
  stat_function(fun = veros1, geom = "line", aes(linetype="Verosimilitud"))+
  stat_function(fun = posterior1, geom = "line", aes(linetype="Posterior"))+
  scale_linetype_manual("", values = c("solid","dashed","dotted")) + ggtitle("n=5")

p2 <- ggplot(data=data.frame(x = c(0.01, 120)), aes(x)) + ylab("")+ xlab("")+
  stat_function(fun = previa, geom = "line",aes(linetype="Previa")) +
  stat_function(fun = veros2, geom = "line", aes(linetype="Verosimilitud"))+
  stat_function(fun = posterior2, geom = "line", aes(linetype="Posterior"))+
  scale_linetype_manual("", values = c("solid","dashed","dotted")) + ggtitle("n=20")

p3 <- ggplot(data=data.frame(x = c(0.01, 120)), aes(x)) + ylab("")+ xlab("")+
  stat_function(fun = previa, geom = "line",aes(linetype="Previa")) +
  stat_function(fun = veros3, geom = "line", aes(linetype="Verosimilitud"))+
  stat_function(fun = posterior3, geom = "line", aes(linetype="Posterior"))+
  scale_linetype_manual("", values = c("solid","dashed","dotted")) + ggtitle("n=50")

p4 <- ggplot(data=data.frame(x = c(0.01, 120)), aes(x)) + ylab("")+ xlab("")+
  stat_function(fun = previa, geom = "line",aes(linetype="Previa")) +
  stat_function(fun = veros4, geom = "line", aes(linetype="Verosimilitud"))+
  stat_function(fun = posterior4, geom = "line", aes(linetype="Posterior"))+
  scale_linetype_manual("", values = c("solid","dashed","dotted")) + ggtitle("n=100")
		
library(gridExtra)
grid.arrange(p1, p2, p3, p4, ncol=2)

# Figura 1.16 Previa de Jeffreys y densidad posterior de sigma2
library(pscl)
n <- 50
sigma2_C <- 10

previa1 <- function(x){
	1/x	
}
previa <- function(x){
	previa1(x)/integrate(previa1,0.01,50)$value
}

verosi1 <- function(x){
	x^(-n/2)*exp(-n*sigma2_C/(2*x))
}
veros1 <- function(x){
	verosi1(x)/integrate(verosi1,0,50)$value
}
posterior1 <- function(x){
	densigamma(x, n/2, n*sigma2_C/2)
}

p <- ggplot(data=data.frame(x = c(0.01, 30)), aes(x)) + ylab("")+ xlab("")+ylim(0,0.22)+
  stat_function(fun = veros1, geom = "line", aes(linetype="Verosimilitud"))+
  stat_function(fun = previa, geom = "line",aes(linetype="Previa")) +
  stat_function(fun = posterior1, geom = "line", aes(linetype="Posterior"))+
  scale_linetype_manual("", values = c("solid","dashed","dotted")) 

# Figura 1.16 Previa de Inv-Gamma(0.001,0.001)  para sigma2
library("MCMCpack")
pre <- function(x){
	dinvgamma(x,0.001,0.001)
}

p <- ggplot(data=data.frame(x = c(0.01, 15)), aes(x)) + ylab("")+ xlab("")+
  stat_function(fun = pre, geom = "line") + 	ggtitle("Densidad Gamma-Inversa(0.001,0.001)")

# Figura 1.17 Predictiva de una nueva Y con theta conocida, sigma2 desconocida
y <- c()
theta <- 10; 
n_0 <- 20; sigma2_0 <- 5
n <- 10; sigma2_c <- 7
for(i in 1:10000){
	alfa <- (n_0+n)/2; beta <- (n_0*sigma2_0+n*sigma2_c)/2
	sigma2 <- rinvgamma(1, alfa, beta)
	y[i] <- rnorm(1, theta, sqrt(sigma2))
}
tn <- function(x){
	v0 <- n_0*sigma2_0+n*sigma2_c
	(pi*v0)^(-0.5)*gamma((n_0+n+1)/2)*(1+(x-theta)^2/v0)^(-(n_0+n+1)/2)/gamma((n_0+n)/2)
}
ggplot()+geom_histogram(data=data.frame(y), aes(x=y,y=..density..),binwidth=0.5,fill="gray",color="black")+
    stat_function(fun=tn,size=1,aes(x=y))  +theme_grey(base_size=25)

# Figura 2.1. Predictiva de una nueva Y con theta y sigma2 independiente en la priori informativa
y.tilde <- c()
nsim <- 10000
theta.pos <- rep(NA,nsim)
sigma2.pos <- rep(NA,nsim)
theta.pos[1] <- 0

mu <- 10; tau2 <- 4
n_0 <- 20; sigma2_0 <- 5
n <- 10; y.bar <- 12; S2 <- 7
a.n <- (n_0+n)/2
b.n <- (n_0*sigma2_0 + (n-1)*S2 + n*(mean(y)-theta.pos[1]))/2
sigma2.pos[1] <- rinvgamma(1, a.n, b.n)
for(i in 2:10000){
	tau2.n <- 1/(n/sigma2.pos[i-1]+1/tau2)
	mu.n <- tau2.n*(n*y.bar/sigma2.pos[i-1]+mu/tau2)
	theta.pos[i] <- rnorm(1, mean=mu.n, sd=sqrt(tau2.n))
	alfa <- (n_0+n)/2; beta <- (n_0*sigma2_0 + (n-1)*S2 + n*(y.bar - theta.pos[i]))/2
	sigma2.pos[i] <- rinvgamma(1, alfa, beta)
	y.tilde[i] <- rnorm(1, theta.pos[i], sqrt(sigma2.pos[i]))
}

ggplot()+geom_histogram(data=data.frame(y.tilde), aes(x=y.tilde,y=..density..),binwidth=0.5,fill="gray",color="black") + theme_grey(base_size=25)



tn <- function(x){
	v0 <- n_0*sigma2_0+n*sigma2_c
	(pi*v0)^(-0.5)*gamma((n_0+n+1)/2)*(1+(x-theta)^2/v0)^(-(n_0+n+1)/2)/gamma((n_0+n)/2)
}
ggplot()+geom_histogram(data=data.frame(y.tilde), aes(x=y.tilde,y=..density..),binwidth=0.5,fill="gray",color="black")+
    stat_function(fun=tn,size=1,aes(x=y))  +theme_grey(base_size=25)

# Figura 2.2. Predictiva de una nueva Y con apriori no informativa
nsim <- 10000
y.tilde <- theta.pos <- c()
n <- 20
y <- rnorm(n,12,3)
y.bar <- mean(y); S2 <- var(y)
a.n <- (n-1)/2; b.n <- (n-1)*S2/2
sigma2.pos <- rinvgamma(nsim,a.n,b.n)
for(i in 1:10000){
	theta.pos[i] <- rnorm(1, y.bar, sqrt(sigma2.pos[i]/n))
	y.tilde[i] <- rnorm(1, theta.pos[i], sqrt(sigma2.pos[i]))
}
tn <- function(x){
	v2 <- (n+1)*S2/n
	(pi*(n-1)*v2)^(-0.5)*gamma(n/2)*(1+(x-y.bar)^2/((n-1)*v2))^(-n/2)/gamma((n-1)/2)
}

ggplot() + geom_histogram(data=data.frame(y.tilde), aes(x=y.tilde,y=..density..),binwidth=0.5,fill="gray",color="black") +     stat_function(fun=tn,size=1,aes(x=y.tilde)) + theme_grey(base_size=20)


# Capitulo 6, GLM

# Graficar la función de vínculo logístico
p <- ggplot(data = data.frame(x = 0), mapping = aes(x = x)) +
ggtitle("Función de vínculo logística") + ylab(expression(paste("g(",theta,")")))+ xlab(expression(theta))
vinc_logit <- function(x) log(x/(1 - x))
p + stat_function(fun = vinc_logit) + xlim(0.001,0.999)

# Graficar la función de vínculo probit
p <- ggplot(data = data.frame(x = 0), mapping = aes(x = x)) +
ggtitle("Función de vínculo probit (línea negra) y el vínculo logit (línea roja)") + ylab(expression(paste("g(",theta,")")))+ xlab(expression(theta))

vinc_probit <- function(x) qnorm(x)
p + stat_function(fun = vinc_probit) + xlim(0.001,0.999) + stat_function(fun = vinc_logit, colour="red")

# Graficar la función de vínculo log-log-complemntario para regresión Beta
p <- ggplot(data = data.frame(x = 0), mapping = aes(x = x)) +
ggtitle("Función de vínculo log-log complementario") + ylab(expression(paste("g(",theta,")")))+ xlab(expression(theta))

log_log <- function(x){
	log(-log(1-x))
}
p + stat_function(fun = log_log) + xlim(0.001,0.999) 


