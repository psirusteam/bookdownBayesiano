library(faraway)
attach(orings)
x <- temp
y <- damage

### Estimacion MV

n.sim <- 200
beta0 <- beta1 <- rep(0, n.sim)
for(i in 2:n.sim){
  eta <- beta0[i-1] + beta1[i-1]*x
  z <- eta + (1+exp(eta))^2/exp(eta)*(y-exp(eta)/(1+exp(eta)))
  sigma2 <- (1+exp(eta))^2/exp(eta)
  M <- lm(z ~ x, weight=1/(sigma2))
  beta0[i] <- M$coef[1]
  beta1[i] <- M$coef[2]
}
mean(beta0)
mean(beta1)

### Bayesiano con aproximación a la normal



### Bayesiano Gibbs MH en \verb'R'
MHlogit <- function(beta0, beta1, x, y, Nsim){
  
  # Creación de la distribución a posteriori.
  # Toma la forma de la verosimilitud pues las apriori son vagues
  post <- function(beta0, beta1, x, y) {
    lp <- beta0 + beta1 * x
    p <- exp(lp)/(1 + exp(lp))
    #modificación para evitar problemas numéricos
    sum(y*log(p)+(6-y)*log(1-p))
  }
  
  ind <-rep(0,Nsim)
  
  betas <- matrix(NA,ncol=2,nrow=Nsim)
  for(k in 1:Nsim){
    betas[k,1] <- beta0
    betas[k,2] <- beta1
    
    #Genera un valor candidato
    beta0.can <- rnorm(1,beta0,0.1)
    beta1.can <- rnorm(1,beta1,0.1)
    
    #Jumping distribution
    q1 <- dnorm(beta0,beta0.can,0.1)*dnorm(beta1,beta1.can,0.1)
    q2 <- dnorm(beta0.can,beta0,0.1)*dnorm(beta1.can,beta1,0.1)
    
    #A posteriori
    p1 <- post(beta0.can,beta1.can, x, y)
    p2 <- post(beta0,beta1, x, y)
    
    #Aceptación beta0
    #modificación para evitar problemas numéricos
    T.val <- min(1,(exp(p1-p2)*(q1/q2)))
    u <- runif(1)
    if (u<= T.val){
      beta0=beta0.can
      beta1=beta1.can
      ind[k]<-1
    }  }
  return(list(ind=ind,betas=betas))
}
beta0=10
beta1=0

res <- MHlogit(beta0,beta1,x,y,Nsim=50000)

mean(res$betas[-(1:20000),1])
mean(res$betas[-(1:20000),2])


sd(res$betas[-(1:20000),1])
sd(res$betas[-(1:20000),2])

#### Bayesiano JAGS

n <- length(y)
Logit.Model <- function(){
  for(i in 1 : n){
    y[i] ~ dbin(theta[i], 6)	
    theta[i] <- exp(beta0+beta1*x[i])/(1+exp(beta0+beta1*x[i]))
  }
  beta0 ~ dnorm(0, 0.0001)
  beta1 ~ dnorm(0, 0.0001)
}

Logit.Model.data <- list("y", "x", "n")
Logit.Model.param <- c("beta0", "beta1")
Logit.Model.inits <- function(){
  list("beta0"=c(0), "beta1"=c(0))
}

Logit.Model.fit <- jags(data=Logit.Model.data, inits=Logit.Model.inits, 
                        Logit.Model.param, n.iter=10000, n.burnin=1000, model.file=Logit.Model)

print(Logit.Model.fit)