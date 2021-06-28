


# Modelos multiparamétricos

En este capítulo, discutimos situaciones donde se requiere estimar simultáneamente más de un parámetro; es decir, los datos que enfrentamos se ajustan a una distribución de probabilidad que involucra múltiples parámetros de forma simultanea. Específicamente, se estudiarán las siguientes distribuciones

* La distribución normal univariada que tiene dos parámetros: la media $\theta$ y la varianza $\sigma^2$.
* La distribución normal mutivariada con vector de medias $\btheta$ y matriz de varianzas-covarianzas $\bSigma$.
* La distribución multinomial cuyo parámetro es un vector de probabilidades $\btheta$.

En el contexto de la estimación bayesiana, es necesario hallar la distribución posterior conjunta de estos parámetros, y encontrar la estimación por alguna de las siguientes dos formas: (1) hallar teóricamente la esperanza de la distribución posterior conjunta, o (2) simular valores de la distribución posterior conjunta, de donde se puede obtener la estimación puntual y por intervalo.

Para abordar esta sección se supone que lector tiene una proficiencia media en términos de simulación estadística, cadenas de markov y métodos de Montecarlo. El apéndice de este libro contiene un resumen exhaustivo de los principales métodos de simulación de distribuciones de probabilidad, los cuales son la base fundamental para realizar una apropiada inferencia bayesiana.

## Modelo Normal con media y varianza desconocida

Supongamos que se dispone de realizaciones de un conjunto de variables independientes e idénticamente distribuidas $Y_1,\cdots,Y_n\sim N(\theta,\sigma^2)$. Cuando se desconoce tanto la media como la varianza de la distribución es necesario plantear diversos enfoques y situarse en el más conveniente, según el contexto del problema. En términos de la asignación de las distribuciones previas para $\theta$ y $\sigma^2$ es posible:

* Suponer que la distribución previa $p(\theta)$ es independiente de la distribución previa $p(\sigma^2)$ y que ambas distribuciones son informativas. 
* Suponer que la distribución previa $p(\theta)$ es independiente de la distribución previa $p(\sigma^2)$ y que ambas distribuciones son no informativas.
* Suponer que la distribución previa para $\theta$ depende de $\sigma^2$ y escribirla como $p(\theta \mid \sigma^2)$, mientras que la distribución previa de $\sigma^2$ no depende de $\theta$ y se puede escribir como $p(\sigma^2)$.

A continuación, analizamos cada uno de estos planteamientos, y desarrollamos los resultados necesarios para la estimación de $\theta$ y $\sigma^2$.

### Parámetros independientes

El primer enfoque que consideraremos para el análisis de los parámetros de interés $\theta$ y $\sigma^2$ en una distribución normal univariada es aquel que supone que las distribuciones previas de cada uno de ellos son independientes, pero al mismo tiempo informativas. @Gelman03 afirman que este supuesto de independencia es atractivo en problemas para los cuales la información previa para $\theta$ no toma la forma de un número fijo de observaciones con varianza $\sigma^2$. Adicionalmente, este supuesto de independencia es coherente con el hecho de que en la teoría clásica de estimación los estimadores insesgados de varianza mínima de $\theta$ y $\sigma^2$ son independientes [@Zhang, Sección 2.4]

En este orden de ideas, y siguiendo la argumentación del capítulo anterior, la distribución previa para el parámetro $\theta$ será

\begin{equation*}
\theta \sim Normal(\mu,\tau^2)
\end{equation*}

Y la distribución previa para el parámetro $\sigma^2$ será
\begin{equation*}
\sigma^2 \sim Inversa-Gamma(n_0/2,n_0\sigma^2_0/2)
\end{equation*}

Asumiendo independencia previa, la distribución previa conjunta estará dada por

\begin{equation}
p(\theta,\sigma^2)\propto (\sigma^2)^{-n_0/2-1}\exp\left\{-\dfrac{n_0\sigma^2_0}{2\sigma^2}\right\}
\exp\left\{-\frac{1}{2\tau^2}(\theta-\mu)^2\right\}
\end{equation}

Una vez que se conoce la forma estructural de la distribución previa conjunta, es posible establecer la distribución posterior conjunta puesto que la verosimilitud de los datos, $p(\mathbf{Y} \mid \theta,\sigma^2)$, está dada por la expresión \@ref(eq:veronormal) y

\begin{equation*}
p(\theta,\sigma^2 \mid \mathbf{Y})\propto p(\mathbf{Y} \mid \theta,\sigma^2)p(\theta,\sigma^2)
\end{equation*}

\BeginKnitrBlock{proposition}
<span class="proposition" id="prp:unnamed-chunk-1"><strong>(\#prp:unnamed-chunk-1) </strong></span>La distribución posterior conjunta de los parámetros de interés está dada por

\begin{align}
p(\theta,\sigma^2 \mid \mathbf{Y})&\propto (\sigma^2)^{-(n+n_0)/2-1} \notag \\
&\times
\exp\left\{-\frac{1}{2\sigma^2}\left[n_0\sigma^2_0+(n-1)S^2+n(\bar{y}-\theta)^2\right]-\frac{1}{2\tau^2}(\theta-\mu)^2\right\}
\end{align}
\EndKnitrBlock{proposition}
<br>
 
\BeginKnitrBlock{proof}
\iffalse{} <span class="proof"><em>Prueba. </em></span>  \fi{}Tenemos que
\begin{align*}
p(\theta,\sigma^2 \mid \mathbf{Y})&\propto p(\mathbf{Y} \mid \theta,\sigma^2)p(\theta,\sigma^2)\\
&\propto(\sigma^2)^{-n/2}\exp\left\{-\frac{1}{2\sigma^2}\sum_{i=1}^n(y_i-\theta)^2\right\}(\sigma^2)^{-n_0/2-1}\exp\left\{-\dfrac{n_0\sigma^2_0}{2\sigma^2}\right\}
\exp\left\{-\frac{1}{2\tau^2}(\theta-\mu)^2\right\} \\
&=(\sigma^2)^{-(n+n_0)/2-1}\exp\left\{-\frac{1}{2\sigma^2}\left[n_0\sigma^2_0+\sum_{i=1}^n(y_i-\theta)^2\right]-\frac{1}{2\tau^2}(\theta-\mu)^2\right\}\\
&\propto (\sigma^2)^{-(n+n_0)/2-1}\exp\left\{-\frac{1}{2\sigma^2}\left[n_0\sigma^2_0+(n-1)S^2+n(\bar{y}-\theta)^2\right]-\frac{1}{2\tau^2}(\theta-\mu)^2\right\}
\end{align*}
donde la última expresión se obtiene al sumar y restar $\bar{y}$ dentro de $(y_i-\theta)^2$.
\EndKnitrBlock{proof}
<br>

Nótese que la distribución posterior conjunta no tiene una forma estructural conocida y por lo tanto no es posible realizar el método de integración analítica para obtener una constante de integración [@Migon]. Sin embargo, sí es posible obtener las distribuciones condicionales posteriores de $\theta$ y de $\sigma^2$, notando que

\begin{align*}
p(\theta \mid \sigma^2,\mathbf{Y})\propto p(\theta,\underbrace{\sigma^2}_{fijo} \mid \mathbf{Y})
\ \ \ \ \ \ \ \ \ \text{y} \ \ \ \ \ \ \ \ \ \
p(\sigma^2 \mid \theta,\mathbf{Y})\propto p(\underbrace{\theta}_{fijo},\sigma^2 \mid \mathbf{Y})
\end{align*}

Es decir, para encontrar la distribución posterior condicional de $\theta | \sigma^2$, se utiliza la distribución posterior conjunta y los términos que no dependan de $\theta$ se incorporan en la constante de proporcionalidad. El mismo razonamiento se aplica para encontrar la distribución posterior condicional de $\sigma^2 | \theta$.

\BeginKnitrBlock{proposition}
<span class="proposition" id="prp:unnamed-chunk-3"><strong>(\#prp:unnamed-chunk-3) </strong></span>La distribución posterior condicional de $\theta$ es

\begin{equation}
(\#eq:PostThetaGibbs)
\theta  \mid  \sigma^2,\mathbf{Y} \sim Normal(\mu_n,\tau_n^2)
\end{equation}
  
En donde las expresiones para $\mu_n$ y $\tau_n^2$ están dadas por \@ref(eq:TauSigman). Por otro lado, la distribución posterior condicional de $\sigma^2$ es 

\begin{equation}
(\#eq:PostSigma2Gibbs)
\sigma^2  \mid  \theta,\mathbf{Y} \sim Inversa-Gamma\left(\dfrac{n_0+n}{2},\dfrac{v_0}{2}\right)
\end{equation}
  
con $v_0=n_0\sigma^2_0+(n-1)S^2+n(\bar{y}-\theta)^2$.
\EndKnitrBlock{proposition}
<br>

\BeginKnitrBlock{proof}
\iffalse{} <span class="proof"><em>Prueba. </em></span>  \fi{}Acudiendo a la distribución posterior conjunta, e incorporando los términos que no dependen de $\theta$ en la constante de proporcionalidad, se tiene que

\begin{align*}
p(\theta \mid \sigma^2,\mathbf{Y})&\propto \exp\left\{-\frac{n}{2\sigma^2}(\bar{y}-\theta)^2-\frac{1}{2\tau^2}(\theta-\mu)^2\right\}
\end{align*}

Completando los cuadrados se encuentra una expresión idéntica a la función de distribución de una variable aleatoria con distribución $Normal(\mu_n, \tau^2_n)$. Similarmente, después de un desarrollo algebraico breve, se tiene la distribución posterior condicional de $\sigma^2$.
\EndKnitrBlock{proof}
<br>

Una vez encontradas las distribuciones posteriores condicionales de $\theta$ y $\sigma^2$, es posible obtener la estimación puntual de estos parámetros usando métodos de Montecarlo, específicamente el muestreo de Gibbs que, en el contexto de este capítulo, se resume en los siguientes pasos:

1. Fijar un valor inicial para $\theta$, denotado por $\theta_{(1)}$.
2. Simular un valor de la distribución de $\sigma^2|\theta,\mathbf{Y}$ en \@ref(eq:PostSigma2Gibbs). Nótese que el parámetro $v_0$, que depende de $\theta$, debe ser reemplazado por $\theta_{(1)}$ del paso anterior. Este valor simulado se denotará como $\sigma^2_{(1)}$.
3. Simular un valor de la distribución de $\theta|\sigma^2,\mathbf{Y}$ en \@ref(eq:PostThetaGibbs). Nótese     que en $\mu_n$ y $\tau^2_n$ se debe reemplazar $\sigma^2$ por $\sigma^2_{(1)}$. Este valor simulado se denotará como $\theta_{(2)}$.
4. Repetir los pasos (2) y (3) hasta completar un número de iteraciones suficientes para alcanzar la convergencia en ambos parámetros.

Después de ejecutar el muestreo de Gibbs, se eliminan los primeros valores simulados para descartar la influencia del valor inicial; además, posiblemente se deba efectuar la fase de *thinning* para eliminar las correlaciones que pueden estar presentes en las cadenas generadas. Posteriormente, se obtendrán los valores finales simulados de las distribuciones de $\theta$ y $\sigma^2$, con los cuales se podrá calcular las estimaciones respectivas, además de estimar los intervalos de credibilidad resultantes con los percentiles muestrales de los valores simulados. 

En cuanto a la distribución predictiva para una nueva observación $\tilde{y}$, esta está dada por la siguiente expresión
\begin{equation*}
p(\tilde{y}\mid\mathbf{Y})=\int_0^\infty\int_{-\infty}^\infty p(\tilde{y}\mid\theta,\ \sigma^2)p(\theta,\ \sigma^2\mid\mathbf{Y})d\theta\ d\sigma^2
\end{equation*}

Hallar esta distribución de forma exacta no es fácil, y podemos optar por conocer el comportamiento probabilístico de $\tilde{y}$ por medio de la simulación estadística. Tal como se explicó en el capítulo anterior, se debe simular en primer lugar valores de $\theta$ y de $\sigma^2$ de la distribución posterior $p(\theta,\ \sigma^2\mid\mathbf{Y})$ - usando el muestreo de Gibbs - y posteriormente, valores de $\tilde{y}$ de la distribución $p(\tilde{y}\mid\theta,\ \sigma^2)$. 

Por otro lado, si se quiere conocer el comportamiento de una nueva muestra aleatoria $Y_1^{*},\cdots,Y_{n^*}^{*}$, es posible hacerlo por medio de la distribución predictiva de la media $\bar{Y}^*$, simulando en primer lugar valores de $\theta$ y de $\sigma^2$ de la distribución posterior $p(\theta,\ \sigma^2\mid\mathbf{Y})$ - usando el muestreo de Gibbs - y posteriormente valores de $\bar{Y}^*$ de la distribución $N(\theta,\frac{\sigma^2}{n^*})$.

\BeginKnitrBlock{example}
<span class="example" id="exm:EjeRenal"><strong>(\#exm:EjeRenal) </strong></span>@Efronims consideró un conjunto de datos que muestran la función renal de 157 individuos que se sometieron a una prueba médica exhaustiva en un hospital. Los resultados de la prueba renal están en un intervalo de -6 puntos a 4 puntos. Entre más alto sea el resultado, se concluye que el riñón del individuo es más sano. Nótese que estas pruebas son importantes para predecir el comportamiento de un riñón donado a un paciente con problemas renales. 

En principio, es de interés para el investigador conocer la media y la dispersión de estos datos, para poder analizar a fondo la situación de los pacientes que esperan un transplante. Dado que se trata de una primera aproximación, se prefiere utilizar distribuciones previas no informativas para los parámetros de la media y varianza. Lo anterior se logra en `STAN` definiendo las distribuciones previas de `theta ~ normal(0, 10000)` y de `sigma2 ~ inv_gamma(0.001, 0.001)`. De esta forma, la distribución previa de $\theta$ está centrada en cero, pero con una varianza muy grande al igual que la distribución de la varianza, los cuales representan distribuciones previas no informativas.

El siguiente código en `STAN` muestra cómo se lleva a cabo la inferencia.
\EndKnitrBlock{example}


```r
NormalMeanVar <- '
data {
  int<lower=0> n;
  real y[n];
}
parameters {
  real sigma;
  real theta;
}
transformed parameters {
  real sigma2;
  sigma2 = pow(sigma, 2);
}
model {
  y ~ normal(theta, sigma);
  theta ~ normal(0, 1000);
  sigma2 ~ inv_gamma(0.001, 0.001);
}
'

y <- c( 1.69045085, -1.41076082, -0.27909483, 
       -0.91387987,  3.21868429, -1.47282460, 
       -0.96524353, -2.45084934,  1.03838153, 
        1.79928679,  0.97826621,  0.67463830, 
       -1.08665864, -0.00509027,  0.43708128)
n <- length(y)

sample_data <- list(y = y,
                    n = n)
NormalMVfit <- stan(model_code = NormalMeanVar,
                   data = sample_data, verbose = FALSE)
```


```r
print(NormalMVfit, digits = 4, 
      pars = c("theta", "sigma", "sigma2"), probs = c(0.025, 0.975))
```

```
## Inference for Stan model: de363179728705eed697c67cf38b3332.
## 4 chains, each with iter=2000; warmup=1000; thin=1; 
## post-warmup draws per chain=1000, total post-warmup draws=4000.
## 
##          mean se_mean     sd    2.5%  97.5% n_eff   Rhat
## theta  0.0657  0.0090 0.4066 -0.7606 0.8497  2033 1.0025
## sigma  1.5419  0.0078 0.3092  1.0756 2.2731  1573 1.0007
## sigma2 2.4729  0.0287 1.0605  1.1568 5.1670  1364 1.0012
## 
## Samples were drawn using NUTS(diag_e) at Mon Jun 21 23:19:07 2021.
## For each parameter, n_eff is a crude measure of effective sample size,
## and Rhat is the potential scale reduction factor on split chains (at 
## convergence, Rhat=1).
```

Después de ejecutar las iteraciones necesarias, la salida del anterior código muestra una estimación puntual para la esperanza de $Y$ de 0.09 con un intervalo de credibilidad del 95\% dado por (-0.71, 0.91). Por otro lado, la estimación puntual de la desviación estándar de $Y$ es de 1.55 con un intervalo de credibilidad del 95\% dado por (1.09, 2.31).

Las figuras \@ref(fig:posNormalMVStan2) muestran la distribución posterior para este ejemplo, junto con la estimación puntual, correspondiente a la desviación estándar. 

\begin{figure}

{\centering \includegraphics{4Multiparametricos_files/figure-latex/posNormalMVStan2-1} 

}

\caption{Distribuciones posteriores.}(\#fig:posNormalMVStan2)
\end{figure}

A continuación se ilustra el uso de `R` en la programación de un algoritmo de Gibbs para los datos del ejemplo anterior. Se recalca que se utiliza la librería `MCMCpack` [@MCMCpack] para generar las realizaciones de la distribución Inversa-Gamma.


```r
library(MCMCpack)

#parametros previos de theta
mu <- 0
tau2 <- 1000
#parametros previos de sigma2
a <- 0.001
b <- 0.001

nsim <- 5000
theta.pos <- rep(NA, nsim)
sigma2.pos <- rep(NA, nsim)

# Valor inicial de theta
theta.pos[1] <- 0

#parametros posteriores de sigma2	
a.n <- a + n/2
b.n <- b + ((n - 1) * var(y) + n * (mean(y) - theta.pos[1]))/2
#simulacion de la distribucion posterior condicional de theta
sigma2.pos[1] <- rinvgamma(1, a.n, b.n)

#####################
# Muestreo de Gibbs #
#####################

for(i in 2:nsim){
  #parametros posteriores de theta	
  tau2.n <- 1 / ((n/sigma2.pos[i - 1]) + (1/tau2))
  mu.n <- tau2.n * (mean(y) * (n/sigma2.pos[i - 1]) + mu/tau2)
  #simulacion de la distribucion posterior condicional de theta
  theta.pos[i] <- rnorm(1, mean=mu.n, sd=sqrt(tau2.n))
  #parametros posteriores de sigma2	
  a.n <- a + n/2
  b.n <- b + ((n - 1) * var(y) + n * (mean(y) - theta.pos[i])) / 2
  #simulacion de la distribucion posterior condicional de theta
  sigma2.pos[i] <- rinvgamma(1, a.n, b.n)
}
```


```r
mean(theta.pos)
```

```
## [1] 0.09059077
```

```r
quantile(theta.pos, c(0.025, 0.975))
```

```
##       2.5%      97.5% 
## -0.6912883  0.9094063
```

```r
mean(sqrt(sigma2.pos))
```

```
## [1] 1.533898
```

```r
quantile(sqrt(sigma2.pos), c(0.025, 0.975))
```

```
##      2.5%     97.5% 
## 0.9983506 2.2915069
```

De donde podemos concluir que los resultados de este algoritmo de Gibbs coinciden con la estimación puntual para la esperanza de $Y$, su desviación estándar, y sus intervalos de credibilidad del 95%, obtenidos con `STAN`. 

\begin{figure}

{\centering \includegraphics{4Multiparametricos_files/figure-latex/plotsNormalMV-1} 

}

\caption{Convergencia de las distribuciones posteriores y diagramas de la función de autocorrelación en las cadenas.}(\#fig:plotsNormalMV)
\end{figure}

Las figuras \@ref(fig:plotsNormalMV) muestran que la convergencia de las cadenas es inmediata y que además no existen correlaciones importantes en los valores simulados de $\theta$ y $\sigma^2$. Por tanto, se concluye que se pueden utilizar directamente estos valores para la obtención de las estimaciones. Finalmente se ilustra la forma de obtener la distribución predictiva para el promedio muestral de 5 nuevos pacientes. 


```r
n.ast <- 5
y.bar <- c()

for(i in 1:(nsim/2)){
  y.bar[i] <- rnorm(1,
                    theta.pos[i + nsim/2],
                    sqrt(sigma2.pos[i + nsim/2]/n.ast))
}

mean(y.bar)
```

```
## [1] 0.1106463
```

```r
sd(y.bar)
```

```
## [1] 0.7960929
```

```r
quantile(y.bar,c(0.025,0.975))
```

```
##      2.5%     97.5% 
## -1.523240  1.607851
```

Podemos ver que se espera que el promedio de las pruebas en 5 nuevos pacientes es de 0.1106463, con un intervalo de credibilidad del 95% que es mucho más ancho que el de $\theta$, pues naturalmente $\bar{Y}$ tiene mayor incertidumbre que los parámetros del modelo; además, el tamaño de nuevos datos es de cinco, el cual es pequeño y hace que el pronóstico para $\bar{Y}^*$ no sea muy preciso. 

### Parámetros dependientes

En algunas situaciones es muy útil asumir una distribución previa conjugada, y para lograr eso no es posible establecer que los parámetros tengan distribuciones previas independientes. Bajo esta situación, la inferencia posterior de los parámetros de interés debe ser llevada a cabo en dos etapas: En la primera, se debe establecer la distribución previa conjunta para ambos parámetros siguiendo la siguiente regla
\begin{equation*}
p(\theta,\sigma^2)=p(\sigma^2)p(\theta \mid \sigma^2)
\end{equation*}

En la segunda etapa ya es posible analizar propiamente cada uno de los parámetros de interés siguiendo otra sencilla regla análoga
\begin{equation*}
p(\theta,\sigma^2 \mid \mathbf{Y})\propto p(\mathbf{Y} \mid \theta,\sigma^2)p(\theta,\sigma^2)
\end{equation*}

La anterior formulación conlleva a asignar una distribución previa para $\theta$ dependiente del parámetro $\sigma^2$. Esto quiere decir que en la distribución $p(\theta \mid \sigma^2)$, el valor de $\sigma^2$ se considera una constante fija y conocida, esta distribución previa está dada por^[La forma como la distribución previa de $\theta$ dependa de $\sigma^2$ es coherente con la información de Fisher sobre $\theta$ dada por $\sigma^{-2}$.]

\begin{equation*}
p(\theta \mid \sigma^2) \sim Normal(\mu,\sigma^2/c_0)
\end{equation*}

donde $c_0$ es una constante. Por otro lado, y siguiendo los argumentos del capítulo anterior, una posible opción para la distribución previa de $\sigma^2$, que no depende de $\theta$, corresponde a
\begin{equation*}
p(\sigma^2)\sim Inversa-Gamma(n_0/2,n_0\sigma^2_0/2)
\end{equation*}

De esta forma, podemos encontrar la distribución conjunta previa de $\theta$ y $\sigma^2$ como sigue:

\BeginKnitrBlock{proposition}
<span class="proposition" id="prp:unnamed-chunk-10"><strong>(\#prp:unnamed-chunk-10) </strong></span>La distribución conjunta previa de los parámetros $\theta$ y $\sigma^2$ está dada por
\begin{equation*}
\theta,\sigma^2 \sim Normal-Inversa-Gamma\left(\mu, c_0, \frac{n_0+1}{2},\frac{n_0\sigma^2_0}{2}\right)
\end{equation*}
\EndKnitrBlock{proposition}
<br>

\BeginKnitrBlock{proof}
\iffalse{} <span class="proof"><em>Prueba. </em></span>  \fi{}\begin{align*}
p(\theta,\sigma^2)&=p(\sigma^2)p(\theta \mid \sigma^2)\\
&\propto (\sigma^2)^{-\frac{n_0}{2}-1}\exp\left\{-\dfrac{n_0\sigma_0^2}{2\sigma^2}\right\}
(\sigma^2)^{-\frac{1}{2}}\exp\left\{-\frac{c_0}{2\sigma^2}(\theta-\mu)^2\right\}\\
&= (\sigma^2)^{-\frac{n_0+1}{2}-1}\exp\left\{-\frac{1}{2\sigma^2}\left[n_0\sigma^2_0+c_0(\theta-\mu)^2\right]\right\}
\end{align*}
la cual corresponde a la forma de la función de densidad de la distribución Normal-Inversa-Gamma.
\EndKnitrBlock{proof}
<br>

Una vez encontrada la distribución conjunta previa, procedemos a encontrar la distribución conjunta posterior, y así poder encontrar las estimaciones de $\theta$ y $\sigma^2$.

\BeginKnitrBlock{proposition}
<span class="proposition" id="prp:unnamed-chunk-12"><strong>(\#prp:unnamed-chunk-12) </strong></span>La distribución posterior conjunta de los parámetros $\theta$ y $\sigma^2$ está dada por
\begin{equation*}
\theta,\sigma^2\mid\mathbf{Y} \sim Normal-Inversa-Gamma\left(\mu_n, c_0+n, \frac{n_0+n+1}{2},\beta \right).
\end{equation*}

con 
\begin{equation*}
\beta=\dfrac{1}{2}\left(n_0\sigma^2_0+(n-1)S^2+\dfrac{c_0n}{c_0+n}(\mu-\bar{y})^2\right)
\end{equation*}

y
\begin{equation*}
\mu_n=\frac{\frac{n}{\sigma^2}\bar{Y}+\frac{c_0}{\sigma^2}\mu}{\frac{n}{\sigma^2}+\frac{c_0}{\sigma^2}}
=\frac{n\bar{Y}+c_0\mu}{n+c_0}
\end{equation*}
\EndKnitrBlock{proposition}
<br>


\BeginKnitrBlock{proof}
\iffalse{} <span class="proof"><em>Prueba. </em></span>  \fi{}En primer lugar, recordamos que la función de verosimilitud de la muestra está dada por
\begin{align}
p(\mathbf{Y} \mid \theta,\sigma^2)= \frac{1}{(2\pi\sigma^2)^{n/2}}
\exp\left\{-\frac{1}{2\sigma^2}\left[(n-1)S^2+n(\bar{y}-\theta)^2\right]\right\}
\end{align}

Por otro lado, se tiene que
\begin{align}
(\#eq:desarro1)
p(\theta,\sigma^2 \mid \mathbf{Y}) & \propto p(\mathbf{Y} \mid \theta,\sigma^2)p(\theta,\sigma^2)\notag\\
&\propto (\sigma^2)^{-\frac{n_0+n+1}{2}-1}
\exp\left\{-\frac{1}{2\sigma^2}\left[n_0\sigma^2_0+c_0(\theta-\mu)^2+(n-1)S^2+n(\bar{y}-\theta)^2\right]\right\}\notag\\
&= (\sigma^2)^{-\frac{n_0+n+1}{2}-1}\\
&\times
\exp\left\{-\frac{1}{2\sigma^2}\left[n_0\sigma^2_0+(n-1)S^2+(c_0+n)(\theta-\mu_n)^2+\frac{c_0n}{c_0+n}(\mu-\bar{y})^2\right]\right\}
\end{align}
puesto que
\begin{align*}
c_0(\theta-\mu)^2+n(\bar{y}-\theta)^2=(c_0+n)(\theta-\mu_n)^2+\frac{c_0n}{c_0+n}(\mu-\bar{y})^2
\end{align*}
\EndKnitrBlock{proof}
<br>

Para encontrar las distribuciones marginales posterior de cada uno de los parámetros se procede de la siguiente forma:
  
1. Para hallar la distribución posterior condicional de $\theta$, dada por $P(\theta \mid \sigma^2,\mathbf{Y})$, se debe considerar que $\sigma^2$ es una constante fija y conocida, tal como se consideró al principio de esta sección. Basado en lo anterior, es posible utilizar la siguiente regla de probabilidad
\begin{align*}
P(\theta \mid \sigma^2,\mathbf{Y})=\frac{p(\theta,\sigma^2 \mid \mathbf{Y})}{p(\sigma^2,\mathbf{Y})}p(\mathbf{Y})\propto p(\theta,\sigma^2 \mid \mathbf{Y})
\end{align*}
Lo anterior sugiere que la distribución marginal posterior de $\theta$, $p(\theta \mid \sigma^2,\mathbf{Y})$, se encuentra utilizando la distribución posterior conjunta, $p(\theta,\sigma^2 \mid \mathbf{Y})$, suponiendo que todas las expresiones que involucren al valor $\sigma^2$ se pueden incluir en la constante de proporcionalidad
2. Dado que $\sigma^2$ no depende de ningún otro parámetro entonces, utilizando la distribución posterior conjunta, es posible encontrar su distribución marginal posterior de la siguiente forma
\begin{align*}
p(\sigma^2 \mid \mathbf{Y})=\int p(\theta,\sigma^2 \mid \mathbf{Y}) d\theta
\end{align*}
3. Un razonamiento similar se puede formular para el parámetro $\theta$; utilizando la distribución posterior conjunta, es posible encontrar su distribución marginal posterior de la siguiente forma
\begin{align*}
p(\theta \mid \mathbf{Y})=\int p(\theta,\sigma^2 \mid \mathbf{Y}) d\sigma^2
\end{align*}

\BeginKnitrBlock{proposition}
<span class="proposition" id="prp:unnamed-chunk-14"><strong>(\#prp:unnamed-chunk-14) </strong></span>La distribución posterior de $\theta$ condicional a $\sigma^2,\mathbf{Y}$ está dada por
\begin{equation*}
\theta \mid \sigma^2,\mathbf{Y} \sim Normal(\mu_n,\sigma^2/(n+c_0))
\end{equation*}
con $\mu_n=\dfrac{n\bar{y}+c_0\mu}{n+c_0}$.
\EndKnitrBlock{proposition}
<br>

\BeginKnitrBlock{proof}
\iffalse{} <span class="proof"><em>Prueba. </em></span>  \fi{}Acudiendo a la distribución posterior conjunta dada en \@ref(eq:desarro1), tenemos que
\begin{align*}
p(\theta \mid \sigma^2,\mathbf{Y})&\propto 
p(\theta,\sigma^2 \mid \mathbf{Y}) \\
&\propto(\sigma^2)^{-\frac{n_0+n+1}{2}-1}\\
&\hspace{1cm}\times
\exp\left\{-\frac{1}{2\sigma^2}\left[n_0\sigma^2_0+(n-1)S^2+(c_0+n)(\theta-\mu_n)^2+\frac{c_0n}{c_0+n}(\mu-\bar{y})^2\right]\right\}\\
&\propto \exp\left\{-\frac{1}{2\sigma^2}(c_0+n)(\theta-\mu_n)^2\right\}
\end{align*}
la cual corresponde a la forma de la función de densidad de la distribución $Normal(\mu_n, \sigma^2/(n+c_0))$.
\EndKnitrBlock{proof}

En el anterior resultado, la media de la distribución condicional posterior $\mu_n$ se puede escribir como $\mu_n=\frac{n}{n+c_0}\bar{y}+\frac{c_0}{n+c_0}\mu$, promedio ponderado entre la estimación clásica $\bar{y}$ y la estimación previa $\mu$. Observando que las ponderaciones $\frac{n}{n+c_0}$ y $\frac{c_0}{n+c_0}$ forman una combinación lineal convexa, se puede definir a $c_0$ como el número de observaciones en la información previa. De esta forma, las ponderaciones de la estimación clásica y la estimación previa dependerán directamente de los tamaños muestrales respectivos.

\BeginKnitrBlock{proposition}
<span class="proposition" id="prp:PosterSigma2IG"><strong>(\#prp:PosterSigma2IG) </strong></span>La distribución marginal posterior del parámetro $\sigma^2$ es
\begin{equation*}
\sigma^2 \mid \mathbf{Y} \sim Inversa-Gamma\left(\frac{n+n_0}{2},\frac{(n+n_0)\sigma^2_n}{2}\right)
\end{equation*}

Donde $(n+n_0)\sigma^2_n=n_0\sigma^2_0+(n-1)S^2+\frac{c_0n}{c_0+n}(\mu-\bar{y})^2$ corresponde a una suma ponderada de la varianza previa, la varianza muestral y la diferencia entre la media muestral y la media previa.
\EndKnitrBlock{proposition}
<br>

\BeginKnitrBlock{proof}
\iffalse{} <span class="proof"><em>Prueba. </em></span>  \fi{}De la distribución posterior conjunta \@ref(eq:desarro1), e integrando con respecto a $\theta$, se tiene que

\begin{align*}
p(\sigma^2 \mid \mathbf{Y})&=\int p(\theta,\sigma^2 \mid \mathbf{Y}) \ d\theta\\
&\propto (\sigma^2)^{-\frac{n_0+n+1}{2}-1}
\exp\left\{-\frac{1}{2\sigma^2}\left[n_0\sigma^2_0+(n-1)S^2+\frac{c_0n}{c_0+n}(\mu-\bar{y})^2\right]\right\}\\
&\hspace{2cm}\times
\int_{-\infty}^{\infty}\exp\left\{-\frac{n+c_0}{2\sigma^2}(\theta-\mu_n)^2\right\} \ d\theta\\
&\propto (\sigma^2)^{-\frac{n_0+n}{2}-1}
\exp\left\{-\frac{1}{2\sigma^2}\left[n_0\sigma^2_0+(n-1)S^2+\frac{c_0n}{c_0+n}(\mu-\bar{y})^2\right]\right\}\\
&\hspace{2cm}\times
\int_{-\infty}^{\infty}\frac{\sqrt{n+c_0}}{\sqrt{2\pi\sigma^2}}
\exp\left\{-\frac{n+c_0}{2\sigma^2}(\theta-\mu_n)^2\right\} \ d\theta\\
&\propto (\sigma^2)^{-\frac{n_0+n}{2}-1}
\exp\left\{-\frac{(n+n_0)\sigma^2_n}{2\sigma^2}\right\}
\end{align*}

la cual corresponde a la forma funcional de la densidad $Inversa-Gamma(\frac{n+n_0}{2},\frac{(n+n_0)\sigma^2_n}{2})$.
\EndKnitrBlock{proof}
<br>

Dadas las distribuciones $p(\sigma^2\mid \mathbf{Y})$ y $p(\theta\mid \sigma^2, \mathbf{Y})$, podemos proceder de la siguiente forma para obtener valores simulados de $\theta$ y $\sigma^2$ y por consiguiente sus estimaciones puntuales. Si el número de iteraciones se fija como $G$, entonces se debe seguir el siguiente algoritmo:

1. Simular $G$ valores de la distribución de $\sigma^2|\mathbf{Y}$; es decir, de la distribución $Inversa-Gamma$ encontrada en el anterior resultado. Estos valores se denotan por $\sigma^2_{(1)},\sigma^2_{(2)},\cdots,\sigma^2_{(G)}$.
2. Para cada valor de $\sigma^2_{(g)}$, con $g=1,\cdots,G$, simular un valor de la distribución de $\theta|\sigma^2,\mathbf{Y}$; es decir, de la distribución $N(\mu_n,\sigma^2/(n+c_0))$, donde $\sigma^2$ se reemplaza por $\sigma^2_{(g)}$. De esta forma, se obtiene los valores $\theta_{(1)},\theta_{(2)},\cdots,\theta_{(G)}$.

Es claro que en el anterior algoritmo, no es necesario fijar algún valor inicial para $\theta$ o para $\sigma^2$. De la misma forma no existirán correlaciones sustantivas entre los valores simulados para ningún parámetro. Por lo tanto, estos valores se pueden usar directamente para el cálculo de las estimaciones, y no es necesario descartar los primeros valores simulados, ni realizar la fase de *thinning*.

Ahora bien, existe otra alternativa para obtener la estimación de $\theta$ y $\sigma^2$, la cual depende directamente de la distribución posterior de cada parámetro. En efecto, la distribución posterior de $\sigma^2$ ya se encontró en el resultado \@ref(prp:PosterSigma2IG), resta encontrar la distribución posterior de $\theta$, la cual se presenta en el siguiente resultado. 

\BeginKnitrBlock{proposition}
<span class="proposition" id="prp:PosThetaTnoestandar"><strong>(\#prp:PosThetaTnoestandar) </strong></span>La distribución posterior del parámetro $\theta$ corresponde a una $t$ no estandarizada con $n_0+n$ grados de libertad, parámetro de localización $\mu_n=\dfrac{n\bar{Y}+c_0\mu}{n+c_0}$ y parámetro de escala $\dfrac{\sigma_n}{\sqrt{c_0+n}}$, con $(n+n_0)\sigma^2_n=n_0\sigma^2_0+(n-1)S^2+\frac{c_0n}{c_0+n}(\mu-\bar{y})^2$. Esto es, 
\begin{equation*}
\theta \mid \mathbf{Y} \sim t_{n+n_0}\left(\mu_n, \frac{\sigma^2_n}{c_0+n}\right)
\end{equation*}
\EndKnitrBlock{proposition}
<br>

\BeginKnitrBlock{proof}
\iffalse{} <span class="proof"><em>Prueba. </em></span>  \fi{}Partiendo de la distribución posterior conjunta e integrando con respecto a $\sigma^2$, se tiene que

\begin{align*}
p(\theta \mid \mathbf{Y})&= \int_0^{\infty} p(\theta,\sigma^2 \mid \mathbf{Y}) \ d\sigma^2 \\
&\propto \int_0^{\infty} \left(\frac{1}{\sigma^2}\right)^{\frac{n_0+n+1}{2}+1}
\exp\left\{-\frac{1}{2\sigma^2}\left[(n_0+n)\sigma^2_n+(c_0+n)(\theta-\mu_n)^2\right]\right\} \ d\sigma^2
\end{align*}

Haciendo un cambio de variable tal que
\begin{equation*}
z=\frac{A}{2\sigma^2}, \ \ \ \ \ \ \ \ \ \ \ \text{donde} \ \ \ A=(n_0+n)\sigma^2_n+(c_0+n)(\theta-\mu_n)^2
\end{equation*}

Por tanto
\begin{equation*}
d\sigma^2=-\frac{A}{2z^2} \ dz
\end{equation*}

Volviendo a la integral en cuestión, se tiene que
\begin{align*}
p(\theta \mid \mathbf{Y})& \propto
\left(\frac{1}{A}\right)^{\frac{n_0+n+1}{2}+1}\int_{\infty}^{0} \frac{-A}{2z^2} (2z)^{\frac{n_0+n+1}{2}+1}e^{-z} \ dz \\
&\propto A^{-\frac{n_0+n+1}{2}}\underbrace{\int_{0}^{\infty} z^{\frac{n_0+n+1}{2}-1}e^{-z}\ dz}_{Gamma\left(\frac{n_0+n+1}{2},1\right)}\\
&\propto A^{-\frac{n_0+n+1}{2}}\\
&= \left[(n_0+n)\sigma^2_n+(c_0+n)(\theta-\mu_n)^2\right]^{-\frac{n_0+n+1}{2}}\\
&\propto \left[1+\frac{(c_0+n)(\theta-\mu_n)^2}{(n_0+n)\sigma^2_n}\right]^{-\frac{n_0+n+1}{2}}\\
&=\left[1+\frac{1}{n_0+n}\left(\frac{\theta-\mu_n}{\sigma_n/\sqrt{c_0+n}}\right)^2\right]^{-\frac{n_0+n+1}{2}}
\end{align*}

la cual corresponde a la forma de la función de densidad de la distribución deseada.
\EndKnitrBlock{proof}
<br>

Las distribuciones encontradas en los resultados \@ref(prp:PosterSigma2IG) y \@ref(prp:PosThetaTnoestandar) permiten estimar directamente los parámetros $\theta$ y $\sigma^2$ usando las esperanzas teóricas de las distribuciones posteriores. Por ende, las estimaciones puntuales son:

\begin{align*}
\hat{\theta}
&=\mu_n=\dfrac{n\bar{Y}+c_0\mu}{n+n_0}\\
\hat{\sigma}^2
&=\dfrac{(n+n_0)\sigma^2_n/2}{(n+n_0)/2-1}=\dfrac{(n+n_0)\sigma^2_n}{n+n_0-2}\approx\sigma^2_n=\dfrac{n_0\sigma^2_0+(n-1)S^2+\frac{c_0n}{c_0+n}(\mu-\bar{y})^2}{n+n_0}
\end{align*}

Los intervalos del $(1-\alpha)\times 100\%$ de credibilidad para $\theta$ y $\sigma^2$ se construyen usando los percentiles $\alpha/2$ y $1-\alpha/2$ de las respectivas distribuciones posteriores dadas en los resultados mencionados anteriormente. Ilustramos el uso de la metodología en el siguiente ejemplo.

\BeginKnitrBlock{example}
<span class="example" id="exm:unnamed-chunk-18"><strong>(\#exm:unnamed-chunk-18) </strong></span>Para los datos de función renal \cite{Efronims} que se muestran en el Ejemplo \@ref(exm:EjeRenal), suponga que la información previa está contenida en la medición de función renal para una muestra de 12 pacientes dadas por: -1.3619, -1.1116, -0.4744, -0.5663, 2.2056, 0.9491, 0.2298, -0.7933, 1.0198, -0.9850, 3.5679 y -1.9504. La media y la varianza muestral de estas 12 observaciones corresponden a $\mu=0.060775$ y $\sigma^2_0=2.598512$; por consiguiente $c_0=n_0=12$. 

Por otro lado, la media y la varianza muestral de los 15 pacientes en la información actual son $\bar{y}=0.08349249$ y $S^2=2.301684$. De esta forma, los parámetros de las distribuciones marginales posteriores de $\theta$ y $\sigma^2$ se pueden calcular como $\mu_n=\frac{15}{15+12}\times 0.08349249+\frac{12}{15+12}\times 0.060775=0.07339583$ y $$\sigma^2_n=\dfrac{12*2.598512+14*2.301684+6.666667*(0.060775-0.08349249)^2}{15+12}=2.348487$$ En conclusión, las distribuciones marginales posterior de $\theta$ y $\sigma^2$ están dadas por
\begin{equation*}
\theta|\mathbf{Y}\sim t_{27}(0.07339583,2.348487/27=0.086981)
\end{equation*}

y
\begin{equation*}
\sigma^2|\mathbf{Y}\sim Inversa-Gamma(27/2=13.5,\ 27*2.348487/2=31.70457)
\end{equation*}

Así, la estimación Bayesiana de $\theta$ es $\mu_n=0.073$ y un intervalo de credibilidad del $95\%$ para $\theta$ se puede calcular como $0.073\pm t_{27,0.975}*\sqrt{0.086981}=(-0.53,\ 0.68)$. Por otro lado, la estimación Bayesiana de $\sigma^2$ está dada por $31.70457/(13.5-1)=2.53$, y un intervalo de credibilidad del $95\%$ para $\sigma^2$ se puede calcular como los percentiles $2.5\%$ y $97.5\%$ de la distribución $Inversa-Gamma(13.5,\ 31.70457)$, dados por $(1.5, 4.4)$. Los anteriores cálculos se ilustran en el siguiente código `R`.
\EndKnitrBlock{example}
<br>


```r
library(pscl)
# Datos de la informacion previa
x <- c(-1.3619, -1.1116, -0.4744, -0.5663, 
        2.2056,  0.9491,  0.2298, -0.7933, 
        1.0198, -0.9850,  3.5679, -1.9504)
# Datos de la informacion actual
y <- c( 1.69045085, -1.41076082, -0.27909483, 
       -0.91387987,  3.21868429, -1.47282460, 
       -0.96524353, -2.45084934,  1.03838153, 
        1.79928679,  0.97826621,  0.67463830, 
       -1.08665864, -0.00509027,  0.43708128)
# Paramatros de la distribucion previa
n0 <- c0 <- 12
mu <- mean(x)
sigma2_0 <- var(x)

# Informacion
n <- length(y)
bar.y <- mean(y)
S2 <- var(y)

# Algunos paramatros de la distribucion posterior
mu.n <- (n * bar.y + c0 * mu)/(n + n0)
sigma2_n <- (n0 * sigma2_0 + (n - 1) * S2 +
               c0 * n * (mu - bar.y)^2/(c0 + n))/(n + n0)
# Estimacion puntual
theta.hat <- mu.n 
sigma2.hat <- (n + n0) * sigma2_n/(n + n0 - 2)
theta.hat
```

```
## [1] 0.07339583
```

```r
sigma2.hat
```

```
## [1] 2.536367
```

```r
# Intervalo de credibilidad de 95% para theta
mu.n + qt(c(0.025, 0.975), df = n + n0) * 
  sqrt(sigma2_n/(n + n0))
```

```
## [1] -0.5317412  0.6785329
```

```r
# Intervalo de credibilidad de 95% para sigma2
qigamma(0.025, alpha = (n + n0)/2, 
        beta = (n + n0) * sigma2_n/2)
```

```
## [1] 1.467991
```

```r
qigamma(0.975, alpha = (n + n0)/2, 
        beta = (n + n0) * sigma2_n/2)
```

```
## [1] 4.351026
```

Otra forma de estimar los parámetros $\theta$ y $\sigma^2$ es utilizando métodos de simulación directa. De esta forma, como se expuso anteriormente, se generan primero los valores de $\sigma^2$ y posteriormente los valores de $\theta$.


```r
n.sim <- 5000
sigma2.res <- rigamma(n.sim, (n + n0)/2, 
                        (n + n0) * sigma2_n/2)
theta.res <- c()
for(i in 1:n.sim){
  theta.res[i] <- rnorm(1, mu.n, sqrt(sigma2.res[i]/(n + c0)))
}

# Estimaciones puntuales
mean(theta.res)
```

```
## [1] 0.06467424
```

```r
mean(sigma2.res)
```

```
## [1] 2.53312
```

```r
# Intervalos de credibilidad del 95%
quantile(theta.res, c(0.025, 0.975))
```

```
##       2.5%      97.5% 
## -0.5560560  0.6429974
```

```r
quantile(sigma2.res, c(0.025, 0.975))
```

```
##    2.5%   97.5% 
## 1.48546 4.34731
```

### Parámetros no informativos

En esta sección consideramos el tratamiento de los datos cuando no tenemos información previa disponible. Suponga que $\mathbf{Y}=\{Y_1,\ldots,Y_n\}$ corresponde a una muestra de variables aleatorias con distribución $Normal(\theta, \sigma^2)$. Luego, la función de distribución conjunta o verosimilitud está dada por 

\begin{equation*}
(\#eq:VeroNormal)
p(\mathbf{Y} \mid \theta, \sigma^2)=\frac{1}{(2\pi\sigma^2)^{n/2}}\exp\left\{-\frac{1}{2\sigma^2}\sum_{i=1}^n(y_i-\theta)^2\right\}
\end{equation*}

En primer lugar suponga que los parámetros tienen distribuciones previas independientes. Por ende, en esta primera etapa se realizará el análisis suponiendo que estas distribuciones son no informativas. Lo anterior implica que la distribución previa conjunta de los parámetros de interés está dada por
\begin{equation}
p(\theta,\sigma^2)=p(\theta)p(\sigma^2)
\end{equation}

Como la distribución previa de $\theta$ es normal, es fácil verificar que ésta empieza a tener las características propias de una distribución no informativa cuando la varianza de la misma se vuelve muy grande, sin importar el valor de la media. Cuando esto sucede, la forma de la distribución previa de $\theta$ se torna plana y es lógico pensar que puede ser acercada mediante una distribución constante, tal que

\begin{equation*}
p(\theta)\propto cte
\end{equation*}

Por otro lado, @Gelman03 afirmaN que la distribución Inversa-Gamma, la cual es la distribución previa para el parámetro $\sigma^2$, se vuelve no informativa cuando los hiper-parámetros toman valores muy cercanos a cero. De esta forma haciendo tender $\alpha \longrightarrow 0$ y $\beta \longrightarrow 0$, entonces la distribución previa de $\sigma^2$ se convierte en

\begin{equation*}
p(\sigma^2)\propto \sigma^{-2}
\end{equation*}

la cual coincide con la distribución previa no informativa de Jeffreys discutida en las secciones anteriores. Por lo anterior, la distribución previa no informativa conjunta está dada por

\begin{equation}
p(\theta,\sigma^2)\propto \sigma^{-2}
\end{equation}

Bajo este marco de referencia se tiene el siguiente resultado sobre la distribución posterior de $\theta$

\BeginKnitrBlock{proposition}
<span class="proposition" id="prp:PosThetaNoInformativa"><strong>(\#prp:PosThetaNoInformativa) </strong></span>La distribución posterior del parámetro $\theta$ sigue una distribución $t$ no estandarizada con $n-1$ grados de libertad, parámetro de localización $\bar{Y}$ y parámetro de escala $\frac{S^2}{n}$; esto es, 

\begin{equation*}
\theta \mid \mathbf{Y}\sim t_{n-1}\left(\bar{y},\frac{S^2}{n}\right).
\end{equation*}

Donde $(n-1)S^2=\sum_{i=1}^n(Y_i-\bar{Y})^2$. Esta distribución también puede expresarse como
\begin{equation*}
\frac{\theta-\bar{y}}{S/\sqrt{n}} \mid \mathbf{Y} \sim t_{n-1}
\end{equation*}

donde $t_{n-1}$ denota la distribución $t$ estandarizada con $n-1$ grados de libertad.
\EndKnitrBlock{proposition}
<br>

\BeginKnitrBlock{proof}
\iffalse{} <span class="proof"><em>Prueba. </em></span>  \fi{}En primer lugar nótese que la distribución posterior conjunta de los parámetros de interés es
\begin{align}
(\#eq:PropThetaSigma2)
p(\theta,\sigma^2 \mid \mathbf{Y})& \propto p(\theta,\sigma^2)p(\mathbf{Y} \mid \theta,\sigma^2) \notag \\
& \propto \frac{1}{\sigma^2}\frac{1}{(2\pi\sigma^2)^{n/2}}\exp\left\{-\frac{1}{2\sigma^2}\sum_{i=1}^n(y_i-\theta)^2\right\} \notag\\
& \propto \left(\frac{1}{\sigma^2}\right)^{n/2+1}
\exp\left\{-\frac{1}{2\sigma^2}\left[\sum_{i=1}^n(y_i-\bar{y})^2+n(\bar{y}-\theta)^2\right]\right\} \notag \\
&= \left(\frac{1}{\sigma^2}\right)^{n/2+1}
\exp\left\{-\frac{1}{2\sigma^2}\left[(n-1)S^2+n(\bar{y}-\theta)^2\right]\right\}
\end{align}
  
Para hallar la distribución marginal posterior de $\theta$ es necesario integrar la anterior expresión con respecto a $\sigma^2$. Con esto, se tiene que

\begin{align*}
p(\theta \mid \mathbf{Y})&= \int_0^{\infty} p(\theta,\sigma^2 \mid \mathbf{Y}) \ d\sigma^2 \\
&\propto \int_0^{\infty} \left(\frac{1}{\sigma^2}\right)^{n/2+1}
\exp\left\{-\frac{1}{2\sigma^2}\left[(n-1)S^2+n(\bar{y}-\theta)^2\right]\right\} \ d\sigma^2
\end{align*}

Haciendo un cambio de variable tal que
\begin{equation*}
z=\frac{A}{2\sigma^2}, \ \ \ \ \ \ \ \ \ \ \ \text{donde} \ \ \ A=(n-1)S^2+n(\bar{y}-\theta)^2
\end{equation*}

Por tanto
\begin{equation*}
d\sigma^2=-\frac{A}{2z^2} \ dz
\end{equation*}

Entonces, volviendo a la integral en cuestión, se tiene que
\begin{align*}
p(\theta \mid \mathbf{Y})& \propto
\left(\frac{1}{A}\right)^{n/2+1}\int_{\infty}^{0} \frac{-A}{2z^2} (2z)^{n/2+1}e^{-z} \ dz \\
&\propto A^{-n/2}\underbrace{\int_{0}^{\infty} z^{n/2-1}e^{-z}\ dz}_{Gamma(n/2)}\\
&\propto A^{-n/2}\\
&= [(n-1)S^2+n(\bar{y}-\theta)^2]^{-n/2}\\
&\propto \left[1+\frac{n(\bar{y}-\theta)^2}{(n-1)S^2}\right]^{-n/2}
=\left[1+\frac{1}{n-1}\left(\frac{\bar{y}-\theta}{S/\sqrt{n}}\right)^2\right]^{-\frac{(n-1)+1}{2}}
\end{align*}

la cual corresponde a la función de densidad de distribución de una variable aleatoria con distribución $t_{n-1}(\bar{y},S^2/n)$.
\EndKnitrBlock{proof}
<br>

\BeginKnitrBlock{proposition}
<span class="proposition" id="prp:PosSigma2NoInformativa"><strong>(\#prp:PosSigma2NoInformativa) </strong></span>La distribución posterior del parámetro $\sigma^2$ sigue una distribución
\begin{equation*}
\sigma^2 \mid \mathbf{Y} \sim Inversa-Gamma((n-1)/2,(n-1)S^2/2).
\end{equation*}
\EndKnitrBlock{proposition}
<br>

\BeginKnitrBlock{proof}
\iffalse{} <span class="proof"><em>Prueba. </em></span>  \fi{}Utilizando el mismo argumento del anterior resultado, se tiene que
\begin{align*}
p(\sigma^2 \mid \mathbf{Y})&= \int_{-\infty}^{\infty} p(\theta,\sigma^2 \mid \mathbf{Y}) \ d\theta \\
& \propto \int_{-\infty}^{\infty} \left(\frac{1}{\sigma^2}\right)^{n/2+1}
\exp\left\{-\frac{1}{2\sigma^2}\left[(n-1)S^2+n(\bar{y}-\theta)^2\right]\right\} \ d\theta \\
& = \left(\frac{1}{\sigma^2}\right)^{n/2+1} \sqrt{2\pi\sigma^2/n}\exp\left\{-\frac{1}{2\sigma^2}(n-1)S^2\right\}\underbrace{\int_{-\infty}^{\infty} \frac{1}{\sqrt{2\pi\sigma^2/n}} \exp\left\{-\frac{n}{2\sigma^2}(\bar{y}-\theta)^2\right\} \ d\theta}_{\text{vale $1$}} \\
& \propto (\sigma^2)^{-n/2-1/2}\exp\left\{-\frac{1}{2\sigma^2}(n-1)S^2\right\}\\
&= (\sigma^2)^{-\frac{n-1}{2}-1}\exp\left\{-\frac{1}{2\sigma^2}(n-1)S^2\right\}
\end{align*}

La cual corresponde a la función de densidad de la distribución $Inversa-Gamma((n-1)/2,(n-1)S^2/2)$.
\EndKnitrBlock{proof}
<br>

De los resultados \ref{PosThetaNoInformativa} y \ref{PosSigma2NoInformativa}, podemos ver que cuando no se dispone de información previa, la estimación bayesiana de $\theta$ y $\sigma^2$ están dadas por

\begin{align*}
\hat{\theta}_B&=E(\theta\mid\mathbf{Y})=\bar{Y}\\
\hat{\sigma}^2_B&=E(\sigma^2\mid\mathbf{Y})=\dfrac{(n-1)S^2/2}{(n-1)/2-1}=\dfrac{n-1}{n-3}S^2\approx S^2
\end{align*}

Por consiguiente, podemos concluir que la estimación bayesiana de $\theta$ cuando no hay información previa es idéntica a la estimación clásica de $\theta$, mientas que la de $\sigma^2$ es muy similar a la estimación clásica. En cuanto a la estimación por intervalo de credibilidad, podemos ver que un intervalo de crediblidad de $(1-\alpha)\times 100\%$ está dado por los percentiles $\alpha/2$ y $1-\alpha/2$ de la distribución $t_{n-1}\left(\bar{Y},\dfrac{S^2}{n}\right)$. Se puede ver que estos corresponden a $\bar{Y}+t_{n-1,\alpha/2}\dfrac{S}{\sqrt{n}}$ y $\bar{Y}+t_{n-1,1-\alpha/2}\dfrac{S}{\sqrt{n}}$. En conclusión, un intervalo de credibildad para $\theta$ está dado por $\bar{Y}\pm t_{n-1,1-\alpha/2}\dfrac{S}{\sqrt{n}}$, el cual es idéntico al intervalo de confianza para $\theta$ en la estadística clásica.

En cuanto al intervalo de crediblidad para $\sigma^2$, este está dado por los percentiles $\alpha/2$ y $1-\alpha/2$ de la distribución $Inversa-Gamma((n-1)/2,\ (n-1)S^2/2)$. En la estadística clásica, el intervalo de confianza para $\sigma^2$ está dada por \begin{equation*}
IC(\sigma^2)=\left(\dfrac{(n-1)S^2}{\chi^2_{n-1,1-\alpha/2}},\ \dfrac{(n-1)S^2}{\chi^2_{n-1,\alpha/2}}\right)
\end{equation*}

Aunque la forma de estos dos intervalos son muy diferentes, resultan ser idénticos. A continuación mostramos el porqué. Suponga que $a$ es el percentil $\alpha/2$ de la distribución $Inversa-Gamma((n-1)/2,\ (n-1)S^2/2)$, esto es, si $X\sim Inversa-Gamma((n-1)/2,\ (n-1)S^2/2)$, entonces $Pr(X<a)=\alpha/2$. Ahora por propiedades de la distribución $Inversa-Gamma$, se tiene que $\dfrac{X}{(n-1)S^2}\sim Inversa-Gamma(\frac{n-1}{2},\ \frac{1}{2})$. Por la relación entre la distribución $Gamma$ y la distribución $Inversa-Gamma$, tenemos que $\dfrac{(n-1)S^2}{X}\sim Gamma(\frac{n-1}{2},\ 2)$, es decir, $\dfrac{(n-1)S^2}{X}\sim\chi^2_{n-1}$, de donde tenemos que
\begin{align*}
\frac{\alpha}{2}&=Pr(X<a)\\
&=Pr\left(\dfrac{(n-1)S^2}{X}>\dfrac{(n-1)S^2}{a}\right)
\end{align*}

Esto es, $\dfrac{(n-1)S^2}{a}$ es el percentil $1-\alpha/2$ de la distribución $\chi^2_{n-1}$, esto es,  $\dfrac{(n-1)S^2}{a}=\chi^2_{n-1,1-\alpha/2}$, de donde $a=\dfrac{(n-1)S^2}{\chi^2_{n-1,1-\alpha/2}}$, así concluimos que el límite inferior del intervalo de credibilidad coincide con el límite inferior del intervalo de confianza. Análogamente se puede ver que también los límites superiores coinciden, y así vemos que el intervalo para $\sigma^2$ coincide en la estadística clásica y la estadística bayesiana sin información previa.

#### Enfoque alterno para estimar $\theta$ y $\sigma^2$ {-}

Existe otra forma de obtener las estimaciones para el parámetro $\theta$. Recordando la expresión \ref{PropThetaSigma2}, podemos afirmar que 
\begin{equation*}
\theta \mid \sigma^2, \mathbf{Y} \sim Normal(\bar{y},\sigma^2/n)
\end{equation*}

puesto que 
\begin{align*}
p(\theta \mid \sigma^2,\mathbf{Y})&\propto p(\theta, \sigma^2 \mid\mathbf{Y})\\
&\propto\exp\left\{-\frac{1}{2\sigma^2}\left[(n-1)S^2+n(\bar{y}-\theta)^2\right]\right\}\\
&=\exp\left\{-\frac{n}{2\sigma^2}(\bar{y}-\theta)^2\right\}
\end{align*}

La cual corresponde a la función de densidad de la distribución $Normal(\bar{y},\sigma^2/n)$. De esta forma, usando las distribución $p(\sigma^2\mid\mathbf{Y})$ y $p(\theta\mid\sigma^2,\mathbf{Y})$, podemos implementar el siguiente procedimiento para obtener valores simulados de $\theta$ y $\sigma^2$. Si el número de iteraciones se fija como $G$, entonces se procede a:

1. Simular $G$ valores de la distribución de $\sigma^2|\mathbf{Y}$ - es decir, de la distribución $Inversa-Gamma((n-1)/2,(n-1)S^2/2)$. Estos valores se denotan por $\sigma^2_{(1)},\sigma^2_{(2)},\cdots,\sigma^2_{(G)}$.
2. Para cada valor de $\sigma^2_{(g)}$, con $g=1,\cdots,G$, simular un valor de la distribución de $\theta|\sigma^2,\mathbf{Y}$ - es decir, de la distribución $N(\bar{y},\sigma^2/n)$, donde $\sigma^2$ se reemplaza por $\sigma^2_{(g)}$. De esta forma, se obtiene los valores $\theta_{(1)},\theta_{(2)},\cdots,\theta_{(G)}$.

Las estimaciones de $\theta$ y $\sigma^2$ se pueden obtener de los valores obtenidos $\theta_{(1)},\theta_{(2)},\cdots,\theta_{(G)}$ y $\sigma^2_{(1)},\sigma^2_{(2)},\cdots,\sigma^2_{(G)}$.

### Distribución predictiva

La distribución predictiva para una nueva observación $\tilde{Y}$ está dada por 
\begin{align*}
p(\tilde{y}\mid\mathbf{Y})
&=\int\int p(\tilde{y}\mid\theta,\sigma^2) p(\theta,\sigma^2\mid\mathbf{Y})\ d\theta\ d\sigma^2\\
&=\int\int p(\tilde{y}\mid \theta,\sigma^2)p(\theta\mid\sigma^2,\mathbf{Y})p(\sigma^2\mid\mathbf{Y})\ d\theta\ d\sigma^2\\
&=\int\left(\int p(\tilde{y}\mid \theta,\sigma^2)p(\theta\mid\sigma^2,\mathbf{Y})\ d\theta\right)p(\sigma^2\mid\mathbf{Y})\ d\sigma^2
\end{align*}

En la integral dentro del paréntesis, el parámetro $\sigma^2$ permanece fijo; por lo cual, dicha integral corresponde a la distribución $N\left(\bar{y},\left(1+\dfrac{1}{n}\right)\sigma^2\right)$. De esta forma, al combinarla con la distribución posterior de $\sigma^2$, tenemos que
\begin{align*}
&\ \ \ \ p(\tilde{y}\mid\mathbf{Y})\\
&=\int_0^\infty \dfrac{1}{\sqrt{2\pi(1+\frac{1}{n})\sigma^2}}\exp\left\{-\dfrac{1}{2\sigma^2(1+\frac{1}{n})}(\tilde{y}-\bar{y})^2\right\}\dfrac{\left(\frac{(n-1)S^2}{2}\right)^{(n-1)/2}}{\Gamma\left(\frac{n-1}{2}\right)}(\sigma^2)^{-\frac{n-1}{2}-1}\exp\left\{-\dfrac{(n-1)S^2}{2\sigma^2}\right\}\ d\sigma^2
\end{align*}

Después de realizar los pasos algebraicos necesarios, se encuentra que 
\begin{equation}
p(\tilde{y}\mid\mathbf{Y})=\dfrac{\Gamma(n/2)}{\Gamma((n-1)/2)}\dfrac{1}{\sqrt{\pi(n-1)}}\left(\left(1+\frac{1}{n}\right)S^2\right)^{-1/2}\left(1+\dfrac{1}{n-1}\dfrac{(\tilde{y}-\bar{y})^2}{\left(1+\frac{1}{n}\right)S^2}\right)^{-n/2}
\end{equation}

La cual corresponde a la distribución $t$ no estandarizada con $n-1$ grados de libertad, parámetro de localización $\bar{y}$ y parámetro de escala $(1+\frac{1}{n})S^2$. De esta forma, podemos ver que los dos primeros momentos de esta distribución están dados por
\begin{align*}
E(\tilde{Y}\mid\mathbf{Y})&=\bar{y}\\
Var(\tilde{Y}\mid\mathbf{Y})&=\dfrac{n-1}{n-3}\left(1+\frac{1}{n}\right)S^2=\dfrac{(n-1)(n+1)}{n(n-3)}S^2
\end{align*}

Otra manera equivalente de conocer el comportamiento probabilístico de $\tilde{y}$ es por medio de la simulación. Se debe simular en primer lugar valores de $\theta$ y de $\sigma^2$ de la distribución posterior $p(\theta,\ \sigma^2\mid\mathbf{Y})$ usando el muestreo de Gibbs y posteriormente simulando valores de $\tilde{y}$ de la distribución $p(\tilde{y}\mid\theta,\ \sigma^2)$. En la figura \ref{PredictivaYprioriNoninformativa} se muestran el histograma de 10 mil valores de $\tilde{Y}$ simulados de esta forma, donde los datos  corresponden a 20 datos simulados de la distribución $N(12, 3^2)$. En la misma gráfica se observa también la función de densidad de la distribución $t$, podemos ver que los valores simulados de $\tilde{Y}$ efectivamente coinciden con la distribución predictiva de $\tilde{Y}$. Por lo anterior, se puede calcular un predictor de $\tilde{Y}$ como el promedio de los 10 mil valores simulados, y calcular el intervalo de predicción usando los percentiles de estos 10 mil valores.


```r
library("MCMCpack")
nsim <- 10000
y.tilde <- theta.pos <- c()
n <- 20
y <- rnorm(n, 12, 3)
y.bar <- mean(y)

S2 <- var(y)
a.n <- (n - 1)/2
b.n <- (n - 1) * S2/2
sigma2.pos <- rinvgamma(nsim, a.n, b.n)

for(i in 1:10000){
	theta.pos[i] <- rnorm(1, y.bar, sqrt(sigma2.pos[i]/n))
	y.tilde[i] <- rnorm(1, theta.pos[i], sqrt(sigma2.pos[i]))
}

tn <- function(x){
	v2 <- (n + 1) * S2/n
	tnpred <- (pi * (n - 1) * v2)^(- 0.5) * gamma(n/2) * 
	  (1 + (x - y.bar)^2/((n - 1) * v2))^(- n/2)/gamma((n - 1)/2)
	return(tnpred)
}

ggplot() + 
  geom_histogram(data = data.frame(y.tilde),
                 aes(x = y.tilde, y=..density..),
                 fill = "gray",
                 color = "black") +
  stat_density(fun = tn, aes(x = y.tilde), alpha = 0, col = 2) 
```

\begin{figure}

{\centering \includegraphics{4Multiparametricos_files/figure-latex/PredictivaYprioriNoninformativa-1} 

}

\caption{10 mil valores simulados de $\tilde{Y}$ y la función de densidad de la distribución predictiva de $\tilde{Y}$ con parámetros no informativos.}(\#fig:PredictivaYprioriNoninformativa)
\end{figure}


Por otro lado, si consideramos parámetros informativos podríamos usar este mismo método de simulación que es realmente efectivo y que se convertirá de ahora en adelante en nuestro medio para realizar las importantes validaciones predictivas en los diferentes modelos bayesianos que se trabajarán en este libro. Por ejemplo, para el mismo conjunto de datos anterior supongamos que en un experimento anterior se recolectaron 20 datos con media 10 y varianza 5. En la figura \@ref(fig:predictivaTsimu) se muestran el histograma de 10 mil valores simulados para una nueva observación. Nótese que convenientemente la forma de la distribución predictiva teórica coincide plenamente con ls distribución predictiva simulada. Esta agradable propiedad nos acompañará en el resto de los capítulos subsecuentes. 



```r
library(MCMCpack)
nsim <- 10000
y.tilde <- theta.pos <- sigma2.pos <- c()
n <- 20
y <- rnorm(n, 12, 3)
y.bar <- mean(y)
S2 <- var(y)

#parametros previos de theta
mu <- 10
tau2 <- 4
#parametros previos de sigma2
n_0 <- 20
sigma2_0 <- 5

# Valor inicial de theta
theta.pos[1] <- 0

#parametros posteriores de sigma2	
a.n <- (n_0 + n)/2
b.n <- (n_0 * sigma2_0 + (n - 1) * S2 + 
          n * (mean(y) - theta.pos[1]))/2
sigma2.pos[1] <- rinvgamma(1, a.n, b.n)

# Valor inicial de y.tilde
y.tilde[1] <- rnorm(1, theta.pos[1], sqrt(sigma2.pos[1]))
#####################
# Muestreo de Gibbs #
#####################

for(i in 2:nsim){
  #parametros posteriores de theta	
  tau2.n <- 1 / ((n/sigma2.pos[i - 1]) + (1/tau2))
  mu.n <- tau2.n * (mean(y) * (n/sigma2.pos[i - 1]) + mu/tau2)
  #simulacion de la distribucion posterior condicional de theta
  theta.pos[i] <- rnorm(1, mean = mu.n, sd = sqrt(tau2.n))
  #parametros posteriores de sigma2	
  b.n <- (n_0 * sigma2_0 + (n - 1) * S2 + 
          n * (mean(y) - theta.pos[i - 1]))/2
  #simulacion de la distribucion posterior condicional de theta
  sigma2.pos[i] <- rinvgamma(1, a.n, b.n)
  #simulacion de la distribucion predictiva
  y.tilde[i] <- rnorm(1, theta.pos[i], sqrt(sigma2.pos[i]))
}

tn <- function(x){
	v0 <- n_0 * sigma2_0 + n * sigma2_c
	tnpred <- (pi * v0)^(- 0.5) * gamma((n_0 + n + 1)/2) * 
	  (1 + (x - theta)^2/v0)^(- (n_0 + n + 1)/2)/gamma((n_0 + n)/2)
	return(tnpred)
}

ggplot() + 
  geom_histogram(data = data.frame(y.tilde),
                 aes(x = y.tilde, y=..density..),
                 fill = "gray",
                 color = "black") +
  stat_density(fun = tn, aes(x = y.tilde), alpha = 0, col = 2) 
```

\begin{figure}

{\centering \includegraphics{4Multiparametricos_files/figure-latex/predictivaTsimu-1} 

}

\caption{10 mil valores simulados de $\tilde{Y}$ y la función de densidad de la distribución predictiva de $\tilde{Y}$ con parámetros informativos.}(\#fig:predictivaTsimu)
\end{figure}



## Modelo Normal multivariante con media desconocida y varianza conocida

Cuando la distribución usada para describir el comportamiento de los datos es una distribución normal multivariante, las técnicas de inferencia no se distancian mucho del caso univariado. Se debe tener en cuenta el manejo matricial de las formas cuadráticas y las propiedades básicas del cálculo de matrices. Los desarrollos y resultados derivados de esta sección redundarán en el análisis de los modelos lineales con el enfoque bayesiano.

Sea $\mathbf{Y}=(Y_1,\ldots,Y_p)'$ un vector aleatorio cuya distribución es normal multivariante dada por

\begin{equation}
p(\mathbf{Y} \mid \btheta,\bSigma)\propto \mid \bSigma \mid ^{-1/2}
\exp\left\{-\frac{1}{2}(\mathbf{y}-\btheta)'\bSigma^{-1}(\mathbf{y}-\btheta)\right\}
\end{equation}
  
En donde $\btheta=(\theta_1,\ldots,\theta_p)'$ es el vector que contiene la media de cada uno de los componentes del vector $\mathbf{Y}$ y $\bSigma$ es la matriz de varianzas y covarianzas de orden $p\times p$, simétrica y definida positiva. La verosimilitud para una muestra de $n$ vectores aleatorios  independientes e idénticamente distribuidos está dada por

\begin{equation*}
  p(\mathbf{Y}_1\ldots,\mathbf{Y}_n \mid \btheta,\bSigma)\propto \mid \bSigma \mid ^{-n/2}
  \exp\left\{-\frac{1}{2}\sum_{i=1}^n(\mathbf{y}_i-\btheta)'\bSigma^{-1}(\mathbf{y}_i-\btheta)\right\}
\end{equation*}
    
Los parámetros que requieren estimación corresponden al vector de medias $\btheta$ y la matriz de varianzas y covarianzas $\bSigma$. Por ahora, se asume que $\bSigma$ es conocida y nos centramos en la estimación del vector de medias $\btheta$. Para la distribución previa, considerando que en general no hay restricción sobre los valores de los componentes de $\btheta$, asumimos que $\btheta$ sigue una distribución previa normal multivariante informativa y parametrizada por los hiper-parámetros $\bmu$ y $\bGamma$

\begin{equation*}
p(\btheta \mid \bmu,\bGamma)\propto \mid \bGamma \mid ^{-1/2}
\exp\left\{-\frac{1}{2}(\btheta-\bmu)'\bGamma^{-1}(\btheta-\bmu)\right\}
\end{equation*}
    
Nótese que esta distribución se hace no informativa cuando $\mid \bGamma^{-1} \mid \longrightarrow 0$, sin importar el valor del vector de medias previa $\bmu$. En el siguiente resultado, encontramos la distribución posterior del parámetro $\btheta$.

\BeginKnitrBlock{proposition}
<span class="proposition" id="prp:unnamed-chunk-23"><strong>(\#prp:unnamed-chunk-23) </strong></span>La distribución posterior del vector $\btheta$ sigue una distribución normal multivariante
\begin{equation*}
\btheta \mid \mathbf{Y},\bSigma \sim N_p (\bmu_n,\bGamma_n).
\end{equation*}

En donde
\begin{align}
(\#eq:Gamman)
\bGamma_n &= \left(\bGamma^{-1}+n\bSigma^{-1}\right)^{-1}\\
(\#eq:mun)
\bmu_n &= \bGamma_n(\bGamma^{-1}\bmu+n \bSigma^{-1}\bar{\mathbf{y}})
\end{align}
\EndKnitrBlock{proposition}
<br>

\BeginKnitrBlock{proof}
\iffalse{} <span class="proof"><em>Prueba. </em></span>  \fi{}En primer lugar, nótese la siguiente identidad
\begin{equation}
\sum_{i=1}^n(\mathbf{Y}_i-\btheta)'\bSigma^{-1}(\mathbf{Y}_i-\btheta)
=\sum_{i=1}^n(\mathbf{Y}_i-\bar{\mathbf{Y}})'\bSigma^{-1}(\mathbf{Y}_i-\bar{\mathbf{Y}})
+n(\bar{\mathbf{Y}}-\btheta)'\bSigma^{-1}(\bar{\mathbf{Y}}-\btheta)
\end{equation}

puesto que
\begin{align*}
&\ \ \ \ \sum_{i=1}^n(\mathbf{Y}_i-\btheta)'\bSigma^{-1}(\mathbf{Y}_i-\btheta)\\
&=\sum_{i=1}^n(\mathbf{Y}_i-\bar{\mathbf{Y}}+\bar{\mathbf{Y}}-\btheta)'
\bSigma^{-1}(\mathbf{Y}_i-\bar{\mathbf{Y}}+\bar{\mathbf{Y}}-\btheta)\\
&=\sum_{i=1}^n(\mathbf{Y}_i-\bar{\mathbf{Y}})'\bSigma^{-1}(\mathbf{Y}_i-\bar{\mathbf{Y}})
+\sum_{i=1}^n(\mathbf{Y}_i-\bar{\mathbf{Y}})'\bSigma^{-1}(\bar{\mathbf{Y}}-\btheta)\\
&\hspace{1.5cm}
+(\bar{\mathbf{Y}}-\btheta)'\bSigma^{-1}\sum_{i=1}^n(\mathbf{Y}_i-\bar{\mathbf{Y}})'
+\sum_{i=1}^n(\bar{\mathbf{Y}}-\btheta)'\bSigma^{-1}(\bar{\mathbf{Y}}-\btheta)\\
&=\sum_{i=1}^n(\mathbf{Y}_i-\bar{\mathbf{Y}})'\bSigma^{-1}(\mathbf{Y}_i-\bar{\mathbf{Y}})
+n(\bar{\mathbf{Y}}-\btheta)'\bSigma^{-1}(\bar{\mathbf{Y}}-\btheta)
\end{align*}

Por otro lado, de la definición de distribución previa, se tiene que
\begin{align*}
p(\btheta \mid \mathbf{Y},\bSigma)
&\propto p(\mathbf{Y} \mid \btheta,\bSigma)p(\btheta,\bSigma)\\
&\propto \exp\left\{-\frac{1}{2}\left[\sum_{i=1}^n(\mathbf{y}_i-\btheta)'\bSigma^{-1}(\mathbf{y}_i-\btheta)+(\btheta-\bmu)'\bGamma^{-1}(\btheta-\bmu)\right]\right\}\\
&\propto \exp\left\{-\frac{1}{2}\left[n(\bar{\mathbf{y}}-\btheta)'\bSigma^{-1}(\bar{\mathbf{y}}-\btheta)+(\btheta-\bmu)'\bGamma^{-1}(\btheta-\bmu)\right]\right\}\\
&\propto \exp\left\{-\frac{1}{2}\left[
  -n\bar{\mathbf{y}}'\bSigma^{-1}\btheta-n\btheta'\bSigma^{-1}\bar{\mathbf{y}}+n\btheta'\bSigma^{-1}\btheta+\btheta'\bGamma^{-1}\btheta-\btheta'\bGamma^{-1}\bmu-\bmu'\bGamma^{-1}\btheta\right]\right\}\\
&= \exp\left\{-\frac{1}{2}\left[
  \btheta'(\bGamma^{-1}+n\bSigma^{-1})\btheta-2\btheta'(\bGamma^{-1}\bmu+n\bSigma^{-1}\bar{\mathbf{y}})\right]\right\}\\
&= \exp\left\{-\frac{1}{2}\left[\btheta'\bGamma^{-1}_n\btheta-2\btheta'\bGamma^{-1}_n\bmu_n\right]\right\}\\
&\propto \exp\left\{-\frac{1}{2}\left[\btheta'\bGamma^{-1}_n\btheta-2\btheta'\bGamma^{-1}_n\bmu_n+\bmu_n\bGamma_n^{-1}\bmu_n\right]\right\}\\
&= \exp\left\{-\frac{1}{2}(\btheta-\bmu_n)'\bGamma_n^{-1}(\btheta-\bmu_n)\right\}
  \end{align*}

La cual corresponde al núcleo de una distribución normal multivariante con vector de medias $\bmu_n$ y matriz de varianzas $\bGamma_n$.
\EndKnitrBlock{proof}
<br>

Observando los parámetros de la distribución posterior, podemos ver que $\bGamma_n^{-1} = \bGamma^{-1}+(\bSigma/n)^{-1}$. Teniendo en cuenta que la matriz de varianzas y covarianzas es una medida de dispersión de la distribución alrededor de su media, la inversa de dicha matriz se puede ver como una medida de precisión de qué tanto se concentra la distribución alrededor de la media. Así, podemos ver que la precisión posterior viene siendo la suma entre la precisión previa y la precisión de la estimación clásica del parámetro $\btheta$.

En cuanto a la media posterior $\bmu_n$, tenemos que
\begin{align*}
\bmu_n&=\left(\bGamma^{-1}+n\bSigma^{-1}\right)^{-1}
(\bGamma^{-1}\bmu+n \bSigma^{-1}\bar{\mathbf{y}})\\
&=\left(\mathbf{I}+n\bGamma\bSigma^{-1}\right)^{-1}\bmu+\left(\frac{1}{n}\bSigma\bGamma^{-1}+\mathbf{I}\right)^{-1}\bar{\mathbf{y}}\\
&=\underbrace{\bSigma\left(\bSigma+n\bGamma\right)^{-1}}_{\mathbf{A}_1}\bmu+\underbrace{n\bGamma\left(\bSigma+n\bGamma\right)^{-1}}_{\mathbf{A}_2}\bar{\mathbf{y}}
\end{align*}

De donde podemos que ver que la media posterior $\bmu_n$ se puede escribir como $\bmu_n=\mathbf{A}_1\bmu+\mathbf{A}_2\bar{\mathbf{y}}$ donde $\mathbf{A}_1+\mathbf{A}_2=\mathbf{I}$. Es claro que en el caso univariado, $\mathbf{A}_1$, $\bmu$, $\mathbf{A}_2$ y $\bar{\mathbf{y}}$ son todos escalares, y $\bmu_n$ es un valor intermedio entre $\bmu$ y $\bar{\mathbf{y}}$. Mientras que en caso multivariado, $\mathbf{A}_1\bmu+\mathbf{A}_2\bar{\mathbf{y}}$ es similar a una combinación convexa entre los vectores $\bmu$ y $\bar{\mathbf{y}}$, pero los coeficientes son matrices en vez de escalares. 

Para ilustrar la relación de $\bmu_n$ con $\bmu$ y $\bar{\mathbf{y}}$, tomamos el caso de $p=2$, y denotamos $\mathbf{A}=\bSigma\left(\bSigma+n\bGamma\right)^{-1}$. Es claro que $\mathbf{A}$ es una matriz simétrica y definida positiva, lo denotaremos con $\mathbf{A}_1=\begin{pmatrix}a_{11}&a_{12}\\ a_{12}&a_{22}\end{pmatrix}$, donde $a_{11}>0$ y $a_{22}>0$. De esta forma
\begin{align*}
\bmu_n&=\begin{pmatrix}a_{11}&a_{12}\\ a_{12}&a_{22}\end{pmatrix}\begin{pmatrix}\mu_{1}\\ \mu_{2}\end{pmatrix}+\begin{pmatrix}1-a_{11}&-a_{12}\\ -a_{12}&1-a_{22}\end{pmatrix}\begin{pmatrix}\bar{y}_{1}\\ \bar{y}_{2}\end{pmatrix}\\
&=\begin{pmatrix}a_{11}\mu_1+(1-a_{11})\bar{y}_{1}+a_{12}(\mu_2-\bar{y}_{2})\\a_{22}\mu_2+(1-a_{22})\bar{y}_{2}+a_{12}(\mu_1-\bar{y}_{1})\end{pmatrix}
\end{align*}

Al observar la primera entrada de $\bmu_n$, podemos ver que este se compone de una combinación convexa entre $\mu_1$ y $\bar{y}_1$ (pues $a_{11}>0$) y una parte que depende de la diferencia $\mu_2-\bar{y}_{2}$; un comportamiento similar se observa en la segunda entrada de $\bmu_n$. Esta observación es interesante, pues ilustra que cada componente de la media posterior $\bmu_n$ no siempre será un promedio ponderado del componentes correspondiente de la media previa y la estimación clásica.

Ilustramos los resultados encontrados suponiendo las dos siguientes situaciones:

1. Suponga que se quiere estimar el vector de medias $\btheta=(\theta_1,\theta_2)'$ con una matriz de varianzas y covarianzas conocida de $\bSigma=\begin{pmatrix}20&8\\ 8&30\end{pmatrix}$. Para esto se recolectan 10 observaciones con vector de promedios muestrales $\bar{\mathbf{y}}=(150,230)'$. Como información previa, suponga que $\bmu=(100, 200)'$ y $\bGamma=\begin{pmatrix}5&3\\ 3&10\end{pmatrix}$. Realizando los cálculos correspondiente, se tiene que $\bGamma_n=\begin{pmatrix}1.42&0.64\\ 0.64&2.31\end{pmatrix}$, $\mathbf{A}=\begin{pmatrix}0.3&-0.026\\ -0.026&0.235\end{pmatrix}$ y $\bmu_n=(136,224)$. Podemos ver que en este caso cada componente de $\bmu_n$ se encuentra entre los componentes correspondientes de $\bmu$ y $\bar{\mathbf{y}}$. Los códigos computacionales se muestran a continuación.


```r
n <- 10
mu <- matrix(c(100, 200))
y.bar <- matrix(c(150, 230))
Gamma <- matrix(c(5, 3, 3, 10), 2, 2)
Sigma <- matrix(c(20, 8, 8, 30), 2, 2)
Gamma.n <- solve(solve(Gamma) + n * solve(Sigma))
A <- Sigma %*% solve(Sigma + n * Gamma) 
mu.n <- Gamma.n %*% 
  (solve(Gamma) %*% mu + n * solve(Sigma) %*% y.bar)
mu.n
```




\begin{tabular}{r}
\hline
135.7889\\
\hline
223.6155\\
\hline
\end{tabular}


2. Tomando los mismos datos del caso anterior, también podemos suponer que $\bar{\mathbf{y}}=(150, 2300)'$, las matrices $\bGamma_n$ y $\mathbf{A}$ no cambian de valor, pero la media de la distribución posterior está dada por $\bmu_n=(190,1808)$. Podemos ver que el primer componente de $\bmu_n$ no está entre 100 y 150 que corresponden a las estimaciones previa y clásica, respectivamente. Esto se debe a que la diferencia entre $\mu_2$ y $\bar{y}_2$ es muy grande.


```r
n <- 10
mu <- matrix(c(100, 200))
y.bar <- matrix(c(150, 2300))
Gamma <- matrix(c(5, 3, 3, 10), 2, 2)
Sigma <- matrix(c(20, 8, 8, 30), 2, 2)
Gamma.n <- solve(solve(Gamma) + n * solve(Sigma))
A <- Sigma %*% solve(Sigma + n * Gamma) 
mu.n <- Gamma.n %*% 
  (solve(Gamma) %*% mu + n * solve(Sigma) %*% y.bar)
mu.n
```




\begin{tabular}{r}
\hline
189.8642\\
\hline
1808.0199\\
\hline
\end{tabular}

### Distribución previa no informativa

Al tener en cuenta que la distribución previa del parámetro $\btheta$ es la distribución normal multivariada, y al observar la forma de la función de densidad, se puede afirmar que cuando $|\bGamma^{-1}|$ es muy pequeño, los parámetros previos $\bmu$ y $\bGamma$ pierden peso en los cálculos de $\bmu_n$ y $\bGamma_n$. En este caso se puede ver que
\begin{align*}
\bGamma_n&\approx n^{-1}\bSigma\\
\bmu_n&\approx \bar{\mathbf{y}}
\end{align*}

De donde podemos concluir que la estimación bayesiana será muy cercana a la estimación clásica $\bar{y}$, más aún, el intervalo de credibilidad también será muy similar al intervalo de confianza del enfoque clásico.

\BeginKnitrBlock{example}
<span class="example" id="exm:EjeStudent"><strong>(\#exm:EjeStudent) </strong></span>@Student introdujo un conjunto de datos clásicos sobre el incremento en horas de sueño producido con 2 medicamentos soporíferos diferentes comparados con grupo control en 10 pacientes. Estos datos se pueden encontrar en `R` con el nombre `sleep` y se pueden definir como realizaciones de vectores aleatorios con distribución normal bivariada. Supongamos que la matriz de varianzas y covarianzas de la distribución es conocida e igual a $\Sigma=\begin{pmatrix}1&0.6\\ 0.6&2\end{pmatrix}$.
  
El parámetro de interés es el vector de medias $\btheta=(\theta_1,\theta_2)'$. Para la distribución previa, suponemos que $\bmu=(0,1)'$, es decir que el primer medicamento no tiene ningún efecto soporífero, mientras que el segundo medicamento tiene un efecto promedio de aumentar 1 hora de sueño, también asumimos que $\Gamma=\begin{pmatrix}2&0\\ 0&2\end{pmatrix}$. Los siguientes códigos de `STAN` ilustran el procedimiento de estimación del parámetro de interés.
\EndKnitrBlock{example}


```r
NormalMultMedia <- '
data {
  int<lower=0> n;
  int<lower=0> P;
  vector[P] y[n];
  vector[P] mu;
  matrix[P, P] Sigma;
  matrix[P, P] Gamma;
}
parameters {
  vector[P] theta;
}
transformed parameters {
  real diftheta;
  diftheta = theta[2] - theta[1];
}
model {
  y ~ multi_normal(theta, Sigma);
  theta ~ multi_normal(mu, Gamma);
}
'

y <- structure(.Data = sleep[,1], .Dim=c(10,2))
n <- nrow(y)
P <- ncol(y)
Sigma  <- matrix(c(1, 0.6, 0.6, 2), 2, 2)
mu <- as.vector(c(0, 1))
Gamma <- matrix(c(2, 0, 0, 2), 2, 2)

sample_data <- list(y = y, n = n, P = P,
                    Sigma = Sigma, mu = mu,
                    Gamma = Gamma)
set.seed(1234)
NormalMultifit <- stan(model_code = NormalMultMedia,
                   data = sample_data, verbose = FALSE)
```

Después de la convergencia del proceso inferencial, la estimación bayesiana de $\btheta$ es $(0.6866, 2.1812)'$ para los dos medicamentos; mientras que los intervalos de crediblidad del 95\% corresponden a (0.0834, 1.2776) y (1.3391, 3.0109).


```r
print(NormalMultifit, digits = 4, 
      pars = "theta", probs = c(0.025, 0.975))
```

```
## Inference for Stan model: 0b2360136b3d61c08a460b98848640c6.
## 4 chains, each with iter=2000; warmup=1000; thin=1; 
## post-warmup draws per chain=1000, total post-warmup draws=4000.
## 
##            mean se_mean     sd   2.5%  97.5% n_eff   Rhat
## theta[1] 0.6874  0.0064 0.3110 0.1038 1.3055  2361 0.9997
## theta[2] 2.2060  0.0090 0.4283 1.3620 3.0268  2245 1.0000
## 
## Samples were drawn using NUTS(diag_e) at Sun Jun 27 21:00:32 2021.
## For each parameter, n_eff is a crude measure of effective sample size,
## and Rhat is the potential scale reduction factor on split chains (at 
## convergence, Rhat=1).
```

Las figuras \@ref(fig:posNormalMultiStan) muestra la distribución posterior para este ejemplo, junto con la estimación puntual, correspondiente a la media. 


```r
bayesplot::mcmc_areas(NormalMultifit, pars = c("theta[1]", "theta[2]"), 
                      prob = 0.95)
```

\begin{figure}

{\centering \includegraphics{4Multiparametricos_files/figure-latex/posNormalMultiStan-1} 

}

\caption{Distribución posterior.}(\#fig:posNormalMultiStan)
\end{figure}

A continuación, mostramos los códigos de `R` para llevar a cabo los cálculos directamente.


```r
y.bar <- colMeans(y)
Gamma.n <- solve(solve(Gamma) + n * solve(Sigma))
mu.n <- Gamma.n %*% 
  (solve(Gamma) %*% mu + n * solve(Sigma) %*% y.bar)
mu.n
```




\begin{tabular}{r}
\hline
0.6802703\\
\hline
2.1905381\\
\hline
\end{tabular}

```r
Gamma.n
```




\begin{tabular}{r|r}
\hline
0.0937527 & 0.0519886\\
\hline
0.0519886 & 0.1804003\\
\hline
\end{tabular}

De los resultados arrojados, vemos que la distribución posterior del parámetro está dada por
\begin{equation*}
\begin{pmatrix}
\theta_1\\
\theta_2
\end{pmatrix}
\sim N_2\left(\begin{pmatrix}
0.68\\
2.19
\end{pmatrix},\begin{pmatrix}
0.094&0.052\\
0.052&0.180
\end{pmatrix}\right)
\end{equation*}

De esta forma, la estimación bayesiana obtenida para los efectos promedios corresponde a 0.68 horas y 2.19 horas, respectivamente; los cuales son similares a los obtenidos por `STAN`. En cuanto a los intervalos de credibilidad del 95\%, estas son dados por los percentiles 2.5\% y 97.5\% de las dos distribuciones posteriores marginales de $\theta_1$ y $\theta_2$. Estos intervalos se pueden obtener así:


```r
qnorm(c(0.025, 0.975), mu.n[1], sqrt(Gamma.n[1, 1]))
```

```
## [1] 0.08014771 1.28039297
```

```r
qnorm(c(0.025, 0.975), mu.n[2], sqrt(Gamma.n[2, 2]))
```

```
## [1] 1.358072 3.023005
```

Ahora, suponga que el objetivo es comparar los medicamentos para concluir si el segundo medicamento es más efectivo que el primero, podemos encontrar la distribución posterior de la diferencia $\theta_2-\theta_1$. Utilizando propiedades de la distribución normal multivariante, podemos encontrar la distribución posterior de $\theta_2-\theta_1$, calcular un intervalo de credibilidad para $\theta_2-\theta_1$ e indagar cuál es la probabilidad de que $\theta_2$ sea mayor a $\theta_1$. Estos cálculos se pueden llevar a cabo de la siguiente forma en `STAN`


```r
print(NormalMultifit, digits = 4, 
      pars = "diftheta", probs = c(0.025, 0.975))
```

```
## Inference for Stan model: 0b2360136b3d61c08a460b98848640c6.
## 4 chains, each with iter=2000; warmup=1000; thin=1; 
## post-warmup draws per chain=1000, total post-warmup draws=4000.
## 
##            mean se_mean     sd   2.5%  97.5% n_eff   Rhat
## diftheta 1.5186  0.0067 0.4136 0.7162 2.3123  3815 0.9998
## 
## Samples were drawn using NUTS(diag_e) at Sun Jun 27 21:00:32 2021.
## For each parameter, n_eff is a crude measure of effective sample size,
## and Rhat is the potential scale reduction factor on split chains (at 
## convergence, Rhat=1).
```

La figuras \@ref(fig:posNormalDifMultiStan) muestra la distribución posterior para este ejemplo, junto con la estimación puntual, correspondiente a la media. 


```r
bayesplot::mcmc_areas(NormalMultifit, pars = "diftheta", 
                      prob = 0.95)
```

\begin{figure}

{\centering \includegraphics{4Multiparametricos_files/figure-latex/posNormalDifMultiStan-1} 

}

\caption{Distribución posterior de los parámetros transformados (diferencia de las medias teóricas).}(\#fig:posNormalDifMultiStan)
\end{figure}

Los mismos cálculos pueden reproducirse en `R`, por medio de los siguientes códigos computacionales.


```r
vec <- matrix(c(-1, 1), 1, 2)
difmedia <- vec %*% mu.n
varianza <- vec %*% Gamma.n %*% t(vec)
qnorm(c(0.025, 0.975), difmedia, sqrt(varianza))
```

```
## [1] 0.7017359 2.3187996
```

```r
1 - pnorm(0, difmedia, sqrt(varianza))
```

```
## [1] 0.9998744
```

Observando los anteriores resultados, vemos que el intervalo de credibildad para $\theta_2-\theta_1$ está dado por $(0.7, 2.3)$, el cual no contiene el valor 0, indicando que el segundo medicamento tiene un efecto mayor que el segundo. Adicionalmente, se observa que con probabilidad muy cercana a uno, el segundo medicamento tiene efecto mayor al primero y por ende tiene un desempeño superior al primero.

Finalmente, ilustramos los resultados obtenidos al usar una distribución previa no informativa, para eso, usaremos $\bGamma=\begin{pmatrix}100&0\\ 0&100\end{pmatrix}$, con $|\bGamma^{-1}|=0.0001$, representando una distribución previa no informativa. Los resultados de estimación arroja la siguiente distribución posterior para el vector de parámetros:
\begin{equation*}
\begin{pmatrix}
\theta_1\\
\theta_2
\end{pmatrix}
\sim N_2\left(\begin{pmatrix}
0.75\\
2.33
\end{pmatrix},\begin{pmatrix}
0.10&0.06\\
0.06&0.20
\end{pmatrix}\right)
\end{equation*}

Los intervalos de credibilidad del 95\% para los parámetros $\theta_1$ y $\theta_2$ están dados por $(0.129, 1.368)$ y $(1.451, 3.202)$, respectivamente. Observamos que estos intervalos de credibilidad son muy similares a los intervalos de confianza del 95\% del enfoque clásico, calculados con la expresión $\bar{y}\pm z_{1-\alpha/2}\frac{\sigma}{\sqrt{n}}$, y que para los datos del ejemplo están dados por $(0.130, 1.370)$ y $(1.453, 3.207)$.

Finalmente, recordamos los dos siguientes resultados relacionados con la distribución normal multivariante que pueden resultar útiles en otros análisis.

\BeginKnitrBlock{proposition}
<span class="proposition" id="prp:unnamed-chunk-33"><strong>(\#prp:unnamed-chunk-33) </strong></span>La distribución posterior marginal de un subconjunto de parámetros, digamos $\btheta^{(1)}$ es también normal multivariante con media igual a la del subvector de medias apropiado, $\bmu_n^{(1)}$ y similar matriz de varianzas $\bGamma_n^{(11)}$.
\EndKnitrBlock{proposition}
  
\BeginKnitrBlock{proposition}
<span class="proposition" id="prp:unnamed-chunk-34"><strong>(\#prp:unnamed-chunk-34) </strong></span>La distribución posterior condicional de un subconjunto de parámetros, digamos $\btheta^{(1)}$, dado $\btheta^{(2)}$ es también normal multivariante dada por
\begin{equation*}
  \btheta^{(1)} \mid \btheta^{(2)} \sim N_p \left(\bmu_n^{(1)}+\bGamma_n^{(12)}\left(\bGamma_n^{(22)}\right)^{-1}
  \left(\theta^{(2)}-\mu_n^{(2)}\right),\bGamma_n^{(1 \mid 2)}\right).
\end{equation*}

En donde
\begin{align}
  \bGamma_n^{(1 \mid 2)}&= \bGamma_n^{(11)}-\bGamma_n^{(12)}\left(\bGamma_n^{(22)}\right)^{-1}\bGamma_n^{(21)}
\end{align}

Con $\bmu^{(1)}$ y $\bmu^{(2)}$ correspondendientes al vector de medias y $\bGamma_n^{(11)}$, $\bGamma_n^{(22)}$ denotan la matriz de varianzas y covarianzas de $\btheta^{(1)}$ y $\btheta^{(2)}$, respectivamente. $\bGamma_n^{(12)}$ es la matriz de covarianzas entre $\btheta^{(1)}$ y $\btheta^{(2)}$, $\bGamma_n^{(21)}$ es la matriz de covarianzas entre $\btheta^{(2)}$ y $\btheta^{(1)}$. 
\EndKnitrBlock{proposition}
  
La prueba de los dos resultados anteriores se sigue inmediatamente de las propiedades de la distribución normal multivariante.

## Modelo normal multivariante con media y varianza desconocida
  
Al igual que en la distribución normal univariada, cuando se desconoce tanto el vector de medias como la matriz de varianzas y covarianzas de la distribución, es necesario plantear diversos enfoques y situarse en el más conveniente. Nótese que en términos de parámetros, existen $p$ parámetros correspondientes al vector de medias $\btheta$ y $\binom{p}{2}=\dfrac{p(p+1)}{2}$ parámetros correspondientes a la matriz de varianzas $\bSigma$. Pensando en la gran cantidad de parámetros que se deben modelar, es necesario tener en cuenta que el número de datos en la muestra aleatoria sea lo suficientemente grande. Suponiendo que el número de observaciones en la muestra aleatoria sea suficiente, existe otra situación que se debe surtir y es la asignación de las distribuciones previas para $\btheta$ y $\bSigma$. En estos términos, es posible
  
1. Suponer que la distribución previa $p(\btheta)$ es independiente de la distribución previa $p(\bSigma)$ y que ambas distribuciones son informativas. Luego, utilizar un análisis de simulación condicional conjunta para extraer muestras provenientes de las respectivas distribuciones posterior.
1. Suponer que la distribución previa para $\btheta$ depende de $\bSigma$ y escribirla como $p(\btheta \mid \bSigma)$, mientras que la distribución previa de $\bSigma$ no depende de $\btheta$ y se puede escribir como $p(\bSigma)$. El análisis posterior de este enfoque encuentra la distribución posterior de $\bSigma \mid \mathbf{Y}$ y con esta se encuentra la distribución posterior de $\btheta \mid \bSigma,\mathbf{Y}$.
1. Suponer que la distribución conjunta previa para $\btheta$ y $\bSigma$ es una distribución no informativa.


### Parámetros independientes con distribuciones previas informativas

En este enfoque se supone que las distribuciones previas para los parámetros de interés son independientes e informativas. Hacemos siguiente observación para lograr que las resultantes distribuciones posterior sean conjugadas. 
\begin{align*}
\sum_{i=1}^n(\mathbf{Y}_i-\theta)'\bSigma^{-1}(\mathbf{Y}_i-\theta)&=traza \left(\sum_{i=1}^n(\mathbf{Y}_i-\theta)'\bSigma^{-1}(\mathbf{Y}_i-\theta)\right)\\
&= \sum_{i=1}^n traza\left((\mathbf{Y}_i-\theta)'\bSigma^{-1}(\mathbf{Y}_i-\theta)\right)\\
&= \sum_{i=1}^n traza\left(\bSigma^{-1}(\mathbf{Y}_i-\theta)(\mathbf{Y}_i-\theta)'\right)\\
&= traza\left(\bSigma^{-1}\sum_{i=1}^n(\mathbf{Y}_i-\theta)(\mathbf{Y}_i-\theta)'\right)\\
&= traza\left(\bSigma^{-1}\mathbf{S}_{\btheta}\right)
\end{align*}

Donde $\mathbf{S}_{\btheta}=\sum_{i=1}^n(\mathbf{Y}_i-\btheta)(\mathbf{Y}_i-\btheta)'$. En cuanto a la asignación de las distribuciones previas, para el vector de medias $\btheta$ es posible usar la distribución normal, esto es,
\begin{equation*}
\btheta \sim Normal_p(\bmu,\bGamma)
\end{equation*}

Por otro lado, la distribución para la matriz de varianzas $\bSigma$ es
\begin{equation*}
\bSigma \sim Inversa-Wishart(\bLambda,v)
\end{equation*}

donde $v$ denota los grados de libertad y $\bLambda$ la matriz de escala. Esto es, la función de densidad está dada por
\begin{equation*}
p(\bSigma)\propto |\bSigma|^{-\frac{v+p+1}{2}}\exp\left\{-\frac{1}{2}traza(\bLambda\bSigma^{-1})\right\}
\end{equation*}

Asumiendo independencia previa, la distribución previa conjunta resulta estar dada por
\begin{align}
p(\btheta,\bSigma)&=p(\btheta)p(\bSigma)\notag\\
&\propto \mid \bSigma \mid ^{-(v+p+1)/2}\notag\\
&\times
\exp\left\{ -\frac{1}{2}\left[traza(\bLambda\bSigma^{-1})+
(\btheta-\bmu)'\bGamma^{-1}(\btheta-\bmu)\right]\right\}
\end{align}


Una vez que se conoce la forma estructural de la distribución previa conjunta, es posible establecer la distribución posterior conjunta teniendo en cuenta la forma de la función de verosimilitud $p(\mathbf{Y} \mid \btheta,\bSigma)$ y la expresión equivalente para $\sum_{i=1}^n(\mathbf{Y}_i-\btheta)'\bSigma^{-1}(\mathbf{Y}_i-\btheta)$ mostrada al inicio de esta sección. Adicionalmente, acudiendo a la simetría de las matrices $\bLambda$, $\bSigma$ y $\mathbf{S}_{\btheta}$, se tiene que
\begin{align}
p(\btheta,\bSigma \mid \mathbf{Y})&\propto p(\btheta,\bSigma)p(\mathbf{Y} \mid \btheta,\bSigma)\notag\\
&\propto \mid \bSigma \mid ^{-(v+n+p+1)/2}\notag\\
&\times
\exp\left\{ -\frac{1}{2}\left[traza(\bLambda\bSigma^{-1}+\bSigma^{-1}\mathbf{S}_{\btheta})+
                                (\btheta-\bmu)'\bGamma^{-1}(\btheta-\bmu)\right]\right\}\notag\\
                              &\propto \mid \bSigma \mid ^{-(v+n+p+1)/2}\notag\\
                              &\times
                              \exp\left\{ -\frac{1}{2}\left[traza(\bSigma^{-1}(\bLambda+\mathbf{S}_{\btheta}))+
                              (\btheta-\bmu)'\bGamma^{-1}(\btheta-\bmu)\right]\right\}
\end{align}

Dado que la distribución posterior conjunta no tiene una forma estructural conocida, no es posible utilizar el método de integración analítica. Sin embargo, es posible obtener las distribuciones condicionales de cada uno de los parámetros suponiendo fijos los restantes y teniendo en cuenta que
\begin{align*}
p(\btheta \mid \bSigma,\mathbf{Y})\propto p(\btheta,\underbrace{\bSigma}_{fijo} \mid \mathbf{Y})
\ \ \ \ \ \ \ \ \ \text{y} \ \ \ \ \ \ \ \ \ \
p(\bSigma \mid \btheta,\mathbf{Y})\propto p(\underbrace{\btheta}_{fijo},\bSigma \mid \mathbf{Y})
\end{align*}

\BeginKnitrBlock{proposition}
<span class="proposition" id="prp:unnamed-chunk-35"><strong>(\#prp:unnamed-chunk-35) </strong></span>La distribución posterior de la matriz de parámetros $\bSigma$ condicional a $\btheta,\mathbf{Y}$ es
\begin{equation*}
(\#eq:PosCondbSigma)
\bSigma \mid \btheta,\mathbf{Y} \sim Inversa-Wishart_{v+n}(\bLambda+\mathbf{S}_{\btheta})
\end{equation*}
\EndKnitrBlock{proposition}
<br>

\BeginKnitrBlock{proof}
\iffalse{} <span class="proof"><em>Prueba. </em></span>  \fi{}La prueba es inmediata notando que
\begin{align*}
\bSigma \mid \btheta,\mathbf{Y}&\propto\mid \bSigma \mid ^{-(v+n+p+1)/2}\notag\\
&\times\exp\left\{ -\frac{1}{2}\left[traza(\bSigma^{-1}(\bLambda+\mathbf{S}_{\btheta}))+(\btheta-\bmu)'\bGamma^{-1}(\btheta-\bmu)\right]\right\}
\end{align*}

Por lo tanto, factorizando convenientemente, se encuentra una expresión idéntica a la función de distribución de una variable aleatoria con distribución $Inversa-Wishart_{v+n}(\bLambda+\mathbf{S}_{\btheta})$.
\EndKnitrBlock{proof}
<br>
                          
\BeginKnitrBlock{proposition}
<span class="proposition" id="prp:unnamed-chunk-37"><strong>(\#prp:unnamed-chunk-37) </strong></span>La distribución posterior del vector de parámetros $\btheta$ condicional a $\bSigma,\mathbf{Y}$ es
\begin{equation}
(\#eq:PosCondbtheta)
\btheta \mid  \bSigma,\mathbf{Y} \sim Normal_p(\bmu_n,\bGamma_n)
\end{equation}
  
donde $\bmu_n$ y $\bGamma_n$ están dadas por las expresiones \@ref(eq:Gamman) y \@ref(eq:mun), respectivamente.
\EndKnitrBlock{proposition}
<br>

Una vez encontradas las distribuciones posteriores condicionales de $\btheta$ y $\bSigma$, se puede obtener la estimación de estos parámetros vía el muestreo de Gibbs, que en este caso se resume en los siguientes pasos:

1. Fijar un valor inicial para $\btheta$; lo denotamos por $\btheta_{(1)}$.
1. Simular un valor de la distribución de $\bSigma|\btheta,\mathbf{Y}$ en \@ref(eq:PosCondbSigma) donde el parámetro $\mathbf{S_\mathbf{\btheta}}$ que depende de $\btheta$, debe ser reemplazado por $\btheta_{(1)}$ del paso anterior; este valor simulado se denotará por $\bSigma_{(1)}$
1. Simular un valor de la distribución de $\btheta|\bSigma,\mathbf{Y}$ en \@ref(eq:PosCondbtheta) donde en $\mathbf{mu}_n$ y $\bGamma_n$ se debe reemplazar $\bSigma$ por $\bSigma$; este valor simulado se denota por $\btheta$.
1. Repetir los pasos (2) y (3) hasta completar un número de iteraciones suficientes para alcanzar la convergencia en ambos parámetros


Una vez tengamos los valores muestreados, se debe garantizar la convergencia y la correlación nula entre estos valores, con el fin de calcular las estimaciones. En el siguiente ejemplo ilustramos la implementación de este muestreo de Gibbs en `R`.

\BeginKnitrBlock{example}
<span class="example" id="exm:EjeStudent2"><strong>(\#exm:EjeStudent2) </strong></span>Retomamos los datos del efecto de dos medicamentos soporíferos introducidos por @Student, los cuales fueron estudiados en el ejemplo \@ref(exm:EjeStudent) asumiendo que la matriz de varianzas y covarianzas era conocida. El vector de medias muestrales de estos datos están dados por $\bar{y}=(0.75, 2.33)'$, y la matriz de varianzas y covarianzas muestrales está dada por $\mathbf{S}=\begin{pmatrix}3.20&2.85\\2.85&4.01\end{pmatrix}$. 

Ahora supongamos que tanto el vector de medias como la matriz de varianzas y covarianzas son desconocidos. Para el vector de medias, asumimos la distribución previa del ejemplo \ref{EjeStudent}, es decir, $\bmu=(0, 1)'$ y $\bGamma=\begin{pmatrix}2&0\\0&2\end{pmatrix}$. Para la matriz de varianzas y covarianzas asumimos la distribución inversa-Wishart con matriz de escala igual a $\bLambda=\begin{pmatrix}20&8\\8&20\end{pmatrix}$ y $v=10$ grados de libertad. De esta forma, la estimación previa de $\bSigma$ viene dada por $\frac{1}{v-2-1}\bLambda=\begin{pmatrix}2.86&1.14\\ 1.14&2.86\end{pmatrix}$. 

Ilustramos los códigos de `STAN` a continuación.
\EndKnitrBlock{example}


```r
NormalMultMediaCov <- '
data {
  int<lower=0> n;
  int<lower=0> P;
  vector[P] y[n];
  vector[P] mu;
  matrix[P, P] Gamma;
  matrix[P, P] Lambda;
  int<lower=0> v;
}
parameters {
  vector[P] theta;
  cov_matrix[P] Sigma;
}
transformed parameters {
  real diftheta;
  diftheta = theta[2] - theta[1];
}
model {
  theta ~ multi_normal(mu, Gamma);
  Sigma ~ inv_wishart(v, Lambda);
  for (i in 1:n)
    y[i] ~ multi_normal(theta, Sigma);
}
'

y <- structure(.Data = sleep[,1], .Dim=c(10,2))
n <- nrow(y)
P <- ncol(y)
Sigma  <- matrix(c(1, 0.6, 0.6, 2), 2, 2)
mu <- as.vector(c(0, 1))
Gamma <- matrix(c(2, 0, 0, 2), 2, 2)
v <- 10
Lambda <- matrix(c(20, 8, 8, 20), 2, 2)

sample_data <- list(y = y, n = n, P = P,
                    Sigma = Sigma, mu = mu,
                    Gamma = Gamma, v = v,
                    Lambda = Lambda)
set.seed(1234)
NormalMultMediaCovfit <- stan(model_code = NormalMultMediaCov,
                   data = sample_data, verbose = FALSE)
```

Con base en los resultados de las cadenas simuladas, se observa que la estimación bayesiana para el número de horas de sueño producidas por los dos medicamentos son 0.5772 y 2.1011, respectivamente. En cuanto a la estimación de la matriz de varianzas y covarianzas, ésta está dada por $\hat{\bSigma}=\begin{pmatrix}3.0132&2.1011\\2.1011&3.5222\end{pmatrix}$.


```r
print(NormalMultMediaCovfit, digits = 4, 
      pars = c("theta", "Sigma"), probs = c(0.025, 0.975))
```

```
## Inference for Stan model: 8ef9836716ef08fd4cd1136cfef43904.
## 4 chains, each with iter=2000; warmup=1000; thin=1; 
## post-warmup draws per chain=1000, total post-warmup draws=4000.
## 
##              mean se_mean     sd    2.5%  97.5% n_eff   Rhat
## theta[1]   0.5482  0.0103 0.5001 -0.4113 1.5300  2337 1.0001
## theta[2]   2.0732  0.0122 0.5480  0.8840 3.0624  2033 1.0008
## Sigma[1,1] 3.0355  0.0241 1.1310  1.5812 5.8057  2206 0.9999
## Sigma[1,2] 2.1079  0.0232 1.0193  0.7354 4.7245  1928 1.0002
## Sigma[2,1] 2.1079  0.0232 1.0193  0.7354 4.7245  1928 1.0002
## Sigma[2,2] 3.5141  0.0280 1.3193  1.8198 6.8809  2214 0.9995
## 
## Samples were drawn using NUTS(diag_e) at Sun Jun 27 21:01:20 2021.
## For each parameter, n_eff is a crude measure of effective sample size,
## and Rhat is the potential scale reduction factor on split chains (at 
## convergence, Rhat=1).
```

Las figuras \@ref(fig:posNormalMVCov1) y \@ref(fig:posNormalMVCov2) muestran las distribuciones posteriores en este ejemplo. 

\begin{figure}

{\centering \includegraphics{4Multiparametricos_files/figure-latex/posNormalMVCov1-1} 

}

\caption{Distribuciones posteriores para el vector de medias y la diferencia de medias.}(\#fig:posNormalMVCov1)
\end{figure}
\begin{figure}

{\centering \includegraphics{4Multiparametricos_files/figure-latex/posNormalMVCov2-1} 

}

\caption{Distribuciones posteriores para los elementos de la matriz de covarinzas.}(\#fig:posNormalMVCov2)
\end{figure}

A continuación se muestran los códigos necesarios para implementar el muestreo de Gibbs de forma manual en `R`.


```r
library(MCMCpack)
library(mvtnorm)

y.bar <- colMeans(y)
n <- nrow(y)
nsim <- 10000

theta.pos <- matrix(NA, nsim, P)
Sigma.pos <- array(NA, c(nsim, P, P))

# Valor inicial de theta
theta.pos[1,] <- c(0, 1)

# Parámetros posteriores de Sigma
v.pos <- v + n
matrix.theta <- kronecker(matrix(rep(1, n)),
                          t(theta.pos[1, ])) 
S.theta <- t(y - matrix.theta) %*% (y-matrix.theta)
Lambda.pos <- Lambda + S.theta
# Simulación de la distribución posterior condicional de Sigma
Sigma.pos[1, , ] <- riwish(v.pos, Lambda.pos)

#####################
# muestreo de Gibbs #
#####################

for(i in 2:nsim){
  # Parámetros posteriores de theta	
  Gamma.n <- solve(solve(Gamma) + 
                     n * solve(Sigma.pos[i - 1, , ]))
  mu.n <- Gamma.n %*%
    (solve(Gamma) %*% mu + 
       n * solve(Sigma.pos[i - 1, , ]) %*% y.bar)
  # Simulación de la distribucion posterior condicional de theta
  theta.pos[i, ] <- rmvnorm(1, mu.n, Gamma.n)
  # Parámetros posteriores de Sigma
  matrix.theta <- kronecker(matrix(rep(1, n)),
                          t(theta.pos[i, ]))  
  S.theta <- t(y - matrix.theta) %*% (y - matrix.theta)
  Lambda.pos <- Lambda + S.theta
  # Simulación de la distribución posterior condicional de Sigma
  Sigma.pos[i, , ] <- riwish(v.pos, Lambda.pos)
}
```

Una vez finalizada la ejecución del muestreo de Gibbs, debemos examinar la calidad de los valores muestreados para asegurar que las estimaciones bayesianas sean obtenidas de una muestra de valores que hayan convergido, y en segundo lugar que no estén correlacionados. Para eso a continuación observamos la gráfica de los valores generados para algunos parámetros (en particular, consideramos los parámetros $\theta_1$, $\theta_2$, $\sigma^2_1$ y $\sigma_{12}$), así como la gráfica de las autocorrelaciones muestrales.

\begin{figure}

{\centering \includegraphics{4Multiparametricos_files/figure-latex/unnamed-chunk-41-1} 

}

\caption{Convergencia de las distribuciones posteriores y diagramas de la función de autocorrelación en las cadenas.}(\#fig:unnamed-chunk-41)
\end{figure}

Con estas gráficas, observamos que los valores muestreados han alcanzado la convergencia; además estos tienen correlaciones cercanas a cero. De esta forma, podemos usar los valores muestreados para calcular las estimaciones y los intervalos de credibilidad.


```r
theta.Bayes <- colMeans(theta.pos)
Sigma.Bayes <- matrix(c(mean(Sigma.pos[,1,1]),
                        mean(Sigma.pos[,2,1]),
                        mean(Sigma.pos[,1,2]),
                        mean(Sigma.pos[,2,2])),
                      nrow = 2, ncol = 2)
theta.Bayes
```

```
## [1] 0.5638342 2.0993759
```

```r
Sigma.Bayes
```




\begin{tabular}{r|r}
\hline
3.042341 & 2.093069\\
\hline
2.093069 & 3.490823\\
\hline
\end{tabular}

El procedimiento inferencial sobre la comparación entre los efectos de los dos medicamentos se puede realizar de la misma manera como ilustró el ejemplo \@ref(exm:EjeStudent).


### Parámetros dependientes

Al igual que en el caso univariado, la inferencia posterior de los parámetros de interés debe ser llevada a cabo en dos etapas: En la primera, se debe establecer la distribución previa conjunta para ambos parámetros mediante
\begin{equation*}
p(\btheta,\bSigma)=p(\bSigma)p(\btheta \mid \bSigma)
\end{equation*}

Luego, en la segunda etapa es posible analizar posterior propiamente cada uno de los parámetros de interés puesto que
\begin{equation*}
p(\btheta,\bSigma \mid \mathbf{Y})\propto p(\mathbf{Y} \mid \btheta,\bSigma)p(\btheta,\bSigma)
\end{equation*}

Al igual que en el caso univariado, la anterior formulación conlleva a asignar una distribución previa para $\btheta$ dependiente de la matriz $\bSigma$. Esto quiere decir que en la distribución $p(\btheta \mid \bSigma)$ el valor de $\bSigma$ se considera una constante fija y conocida. Siguiendo los lineamientos del capítulo anterior, una distribución previa para $\btheta$  condicional a $\bSigma$ es
\begin{equation*}
p(\btheta \mid \bSigma)\sim Normal_p(\bmu,\bSigma/c_0)
\end{equation*}

Donde $c_0$ es una constante. Por otro lado, y siguiendo los argumentos de la sección anterior, una posible opción para la distribución previa de $\bSigma$, corresponde a
\begin{equation*}
p(\bSigma)\sim Inversa-Wishart_{v_0}(\bLambda)
\end{equation*}

\BeginKnitrBlock{proposition}
<span class="proposition" id="prp:unnamed-chunk-43"><strong>(\#prp:unnamed-chunk-43) </strong></span>La distribución previa conjunta de los parámetros $\btheta$ y $\bSigma$ está dada por
\begin{equation*}
p(\btheta,\bSigma) \propto \mid \bSigma \mid ^{-(v_0+p)/2-1}
\exp\left\{ -\frac{1}{2}\left[traza(\bLambda_0\bSigma^{-1})+
c_0(\btheta-\bmu)'\bSigma^{-1}(\btheta-\bmu)\right]\right\}
\end{equation*}
\EndKnitrBlock{proposition}
<br>
                              
\BeginKnitrBlock{proof}
\iffalse{} <span class="proof"><em>Prueba. </em></span>  \fi{}La prueba es inmediata al multiplicar las densidades y asignar los términos que no dependen de los parámetros de interés a la constante de proporcionalidad.
\EndKnitrBlock{proof}
<br>
                              
Para encontrar las distribuciones posteriores de cada uno de los parámetros de interés se utilizan argumentos similares a los del capítulo anterior.
                              
\BeginKnitrBlock{proposition}
<span class="proposition" id="prp:unnamed-chunk-45"><strong>(\#prp:unnamed-chunk-45) </strong></span>La distribución posterior de $\btheta$ condicional a $\bSigma,\mathbf{Y}$ está dada por
\begin{equation*}
\theta \mid \bSigma,\mathbf{Y} \sim Normal_p(\bmu_n,\bSigma/(n+c_0))
\end{equation*}

donde
\begin{equation*}
\mu_n=\frac{n\bar{\mathbf{Y}}+c_0\bmu}{n+c_0}
\end{equation*}
\EndKnitrBlock{proposition}
<br>

\BeginKnitrBlock{proof}
\iffalse{} <span class="proof"><em>Prueba. </em></span>  \fi{}Utilizando propiedades de la distribución condicional, tenemos que
\begin{align*}
p(\btheta|\bSigma,\mathbf{Y})&\propto p(\btheta, \bSigma|\mathbf{Y})\\
&\propto \mid \bSigma \mid ^{-(v_0+p)/2-1}
\exp\left\{ -\frac{1}{2}\left[traza(\bLambda_0\bSigma^{-1})+
c_0(\btheta-\bmu)'\bSigma^{-1}(\btheta-\bmu)\right]\right\}\\
&\ \ \ \ \ \ \ \ \mid \bSigma \mid ^{-n/2}\exp\left\{-\frac{1}{2}\sum_{i=1}^n(\mathbf{y}_i-\btheta)'\bSigma^{-1}(\mathbf{y}_i-\btheta)\right\}\\
&\propto \exp\left\{ -\frac{1}{2}
c_0(\btheta-\bmu)'\bSigma^{-1}(\btheta-\bmu)\right\}\exp\left\{-\frac{1}{2}\sum_{i=1}^n(\mathbf{y}_i-\btheta)'\bSigma^{-1}(\mathbf{y}_i-\btheta)\right\}
\end{align*}

La anterior expresión es la misma que apareción en el capítulo anterior para $p(\btheta\mid\bSigma,\mathbf{Y})$, en donde $\bSigma/c_0$ toma el valor de $\bGamma$. Así, teniendo en cuenta las ecuaciones \@ref(eq:Gamman) y \@ref(eq:mun), podemos afirmar que el vector de medias y la matriz de varianzas y covarianzas posterior están dadas por
\begin{align}
\bGamma_n&=((\bSigma/c_0)^{-1}+n\bSigma^{-1})^{-1}=\frac{\bSigma}{n+c_0}\\
\bmu_n&=\frac{\bSigma}{n+c_0}((\bSigma/c_0)^{-1}\bmu+n\bSigma^{-1}\bar{\mathbf{y}})=\frac{n\bar{\mathbf{Y}}+c_0\bmu}{n+c_0}
\end{align}
\EndKnitrBlock{proof}
<br>

En cuanto a la distribución de $\bSigma$, se tiene el siguiente resultado: 

\BeginKnitrBlock{proposition}
<span class="proposition" id="prp:PosSigma"><strong>(\#prp:PosSigma) </strong></span>La distribución marginal posterior de la matriz de parámetros $\bSigma$ es
\begin{equation*}
\bSigma \mid \mathbf{Y} \sim Inversa-Whishart_{n+v_0}(\bLambda_n)
\end{equation*}

Donde
\begin{equation}\label{bLambda_n}
\bLambda_n=\bLambda+(n-1)\mathbf{S}+\frac{c_0n}{c_0+n}(\bmu-\bar{\mathbf{y}})(\bmu-\bar{\mathbf{y}})'
\end{equation}

con $S$ la matriz de varianzas y covarianzas muestrales.
\EndKnitrBlock{proposition}
<br>

\BeginKnitrBlock{proof}
\iffalse{} <span class="proof"><em>Prueba. </em></span>  \fi{}\begin{align*}
&\ \ \ \ \ p(\bSigma\mid\mathbf{Y})\\
&=\int_{R^p} p(\btheta,\bSigma\mid\mathbf{Y})d\btheta\\
&\propto \int_{R^p}\mid\bSigma\mid^{-(v_0+p+n)/2-1}\exp\left\{-\frac{1}{2}\left[traza(\bLambda\bSigma^{-1})+c_0(\btheta-\bmu)'\bSigma^{-1}(\btheta-\bmu)+\sum_{i=1}^n(\mathbf{y}_i-\btheta)'\bSigma^{-1}(\mathbf{y}_i-\btheta)\right]\right\}d\btheta\\
&\propto \mid\bSigma\mid^{-(v_0+p+n)/2-1} \exp\left\{-\frac{1}{2}\left[traza(\bLambda\bSigma^{-1})\right]\right\}\\
&\ \ \ \ \ \ \ \ \ \ \int_{R^p}\exp\left\{-\frac{1}{2}\left[c_0(\btheta-\bmu)'\bSigma^{-1}(\btheta-\bmu)+\sum_{i=1}^n(\mathbf{y}_i-\btheta)'\bSigma^{-1}(\mathbf{y}_i-\btheta)\right]\right\}d\btheta\\
&\propto \mid\bSigma\mid^{-(v_0+p+n)/2-1} \exp\left\{-\frac{1}{2}\left[traza(\bLambda\bSigma^{-1})+c_0\bmu'\bSigma^{-1}\bmu+\sum_{i=1}^n(\mathbf{y}_i-\bar{\mathbf{y}})'\bSigma^{-1}(\mathbf{y}_i-\bar{\mathbf{y}})+n\bar{\mathbf{y}}'\bSigma^{-1}\bar{\mathbf{y}}\right]\right\}\\
&\ \ \ \ \ \ \ \ \ \ \int_{R^p}\exp\left\{-\frac{1}{2}\left[c_0\btheta'\bSigma^{-1}\btheta-2c_0\bmu'\bSigma^{-1}\btheta-2n\bar{\mathbf{y}}'\bSigma^{-1}\btheta+n\btheta'\bSigma^{-1}\btheta\right]\right\}d\btheta\\
&\propto \mid\bSigma\mid^{-(v_0+p+n)/2-1} \exp\left\{-\frac{1}{2}\left[traza\left((\bLambda+c_0\bmu\bmu'+\sum_{i=1}^n(\mathbf{y}_i-\bar{\mathbf{y}})(\mathbf{y}_i-\bar{\mathbf{y}})'+n\bar{\mathbf{y}}\bar{\mathbf{y}}')\bSigma^{-1}\right)\right]\right\}\\
&\ \ \ \mid\frac{\bSigma}{c_0+n}\mid^{1/2}\exp\left\{\frac{1}{2}\frac{c_0\bmu'+n\bar{\mathbf{y}}'}{c_0+n}\left(\frac{\bSigma}{c_0+n}\right)^{-1}\frac{c_0\bmu+n\bar{\mathbf{y}}}{c_0+n}\right\}\\
&\ \ \ \ \ \ \ \underbrace{\int_{R^p}\mid\frac{\bSigma}{c_0+n}\mid^{-1/2}\exp\left\{-\frac{1}{2}\left(\btheta-\frac{c_0\bmu+n\bar{\mathbf{y}}}{c_0+n}\right)'\left(\frac{\bSigma}{c_0+n}\right)^{-1}\left(\btheta-\frac{c_0\bmu+n\bar{\mathbf{y}}}{c_0+n}\right)\right\}d\btheta}_{\text{Igual a 1}}
\end{align*}

Por otro lado, 
\begin{align*}
&\ \ \ \ \frac{c_0\bmu'+n\bar{\mathbf{y}}'}{c_0+n}\left(\frac{\bSigma}{c_0+n}\right)^{-1}\frac{c_0\bmu+n\bar{\mathbf{y}}}{c_0+n}\\
&=\frac{1}{c_0+n}(c_0\bmu'+n\bar{\mathbf{y}}')\bSigma^{-1}(c_0\bmu+n\bar{\mathbf{y}})\\
&=traza\left(\frac{1}{c_0+n}(c_0\bmu+n\bar{\mathbf{y}})(c_0\bmu'+n\bar{\mathbf{y}}')\bSigma^{-1}\right)\\
&=traza\left(\left(\frac{c_0^2\bmu\bmu'}{c_0+n}+\frac{2c_0n\bar{\mathbf{y}}\bmu'}{c_0+n}+\frac{n^2\bar{\mathbf{y}}\bar{\mathbf{y}}'}{c_0+n}\right)\bSigma^{-1}\right)
\end{align*}

Reemplazando la anterior expresión en $p(\bSigma\mid\mathbf{Y})$, se tiene que
\begin{align*}
&\ \ \ \ \ p(\bSigma\mid\mathbf{Y})\\
&\propto \mid\bSigma\mid^{-(v_0+p+n+1)/2}\exp\left\{-\frac{1}{2}traza\left[\left(\bLambda+(n-1)\mathbf{S}+\frac{c_0n}{c_0+n}(\bmu-\bar{\mathbf{y}})(\bmu-\bar{\mathbf{y}})'\right)\bSigma^{-1}\right]\right\}
\end{align*}

la cual corresponde a la distribución deseada.
\EndKnitrBlock{proof}
<br>

En términos de simulación de densidades, para obtener las estimaciones bayesians de $\btheta$ y $\bSigma$ se debe primero simular valores de $\bSigma$ de la distribución $p(\bSigma \mid \mathbf{Y})$ y luego, se debe utilizar estos valores para simular valores de $\btheta$ de la distribución $p(\btheta \mid \bSigma,\mathbf{Y})$.

Una forma equivalente de obtener las estimaciones es calcular directamente la esperanza teórica de las distribuciones posteriores marginales de $\btheta$ y de $\bSigma$. 

Del resultado \@ref(prp:PosSigma), podemos concluir que la estimación bayesiana de la matriz de varianza y covarianzas $\bSigma$ está dada por
\begin{equation*}
\hat{\bSigma}=\dfrac{\bLambda+(n-1)\mathbf{S}+\frac{c_0n}{c_0+n}(\bmu-\bar{\mathbf{y}})(\bmu-\bar{\mathbf{y}})'}{n+v_0-p-1}
\end{equation*}

Teniendo en cuenta que la estimación previa de $\bSigma$ viene dada por $\hat{\bSigma}_{pre}=\frac{\bLambda}{v_0-p-1}$, podemos ver que la estimación bayesiana de $\bSigma$ está conformada por tres componentes: la estimación previa $\hat{\bSigma}_{pre}$, la estimación clásica $\mathbf{S}$ y una medida de discrepancia entre la estimación previa y la clásica de $\btheta$. Para encontrar correctas formas de escoger los parámetros previas de $\bSigma$, por ahora ignoramos el último componente, y vemos que la estimación previa $\hat{\bSigma}_{pre}$ y la estimación clásica $\mathbf{S}$ entran al cómputo de la estimación bayesiana con los pesos de $v_0-p-1$ y $n-1$, de esta forma, podemos escoger $v_0$ tal que $v_0-p$ represente el número de la información previa, y el valor de $\bLambda$ se puede calcular a partir de $v_0$ y $\hat{\bSigma}_{pre}$. 

El siguiente resultado muestra la distribución posterior marginal de $\btheta$. 
                              
\BeginKnitrBlock{proposition}
<span class="proposition" id="prp:Posbtheta"><strong>(\#prp:Posbtheta) </strong></span>La distribución marginal posterior del parámetro $\btheta$ es la distribución $t$ de Student multivariante tal que
\begin{equation*}
\btheta \mid \mathbf{Y} \sim t_{n+v_0-p+1}\left(\bmu_n, \frac{\bLambda_n}{(c_0+n)(n+v_0-p+1)}\right)
\end{equation*}

con $\bmu_n=\frac{c_0\bmu+n\bar{\mathbf{y}}}{c_0+n}$ y $\bLambda_n$ dado en la ecuación (\ref{bLambda_n}).
\EndKnitrBlock{proposition}
<br>

\BeginKnitrBlock{proof}
\iffalse{} <span class="proof"><em>Prueba. </em></span>  \fi{}\begin{align*}
&\ \ \ \ p(\btheta\mid\mathbf{Y})\\
&=\int_{R^p\times R^p}p(\btheta,\bSigma\mid\mathbf{Y})d\bSigma\\
&=\int_{R^p\times R^p} \mid\bSigma\mid^{-(v_0+p+n)/2-1}\exp\left\{-\frac{1}{2}\left[traza(\bLambda\bSigma^{-1})+c_0(\btheta-\bmu)'\bSigma^{-1}(\btheta-\bmu)+\sum_{i=1}^n(\mathbf{y}_i-\btheta)'\bSigma^{-1}(\mathbf{y}_i-\btheta)\right]\right\}d\bSigma\\
&=\int_{R^p\times R^p} \mid\bSigma\mid^{-(v_0+p+n+2)/2}\exp\left\{-\frac{1}{2}traza\left[\bLambda+c_0(\btheta-\bmu)(\btheta-\bmu)'+\sum_{i=1}^n(\mathbf{y}_i-\btheta)(\mathbf{y}_i-\btheta)'\right]\bSigma^{-1}\right\}d\bSigma\\
&\propto \big|\bLambda+c_0(\btheta-\bmu)(\btheta-\bmu)'+\sum_{i=1}^n(\mathbf{y}_i-\btheta)(\mathbf{y}_i-\btheta)'\big|^{-\frac{v_0+n+1}{2}}\\
&=\big|\bLambda+c_0(\btheta-\bmu)(\btheta-\bmu)'+\sum_{i=1}^n(\mathbf{y}_i-\bar{\mathbf{y}})(\mathbf{y}_i-\bar{\mathbf{y}})'+n(\bar{\mathbf{y}}-\btheta)(\bar{\mathbf{y}}-\btheta)'\big|^{-\frac{v_0+n+1}{2}}\\
&=\big|\bLambda+(n-1)\mathbf{S}+\frac{c_0n}{c_0+n}(\bmu-\bar{\mathbf{y}})(\bmu-\bar{\mathbf{y}})'+(c_0+n)(\btheta-\bmu_n)(\btheta-\bmu_n)'\big|^{-\frac{v_0+n+1}{2}}\\
&=\big|\bLambda_n+(c_0+n)(\btheta-\bmu_n)(\btheta-\bmu_n)'\big|^{-\frac{v_0+n+1}{2}}\\
&\propto\big|\mathbf{I}_p+(c_0+n)\bLambda_n^{-1}(\btheta-\bmu_n)(\btheta-\bmu_n)'\big|^{-\frac{v_0+n+1}{2}}\\
&=\big|1+(c_0+n)(\btheta-\bmu_n)'\bLambda_n^{-1}(\btheta-\bmu_n)\big|^{-\frac{v_0+n+1}{2}}\\
&=\big|1+\frac{1}{n+v_0-p+1}(\btheta-\bmu_n)'\left(\frac{\bLambda_n}{(c_0+n)(n+v_0-p+1)}\right)^{-1}(\btheta-\bmu_n)\big|^{-\frac{v_0+n+1}{2}}
\end{align*}

Esta expresión obtenida corresponde a la forma de la distribución $t$ de Student multivariado. En el desarrollo se utilizó la propiedad $|\mathbf{I}+\mathbf{A}\mathbf{B}|=|\mathbf{I}+\mathbf{B}\mathbf{A}|$ para matrices $\mathbf{A}$ y $\mathbf{B}$ de tamaños compatibles para las multiplicaciones. 
\EndKnitrBlock{proof}
<br>

El anterior resultado indica que la estimación bayesiana del parámetro $\btheta$ está dada por 
\begin{equation*}
\hat{\btheta}=\mathbf{\mu}_n=\dfrac{n\hat{\mathbf{Y}}+c_0\mathbf{\mu}}{n+c_0}=\dfrac{n}{n+c_0}\hat{\mathbf{Y}}+\dfrac{c_0}{n+c_0}\mathbf{\mu}
\end{equation*}

donde se puede observar que $\hat{\btheta}$ se acercará a la estimación clásica $\hat{\mathbf{y}}$ cuando $n$ es grande comparado a $c_0$, de lo contrario se acercará a la estimación previa $\mathbf{\mu}$. La varianza posterior para el $i$-ésimo componente de $\btheta$ está dada por 
\begin{equation*}
var(\theta_i|\mathbf{Y})=\dfrac{\lambda_{ii}}{(c_0+n)(n+v_0-p+1)}\dfrac{n+c_0-p+1}{n+c_0-p-1}\approx\dfrac{\lambda_{ii}}{(c_0+n)(n+v_0-p+1)}
\end{equation*}

donde $\lambda_{ii}$ denota el $i$-ésimo elemento en la diagonal de la matriz $\bLambda_n$.

\BeginKnitrBlock{example}
<span class="example" id="exm:unnamed-chunk-49"><strong>(\#exm:unnamed-chunk-49) </strong></span>@Pena2002 reporta las mediciones de 6 variables indicadoras de desarrollo en 91 paises en los años noventa. Para este ejemplo, utilizamos tres variables: la tasa de natalidad, la tasa de mortalidad y la mortalidad infantil en algunos paises de Suramérica y Asia mostrados en la tabla \@ref(tab:Natalidad). Específicamente, usaremos los datos de los paises de Suramérica como datos muestrales y los de Asia para extraer la información previa.
\EndKnitrBlock{example}

\begin{table}

\caption{(\#tab:Natalidad)Tasa de natalidad, tasa de mortalidad, mortalidad infantil en algunos países.}
\centering
\begin{tabular}[t]{r|r|r}
\hline
Natalidad & Mortalidad & MortalidadInfantil\\
\hline
20.7 & 8.4 & 25.7\\
\hline
46.6 & 18.0 & 111.0\\
\hline
28.6 & 7.9 & 63.0\\
\hline
23.4 & 5.8 & 17.1\\
\hline
27.4 & 6.1 & 40.0\\
\hline
32.9 & 7.4 & 63.0\\
\hline
29.0 & 23.2 & 43.0\\
\hline
34.8 & 6.6 & 42.0\\
\hline
32.9 & 8.3 & 109.9\\
\hline
18.0 & 9.6 & 21.9\\
\hline
27.5 & 4.4 & 23.3\\
\hline
\end{tabular}
\end{table}

A continuación se muestra el proceso necesario para obtener la estimación bayesiana del vector de medias y de la matriz de covarianzas con los datos del ejemplo. 


```r
# Datos muestrales
y.sam <- data.frame(
  Nata = c(20.7, 46.6, 28.6, 23.4, 27.4, 
           32.9, 29, 34.8, 32.9,18,27.5), 
  Mort = c(8.4, 18, 7.9, 5.8, 6.1, 7.4,
           23.2, 6.6, 8.3, 9.6, 4.4),
  Infa = c(25.7, 111, 63, 17.1, 40, 63,
           43, 42, 109.9, 21.9, 23.3))
# Datos de la información previa
y.pre <- data.frame(
  Nata = c(21.2, 30.5, 28.6, 31.6, 
           36.1, 39.6, 17.8),
  Mort = c(6.7, 10.2, 9.4, 5.6,
           8.8, 14.8, 5.2),
  Infa=c(32, 91, 75, 24, 68, 128, 7.5))

n <- nrow(y.sam)
P <- ncol(y.sam)
# Estimación clásica de los parámetros
y.bar <- colMeans(y.sam) 
S <- var(y.sam)

# Estimación previa de los parámetros
mu <- colMeans(y.pre)
c0 <- nrow(y.pre)
v0 <- P + nrow(y.pre)
Lambda <- var(y.pre) * (v0 - P - 1)

# Parámetros de las distribuciones posteriores marginales
mu.n <- (n * y.bar + c0 * mu)/(n + c0)
Lambda.n <- Lambda + (n - 1) * S + 
  matrix(mu - y.bar) %*% 
  t(matrix(mu - y.bar)) * c0 * n
var.theta <- Lambda.n/((c0 + n) * (n + v0 - P + 1))
mu.n
```

```
##      Nata      Mort      Infa 
## 29.288889  9.244444 54.744444
```

```r
var.theta
```




\begin{tabular}{l|r|r|r}
\hline
  & Nata & Mort & Infa\\
\hline
Nata & 2.7991957 & 0.7890556 & 10.675752\\
\hline
Mort & 0.7890556 & 1.3526977 & 2.195668\\
\hline
Infa & 10.6757519 & 2.1956683 & 85.548053\\
\hline
\end{tabular}

```r
Lambda.n
```




\begin{tabular}{l|r|r|r}
\hline
  & Nata & Mort & Infa\\
\hline
Nata & 957.3249 & 269.8570 & 3651.1071\\
\hline
Mort & 269.8570 & 462.6226 & 750.9186\\
\hline
Infa & 3651.1071 & 750.9186 & 29257.4343\\
\hline
\end{tabular}


De los anteriores cálculos, se puede ver que la distribución posterior de $\btheta$ está dada por
\begin{equation*}
\btheta\mid\mathbf{Y}\sim t_{19}\left(\begin{pmatrix}29.3\\9.2\\54.7\end{pmatrix}, \begin{pmatrix}2.80&0.79&10.68\\0.79&1.35&2.20\\10.68&2.20&85.55\end{pmatrix}\right)
\end{equation*}

Usando propiedades de la distribución multivariante $t$ de Student, tenemos que $\theta_1\sim t_{19}(29.29, 2.80)$, $\theta_2\sim t_{19}(9.24, 1.35)$ y $\theta_3\sim t_{19}(0.62, 85.55)$, de allí se puede encontrar fácilmente los intervalos de credibilidad para cada uno de estos tres parámetros.

En cuanto a la distribución posterior de $\bSigma$, ésta está dada por
\begin{equation*}
\bSigma\mid\mathbf{Y}\sim Inversa-Wishart_{21}\left(\begin{pmatrix}957&270&3651\\ 270&463&751\\ 3651&751&29257\end{pmatrix}\right)
\end{equation*}

La estimación bayesiana de $\bSigma$ viene dada por $\hat{\bSigma}=\begin{pmatrix}56.3&15.9&214.8\\15.9&27.2&44.2\\214.8&44.2&1721.0\end{pmatrix}$. Por propiedades de la distribución inversa-Wishart, se puede concluir que los elementos diagonales de $\bSigma$ tienen distribución inversa-Gamma. Por ejemplo, se tiene que $\sigma^2_{1}\sim Inversa-Gamma(21/2, 56.3/2)$ y cualquier inferencia que se desear realizar sobre $\sigma^2_{1}$ es posible a partir de esta distribución. 

Aparte de los análisis anteriores, también podemos realizar ejercicios de comparación y verificar la posible independencia entre parejas de variables. Por ejemplo, si queremos verificar la hipótesis de que la tasa de natalidad es dos veces la tasa de mortalidad, esto es $\theta_1=2 \times \theta_2$. Una forma de confirmar o refutar esta hipótesis es hallar el intervalo de credibilidad de $\theta_1 - 2\theta_2$ que se puede expresar como $(1,-2,0)\begin{pmatrix}\theta_1\\ \theta_2\\ \theta_3\end{pmatrix}$. Por propiedades de la distribución $t$ de Student multivariante, tenemos que $\theta_1 - 2 \theta_2$ tiene distribución $t$ de Student univariada con los mismos grados de libertad que $\btheta$. La esperanza de esta distribución está dada por $(1,-2,0)\begin{pmatrix}29.3\\9.2\\54.7\end{pmatrix}=10.8$ y la escala está dada por $(1,-2,0)\begin{pmatrix}2.80&0.79&10.68\\0.79&1.35&2.20\\10.68&2.20&85.55\end{pmatrix}\begin{pmatrix}1\\-2\\0\end{pmatrix}=5.054$, esto es, 
\begin{equation*}
\theta_1-2\theta_2 \mid \mathbf{Y}\sim t_{19}(10.8, 5.054)
\end{equation*}

De esta forma, un intervalo de credibilidad para $\theta_1-2\theta_2$ viene dado por los percentiles 2.5\% y 97.5\% de la anterior distribución, que a la vez son iguales a los percentiles 2.5\% y 97.5\% de la distribución $t$ de Student estandarizada multiplicado por $\sqrt{5.054}$ y sumando 10.8. Este intervalo es igual a $(6.095, 15.505)$. Al observar que este intervalo no contiene el valor 0, podemos concluir que no es válido afirmar que la tasa de natalidad sea dos veces la tasa de mortalidad.

En el anterior análisis, vemos que el intervalo de credibilidad para $\theta_1-2\theta_2$ contiene solo valores positivos, lo cual es un indicio de que la variable $\theta_1-2\theta_2$ tenga la mayor parte de la función de densidad ubicada en el eje positivo. De hecho podemos indagar por $Pr(\theta_1-2\theta_2>0)$, la cual se puede calcular de la distribución $t_{19}(10.8, 5.054)$ encontrada anteriormente. Esta probabilidad es `1 - pt((0 - 10.8)/sqrt(5.054), 19)` dando como resultado 0.9999383; de donde se muestra una fuerte evidencia de que la tasa de natalidad es superior a dos veces la tasa de mortalidad.

Los anteriores resultados fueron obtenidos directamente de las distribuciones posteriores marginales $p(\btheta|\mathbf{Y})$ y $p(\bSigma|\mathbf{Y})$. De forma equivalente también se puede usar las técnicas de simulación con base en las distribuciones $p(\btheta,\bSigma|\mathbf{Y})$ y $p(\bSigma|\mathbf{Y})$. A continuación se muestran los códigos necesario en `R`:


```r
nsim <- 1000
theta.pos <- matrix(NA, nsim, P)
Sigma.pos <- array(NA, c(nsim, P, P))

for(i in 1:nsim){
  Sigma.pos[i,,] <- riwish(n + v0, Lambda.n)
  theta.pos[i,] <- rmvnorm(1, mu.n, 
                           Sigma.pos[i, , ]/(n + c0))
}

# Estimaciones finales
theta.final <- colMeans(theta.pos)
Sigma.final <- matrix(c(mean(Sigma.pos[, 1, 1]),
                        mean(Sigma.pos[, 1, 2]),
                        mean(Sigma.pos[, 1, 3]),
                        mean(Sigma.pos[, 2, 1]),
                        mean(Sigma.pos[, 2, 2]),
                        mean(Sigma.pos[, 2, 3]),
                        mean(Sigma.pos[, 3, 1]),
                        mean(Sigma.pos[, 3, 2]), 
                        mean(Sigma.pos[, 3, 3])), 
                      nrow = P, ncol = P)
theta.final
```

```
## [1] 29.245516  9.288723 54.650509
```

```r
Sigma.final
```




\begin{tabular}{r|r|r}
\hline
55.88544 & 15.59074 & 211.03550\\
\hline
15.59074 & 26.47268 & 42.52096\\
\hline
211.03550 & 42.52096 & 1687.49974\\
\hline
\end{tabular}

Podemos ver que los resultados obtenidos con los dos métodos son totalmente coincidentes. En cuanto al intervalo de credibilidad para $\theta_1-2\theta_2$, se puede calcular con


```r
quantile(
  theta.pos[,1] - 2 * theta.pos[, 2], 
  c(0.025, 0.975))
```

```
##      2.5%     97.5% 
##  6.120592 15.045185
```

También podemos calcular $Pr(\theta_1-2\theta_2>0)$ como 


```r
sum(theta.pos[, 1] > 2 * theta.pos[, 2])/nsim
```

```
## [1] 1
```
Podemos ver que estos resultados son muy similares a los obtenidos usando $p(\btheta|\mathbf{Y})$.

### Parámetros no informativos

@Gelman03 afirma que la distribución previa no informativa de Jeffreys conjunta para $\btheta,\bSigma$, en este caso está dada por la siguiente expresión
\begin{equation*}
p(\btheta,\bSigma)\propto \mid \bSigma \mid ^{-(p+1)/2}
\end{equation*}

La distribución posterior conjunta para $\btheta,\bSigma$ está dada por

\begin{equation*}
p(\btheta,\bSigma \mid \mathbf{Y})\propto
\mid \bSigma \mid ^{-(p+n+1)/2}
\exp\left\{ -\frac{1}{2}\sum_{i=1}^n
  (\mathbf{Y}_i-\btheta)'\bSigma^{-1}(\mathbf{Y}_i-\btheta)\right\}
  \end{equation*}

De la anterior distribución, podemos encontrar la distribución condicional posterior de $\btheta$ dada en el siguiente resultado.

\BeginKnitrBlock{proposition}
<span class="proposition" id="prp:unnamed-chunk-54"><strong>(\#prp:unnamed-chunk-54) </strong></span>La distribución posterior del vector de parámetros $\btheta$ condicional a $\bSigma,\mathbf{Y}$ es
\begin{equation*}
\btheta \mid \bSigma,\mathbf{Y}\sim Normal_p(\bar{\mathbf{y}},\bSigma/n)
\end{equation*}
\EndKnitrBlock{proposition}
<br>

\BeginKnitrBlock{proof}
\iffalse{} <span class="proof"><em>Prueba. </em></span>  \fi{}Algunas simples operaciones algebráicas muestran que:
\begin{align*}
p(\btheta \mid \bSigma,\mathbf{Y}) &\propto \exp\left\{ -\frac{1}{2}\sum_{i=1}^n (\mathbf{Y}_i-\btheta)'\bSigma^{-1}(\mathbf{Y}_i-\btheta)\right\}\\
&\propto \exp\left\{ -\frac{n}{2}(\btheta-\bar{\mathbf{Y}})'\bSigma^{-1}(\btheta-\bar{\mathbf{Y}})\right\}
\end{align*}
Por lo tanto, factorizando convenientemente, se encuentra una expresión idéntica a la función de distribución de una variable aleatoria con distribución $Normal_p(\bar{y},\bSigma/n)$.
\EndKnitrBlock{proof}
<br>

En cuanto a la estimación de $\bSigma$, en el siguiente resultado encontramos su distribución posterior.

\BeginKnitrBlock{proposition}
<span class="proposition" id="prp:unnamed-chunk-56"><strong>(\#prp:unnamed-chunk-56) </strong></span>La distribución marginal posterior de la matriz de parámetros $\bSigma$ es
\begin{equation*}
\bSigma \mid \mathbf{Y}\sim Inversa-Whishart_{n-1}(\mathbf{S})
\end{equation*}
donde $\mathbf{S}=\sum_{i=1}^n(\mathbf{y}_i-\bar{\mathbf{y}})(\mathbf{y}_i-\bar{\mathbf{y}})'$
\EndKnitrBlock{proposition}
<br>

\BeginKnitrBlock{proof}
\iffalse{} <span class="proof"><em>Prueba. </em></span>  \fi{}En primer lugar recordamos la expresión 
\begin{equation*}
\textbf{S}_{\btheta}=\sum_{i=1}^n(\mathbf{y}_i-\btheta)(\mathbf{y}_i-\btheta)'=\mathbf{S}+n(\btheta-\bar{\mathbf{y}})(\btheta-\bar{\mathbf{y}})'
\end{equation*}

Por otro lado, recurriendo a las propiedades del operador $traza$, e integrando la distribución posterior conjunta con respecto a $\btheta$, se tiene que
\begin{align*}
p(\bSigma \mid \mathbf{Y})&=\int p(\btheta,\bSigma \mid \mathbf{Y}) \ d\btheta\\
&= \mid \bSigma \mid ^{-(p+n+1)/2}\int\exp\left\{ -\frac{1}{2}\sum_{i=1}^n
  (\mathbf{Y}_i-\btheta)'\bSigma^{-1}(\mathbf{Y}_i-\btheta)\right\} \ d\btheta\\
  &= \mid \bSigma \mid ^{-(p+n+1)/2}\int
  \exp\left\{ -\frac{1}{2}traza(\bSigma^{-1}\mathbf{S}_{\btheta})\right\} \ d\btheta\\
  &= \mid \bSigma \mid ^{-(p+n+1)/2}\int
  \exp\left\{ -\frac{1}{2}traza(\bSigma^{-1}
  (\mathbf{S}+n(\btheta-\bar{\mathbf{y}})(\btheta-\bar{\mathbf{y}})'))\right\} \ d\btheta\\
&= \mid \bSigma \mid ^{-(p+n)/2}\exp\left\{ -\frac{1}{2}traza(\bSigma^{-1}\mathbf{S})\right\}\\
&\hspace{2cm}\times
\int \mid \bSigma \mid ^{-1/2}\exp\left\{ -\frac{n}{2}traza(\bSigma^{-1}(\btheta-\bar{\mathbf{y}})(\btheta-\bar{\mathbf{y}})')\right\} \ d\btheta\\
&= \mid \bSigma \mid ^{-(p+n)/2}\exp\left\{ -\frac{1}{2}traza(\bSigma^{-1}\mathbf{S})\right\}\\
&\hspace{2cm}\times\int \mid \bSigma \mid ^{-1/2}\exp\left\{ -\frac{n}{2}traza((\btheta-\bar{\mathbf{y}})'\bSigma^{-1}(\btheta-\bar{\mathbf{y}}))\right\} \ d\btheta\\
&= \mid \bSigma \mid ^{-(p+n)/2}\exp\left\{ -\frac{1}{2}traza(\bSigma^{-1}\mathbf{S})\right\}\\
&\hspace{2cm}\times
\int\underbrace{ \mid \bSigma \mid ^{-1/2}\exp\left\{ -\frac{n}{2}(\btheta-\bar{\mathbf{y}})'\bSigma^{-1}(\btheta-\bar{\mathbf{y}})\right\}}_{Normal_p(\bar{\mathbf{y}},\bSigma/n)} \ d\btheta\\
&= \mid \bSigma \mid ^{-(p+n)/2}\exp\left\{ -\frac{1}{2}traza(\bSigma^{-1}\mathbf{S})\right\}
\end{align*}

Por lo tanto, factorizando convenientemente, se encuentra una expresión idéntica a la función de distribución de una variable aleatoria con distribución $Inversa-Whishart_{n-1}(\mathbf{S})$.
\EndKnitrBlock{proof}
<br>

El anterior resultado indica que la estimación bayesiana de $\bSigma$ cuando se utiliza una previa no informativa está dada por 
\begin{equation*}
\hat{\bSigma}=\frac{\sum_{i=1}^n(\mathbf{y}_i-\bar{\mathbf{y}})(\mathbf{y}_i-\bar{\mathbf{y}})'}{n-p-2}
\end{equation*}

Esta expresión es muy similar a la estimación clásica de la matriz de varianzas y covarianzas dada por $\frac{\sum_{i=1}^n(\mathbf{y}_i-\bar{\mathbf{y}})(\mathbf{y}_i-\bar{\mathbf{y}})'}{n-1}$. Se puede observar que a medida que $n$ se aumente, las dos expresiones darán resultados muy similares, pero siempre la estimación bayesiana será mayor a la estimación clásica, especialmente en situaciones donde el tamaño muestral es pequeño.

Para obtener la estimación de $\bSigma$ junto con la estimación de $\btheta$, podemos proceder de la siguiente forma para obetener valores simulados de $\btheta$ y $\bSigma$ y así, obtener las estimaciones respectivas. Si el número de iteraciones se fija como $G$, entonces se procede a:

1. Simular $G$ valores de la distribución de $\bSigma|\mathbf{Y}$; estos valores se denotan por $\bSigma_{(1)},\bSigma_{(2)},\cdots,\bSigma_{(G)}$.
2. Para cada valor de $\bSigma_{(g)}$, con $g=1,\cdots,G$, simular un valor de la distribución de $\btheta|\bSigma,\mathbf{Y}$; es decir, de la distribución $N_p(\bar{\mathbf{y}}, \bSigma/n)$, donde $\bSigma$ se reemplaza por $\bSigma_{(g)}$. De esta forma, se obtienen los valores $\btheta_{(1)},\btheta_{(2)},\cdots,\btheta_{(G)}$.

El siguiente ejemplo ilustra la forma de obtener las estimaciones siguiendo el anterior procedimiento. 

\BeginKnitrBlock{example}
<span class="example" id="exm:EjeStudent3"><strong>(\#exm:EjeStudent3) </strong></span>Retomamos los datos del efecto de aumento en horas de sueño de dos medicamentos soporíferos utilizados en los ejemplos \@ref(exm:EjeStudent) y \@ref(exm:EjeStudent2). Los siguientes códigos en `R` ilustran el procedimiento computacional para obtener valores de la distribución posterior conjunta de $\btheta$ y $\bSigma$.
\EndKnitrBlock{example}


```r
library(MCMCpack)
library(mvtnorm)

y <- as.matrix(
  data.frame(M1 = sleep[1:10, 1], 
             M2 = sleep[-(1:10), 1]))
n <- nrow(y)

y.bar <- colMeans(y)
S <- var(y) * (n - 1)

nsim <- 1000
theta.pos <- matrix(NA, nsim, 2)
Sigma.pos <- array(NA, c(nsim, 2, 2))

for(i in 1:nsim){
  #simulacion de la distribucion posterior condicional de Sigma
  Sigma.pos[i, , ] <- riwish(n - 1, S)
  #simulacion de la distribucion posterior condicional de theta
  theta.pos[i, ] <- rmvnorm(1, y.bar, Sigma.pos[i, , ]/n)
}
```


Dado que en el cálculo no se hizo uso de valores iniciales y por la forma de las distribuciones posteriores de $p(\btheta|\bSigma,\mathbf{Y})$ y $\bSigma|\mathbf{Y}$, los valores muestrados en la diferentes iteraciones no guardan relación entre sí; por ende, podemos usar directamente todos los valores simulados para realizar el cálculo de las estimaciones bayesianas.


```r
theta.Bayes <- colMeans(theta.pos)
Sigma.Bayes <- matrix(c(mean(Sigma.pos[, 1, 1]),
                        mean(Sigma.pos[, 2, 1]),
                        mean(Sigma.pos[, 1, 2]),
                        mean(Sigma.pos[, 2, 2])), 
                      ncol = 2, nrow = 2)
theta.Bayes
```

```
## [1] 0.7890492 2.3443686
```

```r
Sigma.Bayes
```




\begin{tabular}{r|r}
\hline
4.985367 & 4.337360\\
\hline
4.337360 & 6.002711\\
\hline
\end{tabular}

Por otro lado, la estimación clásica de los parámetros está dada por


```r
y.bar
```

```
##   M1   M2 
## 0.75 2.33
```

```r
var(y)
```




\begin{tabular}{l|r|r}
\hline
  & M1 & M2\\
\hline
M1 & 3.200556 & 2.848333\\
\hline
M2 & 2.848333 & 4.009000\\
\hline
\end{tabular}


Por consiguiente, podemos observar que, en cuanto al parámetro $\btheta$, la estimación bayesiana es igual a la estimación clásica; mientras que el determinante de la estimación bayesiana de $\bSigma$ es mucho mayor que el de la estimación clásica, esto ocurre en situaciones cuando el tamaño muestral es pequeño. En cuanto a la estimación por intervalo de los efectos promedios de los dos medicamentos, tenemos que:


```r
quantile(theta.pos[, 1], c(0.025, 0.975))
```

```
##       2.5%      97.5% 
## -0.6388661  2.1506794
```

```r
t.test(y[, 1])$conf.int
```

```
## [1] -0.5297804  2.0297804
## attr(,"conf.level")
## [1] 0.95
```

```r
quantile(theta.pos[, 2], c(0.025, 0.975))
```

```
##      2.5%     97.5% 
## 0.7742381 3.9578339
```

```r
t.test(y[,2])$conf.int
```

```
## [1] 0.8976775 3.7623225
## attr(,"conf.level")
## [1] 0.95
```

De las anteriores salidas observamos que los resultados obtenidos con el enfoque bayesiano, aunque no son exactamente iguales a los obtenidos con el enfoque clásico, sí son muy similares. En cuanto a la estimación por intervalo de las varianzas y covarianzas. Primero consideramos la varianza del primer medicamento denotada por $\sigma^2_1$. La distribución posterior de la matriz de varianzas y covarianzas está dada por
\begin{equation*}
\bSigma=\begin{pmatrix}\sigma^2_1&\sigma_{12}\\\sigma_{21}&\sigma^2_{2} \end{pmatrix}\sim Inversa-Wishart_9(\mathbf{S})
\end{equation*}

con $\mathbf{S}=\sum_{i=1}^{10}(\mathbf{y}_i-\bar{\mathbf{y}})(\mathbf{y}_i-\bar{\mathbf{y}})'=\begin{pmatrix}28.81&25.64\\25.64&36.08\end{pmatrix}$. Usando las propiedades de la distribución Inversa-Wishart, se puede concluir que la distribución marginal posterior de $\sigma^2_1$ está dada por $Inversa-Gamma(\alpha=\frac{9-1}{2}, \beta=\frac{28.81}{2})$, y su intervalo de credibilidad se puede calcular directamente de dicha distribución, o equivalentemente usando los percentiles muestrales de los valores de $\sigma^2_1$ simulados El intervalo obtenido por estos dos medios son muy similares como se puede ver a continuación.


```r
library(pscl)
qigamma(0.025, alpha = 8/2, beta = 28.81/2)
```

```
## [1] 1.643042
```

```r
qigamma(0.975, alpha = 8/2, beta = 28.81/2)
```

```
## [1] 13.21723
```

```r
quantile(Sigma.pos[, 1, 1], c(0.025, 0.975))
```

```
##     2.5%    97.5% 
##  1.60053 14.00017
```

El intervalo de confianza del 95% se puede obtener con el siguiente código^[Consultar @Zhang[, sección.3.2.1] para mayor información.].


```r
c(9 * var(y[, 1]) / qchisq(0.975, 9), 
  9 * var(y[, 2]) / qchisq(0.025, 9))
```

```
## [1]  1.514238 13.361406
```

En comparación con el intervalo de credibilidad, el intervalo de confianza está ubicado levemente hacia la izquierda del eje real, esto se debe a que la estimación clásica de la varianza siempre será menor a la estimación bayesiana con una previa no informativa.




































