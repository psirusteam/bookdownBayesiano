


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

Las figuras \@ref(fig:posNormalMVStan2) muestra la distribución posterior para este ejemplo, junto con la estimación puntual, correspondiente a la desviación estándar. 

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
\bGamma_n &= \left(\bGamma^{-1}+n\bSigma^{-1}\right)^{-1}\label{Gamma_n}\\
\bmu_n &= \bGamma_n(\bGamma^{-1}\bmu+n \bSigma^{-1}\bar{\mathbf{y}})\label{mu_n}
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
<span class="example" id="exm:unnamed-chunk-27"><strong>(\#exm:unnamed-chunk-27) </strong></span>@Student introdujo un conjunto de datos clásicos sobre el incremento en horas de sueño producido con 2 medicamentos soporíferos diferentes comparados con grupo control en 10 pacientes. Estos datos se pueden encontrar en `R` con el nombre `sleep` y se pueden definir como realizaciones de vectores aleatorios con distribución normal bivariada. Supongamos que la matriz de varianzas y covarianzas de la distribución es conocida e igual a $\Sigma=\begin{pmatrix}1&0.6\\ 0.6&2\end{pmatrix}$.
  
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
## Samples were drawn using NUTS(diag_e) at Tue Jun 22 23:42:17 2021.
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
## Samples were drawn using NUTS(diag_e) at Tue Jun 22 23:42:17 2021.
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
<span class="proposition" id="prp:unnamed-chunk-34"><strong>(\#prp:unnamed-chunk-34) </strong></span>La distribución posterior marginal de un subconjunto de parámetros, digamos $\btheta^{(1)}$ es también normal multivariante con media igual a la del subvector de medias apropiado, $\bmu_n^{(1)}$ y similar matriz de varianzas $\bGamma_n^{(11)}$.
\EndKnitrBlock{proposition}
  
\BeginKnitrBlock{proposition}
<span class="proposition" id="prp:unnamed-chunk-35"><strong>(\#prp:unnamed-chunk-35) </strong></span>La distribución posterior condicional de un subconjunto de parámetros, digamos $\btheta^{(1)}$, dado $\btheta^{(2)}$ es también normal multivariante dada por
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


































