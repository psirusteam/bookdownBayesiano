


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

El primer enfoque que consideraremos para el análisis de los parámetros de interés $\theta$ y $\sigma^2$ en una distribución normal univariada es aquel que supone que las distribuciones previas de cada uno de ellos son independientes, pero al mismo tiempo informativas. \citeasnoun{Gelman03} afirma que este supuesto de independencia es atractivo en problemas para los cuales la información previa para $\theta$ no toma la forma de un número fijo de observaciones con varianza $\sigma^2$. Adicionalmente, este supuesto de independencia es coherente con el hecho de que en la teoría clásica de estimación los estimadores insesgados de varianza mínima de $\theta$ y $\sigma^2$ son independientes [@Zhang, Sección 2.4]

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

\BeginKnitrBlock{proposition}<div class="proposition"><span class="proposition" id="prp:unnamed-chunk-1"><strong>(\#prp:unnamed-chunk-1) </strong></span>La distribución posterior conjunta de los parámetros de interés está dada por

\begin{align}
p(\theta,\sigma^2 \mid \mathbf{Y})&\propto (\sigma^2)^{-(n+n_0)/2-1} \notag \\
&\times
\exp\left\{-\frac{1}{2\sigma^2}\left[n_0\sigma^2_0+(n-1)S^2+n(\bar{y}-\theta)^2\right]-\frac{1}{2\tau^2}(\theta-\mu)^2\right\}
\end{align}</div>\EndKnitrBlock{proposition}
<br>
 
\BeginKnitrBlock{proof}<div class="proof">\iffalse{} <span class="proof"><em>Prueba. </em></span>  \fi{}Tenemos que
\begin{align*}
p(\theta,\sigma^2 \mid \mathbf{Y})&\propto p(\mathbf{Y} \mid \theta,\sigma^2)p(\theta,\sigma^2)\\
&\propto(\sigma^2)^{-n/2}\exp\left\{-\frac{1}{2\sigma^2}\sum_{i=1}^n(y_i-\theta)^2\right\}(\sigma^2)^{-n_0/2-1}\exp\left\{-\dfrac{n_0\sigma^2_0}{2\sigma^2}\right\}
\exp\left\{-\frac{1}{2\tau^2}(\theta-\mu)^2\right\} \\
&=(\sigma^2)^{-(n+n_0)/2-1}\exp\left\{-\frac{1}{2\sigma^2}\left[n_0\sigma^2_0+\sum_{i=1}^n(y_i-\theta)^2\right]-\frac{1}{2\tau^2}(\theta-\mu)^2\right\}\\
&\propto (\sigma^2)^{-(n+n_0)/2-1}\exp\left\{-\frac{1}{2\sigma^2}\left[n_0\sigma^2_0+(n-1)S^2+n(\bar{y}-\theta)^2\right]-\frac{1}{2\tau^2}(\theta-\mu)^2\right\}
\end{align*}
donde la última expresión se obtiene al sumar y restar $\bar{y}$ dentro de $(y_i-\theta)^2$.</div>\EndKnitrBlock{proof}
<br>

Nótese que la distribución posterior conjunta no tiene una forma estructural conocida y por lo tanto no es posible realizar el método de integración analítica para obtener una constante de integración [@Migon]. Sin embargo, sí es posible obtener las distribuciones condicionales posteriores de $\theta$ y de $\sigma^2$, notando que

\begin{align*}
p(\theta \mid \sigma^2,\mathbf{Y})\propto p(\theta,\underbrace{\sigma^2}_{fijo} \mid \mathbf{Y})
\ \ \ \ \ \ \ \ \ \text{y} \ \ \ \ \ \ \ \ \ \
p(\sigma^2 \mid \theta,\mathbf{Y})\propto p(\underbrace{\theta}_{fijo},\sigma^2 \mid \mathbf{Y})
\end{align*}

Es decir, para encontrar la distribución posterior condicional de $\theta | \sigma^2$, se utiliza la distribución posterior conjunta y los términos que no dependan de $\theta$ se incorporan en la constante de proporcionalidad. El mismo razonamiento se aplica para encontrar la distribución posterior condicional de $\sigma^2 | \theta$.

\BeginKnitrBlock{proposition}<div class="proposition"><span class="proposition" id="prp:unnamed-chunk-3"><strong>(\#prp:unnamed-chunk-3) </strong></span>La distribución posterior condicional de $\theta$ es

\begin{equation}
(\#eq:PostThetaGibbs)
\theta  \mid  \sigma^2,\mathbf{Y} \sim Normal(\mu_n,\tau_n^2)
\end{equation}
  
En donde las expresiones para $\mu_n$ y $\tau_n^2$ están dadas por \@ref(eq:TauSigman). Por otro lado, la distribución posterior condicional de $\sigma^2$ es 

\begin{equation}
(\#eq:PostSigma2Gibbs)
\sigma^2  \mid  \theta,\mathbf{Y} \sim Inversa-Gamma\left(\dfrac{n_0+n}{2},\dfrac{v_0}{2}\right)
\end{equation}
  
con $v_0=n_0\sigma^2_0+(n-1)S^2+n(\bar{y}-\theta)^2$.</div>\EndKnitrBlock{proposition}
<br>

\BeginKnitrBlock{proof}<div class="proof">\iffalse{} <span class="proof"><em>Prueba. </em></span>  \fi{}Acudiendo a la distribución posterior conjunta, e incorporando los términos que no dependen de $\theta$ en la constante de proporcionalidad, se tiene que

\begin{align*}
p(\theta \mid \sigma^2,\mathbf{Y})&\propto \exp\left\{-\frac{n}{2\sigma^2}(\bar{y}-\theta)^2-\frac{1}{2\tau^2}(\theta-\mu)^2\right\}
\end{align*}

Completando los cuadrados se encuentra una expresión idéntica a la función de distribución de una variable aleatoria con distribución $Normal(\mu_n, \tau^2_n)$. Similarmente, después de un desarrollo algebraico breve, se tiene la distribución posterior condicional de $\sigma^2$.</div>\EndKnitrBlock{proof}
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

\BeginKnitrBlock{example}<div class="example"><span class="example" id="exm:EjeRenal"><strong>(\#exm:EjeRenal) </strong></span>@Efronims consideró un conjunto de datos que muestran la función renal de 157 individuos que se sometieron a una prueba médica exhaustiva en un hospital. Los resultados de la prueba renal están en un intervalo de -6 puntos a 4 puntos. Entre más alto sea el resultado, se concluye que el riñón del individuo es más sano. Nótese que estas pruebas son importantes para predecir el comportamiento de un riñón donado a un paciente con problemas renales. 

En principio, es de interés para el investigador conocer la media y la dispersión de estos datos, para poder analizar a fondo la situación de los pacientes que esperan un transplante. Dado que se trata de una primera aproximación, se prefiere utilizar distribuciones previas no informativas para los parámetros de la media y varianza. Lo anterior se logra en `STAN` definiendo las distribuciones previas de `theta ~ normal(0, 10000)` y de `sigma2 ~ inv_gamma(0.001, 0.001)`. De esta forma, la distribución previa de $\theta$ está centrada en cero, pero con una varianza muy grande al igual que la distribución de la varianza, los cuales representan distribuciones previas no informativas.

El siguiente código en `STAN` muestra cómo se lleva a cabo la inferencia.</div>\EndKnitrBlock{example}


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
## theta  0.0773  0.0089 0.4111 -0.7173 0.8999  2137 1.0004
## sigma  1.5517  0.0083 0.3157  1.0826 2.3188  1441 1.0021
## sigma2 2.5073  0.0321 1.1161  1.1721 5.3768  1208 1.0025
## 
## Samples were drawn using NUTS(diag_e) at Mon Jun 14 23:55:14 2021.
## For each parameter, n_eff is a crude measure of effective sample size,
## and Rhat is the potential scale reduction factor on split chains (at 
## convergence, Rhat=1).
```

Después de ejecutar las iteraciones necesarias, la salida del anterior código muestra una estimación puntual para la esperanza de $Y$ de 0.09 con un intervalo de credibilidad del 95\% dado por (-0.71, 0.91). Por otro lado, la estimación puntual de la desviación estándar de $Y$ es de 1.55 con un intervalo de credibilidad del 95\% dado por (1.09, 2.31).

Las figuras \@ref(fig:posNormalMVStan2) muestra la distribución posterior para este ejemplo, junto con la estimación puntual, correspondiente a la desviación estándar. 

<div class="figure" style="text-align: center">
<img src="4Multiparametricos_files/figure-html/posNormalMVStan2-1.svg" alt="Distribuciones posteriores." width="576" />
<p class="caption">(\#fig:posNormalMVStan2)Distribuciones posteriores.</p>
</div>

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
## [1] 0.08134119
```

```r
quantile(theta.pos, c(0.025, 0.975))
```

```
##       2.5%      97.5% 
## -0.7351317  0.8825539
```

```r
mean(sqrt(sigma2.pos))
```

```
## [1] 1.54241
```

```r
quantile(sqrt(sigma2.pos), c(0.025, 0.975))
```

```
##     2.5%    97.5% 
## 1.020166 2.347198
```

De donde podemos concluir que los resultados de este algoritmo de Gibbs coinciden con la estimación puntual para la esperanza de $Y$, su desviación estándar, y sus intervalos de credibilidad del 95%, obtenidos con `STAN`. 

<div class="figure" style="text-align: center">
<img src="4Multiparametricos_files/figure-html/plotsNormalMV-1.svg" alt="Convergencia de las distribuciones posteriores y diagramas de la función de autocorrelación en las cadenas." width="768" />
<p class="caption">(\#fig:plotsNormalMV)Convergencia de las distribuciones posteriores y diagramas de la función de autocorrelación en las cadenas.</p>
</div>

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
## [1] 0.1030324
```

```r
sd(y.bar)
```

```
## [1] 0.7967851
```

```r
quantile(y.bar,c(0.025,0.975))
```

```
##      2.5%     97.5% 
## -1.514231  1.553076
```

Podemos ver que se espera que el promedio de las pruebas en 5 nuevos pacientes es de 0.1030324, con un intervalo de credibilidad del 95% que es mucho más ancho que el de $\theta$, pues naturalmente $\bar{Y}$ tiene mayor incertidumbre que los parámetros del modelo; además, el tamaño de nuevos datos es de cinco, el cual es pequeño y hace que el pronóstico para $\bar{Y}^*$ no sea muy preciso.

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

\BeginKnitrBlock{proposition}<div class="proposition"><span class="proposition" id="prp:unnamed-chunk-10"><strong>(\#prp:unnamed-chunk-10) </strong></span>La distribución conjunta previa de los parámetros $\theta$ y $\sigma^2$ está dada por
\begin{equation*}
\theta,\sigma^2 \sim Normal-Inversa-Gamma\left(\mu, c_0, \frac{n_0+1}{2},\frac{n_0\sigma^2_0}{2}\right)
\end{equation*}</div>\EndKnitrBlock{proposition}
<br>

\BeginKnitrBlock{proof}<div class="proof">\iffalse{} <span class="proof"><em>Prueba. </em></span>  \fi{}\begin{align*}
p(\theta,\sigma^2)&=p(\sigma^2)p(\theta \mid \sigma^2)\\
&\propto (\sigma^2)^{-\frac{n_0}{2}-1}\exp\left\{-\dfrac{n_0\sigma_0^2}{2\sigma^2}\right\}
(\sigma^2)^{-\frac{1}{2}}\exp\left\{-\frac{c_0}{2\sigma^2}(\theta-\mu)^2\right\}\\
&= (\sigma^2)^{-\frac{n_0+1}{2}-1}\exp\left\{-\frac{1}{2\sigma^2}\left[n_0\sigma^2_0+c_0(\theta-\mu)^2\right]\right\}
\end{align*}
la cual corresponde a la forma de la función de densidad de la distribución Normal-Inversa-Gamma.</div>\EndKnitrBlock{proof}
<br>

Una vez encontrada la distribución conjunta previa, procedemos a encontrar la distribución conjunta posterior, y así poder encontrar las estimaciones de $\theta$ y $\sigma^2$.

\BeginKnitrBlock{proposition}<div class="proposition"><span class="proposition" id="prp:unnamed-chunk-12"><strong>(\#prp:unnamed-chunk-12) </strong></span>La distribución posterior conjunta de los parámetros $\theta$ y $\sigma^2$ está dada por
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
\end{equation*}</div>\EndKnitrBlock{proposition}
<br>


\BeginKnitrBlock{proof}<div class="proof">\iffalse{} <span class="proof"><em>Prueba. </em></span>  \fi{}En primer lugar, recordamos que la función de verosimilitud de la muestra está dada por
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
\end{align*}</div>\EndKnitrBlock{proof}
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

\BeginKnitrBlock{proposition}<div class="proposition"><span class="proposition" id="prp:unnamed-chunk-14"><strong>(\#prp:unnamed-chunk-14) </strong></span>La distribución posterior de $\theta$ condicional a $\sigma^2,\mathbf{Y}$ está dada por
\begin{equation*}
\theta \mid \sigma^2,\mathbf{Y} \sim Normal(\mu_n,\sigma^2/(n+c_0))
\end{equation*}
con $\mu_n=\dfrac{n\bar{y}+c_0\mu}{n+c_0}$.</div>\EndKnitrBlock{proposition}
<br>

\BeginKnitrBlock{proof}<div class="proof">\iffalse{} <span class="proof"><em>Prueba. </em></span>  \fi{}Acudiendo a la distribución posterior conjunta dada en \@ref(eq:desarro1), tenemos que
\begin{align*}
p(\theta \mid \sigma^2,\mathbf{Y})&\propto 
p(\theta,\sigma^2 \mid \mathbf{Y}) \\
&\propto(\sigma^2)^{-\frac{n_0+n+1}{2}-1}\\
&\hspace{1cm}\times
\exp\left\{-\frac{1}{2\sigma^2}\left[n_0\sigma^2_0+(n-1)S^2+(c_0+n)(\theta-\mu_n)^2+\frac{c_0n}{c_0+n}(\mu-\bar{y})^2\right]\right\}\\
&\propto \exp\left\{-\frac{1}{2\sigma^2}(c_0+n)(\theta-\mu_n)^2\right\}
\end{align*}
la cual corresponde a la forma de la función de densidad de la distribución $Normal(\mu_n, \sigma^2/(n+c_0))$.</div>\EndKnitrBlock{proof}

En el anterior resultado, la media de la distribución condicional posterior $\mu_n$ se puede escribir como $\mu_n=\frac{n}{n+c_0}\bar{y}+\frac{c_0}{n+c_0}\mu$, promedio ponderado entre la estimación clásica $\bar{y}$ y la estimación previa $\mu$. Observando que las ponderaciones $\frac{n}{n+c_0}$ y $\frac{c_0}{n+c_0}$ forman una combinación lineal convexa, se puede definir a $c_0$ como el número de observaciones en la información previa. De esta forma, las ponderaciones de la estimación clásica y la estimación previa dependerán directamente de los tamaños muestrales respectivos.

\BeginKnitrBlock{proposition}<div class="proposition"><span class="proposition" id="prp:PosterSigma2IG"><strong>(\#prp:PosterSigma2IG) </strong></span>La distribución marginal posterior del parámetro $\sigma^2$ es
\begin{equation*}
\sigma^2 \mid \mathbf{Y} \sim Inversa-Gamma\left(\frac{n+n_0}{2},\frac{(n+n_0)\sigma^2_n}{2}\right)
\end{equation*}

Donde $(n+n_0)\sigma^2_n=n_0\sigma^2_0+(n-1)S^2+\frac{c_0n}{c_0+n}(\mu-\bar{y})^2$ corresponde a una suma ponderada de la varianza previa, la varianza muestral y la diferencia entre la media muestral y la media previa.</div>\EndKnitrBlock{proposition}
<br>

\BeginKnitrBlock{proof}<div class="proof">\iffalse{} <span class="proof"><em>Prueba. </em></span>  \fi{}De la distribución posterior conjunta \@ref(eq:desarro1), e integrando con respecto a $\theta$, se tiene que

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

la cual corresponde a la forma funcional de la densidad $Inversa-Gamma(\frac{n+n_0}{2},\frac{(n+n_0)\sigma^2_n}{2})$.</div>\EndKnitrBlock{proof}
<br>

Dadas las distribuciones $p(\sigma^2\mid \mathbf{Y})$ y $p(\theta\mid \sigma^2, \mathbf{Y})$, podemos proceder de la siguiente forma para obtener valores simulados de $\theta$ y $\sigma^2$ y por consiguiente sus estimaciones puntuales. Si el número de iteraciones se fija como $G$, entonces se debe seguir el siguiente algoritmo:

1. Simular $G$ valores de la distribución de $\sigma^2|\mathbf{Y}$; es decir, de la distribución $Inversa-Gamma$ encontrada en el anterior resultado. Estos valores se denotan por $\sigma^2_{(1)},\sigma^2_{(2)},\cdots,\sigma^2_{(G)}$.
2. Para cada valor de $\sigma^2_{(g)}$, con $g=1,\cdots,G$, simular un valor de la distribución de $\theta|\sigma^2,\mathbf{Y}$; es decir, de la distribución $N(\mu_n,\sigma^2/(n+c_0))$, donde $\sigma^2$ se reemplaza por $\sigma^2_{(g)}$. De esta forma, se obtiene los valores $\theta_{(1)},\theta_{(2)},\cdots,\theta_{(G)}$.

Es claro que en el anterior algoritmo, no es necesario fijar algún valor inicial para $\theta$ o para $\sigma^2$. De la misma forma no existirán correlaciones sustantivas entre los valores simulados para ningún parámetro. Por lo tanto, estos valores se pueden usar directamente para el cálculo de las estimaciones, y no es necesario descartar los primeros valores simulados, ni realizar la fase de *thinning*.

Ahora bien, existe otra alternativa para obtener la estimación de $\theta$ y $\sigma^2$, la cual depende directamente de la distribución posterior de cada parámetro. En efecto, la distribución posterior de $\sigma^2$ ya se encontró en el resultado \@ref(prp:PosterSigma2IG), resta encontrar la distribución posterior de $\theta$, la cual se presenta en el siguiente resultado. 

\BeginKnitrBlock{proposition}<div class="proposition"><span class="proposition" id="prp:PosThetaTnoestandar"><strong>(\#prp:PosThetaTnoestandar) </strong></span>La distribución posterior del parámetro $\theta$ corresponde a una $t$ no estandarizada con $n_0+n$ grados de libertad, parámetro de localización $\mu_n=\dfrac{n\bar{Y}+c_0\mu}{n+c_0}$ y parámetro de escala $\dfrac{\sigma_n}{\sqrt{c_0+n}}$, con $(n+n_0)\sigma^2_n=n_0\sigma^2_0+(n-1)S^2+\frac{c_0n}{c_0+n}(\mu-\bar{y})^2$. Esto es, 
\begin{equation*}
\theta \mid \mathbf{Y} \sim t_{n+n_0}\left(\mu_n, \frac{\sigma^2_n}{c_0+n}\right)
\end{equation*}</div>\EndKnitrBlock{proposition}
<br>

\BeginKnitrBlock{proof}<div class="proof">\iffalse{} <span class="proof"><em>Prueba. </em></span>  \fi{}Partiendo de la distribución posterior conjunta e integrando con respecto a $\sigma^2$, se tiene que

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

la cual corresponde a la forma de la función de densidad de la distribución deseada.</div>\EndKnitrBlock{proof}
<br>

Las distribuciones encontradas en los resultados \@ref(prp:PosterSigma2IG) y \@ref(prp:PosThetaTnoestandar) permiten estimar directamente los parámetros $\theta$ y $\sigma^2$ usando las esperanzas teóricas de las distribuciones posteriores. Por ende, las estimaciones puntuales son:

\begin{align*}
\hat{\theta}
&=\mu_n=\dfrac{n\bar{Y}+c_0\mu}{n+n_0}\\
\hat{\sigma}^2
&=\dfrac{(n+n_0)\sigma^2_n/2}{(n+n_0)/2-1}=\dfrac{(n+n_0)\sigma^2_n}{n+n_0-2}\approx\sigma^2_n=\dfrac{n_0\sigma^2_0+(n-1)S^2+\frac{c_0n}{c_0+n}(\mu-\bar{y})^2}{n+n_0}
\end{align*}

Los intervalos del $(1-\alpha)\times 100\%$ de credibilidad para $\theta$ y $\sigma^2$ se construyen usando los percentiles $\alpha/2$ y $1-\alpha/2$ de las respectivas distribuciones posteriores dadas en los resultados mencionados anteriormente. Ilustramos el uso de la metodología en el siguiente ejemplo.

\BeginKnitrBlock{example}<div class="example"><span class="example" id="exm:Eje2.1.2"><strong>(\#exm:Eje2.1.2) </strong></span>Para los datos de función renal \cite{Efronims} que se muestran en el Ejemplo \@ref(exm:EjeRenal), suponga que la información previa está contenida en la medición de función renal para una muestra de 12 pacientes dadas por: -1.3619, -1.1116, -0.4744, -0.5663, 2.2056, 0.9491, 0.2298, -0.7933, 1.0198, -0.9850, 3.5679 y -1.9504. La media y la varianza muestral de estas 12 observaciones corresponden a $\mu=0.060775$ y $\sigma^2_0=2.598512$; por consiguiente $c_0=n_0=12$. 

Por otro lado, la media y la varianza muestral de los 15 pacientes en la información actual son $\bar{y}=0.08349249$ y $S^2=2.301684$. De esta forma, los parámetros de las distribuciones marginales posteriores de $\theta$ y $\sigma^2$ se pueden calcular como $\mu_n=\frac{15}{15+12}\times 0.08349249+\frac{12}{15+12}\times 0.060775=0.07339583$ y $$\sigma^2_n=\dfrac{12*2.598512+14*2.301684+6.666667*(0.060775-0.08349249)^2}{15+12}=2.348487$$ En conclusión, las distribuciones marginales posterior de $\theta$ y $\sigma^2$ están dadas por
\begin{equation*}
\theta|\mathbf{Y}\sim t_{27}(0.07339583,2.348487/27=0.086981)
\end{equation*}

y
\begin{equation*}
\sigma^2|\mathbf{Y}\sim Inversa-Gamma(27/2=13.5,\ 27*2.348487/2=31.70457)
\end{equation*}

Así, la estimación Bayesiana de $\theta$ es $\mu_n=0.073$ y un intervalo de credibilidad del $95\%$ para $\theta$ se puede calcular como $0.073\pm t_{27,0.975}*\sqrt{0.086981}=(-0.53,\ 0.68)$. Por otro lado, la estimación Bayesiana de $\sigma^2$ está dada por $31.70457/(13.5-1)=2.53$, y un intervalo de credibilidad del $95\%$ para $\sigma^2$ se puede calcular como los percentiles $2.5\%$ y $97.5\%$ de la distribución $Inversa-Gamma(13.5,\ 31.70457)$, dados por $(1.5, 4.4)$. Los anteriores cálculos se ilustran en el siguiente código `R`.</div>\EndKnitrBlock{example}
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
## [1] 0.07384158
```

```r
mean(sigma2.res)
```

```
## [1] 2.532906
```

```r
# Intervalos de credibilidad del 95%
quantile(theta.res, c(0.025, 0.975))
```

```
##       2.5%      97.5% 
## -0.5472775  0.6694128
```

```r
quantile(sigma2.res, c(0.025, 0.975))
```

```
##     2.5%    97.5% 
## 1.453781 4.409619
```




































