


# Modelos uniparamétricos

Los modelos que están definidos en términos de un solo parámetro que
pertenece al conjunto de los números reales se definen como modelos
*uniparamétricos*. Este capítulo estudia modelos, discretos y continuos,
que son comunes de implementar en la práctica. Dado que todos ellos son
inducidos por familias de probabilidad conjugadas, entonces las
estimaciones posteriores para los parámetros pueden hallarse sin
necesidad de sofisticaciones computacionales. Es decir, con el uso de
una simple calculadora de bolsillo, es posible realizar inferencia
bayesiana propiamente dicha. Por lo tanto, en este capítulo, será menor
el uso de software estadístico. Sin embargo, para cada modelo se incluye
la sintaxis de programación en `R` y en `STAN` junto con ejemplos
prácticos que permiten la familiarización e interiorización del ambiente
computacional de este software que será indispensable en el desarrollo
de capítulos posteriores.

## Modelo Bernoulli

Suponga que $Y$ es una variable aleatoria con distribución Bernoulli
dada por:

\begin{equation}
p(Y \mid \theta)=\theta^y(1-\theta)^{1-y}I_{\{0,1\}}(y)
\end{equation}

Como el parámetro $\theta$ está restringido al espacio $\Theta=[0,1]$,
entonces es posible formular varias opciones para la distribución previa
del parámetro. En particular, la distribución uniforme restringida al
intervalo $[0,1]$ o la distribución Beta parecen ser buenas opciones.
Puesto que la distribución uniforme es un caso particular de la
distribución Beta, entonces se iniciará con ésta. Por lo tanto la
distribución previa del parámetro $\theta$ estará dada por


\begin{equation}
(\#eq:betadistribution)
p(\theta \mid \alpha,\beta)=
\frac{1}{Beta(\alpha,\beta)}\theta^{\alpha-1}(1-\theta)^{\beta-1}I_{[0,1]}(\theta).
\end{equation}

Bajo este marco de referencia se tienen los siguientes resultados

\BeginKnitrBlock{proposition}<div class="proposition"><span class="proposition" id="prp:unnamed-chunk-1"><strong>(\#prp:unnamed-chunk-1) </strong></span>La distribución posterior del parámetro $\theta$ sigue una distribución
\begin{equation*}
\theta \mid Y \sim Beta(y+\alpha,\beta-y+1)
\end{equation*}</div>\EndKnitrBlock{proposition}

\BeginKnitrBlock{proof}<div class="proof">\iffalse{} <span class="proof"><em>Prueba. </em></span>  \fi{}\begin{align*}
p(\theta \mid Y)&\propto 
p(Y \mid \theta)p(\theta \mid \alpha,\beta)\\
&=\frac{I_{\{0,1\}}(y)}{Beta(\alpha,\beta)}\theta^y\theta^{\alpha-1}(1-\theta)^{\beta-1}(1-\theta)^{1-y}I_{[0,1]}(\theta)\\
&\propto \theta^{y+\alpha-1}(1-\theta)^{\beta-y+1-1}I_{[0,1]}(\theta)
\end{align*}
    
Por lo tanto, factorizando convenientemente, se encuentra una expresión idéntica a la función de distribución de una variable aleatoria con distribución $Beta(y+\alpha,\beta-y+1)$.</div>\EndKnitrBlock{proof}
<br>

Del anterior resultado, podemos ver que la familia de distribuciones
Beta es conjugada con respecto a la familia de distribuciones Bernoulli.
Ahora consideremos cuál sería la distribución previa no informativa de
Jeffreys para el parámetro $\theta$. De acuerdo a la definición
\@ref(def:Jeffreys), se tiene que


\begin{equation*}
p(\theta)\propto I(\theta) ^{1/2}
\end{equation*}

En donde $I(\theta)$ es la información de Fisher del parámetro $\theta$,
que en este caso está dada por

\begin{align*}
I(\theta)&=
-E\left\{\dfrac{\partial^2}{\partial\theta^2}\log{p(\mathbf{Y}\mid\theta)}\right\}\\
&=-E\left\{\dfrac{\partial^2}{\partial\theta^2}\{Y\log\theta+(1-Y)\log(1-\theta)\}\right\}\\
&=E\left\{\frac{Y}{\theta^2}+\frac{1-Y}{(1-\theta)^2}\right\}\\
&=\frac{1}{\theta(1-\theta)}
\end{align*}

De esta forma, la distribución previa no informativa de Jeffreys debe
ser proporcional a $\theta^{-1/2}(1-\theta)^{-1/2}$, que asimismo
corresponde a la distribución $Beta(1/2,1/2)$, cuya función de densidad
se muestra en la figura \@ref(fig:jefber1) la cual asigna iguales pesos
a los valores extremos del parámetro de interés y su característica de
ser no informativa se representa en la simetría de la función alrededor
del valor 0.5.

<div class="figure" style="text-align: center">
<img src="3Uniparametricos_files/figure-html/jefber1-1.svg" alt="Distribución previa no informativa de Jeffreys para el parámetro de una distribución Bernoulli" width="576" />
<p class="caption">(\#fig:jefber1)Distribución previa no informativa de Jeffreys para el parámetro de una distribución Bernoulli</p>
</div>

\BeginKnitrBlock{proposition}<div class="proposition"><span class="proposition" id="prp:unnamed-chunk-3"><strong>(\#prp:unnamed-chunk-3) </strong></span>La distribución predictiva previa para una observación $y$ está dada por
\begin{equation}
(\#eq:Predipreviabernou)
p(Y)=
\frac{Beta(y+\alpha,\beta-y+1)}{Beta(\alpha,\beta)}I_{\{0,1\}}(y)
\end{equation}
La cual define una auténtica función de densidad de probabilidad continua.</div>\EndKnitrBlock{proposition}
<br>

\BeginKnitrBlock{proof}<div class="proof">\iffalse{} <span class="proof"><em>Prueba. </em></span>  \fi{}De la definición de función de distribución predictiva se tiene que

\begin{align*}
p(Y)&=\int p(Y \mid \theta)p(\theta \mid \alpha,\beta)\ d\theta\\
&=\int_0^1 \theta^y(1-\theta)^{1-y}I_{\{0,1\}}(y)\frac{1}{Beta(\alpha,\beta)}\theta^{\alpha-1}(1-\theta)^{\beta-1}\ d\theta\\
&=\frac{Beta(y+\alpha,\beta-y+1)}{Beta(\alpha,\beta)}I_{\{0,1\}}(y)
    \int_0^1\frac{\theta^{y+\alpha-1}(1-\theta)^{\beta-y+1-1}}{Beta(y+\alpha,\beta-y+1)}\ d\theta\\
&=\frac{Beta(y+\alpha,\beta-y+1)}{Beta(\alpha,\beta)}I_{\{0,1\}}(y)
\end{align*}

Nótese que en la anterior expresión, la integral al lado derecho de la tercera igualdad es igual a la unidad, puesto que la expresión matemática dentro de la integral corresponde a la función de densidad de una variable aleatoria con distribucion $Beta$, que tiene rango en el intervalo $(0,1)$. Por otro lado se deben verificar las dos condiciones de función de densidad. Es decir

1. $p(Y)>0 \ \forall y \in Y)$. Esta condición se tiene trivialmente puesto que la función matemática Beta siempre toma valores positivos.
2. $\int p(y)\ dx=1$. En este caso, esta función es discreta definida en el conjunto $\{0,1\}$. Por lo tanto esta condición es equivalente a
\begin{equation*}
\sum_{y\in{\{0,1\}}}P(Y=y)=
  \sum_{y\in{\{0,1\}}}\frac{Beta(y+\alpha,\beta-y+1)}{Beta(\alpha,\beta)}
=1
\end{equation*}

Lo cual se verifica fácilmente teniendo en cuenta las propiedades de la función matemática Beta y de la función matemática Gamma.</div>\EndKnitrBlock{proof}
<br>

La distribución predictiva dada en \@ref(eq:Predipreviabernou) está
basada únicamente en la distribución previa del parámetro $\theta$. Una
vez observada la variable $Y$ se puede pensar en actualizar la
distribución predictiva basando la inferencia en la distribución
posterior del parámetro; esta distribución se da en el siguiente
resultado.

\BeginKnitrBlock{proposition}<div class="proposition"><span class="proposition" id="prp:unnamed-chunk-5"><strong>(\#prp:unnamed-chunk-5) </strong></span>Después de la recolección de los datos, la distribución predictiva posterior para una nueva observación $\tilde{y}$ está dada por

\begin{equation}
p(\tilde{y} \mid Y)=
  \frac{Beta(\tilde{y}+y+\alpha,\beta-\tilde{y}-y+2)}{Beta(y+\alpha,\beta-y+1)}I_{\{0,1\}}(\tilde{y}),
\end{equation}</div>\EndKnitrBlock{proposition}
<br>

\BeginKnitrBlock{proof}<div class="proof">\iffalse{} <span class="proof"><em>Prueba. </em></span>  \fi{}De la definición de función de distribución predictiva se tiene que

\begin{align*}
p(\tilde{y} \mid Y)&=\int p(\tilde{y} \mid \theta)p(\theta \mid Y)\ d\theta\\
&=\int_0^1\theta^{\tilde{y}}(1-\theta)^{1-\tilde{y}}I_{\{0,1\}}(\tilde{y})\frac{\theta^{y+\alpha-1}(1-\theta)^{\beta-y+1-1}}{Beta(y+\alpha,\beta-y+1)}\ d\theta\\
&=\frac{Beta(\tilde{y}+y+\alpha,\beta-\tilde{y}-y+2)}{Beta(y+\alpha,\beta-y+1)}I_{\{0,1\}}(\tilde{y})\\
&\hspace{2cm}\times \int_0^1\frac{\theta^{\tilde{y}+y+\alpha-1}(1-\theta)^{\beta-\tilde{y}-y+2-1}}{Beta(\tilde{y}+y+\alpha,\beta-\tilde{y}-y+2)}\ d\theta\\
&=\frac{Beta(\tilde{y}+y+\alpha,\beta-\tilde{y}-y+2)}{Beta(y+\alpha,\beta-y+1)}I_{\{0,1\}}(\tilde{y})
\end{align*}</div>\EndKnitrBlock{proof}
<br>

En la práctica rara vez se observa la realización de una única variable
aleatoria Bernoulli $Y$, sino una muestra de variables aleatorias $Y_1$,
$\cdots$, $Y_n$. En este caso, la distribución posterior del parámetro
$\theta$ está dada en el siguiente resultado.

\BeginKnitrBlock{proposition}<div class="proposition"><span class="proposition" id="prp:unnamed-chunk-7"><strong>(\#prp:unnamed-chunk-7) </strong></span>Cuando se tiene una muestra aleatoria $Y_1,\ldots,Y_n$ de variables con distribución Bernoulli de parámetro $\theta$, entonces la distribución posterior del parámetro de interés es

\begin{equation*}
\theta \mid Y_1,\ldots,Y_n \sim Beta\left(\sum_{i=1}^ny_i+\alpha,\beta-\sum_{i=1}^ny_i+n\right)
\end{equation*}</div>\EndKnitrBlock{proposition}

\BeginKnitrBlock{example}<div class="example"><span class="example" id="exm:bernoelectoral"><strong>(\#exm:bernoelectoral) </strong></span>Es común en muchos países del mundo que se presenten encuestas de opinión electoral unas semanas antes de las elecciones presidenciales. Dentro de este tipo de encuestas se acostumbra a indagar acerca del favoritismo de los candidatos involucrados en la contienda electoral. Suponga que el candidato presidencial A está interesado en conocer su intención de voto previa a las elecciones. Para esto, él contrata a una firma encuestadora para la realización de una encuesta entre la población votante. El resultado de este estudio puede hacer cambiar o afirmar las estrategias publicitarias y la redefinición de la campaña electoral. La firma encuestadora decide implementar una estrategia de muestreo con un tamaño de muestra de doce mil personas. A cada respondiente se le realiza la siguiente pregunta: 
  
  > Si las elecciones presidenciales fueran mañana. ¿Usted votaría por el candidato A?
    
Las respuestas a esta pregunta son realizaciones de una muestra aleatoria de doce mil variables con densidad Bernoulli. Los resultados del estudio arrojan que 6360 personas de las personas entrevistadas, es decir un 53%, votarían por el suscrito candidato. Técnicamente se debe analizar esta cifra puesto que las implicaciones de ganar en una primera vuelta son grandes en el sentido económico, logístico y administrativo. Claramente, el dato 53% asegura una ventaja dentro de la muestra de doce mil personas. Sin embargo, es necesario realizar un estudio más profundo acerca de la caracterización estructural de la intención de voto del candidato en el electorado.
    
Con base en lo anteriormente expuesto, se decide utilizar la inferencia bayesiana puesto que existe información previa de un estudio anterior, contratado por el mismo candidato unos meses atrás en donde se entrevistaron a mil personas, con un favoritismo que estaba alrededor del 35 por ciento. Esta situación conlleva a la utilización de la metodología bayesiana que incorpora la información pasada acerca del mismo fenómeno.
    
El estadístico de la firma encuestadora decide utilizar una distribución previa Beta, definiendo los parámetros de la distribución previa como $\alpha$ igual al número de votantes a favor y $\beta$ igual al número de votantes en contra. Es decir, $Beta(\alpha=350, \beta=650)$. Por lo anterior, la distribución posterior del parámetro de interés, que representa la probabilidad de éxito en las elecciones presidenciales, es $Beta(6360+350, 650-6360+12000)=Beta(6710, 6290)$. Por lo tanto, utilizando la distribución posterior, se estima que la intención de voto por el candidato es de $\frac{6710}{6710+6290}=\frac{6710}{13000}=0.516$ y este valor equivale a la media de la distribución posterior. 
    
Sin embargo, si no se tuviese información previa como la suministrada por el estudio de meses anteriores, el análisis bayesiano sugeriría trabajar con una distribución previa no informativa, que en este caso, correspondería a una $Beta(\alpha=0.5, \beta=0.5)$. siguiendo el mismo análisis, se tiene que la distribución posterior es $Beta(6360.5, 5640.5)$. Finalmente, se estimaría que la intención de voto por el candidato es de $\frac{6350.5}{12001}=0.529$. 

La figuras \@ref(fig:BernoEj1) muestra el comportamiento de las distribuciones previas y posteriores en ambos escenarios. Nótese que la distribución no informativa influye muy poco en el comportamiento de la distribución posterior.</div>\EndKnitrBlock{example}

<div class="figure" style="text-align: center">
<img src="3Uniparametricos_files/figure-html/BernoEj1-1.svg" alt="Distribuciones previas (línea punteada) y posteriores (línea sólida) para el ejemplo de las encuestas electorales." width="576" />
<p class="caption">(\#fig:BernoEj1)Distribuciones previas (línea punteada) y posteriores (línea sólida) para el ejemplo de las encuestas electorales.</p>
</div>

Utilizando el siguiente código en R, es posible conocer los intervalos
de credibilidad para las dos distribuciones posteriores. Además, es
posible concluir que, en ambos escenarios, el candidato aventaja
significativamente a sus contrincantes y, salvo algún cambio drástico en
el comportamiento del electorado, ganará las elecciones. Lo anterior se
deduce puesto que el intervalo de credibilidad al 95 % no contiene
ningún valor menor a 0.5


```r
qbeta(c(0.025, 0.975), 6710, 6290)
```

```
## [1] 0.5075614 0.5247415
```

```r
qbeta(c(0.025, 0.975), 6350.5, 5640.5)
```

```
## [1] 0.5206678 0.5385340
```

Por otro lado, el siguiente código en `STAN` permite obtener el mismo
tipo de inferencia creando cuatro cadenas cuya distribución de
probabilidad coincide con la distribución posterior del ejemplo.


```r
Bernoulli <- "
data {
  int<lower=0> n;
  int y[n];
}
parameters {
  real<lower=0, upper=1> theta;
}
model {
  y ~ bernoulli(theta);
  theta ~ beta(350, 650);
}

"

library(rstan)
options(mc.cores = parallel::detectCores())

n <- 12000
s <- 6350
y <- c(rep(1, s), rep(0, n - s))
sample_data <- list(n = n, y = y)

Berfit <- stan(model_code = Bernoulli, 
               data = sample_data, verbose = FALSE)
```

La siguiente salida de `STAN` permite conocer la estimación bayesiana posterior y los límites del intervalo de credibilidad al 95%.


```r
print(Berfit, pars = "theta", 
      digits = 4, probs = c(0.025, 0.975))
```

```
## Inference for Stan model: 62d7c91114fcc6227deb68f6059b2b09.
## 4 chains, each with iter=2000; warmup=1000; thin=1; 
## post-warmup draws per chain=1000, total post-warmup draws=4000.
## 
##         mean se_mean     sd   2.5%  97.5% n_eff   Rhat
## theta 0.5154   1e-04 0.0043 0.5065 0.5233  1478 1.0014
## 
## Samples were drawn using NUTS(diag_e) at Sun Jun  6 23:41:45 2021.
## For each parameter, n_eff is a crude measure of effective sample size,
## and Rhat is the potential scale reduction factor on split chains (at 
## convergence, Rhat=1).
```

La figura \@ref(fig:posBernoulliStan) muestra la distribución posterior para este ejemplo, junto con la estimación puntual, correspondiente a la media. 


```r
bayesplot::mcmc_areas(Berfit, pars = "theta", 
                      prob = 0.95)
```

<div class="figure" style="text-align: center">
<img src="3Uniparametricos_files/figure-html/posBernoulliStan-1.svg" alt="Distribución posterior." width="576" />
<p class="caption">(\#fig:posBernoulliStan)Distribución posterior.</p>
</div>

## Modelo Binomial

Cuando se dispone de una muestra aleatoria de variables con distribución
Bernoulli $Y_1,\ldots,Y_n$, la inferencia bayesiana se puede llevar a
cabo usando la distribución Binomial, puesto que es bien sabido que la
suma de variables aleatorias Bernoulli

\begin{equation*}
S=\sum_{i=1}^nY_i
\end{equation*}

sigue una distribución Binomial. Es decir:

\begin{equation}
p(S \mid \theta)=\binom{n}{s}\theta^s(1-\theta)^{n-s}I_{\{0,1,\ldots,n\}}(s),
\end{equation}

Nótese que la distribución binomial es un caso general para la
distribución Bernoulli, cuando $n=1$. Entonces, así como en la
distribución Bernoulli, el parámetro $\theta$ está restringido al
espacio $\Theta=[0,1]$. Luego, es admisible proponer que $\theta$ siga
una distribución Beta. Por tanto la distribución previa del parámetro
$\theta$ está dada por la expresión \@ref(eq:betadistribution). Bajo
este marco de referencia se tienen los siguientes resultados


\BeginKnitrBlock{proposition}<div class="proposition"><span class="proposition" id="prp:unnamed-chunk-11"><strong>(\#prp:unnamed-chunk-11) </strong></span>La distribución posterior del parámetro $\theta$ sigue una distribución

\begin{equation*}
\theta \mid S \sim Beta(s+\alpha,\beta-s+n)
\end{equation*}</div>\EndKnitrBlock{proposition}
<br>

\BeginKnitrBlock{proof}<div class="proof">\iffalse{} <span class="proof"><em>Prueba. </em></span>  \fi{}\begin{align*}
p(\theta \mid S)&\propto p(S \mid \theta)p(\theta \mid \alpha,\beta)\\
&=\frac{\binom{n}{s}I_{\{0,1,\ldots,n\}}(s)}{Beta(\alpha,\beta)}
\theta^s\theta^{\alpha-1} (1-\theta)^{\beta-1}(1-\theta)^{n-s}I_{[0,1]}(\theta)\\
&\propto \theta^{s+\alpha-1} (1-\theta)^{\beta-s+n-1}I_{[0,1]}(\theta)
\end{align*}

Por lo tanto, factorizando convenientemente, se llega a una expresión idéntica a la función de distribución de una variable aleatoria con distribución $Beta(s+\alpha,\beta-s+n)$.</div>\EndKnitrBlock{proof}
<br>

Del resultado anterior podemos ver que el estimador bayesiano de
$\theta$ está dada por la esperanza de la distribución posterior, dada
por

\begin{equation}
(\#eq:estibeta)
\hat{\theta}_{B}=\frac{s+\alpha}{n+\alpha+\beta}
\end{equation}

En la práctica, se acostumbra a escoger los hiperparámetros $\alpha$ y
$\beta$ de tal forma que correspondan respectivamente al número de
éxitos y fracasos obtenidos en datos que pudieron ser recolectados
previamente. De esta forma, $\hat{\theta}_{P}=\alpha/(\alpha+\beta)$
corresponde a la estimación previa del parámetro $\theta$. Por otro
lado, el estimador clásico de $\theta$ está dado por
$\hat{\theta}_{C}=s/n$. Entonces es posible notar que el estimador
bayesiano de $\theta$ en \@ref(eq:estibeta) de alguna forma combina el
estimador clásico con el estimador previo. Más aún, se puede ver que
$\hat{\theta}_{B}$ se puede escribir como un promedio ponderado entre la
estimación clásica y la estimación previa. Puesto que

\begin{align*}
\hat{\theta}_{B}=\frac{s+\alpha}{n+\alpha+\beta}&=\frac{s}{n+\alpha+\beta}+\frac{\alpha}{n+\alpha+\beta}\\
&=\frac{n}{n+\alpha+\beta}\frac{s}{n}+\frac{\alpha+\beta}{n+\alpha+\beta}\frac{\alpha}{\alpha+\beta}\\
&=\frac{n}{n+\alpha+\beta}\hat{\theta}_{C}+\frac{\alpha+\beta}{n+\alpha+\beta}\hat{\theta}_{P}
\end{align*}

De esta forma, queda en evidencia que la estimación bayesiana de
$\theta$ siempre será un valor intermedio entre la estimación clásica y
la estimación previa. La figura \@ref(fig:beta3) da una ilustración
acerca de la anterior afirmación, en donde se puede observar que para
una distribución previa concentrada en $2/7$ y una función de
verosimilitud^[La función de verosimilitud es una función del parámetro
y sólo se puede graficar una vez se hayan observado las realizaciones de
la variable aleatoria.] con máximo en $8/10$, entonces la distribución
posterior estará centrada en $10/17$; es decir, la estimación bayesiana
se encuentra situada entre la estimación previa y la estimación clásica.

<div class="figure" style="text-align: center">
<img src="3Uniparametricos_files/figure-html/beta3-1.svg" alt="Funciones de verosimilitud, previa y posterior para $\alpha=2$, $\beta=5$, $s=8$ y $n=10$." width="576" />
<p class="caption">(\#fig:beta3)Funciones de verosimilitud, previa y posterior para $\alpha=2$, $\beta=5$, $s=8$ y $n=10$.</p>
</div>

Por otro lado, entre más grande sea el tamaño muestral $n$, más
cercano estará $\hat{\theta}_{B}$ de $\hat{\theta}_{C}$ o
equivalentemente la función de densidad posterior de $\theta$ estará más
concentrada en $s/n$; mientras que entre mayor número de datos tenga la
muestra de la distribución previa ($\alpha+\beta$ = número de datos), más
cercano estará $\hat{\theta}_{B}$ de $\hat{\theta}_{P}$ y la densidad
posterior de $\theta$ estará más concentrada en $\alpha/(\alpha+\beta)$.

Para ilustrar lo anterior, suponga que la distribución previa de
$\theta$ está parametrizada con $\alpha=\beta=5$, es decir la estimacion previa es 0.5, y suponga además que la estimación clásica es 0.33, pero el tamaño muestral $n$ incrementa manteniendo constante la estimacion
clásica. En la figura \@ref(fig:betan) se muestra la estimación posterior de
$\theta$, es evidente que a medida que el tamaño muestral $n$ aumenta,
la estimación posterior se acerca más a la estimación clásica.

<div class="figure" style="text-align: center">
<img src="3Uniparametricos_files/figure-html/betan-1.svg" alt="Estimación posterior de $\theta$ para diferentes valores de $n$ y $s$ con $\alpha=\beta=5$." width="576" />
<p class="caption">(\#fig:betan)Estimación posterior de $\theta$ para diferentes valores de $n$ y $s$ con $\alpha=\beta=5$.</p>
</div>

Anteriormente, se comentó que se acostumbra a escoger los parámetros
$\alpha$ y $\beta$ que correspondan al número de éxitos y fracasos en la
información previa. Sin embargo, la información previa puede no
presentarse de esta forma. Por ejemplo, en algunas situaciones, la
información previa puede proveer el valor de $\theta$, es decir, el
valor de $\hat{\theta}_P$, y el valor de la desviación estándar de la
estimación (comúnmente conocido como el error estándar). Por ejemplo,
suponga que $\hat{\theta}_P=0.5$ con un error estándar de 0.1, entonces
podemos encontrar los valores de $\alpha$ y $\beta$ de las expresiones
$\frac{\alpha}{\alpha+\beta}=0.5$ y
$\sqrt{\frac{\alpha\beta}{(\alpha+\beta)^2(\alpha+\beta+1)}}=0.1$, de
donde se tiene que $\alpha=12$ y $\beta=12$, y la distribución previa
correspondiente $Beta(12, 12)$ tiene una esperanza de 0.05 y una
desviación estándar de 0.1. Se puede ver que entre mayor sea la
desviación estándar, menores resultan los valores de $\alpha$ y $\beta$,
que conducen a una distribución previa menos informativa.

Ahora, se vio anteriomente que la distribución previa no informativa de
Jeffreys corresponde a la distribución $Beta(1/2, 1/2)$, la cual conduce
a la distribución posterior $Beta(s+1/2, n-s+1/2)$, que a su vez nos
lleva al estimador 

\begin{equation}
(\#eq:EstithetaJeffreys)
\hat{\theta}_B=\frac{s+1/2}{n+1}
\end{equation}

La anterior expresión es comparable con el estimador clásico
$\hat{\theta}_C=\frac{s}{n}$, en el sentido de que los dos son
aplicables cuando no se dispone de ninguna información previa. Podemos
observar que, aparte del alto grado de similitud que tienen los dos
estimadores, es preferible usar el estimador \@ref(eq:EstithetaJeffreys) en situaciones donde el valor teórico de $\theta$ es muy pequeño, y
como consecuencia $s=0$ en la muestra. Por ejemplo, cuando $\theta$
representa el porcentaje de personas que están infectados con algún
virus poco común. En estos casos, el estimador clásico
$\hat{\theta}_C=0$ sugiriendo que ningún porcentaje de la población está
infectado, conclusión que puede ser errónea. Por otro lado, el estimador
bayesiano $\hat{\theta}_B=\frac{0.5}{n+1}$ tiende a un porcentaje muy pequeño a medida que aumente el tamaño muestral $n$,
pero nunca llega a ser nulo.

En el siguiente resultado, se encuentra la distribución predictiva
previa para una variable binomial $S$.


\BeginKnitrBlock{proposition}<div class="proposition"><span class="proposition" id="prp:unnamed-chunk-13"><strong>(\#prp:unnamed-chunk-13) </strong></span>La distribución predictiva previa para la observación particular de la suma de variables aleatorias Bernoulli, $s$, está dada por una distribución Beta-Binomial

\begin{equation}
p(S)=\binom{n}{s}\frac{Beta(s+\alpha,\beta-s+n)}{Beta(\alpha,\beta)}I_{\{0,1,\ldots,n\}}(s)
\end{equation}</div>\EndKnitrBlock{proposition}
<br>

\BeginKnitrBlock{proof}<div class="proof">\iffalse{} <span class="proof"><em>Prueba. </em></span>  \fi{}De la definición de función de distribución predictiva previa se tiene que

\begin{align*}
p(S)
&=\int p(S \mid \theta)p(\theta \mid \alpha,\beta)\ d\theta\\
&=\int_0^1\binom{n}{s}\theta^s(1-\theta)^{n-s}I_{\{0,1,\ldots,n\}}(s)
\frac{1}{Beta(\alpha,\beta)}\theta^{\alpha-1}(1-\theta)^{\beta-1}\ d\theta\\
&=\binom{n}{s}\frac{Beta(s+\alpha,\beta-s+n)}{Beta(\alpha,\beta)}I_{\{0,1,\ldots,n\}}(s)\\
&\hspace{2cm}\times\int_0^1\frac{\theta^{s+\alpha-1}(1-\theta)^{\beta-s+n-1}}{Beta(s+\alpha,\beta-s+n)}\ d\theta\\
&=\binom{n}{s}\frac{Beta(s+\alpha,\beta-s+n)}{Beta(\alpha,\beta)}I_{\{0,1,\ldots,n\}}(s)
\end{align*}</div>\EndKnitrBlock{proof}
<br>

Una vez observados los valores muestrales, podemos encontrar la
distribución predictiva posterior para una nueva variable binomial
$\tilde{S}$ en una muestra de tamaño $\tilde{n}$. Esta distribución se
encuentra en el siguiente resultado.

\BeginKnitrBlock{proposition}<div class="proposition"><span class="proposition" id="prp:ResPredBinom"><strong>(\#prp:ResPredBinom) </strong></span>Después de la recolección de los datos $y_1$, $\cdots$, $y_n$, la distribución predictiva posterior para una nueva variable $\tilde{S}$ en una muestra del tamaño $\tilde{n}$ está dada por

\begin{equation}
(\#eq:Binompredict)
p(\tilde{s} \mid S)=\binom{\tilde{n}}{\tilde{s}}\frac{Beta(\tilde{s}+s+\alpha,\beta-\tilde{s}-s+n+\tilde{n})}{Beta(s+\alpha,\beta-s+n)}I_{\{0,1,\ldots,\tilde{n}\}}(\tilde{s}),
\end{equation}</div>\EndKnitrBlock{proposition}
<br>

\BeginKnitrBlock{proof}<div class="proof">\iffalse{} <span class="proof"><em>Prueba. </em></span>  \fi{}De la definición de función de distribución predictiva se tiene que

\begin{align*}
p(\tilde{s} \mid S)&=\int p(\tilde{s} \mid \theta)p(\theta \mid S)\ d\theta\\
&=\int_0^1 \binom{\tilde{n}}{\tilde{s}} \theta^{\tilde{s}}(1-\theta)^{\tilde{n}-\tilde{s}}I_{\{0,1,\ldots,\tilde{n}\}}(\tilde{s})
\frac{\theta^{s+\alpha-1}(1-\theta)^{\beta-s+n-1}}{Beta(s+\alpha,\beta-s+n)}\ d\theta\\
&=\binom{\tilde{n}}{\tilde{s}}\frac{Beta(\tilde{s}+s+\alpha,\beta-\tilde{s}-s+n+\tilde{n})}{Beta(s+\alpha,\beta-s+n)}I_{\{0,1,\ldots,\tilde{n}\}}(\tilde{s})\\
& \hspace{2cm}\times
\int_0^1\frac{\theta^{\tilde{s}+s+\alpha-1}(1-\theta)^{\beta-\tilde{s}-s+n+\tilde{n}-1}}
{Beta(\tilde{s}+s+\alpha,\beta-\tilde{s}-s+n+\tilde{n})}\ d\theta\\
&=\binom{\tilde{n}}{\tilde{s}}\frac{Beta(\tilde{s}+s+\alpha,\beta-\tilde{s}-s+n+\tilde{n})}{Beta(s+\alpha,\beta-s+n)}I_{\{0,1,\ldots,\tilde{n}\}}(\tilde{s})
\end{align*}</div>\EndKnitrBlock{proof}
<br>

En la anterior distribución predictiva, se necesita calcular funciones
Beta. Cuando los tamaños muestrales $n$, $\tilde{n}$ y/o los
parámetros de la distribución previa $\alpha$ y $\beta$ son muy grandes,
`R` puede presentar problemas numéricos al momento de calcular directamente estas funciones. Por ejemplo, supongamos que $n=1000$, $s=650$, $\alpha=200$, $\beta=300$ y $\tilde{n}=800$, de esta forma, los posibles
valores para $\tilde{s}$ son $0,1,\cdots,800$, y se tiene que la
probabilidad de que $\tilde{s}$ tome el valor 500 está dada por

\begin{equation}
(\#eq:Ejebinom)
Pr(\tilde{s}=500|S)=
\binom{800}{500}\frac{Beta(1350,950)}{Beta(850,650)}
\end{equation}

y desafortunadamente, en \textsf{R} se presenta error al intentar
ejecutar \texttt{beta(1350,950)} o \texttt{beta(850,650)}. 


```r
beta(1350, 950)
```

```
## [1] 0
```

```r
beta(850, 650)
```

```
## [1] 0
```

Por ende, es posible plantear la
siguiente solución numérica cuando se quiere calcular la función
predictiva \@ref(eq:Binompredict) en muestras grandes. El problema
central es el cómputo de $\frac{Beta(a,b)}{Beta(c,d)}$ con $a \geq c$ y
$b \geq d$, valores enteros. Podemos ver que

\begin{align*}
&\ \ \ \ \frac{Beta(a,b)}{Beta(c,d)}\\
&=\frac{(a-1)!(b-1)!(c+d-1)!}{(c-1)!(d-1)!(a+b-1)!}\\
&=\frac{(a-1)(a-2)\cdots(a-(a-c))(b-1)(b-2)\cdots(b-(b-d))}{(a+b-1)(a+b-2)\cdots(a+b-(a+b-c-d))}\\
&=\frac{a^{a-c}(1-\frac{1}{a})(1-\frac{2}{a})\cdots(1-\frac{a-c}{a})b^{b-d}(1-\frac{1}{b})(1-\frac{2}{b})\cdots(1-\frac{b-d}{b})}{(a+b)^{a+b-c-d}(1-\frac{1}{a+b})(1-\frac{2}{a+b})\cdots(1-\frac{a+b-c-d}{a+b})}\\
&=\underbrace{\left(\frac{a}{a+b}\right)^{a-c}}_{t_1}\underbrace{\left(\frac{b}{a+b}\right)^{b-d}}_{t_2}\underbrace{(1-\frac{1}{a})(1-\frac{2}{a})\cdots(1-\frac{a-c}{a})}_{t_3}\\
&\ \ \ \ \ \ \underbrace{(1-\frac{1}{b})(1-\frac{2}{b})\cdots(1-\frac{b-d}{b})}_{t_4}\underbrace{(1-\frac{1}{a+b})(1-\frac{2}{a+b})\cdots(1-\frac{a+b-c-d}{a+b})}_{t_5}
\end{align*}

Calculando separadamente los términos $t_1$, $t_2$, $t_3$, $t_4$ y $t_5$
podemos calcular $\frac{Beta(a,b)}{Beta(c,d)}$ para valores grandes de
$a$, $b$, $c$ y $d$. La siguiente función `prob` calcula la
densidad \@ref(eq:Binompredict) para un valor particular de $\tilde{s}$
usando la anterior técnica.


```r
prob <- function(s.mono, n.mono, s, n, alfa, beta){
  a <- s.mono + s + alfa
  b <- n.mono - s.mono + n - s + beta 
  c <- s + alfa
  d <- n - s + beta
  t1 <- (a/(a + b))^(a - c)
  t2 <- (b/(a + b))^(b - d) 
  t3 <- prod(1 - c(1:(a - c))/a)
  t4 <- prod(1 - c(1:(b - d))/b) 
  t5 <- prod(1 - c(1:(a + b - c - d))/(a + b))
  if (a==c) 
    resul <- t2 * t4/t5 
  if (b==d)
    resul <- t1 * t3/t5
  if (a > c & b > d)
    resul <- choose(n.mono, s.mono) * t1 * t2 * t3 * t4/t5
  return(resul) 
  }
```


Si queremos examinar la distribución predictiva para todos valores de la
variable $\tilde{S}$, podemos usar los siguientes códigos 

```r
n <- 1000
s <- 650 
alfa <- 200
beta <- 300 
n.mono <-800
res <- rep(NA, (1 + n.mono)) 
for(i in 1:length(res)){
  res[i] <- prob(i - 1, n.mono, s, n, alfa, beta)
}
```

Y como resultado, el objeto `res` contiene las 801 probabilidades asociadas a todos los posibles valores de $\tilde{s}$. Los resultados obtenidos con la anterior técnica son equivalentes a lo
obtenido usando la función `lbeta` que computa el logaritmo natural
de la función beta. Así, para calcular la probabilidad en
\@ref(eq:Ejebinom), simplemente usamos el siguiente código


```r
choose(800, 500) * exp(lbeta(1350, 950) - lbeta(850, 650)) 
```

```
## [1] 0.0005969157
```

Nótese que esta probabilidad es la misma contenida en la posición 501 del objeto `res` igual a 5.969157\times 10^{-4}. Finalmente, se observa que la distribución predictiva \@ref(eq:Binompredict) corresponde a una distribución
Beta-binomial con parámetros $s+\alpha$ y $\beta-s+n$. El paquete
`VGAM` [@VGAM] en `R`  contiene funciones que calculan la
función de densidad, función de distribución, percentiles, además de
generar números aleatorios para la distribución Beta-binomial. Las
probabilidades puntuales de $\tilde{s}$ se puede calcular con la función
`dbetabinom`, teniendo en cuenta que los parámetros utilizados son
$\mu=(s+\alpha)/(n+\alpha+\beta)$ y $\rho=1/(1+n+\alpha+\beta)$. Con el
siguiente código, podemos calcular las probabilidades para todos los
posibles valores de $\tilde{s}$.


```r
library(VGAM) 
mu <- (s + alfa)/(n + alfa + beta)
rho <- 1/(1 + n + alfa + beta) 
res2 <- rep(NA, (1 + n.mono)) 
for(i in 1:length(res2)){
  res2[i] <- dbetabinom(i - 1,
                        size = n.mono,
                        prob = mu,
                        rho = rho)
}
```

Podemos observar que la posición 501 del objeto `res2` es igual a 5.969157\times 10^{-4}, el cual es idéntico a lo obtenido en `res`. Adicionalmente, al escribir la distribución predictiva de
\@ref(eq:Binompredict) como la función de densidad de una distribución
Beta-binomial, se puede encontrar la esperanza de esta distribución, la
cual está dada por

\begin{equation*}
E(\tilde{S}|S)=\tilde{n}\frac{s+\alpha}{n+\alpha+\beta}
\end{equation*}

Nótese que la esperanza en la anterior expresión corresponde simplemente
al tamaño $\tilde{n}$ de la nueva muestra multiplicado por la
estimación bayesiana del parámetro $\theta$. Adicionalmente, la
esperanza de $\tilde{S}$ también se puede obtener multiplicando todos los
posibles valores de $\tilde{S}$ con su respectiva probabilidad, y sumando
al final, como se muestra a continuación.


```r
sum(res * c(0:n.mono)) 
```

```
## [1] 453.3333
```

```r
n.mono * (s + alfa)/(n + alfa + beta)
```

```
## [1] 453.3333
```

Retomando el ejemplo \@ref(exm:bernoelectoral), suponga que la encuesta de opinión electoral
se lleva a cabo en diferentes ciudades de un determinado país, en este
caso, para cada ciudad se tiene una muestra de variables con
distribución Bernoulli o equivalentemente una variable binomial; de esta
forma, se dispone de una muestra de variables con distribución Binomial.
La distribución posterior del parámetro $\theta$ para estos casos se
encuentra en el siguiente resultado.

\BeginKnitrBlock{proposition}<div class="proposition"><span class="proposition" id="prp:postbinom"><strong>(\#prp:postbinom) </strong></span>Cuando se tiene una sucesión de variables aleatorias $S_1,\ldots,S_i, \ldots,S_k$ independientes y con distribución $Binomial(n_i,\theta)$ para $i=1,\ldots,k$, entonces la distribución posterior del parámetro de interés $\theta$ es

\begin{equation*}
\theta \mid S_1,\ldots,S_k \sim Beta\left(\sum_{i=1}^ks_i+\alpha,\beta+\sum_{i=1}^k n_i-\sum_{i=1}^k s_i\right)
\end{equation*}</div>\EndKnitrBlock{proposition}
<br>

\BeginKnitrBlock{proof}<div class="proof">\iffalse{} <span class="proof"><em>Prueba. </em></span>  \fi{}\begin{align*}
p(\theta \mid S_1,\ldots,S_k)&\propto \prod_{i=1}^kp(S_i \mid \theta)p(\theta \mid \alpha,\beta)\\
&\propto \prod_{i=1}^k\theta^{\sum_{i=1}s_i}\theta^{\alpha-1}(1-\theta)^{\beta-1}
(1-\theta)^{\sum_{i=1}^kn_i-\sum_{i=1}^ks_i}I_{[0,1]}(\theta)\\
&= \theta^{\sum_{i=1}^ks_i+\alpha-1}(1-\theta)^{\sum_{i=1}^kn_i-\sum_{i=1}^ks_i+\beta}I_{[0,1]}(\theta)
\end{align*}

Por lo tanto, factorizando convenientemente, se encuentra una expresión idéntica a la función de densidad de la distribución $Beta\left(\sum_{i=1}^ks_i+\alpha,\beta+\sum_{i=1}^k n_i-\sum_{i=1}^n s_i\right)$.</div>\EndKnitrBlock{proof}
<br>

\BeginKnitrBlock{example}<div class="example"><span class="example" id="exm:unnamed-chunk-23"><strong>(\#exm:unnamed-chunk-23) </strong></span>El siguiente conjunto de datos fue estudiado inicialmente por @Efron75 y se ha convertido en uno de los ejemplos prácticos más citados en la historia de la estadística moderna. Se trata de los porcentajes de bateo en una muestra de 18 jugadores profesionales en la temporada regular de béisbol en Estados Unidos en el año 1970. @wikiBat establece que, en términos generales, este valor representa la razón entre la cantidad de *hits* y el número de turnos al bate^[Un *hit* es la conexión efectuada por el bateador que coloca la pelota dentro del terreno de juego, permitiéndole alcanzar al menos una base, sin que se produzca un error de defensa del equipo contrario. Para lograr un hit, el bateador debe llegar a primera base antes de que ningún jugador defensivo lo toque con la bola en el trayecto del home a la inicial, o que el jugador de la defensa que tenga la bola pise la primera base antes que el bateador llegue a la misma.]. La fórmula para calcular esta estadística es $s/n$, donde $s$ es el número de *hits* y $n$ es el total de turnos. Este conjunto de datos está disponible en el paquete `pscl` de `R` y se puede cargar mediante el siguiente código computacional.</div>\EndKnitrBlock{example}



```r
library(pscl)
data(EfronMorris)
```


|name             |team  |league |  r|     y|   n|     p|
|:----------------|:-----|:------|--:|-----:|---:|-----:|
|Roberto Clemente |Pitts |NL     | 18| 0.400| 367| 0.346|
|Frank Robinson   |Balt  |AL     | 17| 0.378| 426| 0.298|
|Frank Howard     |Wash  |AL     | 16| 0.356| 521| 0.276|
|Jay Johnstone    |Cal   |AL     | 15| 0.333| 275| 0.222|
|Ken Berry        |Chi   |AL     | 14| 0.311| 418| 0.273|
|Jim Spencer      |Cal   |AL     | 14| 0.311| 466| 0.270|
|Don Kessinger    |Chi   |NL     | 13| 0.289| 586| 0.263|
|Luis Alvarado    |Bos   |AL     | 12| 0.267| 138| 0.210|
|Ron Santo        |Chi   |NL     | 11| 0.244| 510| 0.269|
|Ron Swoboda      |NY    |NL     | 11| 0.244| 200| 0.230|
|Del Unser        |Wash  |AL     | 10| 0.222| 277| 0.264|
|Billy Williams   |Chi   |AL     | 10| 0.222| 270| 0.256|
|George Scott     |Bos   |AL     | 10| 0.222| 435| 0.303|
|Rico Petrocelli  |Bos   |AL     | 10| 0.222| 538| 0.264|
|Ellie Rodriguez  |KC    |AL     | 10| 0.222| 186| 0.226|
|Bert Campaneris  |Oak   |AL     |  9| 0.200| 558| 0.285|
|Thurman Munson   |NY    |AL     |  8| 0.178| 408| 0.316|
|Max Alvis        |Mil   |NL     |  7| 0.156|  70| 0.200|

En la primera columna se tiene el número del jugador, la segunda columna proporciona el nombre del jugador, la cuarta columna representan el número de *hits* en los primeros 45 turnos al bate. La sexta columna representa el número de turnos al bate al final de la temporada regular y la última columna representa el promedio de bateo en la temporada.

Suponga que, partiendo de la muestra de los 18 jugadores, el objetivo es estimar el porcentaje de bateo, notado como $\theta$, en toda la liga en el año de 1970. En primera instancia es plausible considerar que cada uno de los jugadores se comporta de manera independiente y que el porcentaje de bateo es común a todos, puesto que pertenecen a la misma liga profesional. Por lo tanto, es posible establecer que el número de *hits* $s_i$ ($i=1,\ldots,18$) para cada jugador tiene la siguiente distribución

\begin{equation*}
S_i\sim Binomial (n_i,\theta) \ \ \ \ \ \ \ \ \ i=1,\ldots,18.
\end{equation*}

Utilizando un enfoque bayesiano, es posible sacar provecho de la información recolectada al principio de la temporada, constituida por la tercera y cuarta columna del archivo de datos. En esta instancia, se tuvieron $18+17+\cdots+8+7=215$ hits para un total de $45\times 18= 810$ turnos al bate. Con esta información, se define la caracterización estructural de la distribución previa que, siguiendo las recomendaciones anteriores, está dada por una $Beta(\alpha=215, \beta=810-215)=Beta(\alpha=215, \beta=595)$. Del resultado \@ref(prp:postbinom), y teniendo en cuenta que al final de la temporada se obtuvieron $\sum S_i = 1825$ *hits* para un total de $\sum n_i =6649$ turnos al bate, se tiene que la distribución posterior para este ejemplo es una $Beta(1825+215,6649-1825+595)=Beta(2040, 5419)$. Por lo tanto, utilizando la distribución posterior, se estima que el porcentaje de bateo en la liga profesional en el año de 1970 es de $\frac{2040}{2040+5419}=\frac{2040}{7459}=0.273$. Este valor corresponde a la media de la distribución posterior.

Nótese que los mismos resultados se encuentran cuando se analiza este conjunto de datos en `STAN`, mediante el siguiente código computacional.


```r
Binomial <- 'data {
  int<lower=0> n;
  int<lower=0> m[n];
  int<lower=0> s[n];
}
parameters {
  real<lower=0, upper=1> theta;
}
model {
  for(i in 1:n) {
  s[i] ~ binomial(m[i], theta);
  }
  theta ~ beta(215, 595);
}
'

library(rstan)
options(mc.cores = parallel::detectCores())

s <- round(EfronMorris$n * EfronMorris$p)
sample_data <- list(s = s, 
                    n = nrow(EfronMorris),
                    m = EfronMorris$n)

Binfit <- stan(model_code = Binomial, 
               data = sample_data, verbose = FALSE)
```

La siguiente salida de `STAN` permite conocer la estimación bayesiana posterior y los límites del intervalo de credibilidad al 95%.


```r
print(Binfit, pars = "theta", 
      digits = 4, probs = c(0.025, 0.975))
```

```
## Inference for Stan model: b5a37600d5b0f80332bf311eb740e4c6.
## 4 chains, each with iter=2000; warmup=1000; thin=1; 
## post-warmup draws per chain=1000, total post-warmup draws=4000.
## 
##         mean se_mean     sd   2.5%  97.5% n_eff  Rhat
## theta 0.2736   1e-04 0.0052 0.2633 0.2839  1611 1.001
## 
## Samples were drawn using NUTS(diag_e) at Sun Jun  6 23:42:42 2021.
## For each parameter, n_eff is a crude measure of effective sample size,
## and Rhat is the potential scale reduction factor on split chains (at 
## convergence, Rhat=1).
```

La figura \@ref(fig:posBinomialStan) muestra la distribución posterior para este ejemplo, junto con la estimación puntual, correspondiente a la media. 


```r
bayesplot::mcmc_areas(Binfit, pars = "theta", 
                      prob = 0.95)
```

<div class="figure" style="text-align: center">
<img src="3Uniparametricos_files/figure-html/posBinomialStan-1.svg" alt="Distribución posterior." width="576" />
<p class="caption">(\#fig:posBinomialStan)Distribución posterior.</p>
</div>

Por otro lado, el mismo intervalo de credibilidad del 95% correspondiente se puede hallar mediante el siguiente código computacional de `R`.


```r
qbeta(c(0.025, 0.975), 2040, 5419)
```

```
## [1] 0.2634379 0.2836674
```

La figura \@ref(fig:BinomEj1) muestra el comportamiento de las distribuciones previa y posterior para este ejemplo. Nótese que, con un análisis frecuentista, se hubiese llegado a una estimación cercana de $\frac{1825}{6649}=0.274$. 


<div class="figure" style="text-align: center">
<img src="3Uniparametricos_files/figure-html/BinomEj1-1.svg" alt="Función de densidad previa y función de densidad posterior para el ejemplo de bateo." width="576" />
<p class="caption">(\#fig:BinomEj1)Función de densidad previa y función de densidad posterior para el ejemplo de bateo.</p>
</div>

Es posible analizar este conjunto de datos desde otra perspectiva al suponer que los jugadores no constituyen una muestra aleatoria y cada uno de ellos tiene un promedio de bateo diferente. Sin embargo, este análisis se deja como ejercicio en un capítulo posterior.


\BeginKnitrBlock{example}<div class="example"><span class="example" id="exm:unnamed-chunk-29"><strong>(\#exm:unnamed-chunk-29) </strong></span>Continuando con el conjunto de datos de Efron y Morris, suponga que el entrenador de un equipo de las ligas inferiores está interesado en adquirir los servicios de Max Alvis. Este jugador no tuvo un buen promedio de bateo en la temporada y no tuvo muchos turnos al bate. El entrenador quiere conocer cuál será el número más probable de *hits* que anotará en la siguiente temporada. Teniendo en cuenta que es un jugador que viene de la liga profesional, lo más conveniente es que tenga muchos turnos al bate, digamos 400.

Para resolver este cuestionamiento, es conveniente recurrir a la función predictiva posterior, dada en el resultado \@ref(prp:ResPredBinom). Para este análisis, se define la caracterización estructural de la distribución previa del jugador que está dada por una $Beta(\alpha=7, \beta=38)$. La siguiente función en `R` permite obtener la distribución predictiva para este jugador, que se muestra en la figura \@ref(fig:PredBinom).</div>\EndKnitrBlock{example}


```r
n <- 70
s <- 14
alp <-7
bet <- 38
n.ast <- 400
predictiva <- rep(NA, n.ast + 1)
for(k in 0:n.ast){
  predictiva[k + 1] <-
  choose(n.ast,k) *
    beta(k+s+alp,bet-k-s+n.ast+n)/beta(s+alp,bet-s+n)
}

sum(predictiva)
```

```
## [1] 1
```

```r
which(predictiva==max(predictiva))
```

```
## [1] 71
```

La última línea del código computacional permite concluir que lo más probable es que el jugador realice 71 hits en 400 turnos al bate, cifra que no convence al entrenador para adquirir los servicios del jugador.

<div class="figure" style="text-align: center">
<img src="3Uniparametricos_files/figure-html/PredBinom-1.svg" alt="Función de densidad predictiva posterior para el jugador Max Alvis." width="576" />
<p class="caption">(\#fig:PredBinom)Función de densidad predictiva posterior para el jugador Max Alvis.</p>
</div>

## Modelo Binomial negativo

La distribución binomial negativa describe el número de ensayos necesarios para alcanzar un número determinado y fijo de éxitos $k$ en una secuencia independiente de experimentos tipo Bernoulli. Esta distribución es particularmente útil cuando el parámetro $\theta$ que se quiere estimar es muy pequeño, como la proporción de una población que padece de alguna enfermedad rara. La razón por la que no se utiliza la distribución binomial es que al fijar el número de ensayos $n$, con un aprobabilidad $\theta$ muy pequeña, es muy probable que en en la muestra de tamaño $n$ no se encuetre ningún paciente con la enfermedad; mientras que al utilizar la distribución binomial negativa, de antemano se garantiza que se obtendrá $k$ pacientes con la enfermedad en la muestra.

Suponga que $Y$ es una variable aleatoria cuya distribución es Binomial negativa, y que representa el número de ensayos necesarios $y$ para alcanzar un número determinado y fijo de éxitos $k$ en un experimento. La forma funcional de esta distribución es la siguiente
\begin{equation}
p(Y \mid \theta)=\binom{y-1}{k-1}\theta^k(1-\theta)^{y-k}I_{\{k,k+1,\ldots,\}}(y),
\end{equation}

Así como en la distribución Bernoulli y Binomial, el parámetro $\theta$ está restringido al espacio $\Theta=[0,1]$. Luego, es admisible proponer que $\theta$ siga una distribución Beta. Por tanto, la distribución previa del parámetro $\theta$ está dada por la expresión \@ref(eq:betadistribution). Bajo este marco de referencia se tienen los siguientes resultados

\BeginKnitrBlock{proposition}<div class="proposition"><span class="proposition" id="prp:unnamed-chunk-31"><strong>(\#prp:unnamed-chunk-31) </strong></span>La distribución posterior del parámetro $\theta$ sigue una distribución
\begin{equation*}
\theta \mid Y \sim Beta(\alpha+k,\beta+y-k)
\end{equation*}</div>\EndKnitrBlock{proposition}
<br>

\BeginKnitrBlock{proof}<div class="proof">\iffalse{} <span class="proof"><em>Prueba. </em></span>  \fi{}\begin{align*}
p(\theta \mid Y)&\propto p(Y \mid \theta)p(\theta \mid \alpha,\beta)\\
&\propto \theta^{\alpha+k-1} (1-\theta)^{y+\beta-k-1}I_{[0,1]}(\theta)
\end{align*}
Por lo tanto, factorizando convenientemente, se llega a una expresión idéntica a la función de distribución de una variable aleatoria con distribución $Beta(\alpha+k,\beta+y-k)$.</div>\EndKnitrBlock{proof}
<br>

En algunas situaciones se puede encontrar una muestra de variables con distribución binomial negativa. Por ejemplo, la entrevista de pacientes para encontrar cierta enfermedad puede llevarse a cabo en diferentes puntos de atención médica o en diferentes ciudades del país. Así en cada punto de atención, se tendrá el dato correspondiente a una variable con distribución binomial negativa. El procedimiento inferencial bayesiano para estas situaciones se describe a continuación:

\BeginKnitrBlock{proposition}<div class="proposition"><span class="proposition" id="prp:unnamed-chunk-33"><strong>(\#prp:unnamed-chunk-33) </strong></span>Cuando se tiene una sucesión de variables aleatorias $Y_1,\ldots, Y_n$ independientes y con distribución $BinomialNegativa(k_i,\theta)$ $(i=1,\ldots,n)$, entonces la distribución posterior del parámetro de interés es
\begin{equation}
\theta \mid Y_1,\ldots, Y_n \sim Beta(\alpha+\sum_{i=1}^n k_i,\beta+\sum_{i=1}^n y_i-\sum_{i=1}^n k_i)
\end{equation}</div>\EndKnitrBlock{proposition}
<br>

\BeginKnitrBlock{example}<div class="example"><span class="example" id="exm:unnamed-chunk-34"><strong>(\#exm:unnamed-chunk-34) </strong></span>Una franquicia de investigación farmacéutica ha desarrollado un nuevo tratamiento farmacológico sobre pacientes diabéticos que padezcan, a su vez, de enfermedades cardíacas o cardiopatías (angina de pecho, infarto de miocardio,  insuficiencia mitral, estenosis mitral, entre otras). Para evaluar el nuevo tratamiento, es necesario seleccionar una muestra, mediante el diseño de un experimento clínico, de pacientes que tienen estas características.

Por otro lado, se sabe que la proporción de personas que padecen de diabetes y que además tienen algún tipo de condición cardiaca es muy baja y es necesario obtener una estimación precisa de la proporción de personas con estas condiciones. Con base en lo anteriormente expuesto, se puede pensar en seleccionar una muestra grande de personas y utilizar un acercamiento binomial para estimar esta proporción. Sin embargo, dado que la prevalencia de esta condición es bastante baja, es posible que el número de personas en la muestra que presenten estas enfermedades sea nulo; por consiguiente, la estimación binomial no será, de ninguna forma, precisa.

Por lo tanto, el diseño clínico está supeditado al uso de la distribución Binomial Negativa, en donde se entrevistarán pacientes, de una base de datos de un hospital de la ciudad asociado con la franquicia, hasta conseguir una muestra de cinco pacientes que padezcan de estas condiciones. Después de varios meses de entrevistas, se encontró el quinto paciente en la entrevista número 1106.

Mediante el análisis bayesiano, suponiendo una distribución previa $Beta(0.5, 0.5)$, se llega a que la distribución posterior del parámetros $\theta$ es $Beta(0.5+5, 0.5+1106-5)=Beta(5.5, 1101.5)$. Por lo tanto, la estimación puntual del parámetro de interés, que corresponde a la media de la distribución posterior, es 0.0049, que equivale una proporción de 0.49% de personas con estas enfermedades.

El siguiente código computacional muestra cómo se puede llegar a las mismas conclusiones con `STAN`, haciendo la salvedad de que `STAN` define esta distribución en términos del número de fracasos $m = y - k$ necesarios para obtener $k$ éxitos.</div>\EndKnitrBlock{example}


```r
BinNegativa <- 'data {
  int<lower=0> k;
  int<lower=0> y;
}
transformed data {
  int<lower=0> m;
  m = y - k;
}
parameters {
  real<lower=0> beta;
}
transformed parameters {
  real<lower=0> theta;
  theta = beta/(beta + 1);
}
model {
  m ~ neg_binomial(k, beta);
  theta ~ beta(0.5, 0.5);
}
'

y <- 1106
k <- 5
sample_data <- list(k = k, y = y)

BNfit <- stan(model_code = BinNegativa,
              data = sample_data, verbose = FALSE)
```

La siguiente salida de `STAN` permite conocer la estimación bayesiana posterior y los límites del intervalo de credibilidad al 95%.


```r
print(BNfit, pars = "theta", 
      digits = 4, probs = c(0.025, 0.975))
```

```
## Inference for Stan model: ca45da8b0e94b143378403f155105e09.
## 4 chains, each with iter=2000; warmup=1000; thin=1; 
## post-warmup draws per chain=1000, total post-warmup draws=4000.
## 
##         mean se_mean     sd   2.5%  97.5% n_eff   Rhat
## theta 0.0049   1e-04 0.0021 0.0016 0.0099  1583 1.0009
## 
## Samples were drawn using NUTS(diag_e) at Sun Jun  6 23:43:20 2021.
## For each parameter, n_eff is a crude measure of effective sample size,
## and Rhat is the potential scale reduction factor on split chains (at 
## convergence, Rhat=1).
```

Después de las iteraciones necesarias, la salida del anterior código muestra la estimación puntual dada por 0.00498 y un intervalo de credibilidad al 95\%, dado por $(0.00174, 0.01013)$.

\BeginKnitrBlock{example}<div class="example"><span class="example" id="exm:binnegex"><strong>(\#exm:binnegex) </strong></span>Continuando con la temática del ejemplo anterior, suponga que la franquicia llevó a cabo la misma investigación en las 31 ciudades con mayor densidad poblacional de país. En total, se tuvieron 29620 entrevistas para un total de éxitos de 152, tal como se muestra a continuación.

| Ciudad          | y    | k |
|-----------------|------|---|
| BOGOTA          | 1001 | 4 |
| MEDELLIN        | 978  | 6 |
| CALI            | 999  | 5 |
| BARRANQUILLA    | 860  | 4 |
| CARTAGENA       | 1155 | 4 |
| CUCUTA          | 585  | 6 |
| BUCARAMANGA     | 1030 | 3 |
| IBAGUE          | 960  | 5 |
| SOLEDAD         | 1002 | 6 |
| SANTA MARTA     | 763  | 7 |
| SOACHA          | 1036 | 5 |
| PASTO           | 779  | 5 |
| MONTERIA        | 1158 | 4 |
| VILLAVICENCIO   | 1017 | 5 |
| BELLO           | 888  | 6 |
| MANIZALES       | 977  | 4 |
| VALLEDUPAR      | 1256 | 6 |
| BUENAVENTURA    | 1349 | 6 |
| NEIVA           | 1047 | 5 |
| PALMIRA         | 1088 | 5 |
| ARMENIA         | 649  | 3 |
| POPAYAN         | 765  | 4 |
| FLORIDABLANCA   | 699  | 5 |
| SINCELEJO       | 1042 | 4 |
| ITAGUI          | 1212 | 5 |
| BARRANCABERMEJA | 660  | 5 |
| TULUA           | 671  | 5 |
| ENVIGADO        | 835  | 6 |
| DOSQUEBRADAS    | 997  | 5 |
| RIOHACHA        | 1146 | 4 |
| SINCELEJO       | 1016 | 5 |
  
Mediante el análisis bayesiano, suponiendo una distribución previa^[Nótese que es posible también asignar una previa informativa $Beta(5.5, 1101.5)$, que da cuenta de la información del estudio del ejemplo anterior.] no informativa $Beta(0.5, 0.5)$, se llega a que la distribución posterior del parámetros $\theta$ es $Beta(0.5+152, 0.5+29620-152)=Beta(152.5, 29468.5)$. Por lo tanto, la estimación puntual del parámetro de interés, que corresponde a la media de la distribución posterior, es 0.0051, que equivale a una proporción de 0.51% de personas con estas enfermedades. El siguiente código computacional muestra cómo se puede llegar a las mismas conclusiones con `STAN`</div>\EndKnitrBlock{example}


```r
BinNegativa2 <- 'data {
  int<lower=0> n;
  int<lower=0> k[n];
  int<lower=0> y[n];
}
transformed data {
  int<lower=0> m[n];
  for(i in 1:n){
    m[i] = y[i] - k[i];
  }
}
parameters {
  real<lower=0> b;
}
transformed parameters {
  real<lower=0> theta;
  theta = b/(b + 1);
}
model {
  for(i in 1:n){
    m[i] ~ neg_binomial(k[i], b);
  }
  theta ~ beta(0.5, 0.5);
}
'

y <- c(1001, 978, 999, 860, 1155, 585, 1030, 
       960, 1002, 763, 1036, 779, 1158, 1017, 
       888, 977, 1256, 1349, 1047, 1088, 649, 
       765, 699, 1042, 1212, 660, 671, 835, 
       997, 1146, 1016)
k <- c(4, 6, 5, 4, 4, 6, 3, 5, 6, 7, 5, 5, 4, 
       5, 6, 4, 6, 6, 5, 5, 3, 4, 5, 4, 5, 5, 
       5, 6, 5, 4, 5)
sample_data <- list(k = k, y = y, n = length(y))

BNfit2 <- stan(model_code = BinNegativa2,
               data = sample_data, verbose = FALSE)
```

Después de cinco mil iteraciones, la salida del anterior código muestra la estimación puntual dada por 0.00515 y un intervalo de credibilidad al 95%, dado por $(0.00439, 0.00603)$, mucho más estrecho que el intervalo de credibilidad del anterior ejemplo


```r
print(BNfit2, pars = "theta", 
      digits = 4, probs = c(0.025, 0.975))
```

```
## Inference for Stan model: 756b8df22716f40c50aa8045790299a4.
## 4 chains, each with iter=2000; warmup=1000; thin=1; 
## post-warmup draws per chain=1000, total post-warmup draws=4000.
## 
##         mean se_mean    sd   2.5% 97.5% n_eff   Rhat
## theta 0.0052       0 4e-04 0.0044 0.006  1302 1.0032
## 
## Samples were drawn using NUTS(diag_e) at Sun Jun  6 23:43:53 2021.
## For each parameter, n_eff is a crude measure of effective sample size,
## and Rhat is the potential scale reduction factor on split chains (at 
## convergence, Rhat=1).
```

La figura \@ref(fig:posBNStan) muestra la distribución posterior para este ejemplo, junto con la estimación puntual, correspondiente a la media. 


```r
bayesplot::mcmc_areas(BNfit2, pars = "theta", 
                      prob = 0.95)
```

<div class="figure" style="text-align: center">
<img src="3Uniparametricos_files/figure-html/posBNStan-1.svg" alt="Distribución posterior." width="576" />
<p class="caption">(\#fig:posBNStan)Distribución posterior.</p>
</div>

Una vez observados los datos actuales y encontrada la distribución posterior, se puede encontrar la distribución predictiva posterior de una nueva variable con distribución binomial negativa. Es decir, se puede definir el mecanismo probabilístico para el número de ensayos necesarios para encontrar $\tilde{k}$ éxitos.

\BeginKnitrBlock{proposition}<div class="proposition"><span class="proposition" id="prp:unnamed-chunk-39"><strong>(\#prp:unnamed-chunk-39) </strong></span>Después de la recolección de datos, la distribución predictiva posterior para una nueva variable $\tilde{Y}$ está dada por

\begin{equation*}
p(\tilde{Y}|Y_1,\cdots,Y_n)=\binom{\tilde{y}-1}{\tilde{k}-1}\frac{Beta(\alpha+\tilde{k}+\sum k_i,\beta+\tilde{y}-\tilde{k}+\sum y_i-\sum k_i)}{Beta(\alpha+\sum k_i,\beta+\sum y_i-\sum k_i)}I_{\{\tilde{k},\tilde{k}+1,\cdots\}}(\tilde{y})
\end{equation*}</div>\EndKnitrBlock{proposition}
<br>

\BeginKnitrBlock{proof}<div class="proof">\iffalse{} <span class="proof"><em>Prueba. </em></span>  \fi{}\begin{align*}
&\ \ \ p(\tilde{Y}|Y_1,\cdots,Y_n)\\
&=\int p(\tilde{Y}|\theta)p(\theta|Y_1,\cdots,Y_n)d\theta\\
&=\int_{0}^1\binom{\tilde{y}-1}{\tilde{k}-1}\theta^{\alpha+\tilde{k}}(1-\theta)^{\beta+\tilde{y}-\tilde{k}}I_{\{\tilde{k},\tilde{k}+1,\cdots\}}(\tilde{y})\frac{\theta^{\sum k_i-1}(1-\theta)^{\sum y_i-\sum k_i-1}}{Beta(\alpha+\sum k_i,\beta+\sum y_i-\sum k_i)}d\theta\\
&=\binom{\tilde{y}-1}{\tilde{k}-1}\frac{I_{\{\tilde{k},\tilde{k}+1,\cdots\}}(\tilde{y})}{Beta(\alpha+\sum k_i,\beta+\sum y_i-\sum k_i)}\int_{0}^1\theta^{\alpha+\tilde{k}+\sum k_i-1}(1-\theta)^{\beta+\tilde{y}-\tilde{k}+\sum y_i-\sum k_i-1}d\theta\\
&=\binom{\tilde{y}-1}{\tilde{k}-1}\frac{Beta(\alpha+\tilde{k}+\sum k_i,\beta+\tilde{y}-\tilde{k}+\sum y_i-\sum k_i)}{Beta(\alpha+\sum k_i,\beta+\sum y_i-\sum k_i)}I_{\{\tilde{k},\tilde{k}+1,\cdots\}}(\tilde{y})
\end{align*}</div>\EndKnitrBlock{proof}
<br>

\BeginKnitrBlock{example}<div class="example"><span class="example" id="exm:unnamed-chunk-41"><strong>(\#exm:unnamed-chunk-41) </strong></span>Siguiendo con los datos del ejemplo \@ref(exm:binnegex), suponga que se quiere recolectar información de tres pacientes con cardiopatía en cierta ciudad, y se quiere conocer acerca del número de entrevistas necesarias para . Utilizando la distribución previa $Beta(0.5,0.5)$ y los datos de las 31 ciudades del ejemplo, se tiene que la distribución predictiva para el número de entrevistas necesarias para encontrar 3 pacientes está dada por

\begin{align*}
&\ \ \ \ p(\tilde{Y}|Y_1,\cdots,Y_n)\\
&=\binom{\tilde{y}-1}{4}\frac{Beta(0.5+5+152,0.5+\tilde{y}-5+29620-152)}{Beta(0.5+152,0.5+29620-152)}I_{\{5,6,\cdots\}}(\tilde{y})\\
&=\binom{\tilde{y}-1}{4}\frac{Beta(157.5,\tilde{y}+29463.5)}{Beta(152.5,29468.5)}I_{\{5,6,\cdots\}}(\tilde{y})
\end{align*}

Con los siguientes códigos se puede calcular la anterior función predictiva.</div>\EndKnitrBlock{example}


```r
BNpred <- function(y, alfa, beta, s, n, k){
  choose(y - 1, k - 1) *
    exp(lbeta(alfa + k + s, beta + y - k + n - s) -
          lbeta(alfa+s,beta+n-s))
}

alfa <- beta <- 0.5
s <- sum(k)
n <- sum(y)
k <- 5

fun <- rep(NA)
for(y in 5:5000){
  fun[y - 4] <- BNpred(y, alfa, beta, s, n, k)
}

sum(fun)
```

```
## [1] 0.9999994
```

<div class="figure" style="text-align: center">
<img src="3Uniparametricos_files/figure-html/jefber-1.svg" alt="Distribución predictiva posterior para el número de entrevistas necesarias para encontrar 5 pacientes." width="576" />
<p class="caption">(\#fig:jefber)Distribución predictiva posterior para el número de entrevistas necesarias para encontrar 5 pacientes.</p>
</div>

Se puede ver que el número de entrevistas que tiene mayor probabilidad asociadas es el valor 768, usando el comando `which(fun == max(fun))`. También, se puede observar que la probabilidad de que en menos de 500 entrevistas se encuentren los 5 pacientes es de solo el 0.1200985 usando el comando `sum(fun[1:(500 - 4)])`

## Modelo Poisson

Suponga que $\mathbf{Y}=\{Y_1,\ldots,Y_n\}$ es una muestra aleatoria de variables con distribución Poisson con parámetro $\theta$, la función de distribución conjunta o la función de verosimilitud está dada por

\begin{align*}
p(\mathbf{Y} \mid \theta)&=\prod_{i=1}^n\frac{e^{-\theta}\theta^{y_i}}{y_i!}I_{\{0,1,\ldots\}}(y_i)\\
&=\frac{e^{-n\theta}\theta^{\sum_{i=1}^ny_i}}{\prod_{i=1}^ny_i!}I_{\{0,1,\ldots\}^n}(y_1,\ldots,y_n)
\end{align*}

donde $\{0,1\ldots\}^n$ denota el producto cartesiano $n$ veces sobre el conjunto $\{0,1\ldots\}$. Por otro lado, como el parámetro $\theta$ está restringido al espacio $\Theta=(0,\infty)$, entonces es posible formular varias opciones para la distribución previa del parámetro. Algunas de estas se encuentran considerando la distribución exponencial, la distribución Ji-cuadrado o la distribución Gamma. Nótese que las dos primeras son casos particulares de la última. Por lo tanto, la distribución previa del parámetro $\theta$ está dada por

\begin{equation}
(\#eq:PreviaGamma)
p(\theta \mid \alpha,\beta)=\frac{\beta^\alpha}{\Gamma(\alpha)}\theta^{\alpha-1} e^{-\beta\theta}I_{(0,\infty)}(\theta).
\end{equation}

Bajo este marco de referencia se tienen el siguiente resultado con respecto a la distribución posterior del parámetro de interés $\theta$.

\BeginKnitrBlock{proposition}<div class="proposition"><span class="proposition" id="prp:ResPoissonPost"><strong>(\#prp:ResPoissonPost) </strong></span>La distribución posterior del parámetro $\theta$ está dada por

\begin{equation*}
\theta \mid \mathbf{Y} \sim Gamma\left(\sum_{i=1}^ny_i+\alpha,n+\beta\right)
\end{equation*}</div>\EndKnitrBlock{proposition}
<br>

\BeginKnitrBlock{proof}<div class="proof">\iffalse{} <span class="proof"><em>Prueba. </em></span>  \fi{}\begin{align*}
p(\theta \mid \mathbf{Y})&\propto p(\mathbf{Y} \mid \theta)p(\theta \mid \alpha,\beta)\\
&=\frac{I_{\{0,1,\ldots\}^n}(y_1,\ldots,y_n)}{\prod_{i=1}^ny_i!}\frac{\beta^\alpha}{\Gamma(\alpha)}
\theta^{\alpha-1}\theta^{\sum_{i=1}^ny_i}e^{-\beta\theta}e^{-n\theta}I_{(0,\infty)}(\theta)\\
&\propto \theta^{\sum_{i=1}^ny_i+\alpha-1}e^{-(\beta+n)\theta}I_{(0,\infty)}(\theta)
\end{align*}

Por lo tanto, factorizando convenientemente, se encuentra una expresión idéntica a la función de distribución de una variable aleatoria con distribución $Gamma(\sum_{i=1}^ny_i+\alpha,n+\beta)$.</div>\EndKnitrBlock{proof}
<br>

Utilizando el resultado anterior, se tiene que la estimación Bayesiana del parámetro $\theta$ está dada por
\begin{equation*}
\hat{\theta}=\frac{\sum_{i=1}^ny_i+\alpha}{n+\beta}.
\end{equation*}

La anterior expresión sugiere tomar los parámetros de la distribución previa $\alpha$ y $\beta$ de la siguiente manera: $\beta$ representa el número de observaciones en la información previa, mientras que $\alpha$ representa la suma de los datos de la información previa. De esta forma, $\alpha/\beta$ representa la estimación previa del parámetro $\theta$. Y la estimación Bayesiana de $\theta$ se puede escribir como

\begin{align*}
\hat{\theta}&=\frac{\sum_{i=1}^ny_i+\alpha}{\beta+n}\\
&=\frac{n}{n+\beta}*\frac{\sum y_i}{n}+\frac{\beta}{n+\beta}*\frac{\alpha}{\beta}\\
&=\frac{n}{n+\beta}*\hat{\theta}_C+\frac{\beta}{n+\beta}*\hat{\theta}_P
\end{align*}

Es decir, la estimación Bayesiana de $\theta$ es un promedio ponderado entre la estimación clásica y la estimación previa del parámetro $\theta$, donde los pesos dependen directamente del tamaño muestral de la información actual y de la información previa. 

A continuación estudiamos las distribuciones predictivas previa y posterior para una nueva observación 

\BeginKnitrBlock{proposition}<div class="proposition"><span class="proposition" id="prp:unnamed-chunk-44"><strong>(\#prp:unnamed-chunk-44) </strong></span>La distribución predictiva previa para una observación $\mathbf{y}=\{y_1,\ldots,y_n\}$ de la muestra aleatoria está dada por

\begin{equation}
(\#eq:PrePriorPoisson)
p(\mathbf{Y})=\frac{\Gamma(\sum_{i=1}^ny_i+\alpha)}{\Gamma(\alpha)}\frac{\beta^\alpha}{(n+\beta)^{\sum_{i=1}^ny_i+\alpha}}
\frac{I_{\{0,1,\ldots\}^n}(y_1,\ldots,y_n)}{\prod_{i=1}^ny_i!}
\end{equation}
  
y define una auténtica función de densidad de probabilidad continua.</div>\EndKnitrBlock{proposition}
<br>

\BeginKnitrBlock{proof}<div class="proof">\iffalse{} <span class="proof"><em>Prueba. </em></span>  \fi{}De la definición de función de distribución predictiva se tiene que 

\begin{align*}
p(\mathbf{Y})&=\int p(\mathbf{Y} \mid \theta)p(\theta \mid \alpha,\beta)\ d\theta\\
&=\int_0^{\infty} \frac{e^{-n\theta}\theta^{\sum_{i=1}^ny_i}}{\prod_{i=1}^ny_i!}I_{\{0,1,\ldots\}^n}(y_1,\ldots,y_n)
\frac{\beta^\alpha \theta^{\alpha-1} e^{-\beta\theta}}{\Gamma(\alpha)}\ d\theta\\
&=\frac{\Gamma(\sum_{i=1}^ny_i+\alpha)}{\Gamma(\alpha)}\frac{\beta^\alpha}{(n+\beta)^{\sum_{i=1}^ny_i+\alpha}}
\frac{I_{\{0,1,\ldots\}^n}(y_1,\ldots,y_n)}{\prod_{i=1}^ny_i!}\\
&\hspace{2cm}\times
\int_0^{\infty} \frac{(n+\beta)^{\sum_{i=1}^ny_i+\alpha}}{\Gamma(\sum_{i=1}^ny_i+\alpha)}
\theta^{\sum_{i=1}^ny_i+\alpha-1}e^{-(\beta+n)\theta} \ d\theta\\
&=\frac{\Gamma(\sum_{i=1}^ny_i+\alpha)}{\Gamma(\alpha)}\frac{\beta^\alpha}{(n+\beta)^{\sum_{i=1}^ny_i+\alpha}}
\frac{I_{\{0,1,\ldots\}^n}(y_1,\ldots,y_n)}{\prod_{i=1}^ny_i!}
\end{align*}</div>\EndKnitrBlock{proof}
<br>

En el caso en el que la muestra aleatoria estuviera constituida por una sola variable aleatoria, entonces $n=1$ y si, en particular, los hiper-parámetros de la distribución previa fuesen $\alpha=\beta=1$, entonces no es difícil ver, utilizando la definición de la función matemática Gamma, que la función de distribución predictiva \@ref(eq:PrePriorPoisson) estaría dada por

\begin{align}
(\#eq:PrePriorPoisson1)
p(Y)&=\frac{\Gamma(y+1)}{\Gamma(1)}\frac{1}{2^{y+1}}\frac{I_{\{0,1,\ldots\}}(y)}{y!} \notag\\
&=\frac{1}{2^{y+1}}I_{\{0,1,\ldots\}}(y)
\end{align}

Para chequear la convergencia de la anterior distribución es necesario recurrir a los resultados del análisis matemático @Apostol[pág. 361]. Dado que el espacio de muestreo de la variable aleatoria $Y$ es $\{0,1,\ldots\}$, entonces la suma infinita converge a uno, lo que conlleva a que en este caso particular $P(Y)$ sea una auténtica función de densidad de probabilidad.

\begin{align*}
\sum_{y=0}^{\infty}p(Y=y)=\sum_{y=0}^{\infty}\left(\frac{1}{2}\right)^{y+1}=\frac{1}{2}\sum_{y=0}^{\infty}\left(\frac{1}{2}\right)^{y}
=\frac{1}{2}\frac{1}{1-1/2}=1
\end{align*}

Podemos afirmar que la expresion \@ref(eq:PrePriorPoisson1) sí representa una función de densidad de una variable discreta. Ahora, consideramos la distribución predictiva poseterior de una muestra aleatoria, esta distribución se presenta en el siguiente resultado.

\BeginKnitrBlock{proposition}<div class="proposition"><span class="proposition" id="prp:ResPoissonPred"><strong>(\#prp:ResPoissonPred) </strong></span>Después de la recolección de los datos, la distribución predictiva posterior para una nueva posible observación $\tilde{\mathbf{y}}=\{\tilde{y}_1,\ldots,\tilde{y}_{n^*}\}$, de tamaño $n^*$, está dada por

\begin{align}
p(\tilde{\mathbf{y}} \mid \mathbf{Y})&=\frac{\Gamma(\sum_{i=1}^{n^*}\tilde{y}_i+\sum_{i=1}^ny_i+\alpha)}{\Gamma(\sum_{i=1}^ny_i+\alpha)}
\frac{(\beta+n)^{\sum_{i=1}^ny_i+\alpha}}{({n^*}+\beta+n)^{\sum_{i=1}^{n^*}\tilde{y}_i+\sum_{i=1}^ny_i+\alpha}}\notag\\
&\hspace{5cm}\times
\frac{I_{\{0,1,\ldots\}^{n^*}}(\tilde{y}_1,\ldots,\tilde{y}_{n^*})}{\prod_{i=1}^{n^*}\tilde{y}_i!}
\end{align}</div>\EndKnitrBlock{proposition}
<br>

La anterior distribución corresponde a una distribucion multivariada que nos permite calcular probabilidades predictivas para cualesquiera valores de $\tilde{y}_1$, $\cdots$, $\tilde{y}_{n^*}$; sin embargo, en algunas situaciones, como por ejemplo cuando $\theta$ representa el número promedio de algún suceso en una región geográfica, al momento de la predicción, podemos estar interesados en predecir el número total o el número promedio de sucesos en la nueva muestra aleatoria de regiones geográficas. Es decir, podemos estar más interesados en la distribución de $\sum_{y=1}^{n^*} \tilde{y}_i$ o de $\sum_{y=1}^{n^*} \tilde{y}_i/n^*$ en vez de la distribución conjunta de $\tilde{y}_1$, $\cdots$, $\tilde{y}_{n^*}$. La distribución predictiva de $\sum_{y=1}^{n^*} \tilde{y}_i$ se presenta en el siguiente resultado, y con esta se pueden obtener fácilmente probabilidades predictivas para $\sum_{y=1}^{n^*} \tilde{y}_i/n^*$.

\BeginKnitrBlock{proposition}<div class="proposition"><span class="proposition" id="prp:unnamed-chunk-46"><strong>(\#prp:unnamed-chunk-46) </strong></span>Después de la recolección de los datos, la distribución predictiva posterior para la suma de un vector de observaciones nuevas $\left(\tilde{y}_1,\ldots,\tilde{y}_{n^*}\right)$, $\tilde{s} = \sum_{y=1}^{n^*} \tilde{y}_i$, está dada por:

\begin{align}
(\#eq:PrePosPoissonSum)
p(\tilde{\mathbf{s}} \mid \mathbf{Y})&=\frac{\Gamma(\tilde{s}+\sum_{i=1}^ny_i+\alpha)}{\Gamma(\sum_{i=1}^ny_i+\alpha)}
\frac{(n+\beta)^{\sum_{i=1}^ny_i+\alpha}}{({n^*}+n+\beta)^{\tilde{s}+\sum_{i=1}^ny_i+\alpha}}\frac{(n^*)^{\tilde{s}}I_{\{0,1,\ldots\}}(\tilde{s})}{\tilde{s}!}
\end{align}</div>\EndKnitrBlock{proposition}
<br>

\BeginKnitrBlock{proof}<div class="proof">\iffalse{} <span class="proof"><em>Prueba. </em></span>  \fi{}Usando el hecho de que $\theta|\mathbf{Y}\sim Gamma(\sum_{i=1}^{n}y_i+\alpha,n+\beta)$ y $\tilde{s}|\theta\sim Poisson(n^*\theta)$ se procede a calcular $\tilde{s}/p(\mathbf{y})$,
así:

\begin{align*}
&\ \ \ \ p(\tilde{s}|\mathbf{y}) \\
&= \int_{\Omega} p(\tilde{s}|\theta)p(\theta|\mathbf{y})d\theta\\
& = \int_{\Omega} \frac{(n^{*}\theta)^{\tilde{s}}e^{-n^*\theta}}{\tilde{s}!} I_{\{0,1,\ldots\}}(\tilde{s}) (\beta+n)^{\sum_{i=1}^{n}y_i+\alpha}\frac{\theta^{\tilde{s}+\sum_{i=1}^{n}y_i+\alpha-1}}{\Gamma(\sum_{i=1}^{n}y_i+\alpha)}e^{-(\beta+n)\theta}I_{(0,\infty)}(\theta) d\theta\\
&= \frac{(n^*)^{\tilde{s}}(\beta+n)^{\sum_{i=1}^{n}y_i+\alpha}}{\tilde{s}!\Gamma(\sum_{i=1}^{n}y_i+\alpha)}I_{\{0,1,\ldots\}}(\tilde{s})\int_{0}^{\infty}\theta^{\sum_{i=1}^{n}y_i+\alpha-1}e^{-(n^*+\beta+n)\theta}d\theta
\end{align*}

Agrupando las constantes para obtener la integral de una distribución gamma con $\alpha=\tilde{s}+\sum_{i=1}^{n}y_i+\alpha$ y $\beta=n^*+n+\beta$ se obtiene el resultado.</div>\EndKnitrBlock{proof}
<br>

En la práctica, evaluar directamente la expresión \@ref(eq:PrePosPoissonSum) puede ocasionar problemas numéricas, por la presencia de la función Gamma y las potencias. Para evitar dicha dificultad, podemos usar la siguiente expresión equivalente cuando $\tilde{s}=1,2,\cdots$:

\begin{equation*}
p(\tilde{\mathbf{s}} \mid \mathbf{Y})=\frac{\Gamma(\tilde{s})}{Beta(\tilde{s},\sum_{i=1}^ny_i+\alpha)}
\left(\frac{n+\beta}{n^*+n+\beta}\right)^{\sum_{i=1}^ny_i+\alpha}\frac{(n^*)^{\tilde{s}}}{(n^*+n+\beta)^{\tilde{s}}}\tilde{s}!
\end{equation*}

Cuando $\tilde{s}=0$, la distribución predictiva es simplemente:

\begin{equation*}
p(\tilde{\mathbf{s}} \mid \mathbf{Y})=
\left(\frac{n+\beta}{n^*+n+\beta}\right)^{\sum_{i=1}^ny_i+\alpha}
\end{equation*}

Ahora, debido a la complejidad de la expresión en \@ref(eq:PrePosPoissonSum), es prácticamente imposible comprobar analíticamente $\sum_{i=0}^\infty p(\tilde{s}=i)=1$, y también muy difícil encontrar una expresión matemática cerrada de la esperanza de la variable $\mathbf{\tilde{s}}$. Sin embargo, en situaciones prácticas, se puede usar aproximaciones numéricas tal como se verá en el ejemplo al final de esta sección.

Anteriormente en el ejemplo \@ref(exm:EjemPoisson) se consideró la situación en la cual no se tenía ninguna consideración para formular una distribución previa y se consluyó que la distribución previa no informativa de Jeffreys para este caso es
\begin{equation*}
p(\theta)\propto\theta^{-1/2},
\end{equation*}

Esta es una distribución previa impropia, puesto que $\int_{0}^\infty \theta^{-1/2}=\infty$. Sin embargo, este hecho no afecta que la inferencia posterior se pueda llevar a cabo, puesto que la distribución posterior está dada por

\begin{equation*}
\theta|\mathbf{Y}\sim Gamma(\sum y_i+1/2,n)
\end{equation*}

Por consiguiente, la estimación Bayesiana del parámetro $\theta$ viene dada por 

\begin{equation*}
\hat{\theta}=\frac{\sum y_i+1/2}{n}.
\end{equation*}

la cual es muy similar a la estimación clásica de $\theta$ dada por $\bar{Y}$. Cuando se utiliza la distribución previa no informativa de Jeffreys, la distribución predictiva para nuevas observaciones $\tilde{y}={\tilde{y}_1,\cdots,\tilde{y}_{n^*}}$ y $\tilde{s}=\sum_{i=1}^{n^*}\tilde{y_i}$ están dadas por

\begin{equation}
(\#eq:PredPoissonJeffreys)
p(\tilde{\mathbf{y}} \mid \mathbf{Y})=\frac{\Gamma(\sum_{i=1}^{n^*}\tilde{y}_i+\sum_{i=1}^ny_i+0.5)}{\Gamma(\sum_{i=1}^ny_i+0.5)}
\frac{n^{\sum_{i=1}^ny_i+0.5}}{({n^*}+n)^{\sum_{i=1}^n\tilde{y}_i+\sum_{i=1}^ny_i+0.5}}
\frac{I_{\{0,1,\ldots\}^{n^*}}(\tilde{y}_1,\ldots,\tilde{y}_{n^*})}{\prod_{i=1}^{n^*}\tilde{y}_i!}
\end{equation}

y

\begin{equation}
(\#eq:Pred1PoissonJeffreys)
p(\tilde{\mathbf{s}} \mid \mathbf{Y})=\frac{\Gamma(\tilde{s}+\sum_{i=1}^ny_i+0.5)}{\Gamma(\sum_{i=1}^ny_i+0.5)}
\frac{n^{\sum_{i=1}^ny_i+0.5}}{({n^*}+n)^{\tilde{s}+\sum_{i=1}^ny_i+0.5}}\frac{I_{\{0,1,\ldots\}}(\tilde{s})}{\tilde{s}!}
\end{equation}

\BeginKnitrBlock{example}<div class="example"><span class="example" id="exm:unnamed-chunk-48"><strong>(\#exm:unnamed-chunk-48) </strong></span>Por políticas gubernamentales, las autoridades municipales están obligados a realizar un seguimiento exhaustivo al comportamiento de la accidentalidad en las vías urbanas y medirlo en términos del número de accidentes de tránsito. Lo anterior es necesario para evaluar la gestión de la administración pública y evaluar las políticas que el gobierno de la ciudad ha implementado para disminuir esta cifra.

Suponga que en una ciudad se quiere implementar una estrategia educativa para disminuir el número de accidentes de tránsito generados por manejar en estado de embriaguez. Para esto, se registraron durante diez días 30 días el número de accidentes de tránsito por ebriedad del conductor. Los datos para cada uno de los días son 22, 9, 9, 20, 10, 14, 11, 14, 11, 11, 19, 12, 8, 9, 16, 8, 13, 8, 14, 12, 14, 11, 14, 13, 11, 14, 13, 11, 7, 12.

Es posible modelar la variable aleatoria número de accidentes de tránsito en un día mediante una distribución de Poisson puesto que el promedio muestral y la varianza muestral de los datos son semejantes. Para este conjunto de datos, el promedio equivale a 12.33, mientras que la varianza es de 12.51. El histograma de los valores observados se puede ver en la figura \@ref(fig:EjemPoisson1).</div>\EndKnitrBlock{example}


<div class="figure" style="text-align: center">
<img src="3Uniparametricos_files/figure-html/EjemPoisson1-1.svg" alt="Histograma para los datos de accidentes de tránsito." width="576" />
<p class="caption">(\#fig:EjemPoisson1)Histograma para los datos de accidentes de tránsito.</p>
</div>

En primera instancia, es posible realizar un análisis no informativo, al formular una distribución previa de Jeffreys proporcional a $\theta^{-1/2}$, para lo cual la distribución posterior será $Gamma(\sum_{i=1}^n y_i+1/2, n)=Gamma(370.5, 30)$. Por lo tanto, un estimador de $\theta$ está dado por la media de la distribución posterior que es $370.5/30=12.35$, muy cercano al valor del estimador de máxima verosimilitud correspondiente al promedio muestral. La figura \@ref(fig:EjemPoisson2) (lado izquierdo) muestra el comportamiento de las distribuciones de Jeffreys y posterior para este ejemplo.

Por otro lado, basándose en datos históricos, la alcaldía observó que, en el mismo periodo del año anterior, ocurrieron 37 accidentes en 9 días de observación. Luego, una distribución previa informativa^[En la práctica, se recomienda que los valores de los hiperparámetros $\alpha$ y $\beta$ correspondan a la suma del número de eventos más uno y número de observaciones, respectivamente.] está dada por $Gamma(\alpha=38,\beta=9)$. Luego, apelando al resultado \@ref(prp:ResPoissonPost), la distribución posterior corresponde a una $Gamma(370+38, 30+9)=Gamma(408, 39)$. Para este caso, un estimador de $\theta$ está dado por la media de la distribución posterior que es $480/39=12.31$. La figura \@ref(fig:EjemPoisson2) (lado derecho) muestra el comportamiento de las distribuciones previa (informativa) y posterior para este ejemplo.

<div class="figure" style="text-align: center">
<img src="3Uniparametricos_files/figure-html/EjemPoisson2-1.svg" alt="Distribución previa y distribución posterior para el ejemplo del tránsito con dos distribuciones previas diferentes (el lado izquierdo representa el caso cuando se usa la previa no informativa, el lado derecho la previa informativa)." width="576" />
<p class="caption">(\#fig:EjemPoisson2)Distribución previa y distribución posterior para el ejemplo del tránsito con dos distribuciones previas diferentes (el lado izquierdo representa el caso cuando se usa la previa no informativa, el lado derecho la previa informativa).</p>
</div>

A continuación se examina la distribución predictiva. En la figura \@ref(fig:PredPostPoisson) se grafica la distribución predictiva para una nueva observación cuando se usa la previa no informativa y la previa informativa. Los códigos para el cálculo cuando se usan ambas distribuciones previas es como sigue:


```r
n <- length(Trans)

pre.Transito.NoInf <- function(s){
  if(s > 0){
    val <- gamma(s) * 
      (n/(n + 1)) ^ (sum(Trans) + 0.5)/
      (beta(s, sum(Trans) + 0.5) * 
         prod(1:s) * (n + 1) ^ s)
  }
  if( s == 0){
    val <- (n/(n + 1)) ^ (sum(Trans) + 0.5)
  }
  return(val)
}

pre.Transito <- function(s){
  if(s > 0){
    val <- gamma(s) * 
      ((n + beta)/(n + beta + 1))^(sum(Trans) + alfa)/
      (beta(s, sum(Trans) + alfa) * (1+n+beta) ^ s * 
         prod(1:s))
  }
  if(s == 0){
    val <- ((n + beta)/(n + beta + 1)) ^ (sum(Trans) + alfa)
    }
  return(val)
}

s.max <- 40 
s.val <- 0:s.max

pre.NoInf.val <- pre.Inf.val <- c()
for(i in 1:length(s.val)){
  pre.NoInf.val[i] <- pre.Transito.NoInf(s.val[i])
  pre.Inf.val[i] <- pre.Transito(s.val[i])
}
sum(pre.NoInf.val)
```

```
## [1] 1
```

```r
sum(pre.Inf.val)
```

```
## [1] 1
```

Nótese que en los anteriores códigos se usó como valor máximo 40 para la variable $\mathbf{\tilde{s}}$, a pesar de que esta toma valores infinitos; pero al ver que la suma de las probabilidades desde el valor 0 hasta el 40 es igual a 1, podemos concluir que la probabilidad de que $\mathbf{\tilde{s}} > 40$ es prácticamente nula.

Finalmente, podemos tener una aproximación de la esperanza de la variable $\mathbf{\tilde{s}}$ como


```r
sum(pre.NoInf.val * s.val)
```

```
## [1] 12.35
```
<div class="figure" style="text-align: center">
<img src="3Uniparametricos_files/figure-html/PredPostPoisson-1.svg" alt="Distribución predictiva posterior para $n^*=1$ para el ejemplo del tránsito." width="576" />
<p class="caption">(\#fig:PredPostPoisson)Distribución predictiva posterior para $n^*=1$ para el ejemplo del tránsito.</p>
</div>

A continuación se presenta el código computacional para realizar la inferencia bayesiana en `STAN` utilizando la distribución previa predictiva.


```r
Poisson <- '
data {
  int<lower=0> n;
  int<lower=0> y[n];
}
parameters {
  real<lower=0> theta;
}
model {
  y ~ poisson(theta);
  theta ~ gamma(38, 9);
}
'
sample_data <- list(y = Trans, n = length(Trans))
Poissonfit <- stan(model_code = Poisson,
               data = sample_data, verbose = FALSE)
```

Después de converger, la salida del anterior código muestra la estimación puntual dada por 10.482 y un intervalo de credibilidad al 95%, dado por $(9.470, 11.500)$, mucho más estrecho que el intervalo de credibilidad del anterior ejemplo


```r
print(Poissonfit, pars = "theta", 
      digits = 4, probs = c(0.025, 0.975))
```

```
## Inference for Stan model: 1c8218d83e89323ee40322eb3f4fd32d.
## 4 chains, each with iter=2000; warmup=1000; thin=1; 
## post-warmup draws per chain=1000, total post-warmup draws=4000.
## 
##          mean se_mean     sd  2.5%   97.5% n_eff   Rhat
## theta 10.4414  0.0121 0.4941 9.471 11.4121  1654 1.0005
## 
## Samples were drawn using NUTS(diag_e) at Sun Jun  6 23:44:44 2021.
## For each parameter, n_eff is a crude measure of effective sample size,
## and Rhat is the potential scale reduction factor on split chains (at 
## convergence, Rhat=1).
```

La figura \@ref(fig:posPoissonStan) muestra la distribución posterior para este ejemplo, junto con la estimación puntual, correspondiente a la media. 


```r
bayesplot::mcmc_areas(Poissonfit, pars = "theta", 
                      prob = 0.95)
```

<div class="figure" style="text-align: center">
<img src="3Uniparametricos_files/figure-html/posPoissonStan-1.svg" alt="Distribución posterior." width="576" />
<p class="caption">(\#fig:posPoissonStan)Distribución posterior.</p>
</div>

## Modelo Exponencial

Suponga que $\mathbf{Y}=\{Y_1,\ldots,Y_n\}$ corresponde a una muestra de variables aleatorias con distribución Exponencial. Luego, la función de distribución conjunta o verosimilitud está dada por

\begin{align}
p(\mathbf{Y} \mid \theta)&=\prod_{i=1}^n\theta e^{(-\theta y)}I_{(0,\infty)}(y_i) \notag \\
&=\theta^n e^{(-\theta \sum_{i=1}^ny_i)}I_{(0,\infty)^n}(y_1,\ldots,y_n)
\end{align}

Donde $\{0,1\ldots\}^n$ denota el producto cartesiano $n$ veces sobre el intervalo $(0,\infty)$. Por otro lado, como el parámetro $\theta$ está restringido al espacio $\Theta=(0,\infty)$, entonces es posible formular varias opciones para la distribución previa del parámetro, al igual que en la distribución Poisson. Así mismo, suponga que la distribución previa para el parámetro de interés es la distribución Gamma tal como aparece en la expresión \@ref(eq:PreviaGamma). Bajo este marco de referencia se tienen los siguientes resultados

\BeginKnitrBlock{proposition}<div class="proposition"><span class="proposition" id="prp:unnamed-chunk-53"><strong>(\#prp:unnamed-chunk-53) </strong></span>La distribución posterior del parámetro $\theta$ sigue una distribución
\begin{equation*}
\theta \mid \mathbf{Y} \sim Gamma\left(\alpha+n,\beta+\sum_{i=1}^ny_i\right)
\end{equation*}</div>\EndKnitrBlock{proposition}
<br>

\BeginKnitrBlock{proof}<div class="proof">\iffalse{} <span class="proof"><em>Prueba. </em></span>  \fi{}\begin{align*}
p(\theta \mid \mathbf{Y})&\propto p(\mathbf{Y} \mid \theta)p(\theta \mid \alpha,\beta)\\
&=\theta^n e^{(-\theta \sum_{i=1}^ny_i)}I_{(0,\infty)^n}(y_1,\ldots,y_n)\frac{\beta^\alpha \theta^{\alpha-1} e^{-\beta\theta}}{\Gamma(\alpha)}I_{(0,\infty)}(\theta)\\
&\propto \theta^{\alpha+n-1}e^{-(\beta+\sum_{i=1}^ny_i)}I_{(0,\infty)}(\theta)
\end{align*}

Por lo tanto, factorizando convenientemente, se encuentra una expresión idéntica a la función de distribución de una variable aleatoria con distribución $Gamma(\alpha+n,\beta+\sum_{i=1}^ny_i)$.</div>\EndKnitrBlock{proof}
<br>


\BeginKnitrBlock{proposition}<div class="proposition"><span class="proposition" id="prp:unnamed-chunk-55"><strong>(\#prp:unnamed-chunk-55) </strong></span>La distribución predictiva previa para una observación $\mathbf{y}=\{y_1,\ldots,y_n\}$ de la muestra aleatoria está dada por

\begin{equation}
p(\mathbf{Y})=\frac{\Gamma(\alpha+n)}{\Gamma(\alpha)}\frac{\beta^\alpha}{(\beta+\sum_{i=1}^ny_i)^{\alpha+n}}
I_{(0,\infty)^n}(y_1,\ldots,y_n)
\end{equation}

y define una auténtica función de densidad de probabilidad continua.</div>\EndKnitrBlock{proposition}
<br>

\BeginKnitrBlock{proof}<div class="proof">\iffalse{} <span class="proof"><em>Prueba. </em></span>  \fi{}De la definición de función de distribución predictiva se tiene que

\begin{align*}
p(\mathbf{Y})&=\int p(\mathbf{Y} \mid \theta)p(\theta \mid \alpha,\beta)\ d\theta\\
&=\int_0^{\infty}\theta^n e^{(-\theta \sum_{i=1}^ny_i)}I_{(0,\infty)^n}(y_1,\ldots,y_n)\frac{\beta^\alpha \theta^{\alpha-1} e^{-\beta\theta}}{\Gamma(\alpha)} \ d\theta\\
&=\frac{\Gamma(n+\alpha)}{\Gamma(\alpha)}\frac{\beta^\alpha}{(\beta+\sum_{i=1}^ny_i)^{\alpha+n}}I_{(0,\infty)^n}(y_1,\ldots,y_n)\\
&\hspace{2cm}\times
\int_0^{\infty} \frac{(\beta+\sum_{i=1}^ny_i)^{\alpha+n}}{\Gamma(n+\alpha)} \theta^{\alpha+n-1}e^{-(\beta+\sum_{i=1}^ny_i)\theta}
\ d\theta\\
&=\frac{\Gamma(\alpha+n)}{\Gamma(\alpha)}\frac{\beta^\alpha}{(\beta+\sum_{i=1}^ny_i)^{\alpha+n}}I_{(0,\infty)^n}(y_1,\ldots,y_n)
\end{align*}</div>\EndKnitrBlock{proof}
<br>

Por ejemplo, en el caso en que la muestra aleatoria estuviera constituida por una sola variable aleatoria, entonces no es difícil ver, utilizando la definición de la función matemática Gamma, que la función de distribución predictiva estaría dada por

\begin{align*}
p(Y)&=\frac{\Gamma(\alpha+1)}{\Gamma(\alpha)}\frac{\beta^\alpha}{(\beta+y)^{\alpha+1}}
I_{(0,\infty)}(y)\notag \\
&=\frac{\alpha \beta^\alpha}{(\beta+y)^{\alpha+1}}
I_{(0,\infty)}(y)
\end{align*}

Para chequear la convergencia de la anterior distribución es necesario recurrir a los resultados del cálculo integral. Dado que el espacio de muestreo de la variable aleatoria $Y$ es el intervalo $(0,\infty)$, entonces la integral definida es igual a uno, lo que conlleva a que en este caso particular $P(Y)$ sea una auténtica función de densidad de probabilidad.

\begin{align*}
\int_0^{\infty}p(Y)\ dy=\int_0^{\infty}\frac{\alpha \beta^\alpha}{(\beta+y)^{\alpha+1}} \ dy
=  \beta^\alpha \left[\frac{(\beta+y)^{-\alpha}}{-\alpha} \right]_0^{\infty}
=1
\end{align*}

Volviendo al caso general en donde se tiene una muestra aleatoria, se tiene el siguiente resultado.

\BeginKnitrBlock{proposition}<div class="proposition"><span class="proposition" id="prp:unnamed-chunk-57"><strong>(\#prp:unnamed-chunk-57) </strong></span>Después de la recolección de los datos, la distribución predictiva posterior para una conjunto de nuevas variables aleatorias $\tilde{\mathbf{y}}=\{\tilde{y}_1,\ldots,\tilde{y}_{n^*}\}$, de tamaño $n^*$, está dada por

\begin{align}
p(\tilde{\mathbf{y}} \mid \mathbf{Y})&=
\frac{\Gamma(n+\alpha+n^{*})}{\Gamma(n+\alpha)}
\frac{(\beta+\sum_{i=1}^ny_i)^{n+\alpha}}{(\sum_{i=1}^{n^*}\tilde{y}_i+\beta+\sum_{i=1}^ny_i)^{n^*+\alpha+n}}\notag\\
&\hspace{4cm}\times
I_{(0,\infty)^{n^*}}(\tilde{y}_1,\ldots,\tilde{y}_n)
\end{align}</div>\EndKnitrBlock{proposition}
<br>

El anterior resultado permite calcular la distribución predictiva conjunta de variables aleatorias por observar. En algunas situaciones lo que se quiere pronosticar es el comportamiento probabilístico de promedio muestral de este conjunto de variables aleatorias; es decir, $\bar{Y}^*=\sum_{i=1}^{n^*}\tilde{Y}_i$. En el siguiente resultado se presenta la distribución predictiva de esta variable aleatoria.

\BeginKnitrBlock{proposition}<div class="proposition"><span class="proposition" id="prp:unnamed-chunk-58"><strong>(\#prp:unnamed-chunk-58) </strong></span>Después de la recolección de los datos, la distribución predictiva posterior para el promedio muestral de un nuevo conjunto de variables aleatorias $\bar{Y}^*=\sum_{i=1}^{n^*}\tilde{Y}_i$ está dada por

\begin{equation*}
p(\bar{Y}^*)=\frac{n^*\Gamma(n^*+\alpha+n)}{\Gamma(n^*)\Gamma(\alpha+n)}\frac{(\beta+\sum_{i=1}^ny_i)^{\alpha+n}}{(n^*\bar{Y}^*+\beta+\sum y_i)^{n^*+\alpha+n}}(n^*\bar{Y}^*)^{n^*-1}I_{(0,\infty)}(\bar{Y}^*)
\end{equation*}</div>\EndKnitrBlock{proposition}
<br>

\BeginKnitrBlock{proof}<div class="proof">\iffalse{} <span class="proof"><em>Prueba. </em></span>  \fi{}En primer lugar se halla la distribución predictiva posterior de la variable $\tilde{S}=\sum_{i=1}^{n^*}\tilde{Y}_i$, teniendo en cuenta que $\tilde{S}|\theta\sim Gamma(n^*,\theta)$, de esta forma

\begin{align*}
p(\tilde{S}|\mathbf{Y})&=\int p(\tilde{S}|\theta)p(\theta|\mathbf{Y})\ d\theta\\
&=\int_0^{\infty} \frac{\theta^{n^*}}{\Gamma(n^*)}\tilde{S}^{n^*-1}e^{-\theta\tilde{S}}I_{(0,\infty)}(\tilde{S})\frac{(\beta+\sum_{i=1}^ny_i)^{\alpha+n}}{\Gamma(\alpha+n)}\theta^{{\alpha+n-1}}e^{-(\beta+\sum y_i)\theta}d\theta\\
&=\frac{\tilde{S}^{n^*-1}(\beta+\sum_{i=1}^ny_i)^{\alpha+n}}{\Gamma(n^*)\Gamma(\alpha+n)}I_{(0,\infty)}(\tilde{S})\int_0^{\infty} \theta^{n^*+\alpha+n-1}e^{-(\tilde{S}+\beta+\sum y_i)\theta}\ d\theta\\
&=\frac{\tilde{S}^{n^*-1}(\beta+\sum_{i=1}^ny_i)^{\alpha+n}}{\Gamma(n^*)\Gamma(\alpha+n)}\frac{\Gamma(n^*+\alpha+n)}{(\tilde{S}+\beta+\sum y_i)^{n^*+\alpha+n}}I_{(0,\infty)}(\tilde{S})
\end{align*}

Al aplicar el teorema de transformación a la distribución predictiva, se puede hallar la distribución de $\bar{Y}^*$, dada por

\begin{align*}
p(\bar{Y}^*|\mathbf{Y})=\frac{n^*\Gamma(n^*+\alpha+n)}{\Gamma(n^*)\Gamma(\alpha+n)}\frac{(\beta+\sum_{i=1}^ny_i)^{\alpha+n}}{(n^*\bar{Y}^*+\beta+\sum y_i)^{n^*+\alpha+n}}(n^*\bar{Y}^*)^{n^*-1}I_{(0,\infty)}(\bar{Y}^*)
\end{align*}</div>\EndKnitrBlock{proof}
<br>

En la práctica puede ocurrir que algunos de los valores de $n$, $n^*$, $\sum_{i=1}^ny_i$ y $n^*\bar{Y}^*$ sean muy grandes; por consiguiente, evaluar directamente la expresión anterior puede ocasionar problemas numéricos. Realizando algunas operaciones algebraicas, se encuentra la siguiente expresión equivalente para la distribución predictiva posterior de $\bar{Y}^*$ que evita problemas numéricas:

\begin{equation}
(\#eq:PredExpoInforma2)
p(\bar{Y}^*|\mathbf{Y})=\frac{1}{\bar{Y}^*Beta(n,n^*)}\left(\frac{\beta+\sum_{i=1}^ny_i}{\beta+\sum_{i=1}^ny_i+n^*\bar{Y}^*}\right)^{\alpha+n}\left(\frac{n^*\bar{Y}^*}{\beta+\sum_{i=1}^ny_i+n^*\bar{Y}^*}\right)^{n^*}I_{(0,\infty)}(\bar{Y}^*)
\end{equation}

Por otro lado, se puede ver que al utilizar la distribución previa no informativa de Jeffreys, la distribución predictiva posterior de $\bar{Y}^*$ está dada por

\begin{equation}
(\#eq:PredExpoJeffreys1)
p(\bar{Y}^*|\mathbf{Y})=\frac{n^*\Gamma(n^*+n)}{\Gamma(n^*)\Gamma(n)}\frac{(\sum_{i=1}^ny_i)^n}{(n^*\bar{Y}^*+\sum y_i)^{n^*+n}}(n^*\bar{Y}^*)^{n^*-1}I_{(0,\infty)}(\bar{Y}^*)
\end{equation}

La cual es equivalente a la siguiente expresión que en ocasiones puede ser útil para evitar problemas numéricos

\begin{equation}
(\#eq:PredExpoJeffreys2)
p(\bar{Y}^*|\mathbf{Y})=\frac{1}{\bar{Y}^*Beta(n,n^*)}\left(\frac{\sum_{i=1}^ny_i}{\sum_{i=1}^ny_i+n^*\bar{Y}^*}\right)^n\left(\frac{n^*\bar{Y}^*}{\sum_{i=1}^ny_i+n^*\bar{Y}^*}\right)^{n^*}I_{(0,\infty)}(\bar{Y}^*)
\end{equation}

\BeginKnitrBlock{example}<div class="example"><span class="example" id="exm:unnamed-chunk-60"><strong>(\#exm:unnamed-chunk-60) </strong></span>@survi reportan un conjunto de datos que da cuenta de los tiempos de sobrevivencia de $n=69$ miembros del programa de transplante de corazón de Stanford (los tiempos se reportan en días después del transplante). Los datos pueden ser encontrados en el paquete `survival` @survival de `R`, mediante la implementación del siguiente código computacional.</div>\EndKnitrBlock{example}


```r
library(survival)
library(dplyr)
data(heart)

sobrevida <- heart %>%
  filter(transplant == 1) %>%
  mutate(tiempo = stop - start)
```

A continuación, se muestran los primeros y últimos datos de este estudio. Se recuerda que el total de pacientes atendidos en este estudio fue de $n=69$ y la suma de los tiempos de sobrevida es de $\sum_{i=1}^ny_i=25998.5$.


```r
sobrevida %>% {
  rbind(head(., 5), tail(., 5))
}
```



|   | start| stop| event|        age|      year| surgery|transplant |  id| tiempo|
|:--|-----:|----:|-----:|----------:|---------:|-------:|:----------|---:|------:|
|1  |     1|   16|     1|   6.297057| 0.2655715|       0|1          |   3|     15|
|2  |    36|   39|     1|  -7.737166| 0.4900753|       0|1          |   4|      3|
|3  |    51|  675|     1|   2.869268| 0.7802875|       0|1          |   7|    624|
|4  |    12|   58|     1|  -5.497604| 0.8624230|       0|1          |  10|     46|
|5  |    26|  153|     1|  -0.019165| 0.8733744|       0|1          |  11|    127|
|65 |     2|   16|     1|  -7.718001| 5.9767283|       0|1          |  95|     14|
|66 |    13|  180|     0| -21.349760| 6.0095825|       0|1          |  96|    167|
|67 |    21|  131|     0| -24.383299| 6.1437372|       0|1          |  97|    110|
|68 |    96|  109|     0| -19.370294| 6.2039699|       0|1          |  98|     13|
|69 |    38|   39|     0| -12.939083| 6.3956194|       1|1          | 100|      1|

Estos tiempos de sobrevivencia pueden ser modelados mediante una distribución exponencial. Además de inferir acerca del parámetro de esta distribución, también es posible inferir acerca del tiempo promedio de sobrevivencia de un individuo sometido a este tipo de transplantes. Luego, dadas las implicaciones del estudio, se debe ser muy cuidadoso en la asignación de los parámetros de la distribución previa. Una forma de hacerlo es asignar valores muy pequeños a estos parámetros. Otra forma de hacerlo es utilizando la distribución previa de Jeffreys, que corresponde a una distribución impropia y que conduce a resultados muy cercanos a los del enfoque anterior.

Utilizando parámetros previos muy cercanos a cero, la distribución posterior del parámetro de interés es $Gamma(69, 25998.5)$. Como es bien sabido, una estimación bayesiana para el parámetro $\theta$ está dada por la media de esta distribución posterior, la cual equivale a $69/25998.5=0.0026$. Ahora, como la esperanza de la distribución exponencial es $1/\theta$, entonces el tiempo promedio de sobrevivencia es de $1/0.0026=376.78$ días. Sin embargo, en este tipo de estudios es común que se presenten muchos datos atípicos; por ende el promedio no es una medida de escala válida en este tipo de análisis, puesto que no es una medida robusta y se prefiere la utilización de la mediana. El siguiente código computacional en `STAN` puede ser usado para realizar inferencias sobre el parámetro $\theta$, sobre el tiempo promedio y el tiempo mediano. De la misma forma, es posible obtener intervalos de credibilidad para estos parámetros.


```r
Exponencial <- '
data {
  int<lower=0> n;
  vector[n] y;
}
parameters {
  real<lower=0> theta;
}
transformed parameters {
  real<lower=0> invtheta = 1/theta;
}
model {
  y ~ exponential(theta);
  theta ~ gamma(0.1, 0.1);
}
'
sample_data <- list(y = sobrevida$tiempo, 
                    n = nrow(sobrevida))
Expofit <- stan(model_code = Exponencial,
                data = sample_data, verbose = FALSE)
```

Después de las iteraciones, los resultados de este código muestran una estimación para	 $\theta$ de 0.0026 con un intervalo de credibilidad de (0.00208, 0.00332). Para la media $1/\theta$, se tiene una estimación puntual de 381.4 con un intervalo de credibilidad de (302.2, 485.8). La mediana se estimó en 376.8 días de sobrevivencia.


```r
print(Expofit, digits = 4, 
      pars = c("theta", "invtheta"),
      probs = c(0.025, 0.975))
```

```
## Inference for Stan model: 1e6dccb88d5f711b516d61c512b5ced4.
## 4 chains, each with iter=2000; warmup=1000; thin=1; 
## post-warmup draws per chain=1000, total post-warmup draws=4000.
## 
##              mean se_mean      sd     2.5%    97.5% n_eff   Rhat
## theta      0.0027  0.0000  0.0003   0.0021   0.0033  1634 0.9999
## invtheta 381.9284  1.1482 46.7530 301.9569 483.5668  1658 0.9999
## 
## Samples were drawn using NUTS(diag_e) at Sun Jun  6 23:45:23 2021.
## For each parameter, n_eff is a crude measure of effective sample size,
## and Rhat is the potential scale reduction factor on split chains (at 
## convergence, Rhat=1).
```

Las figuras \@ref(fig:posExponencialStan1) y \@ref(fig:posExponencialStan2) muestra la distribución posterior para este ejemplo, junto con la estimación puntual, correspondiente a la media. 


```r
bayesplot::mcmc_areas(Expofit, pars = "theta", 
                      prob = 0.95)
```

<div class="figure" style="text-align: center">
<img src="3Uniparametricos_files/figure-html/posExponencialStan1-1.svg" alt="Distribución posterior." width="576" />
<p class="caption">(\#fig:posExponencialStan1)Distribución posterior.</p>
</div>


```r
bayesplot::mcmc_areas(Expofit, pars = "invtheta", 
                      prob = 0.95)
```

<div class="figure" style="text-align: center">
<img src="3Uniparametricos_files/figure-html/posExponencialStan2-1.svg" alt="Distribución posterior." width="576" />
<p class="caption">(\#fig:posExponencialStan2)Distribución posterior.</p>
</div>

\BeginKnitrBlock{example}<div class="example"><span class="example" id="exm:unnamed-chunk-65"><strong>(\#exm:unnamed-chunk-65) </strong></span>Suponga ahora que se va a realizar el trasplante de corazón a 5 pacientes, y se quiere conocer el comportamiento probabilístico del tiempo promedio de sobrevida en estos 5 pacientes. Aplicando la distribución predictiva y definiendo una distribución previa no informativa de Jeffreys, se tiene que

\begin{align*}
p(\bar{Y}^*|\mathbf{Y})&=\frac{5\Gamma(5+69)}{\Gamma(5)\Gamma(69)}\frac{25998.5^{69}}{(5\bar{Y}^*+25998.5)^{5+69}}(5\bar{Y}^*)^4\\
&=\frac{1}{\bar{Y}^*Beta(5,69)}\left(\frac{25998.5}{5\bar{Y}^*+25998.5}\right)^{69}\left(\frac{5\bar{Y}^*}{5\bar{Y}^*+25998.5}\right)^5
\end{align*}

El cálculo de esta función predictiva se puede llevar a cabo con el siguiente código en `R`, además de comprobar que la integral de la función es 1.</div>\EndKnitrBlock{example}


```r
pred_exp <- function(x){
  ((s/(s+x*n.mono))^n) * 
    ((x*n.mono/(s+x*n.mono))^n.mono) / 
    (x*beta(n,n.mono))
}


alfa <- beta <- 0
s <- 25998.5
n <- 69
n.mono <- 5
integrate(pred_exp, 0.0001, 10000)
```

```
## 1 with absolute error < 3.2e-10
```

La distribución predictiva de esta función se puede visualizar en la Figura \@ref(fig:PredPostExpoEje), donde se puede ver que la mayor masa de la función se acumula alrededor de los 300 días. Usando el comando `integrate(pred_exp, 800, 10000)`, también se puede observar que la probabilidad de que en promedio los cinco pacientes sobrevivan más de 800 días es de 0.0264485.

<div class="figure" style="text-align: center">
<img src="3Uniparametricos_files/figure-html/PredPostExpoEje-1.svg" alt="Distribución predictiva posterior para el tiempo promedio de sobrevivencia de transplante de corazón." width="576" />
<p class="caption">(\#fig:PredPostExpoEje)Distribución predictiva posterior para el tiempo promedio de sobrevivencia de transplante de corazón.</p>
</div>

## Modelo Normal con media desconocida

En esta sección se consideran datos que pueden ser descritos adecuadamente por medio de la distribución normal la cual, a diferencia de las anteriores distribuciones consideradas, tiene dos parámetros. En esta parte, se asume que la varianza teórica es conocida y el objetivo es estimar la media teórica. En el siguiente capítulo se considerará el caso general cuando ambos parámetros son desconocidos.

Suponga que $Y_1,\cdots,Y_n$ son variables independientes e idénticamente distribuidos con distribución $Normal(\theta,\sigma^2)$ con $\theta$ desconocido pero $\sigma^2$ conocido. De esta forma, la función de verosimilitud de los datos está dada por

\begin{align}
(\#eq:veronormal)
p(\mathbf{Y} \mid \theta)&=\prod_{i=1}^n\frac{1}{\sqrt{2\pi\sigma^2}}\exp\left\{-\frac{1}{2\sigma^2}(y_i-\theta)^2\right\}I_\mathbb{R}(y)\\
&=(2\pi\sigma^2)^{-n/2}\exp\left\{-\frac{1}{2\sigma^2}\sum_{i=1}^n(y_i-\theta)^2\right\}
\end{align}

Como el parámetro $\theta$ puede tomar cualquier valor en los reales, es posible asignarle una distribución previa $\theta \sim Normal(\mu,\tau^2)$. Bajo este marco de referencia se tienen los siguientes resultados

\BeginKnitrBlock{proposition}<div class="proposition"><span class="proposition" id="prp:unnamed-chunk-67"><strong>(\#prp:unnamed-chunk-67) </strong></span>La distribución posterior del parámetro de interés $\theta$ sigue una distribución

\begin{equation*}
\theta|\mathbf{Y} \sim Normal(\mu_n,\tau^2_n).
\end{equation*}

En donde

\begin{equation}
(\#eq:TauSigman)
\mu_n=\frac{\frac{n}{\sigma^2}\bar{Y}+\frac{1}{\tau^2}\mu}{\frac{n}{\sigma^2}+\frac{1}{\tau^2}}
\ \ \ \ \ \ \ \text{y} \ \ \ \ \ \ \
\tau_n^2=\left(\frac{n}{\sigma^2}+\frac{1}{\tau^2}\right)^{-1}
\end{equation}</div>\EndKnitrBlock{proposition}
<br>

\BeginKnitrBlock{proof}<div class="proof">\iffalse{} <span class="proof"><em>Prueba. </em></span>  \fi{}\begin{align*}
p(\theta \mid \mathbf{Y})&\propto p(\mathbf{Y} \mid \theta)p(\theta \mid \mu,\tau^2)\\
&\propto \exp\left\{-\frac{1}{2\sigma^2}\sum_{i=1}^n(y_i-\theta)^2-\frac{1}{2\tau^2}(\theta-\mu)^2\right\}\\
&= \exp\left\{-\frac{1}{2}\left[\frac{\sum_{i=1}^n(y_i-\theta)^2}{\sigma^2}+\frac{(\theta-\mu)^2}{\tau^2}\right]\right\}\\
&\propto \exp\left\{-\frac{1}{2}\left[\frac{n\theta^2}{\sigma^2}-\frac{2\theta\sum_{i=1}^ny_i}{\sigma^2}+\frac{\theta^2}{\tau^2}-\frac{2\theta\mu}{\tau^2}\right]\right\}\\
&= \exp\left\{-\frac{\theta^2}{2}\left[\frac{n}{\sigma^2}+\frac{1}{\tau^2}\right]+
\theta\left[\frac{n\bar{y}}{\sigma^2}+\frac{\mu}{\tau^2}\right]\right\}\\
&= \exp\left\{-\frac{\theta^2}{2\tau^2_n}+\frac{\theta\mu_n}{\tau_n^2}\right\}\\
&= \exp\left\{-\frac{1}{2\tau^2_n}(\theta^2-2\theta\mu_n)\right\}\\
&\propto \exp\left\{-\frac{1}{2\tau^2_n}(\theta^2-2\theta\mu_n+\mu_n^2)\right\}\\
&= \exp\left\{-\frac{1}{2\tau^2_n}(\theta-\mu_n)^2\right\}
\end{align*}

Por lo tanto, se encuentra una expresión idéntica a la función de distribución de una variable aleatoria con distribución $Normal(\mu_n,\tau_n^2)$.</div>\EndKnitrBlock{proof}
<br>

Observando la forma de $\mu_n$, que corresponde a la estimación bayesiana del parámetro $\theta$, podemos concluir que este es una combinación convexa entre el estimador clásico de máxima verosimlitud $\hat{\theta}_C=\bar{y}$ y el estimador previo $\hat{\theta}_P=\mu$, puesto que:

\begin{align*}
\hat{\theta}_B=\mu_n&=\frac{\frac{n}{\sigma^2}\bar{Y}+\frac{1}{\tau^2}\mu}{\frac{n}{\sigma^2}+\frac{1}{\tau^2}}\\
&=\frac{\frac{n}{\sigma^2}}{\frac{n}{\sigma^2}+\frac{1}{\tau^2}}\bar{Y}+\frac{\frac{1}{\tau^2}}{\frac{n}{\sigma^2}+\frac{1}{\tau^2}}\mu\\
&=\frac{\frac{n}{\sigma^2}}{\frac{n}{\sigma^2}+\frac{1}{\tau^2}}\hat{\theta}_C+\frac{\frac{1}{\tau^2}}{\frac{n}{\sigma^2}+\frac{1}{\tau^2}}\hat{\theta}_P
\end{align*}

De donde se puede concluir que, para una distribución previa fija, entre mayor sea el tamaño muestral $n$, más peso tendrá el estimador clásico $\hat{\theta}_C$ en la estimación bayesiana. De la misma forma, para un conjunto fijo de datos $\mathbf{Y}$, entre menor sea la varianza previa, $\tau^2$, más certeza tenemos sobre la información previa y por consiguiente la estimación bayesiana $\mu_n$ se acercará más a la estimación previa. En la Figura \@ref(fig:ComparaNormal) se observa la función de densidad previa, función de verosimilitud y función de densidad posterior con $\mu=5$, $\tau^2=0.01$, $\bar{y}=2$, $\sigma^2=1$ y $n=5,10,50,200$. Podemos observar que a medida que el tamaño muestral $n$ aumente, la función de verosimilitud (vista como la función del parámetro $\theta$) se vuelve más concentrada alrededor del valor de $\bar{y}$, y a consecuencia, la función de densidad posterior de $\theta$ se sitúa más cercana a la función de verosimilitud, y la estimación bayesiana se acerca más a la estimación clásica $\bar{y}$.

<div class="figure" style="text-align: center">
<img src="3Uniparametricos_files/figure-html/ComparaNormal-1.svg" alt="Distribución previa, función de verosimilitud y distribución posterior del parámetro $\theta$ con $\mu=5$, $\tau^2=0.01$, $\bar{y}=2$, $\sigma^2=1$ y $n=5,10,50,200$." width="576" />
<p class="caption">(\#fig:ComparaNormal)Distribución previa, función de verosimilitud y distribución posterior del parámetro $\theta$ con $\mu=5$, $\tau^2=0.01$, $\bar{y}=2$, $\sigma^2=1$ y $n=5,10,50,200$.</p>
</div>


### Distribución previa no informativa para $\theta$

Por otro lado, nótese que en el caso en donde se desconozca el comportamiento estructural de $\theta$, es posible definir su distribución previa tan plana y vaga como sea posible. Para esto, basta con hacer tender al parámetro de precisión de la distribución previa hacia infinito. Es decir $\tau^2 \longrightarrow \infty$, en este caso, la distribución previa de $\theta$ corresponde a una distribución impropia, $p(\theta)\propto cte$. Se puede ver que bajo esta distribución previa, la distribución posterior tendería a una $Normal(\bar{y},\sigma^2/n)$.

La anterior idea intuitiva de usar la distribución previa $p(\theta)\propto cte$ para representar la falta de la información *a priori* corresponde a la distribución previa no informativa de Jeffreys, puesto que la información de Fisher del parámetro $\theta$ en una variable con distribución normal está dada por 

\begin{equation*}
I(\theta)=1/\sigma^2
\end{equation*}

De donde se puede concluir que la previa no informativa de Jeffreys está dada por 
\begin{equation*}
p(\theta)\propto 1/\sigma\propto cte
\end{equation*}

Finalmente, es posible comparar los resultados inferenciales obtenidos con la previa no informativa de Jeffreys y el enfoque inferencial clásico en términos de la estimación puntual y el intervalo de credibilidad y de confianza. En cuanto a la estimación puntual, es claro que ambos enfoques conducen al mismo estimador $\hat{\theta}=\bar{Y}$. Con respecto al intervalo para el parámetro $\theta$, al usar el enfoque bayesiano con la previa no informativa de Jeffreys, el intervalo de credibilidad de $(1-\alpha)\times 100\%$ está dado por los percentiles $\alpha/2$ y $1-\alpha/2$ de la distribución posterior de $\theta$: $Normal(\bar{y},\sigma^2/n)$. Al denotar estos percentiles como $a$ y $b$, respectivamente. 
Por definición tenemos que, si $X\sim N(\bar{y},\sigma^2/n)$

\begin{align*}
\alpha/2&=Pr(X<a)\\
&=Pr(\frac{X-\bar{y}}{\sigma/\sqrt{n}}<\frac{a-\bar{y}}{\sigma/\sqrt{n}})\\
&=Pr(Z < \frac{a-\bar{y}}{\sigma/\sqrt{n}})
\end{align*}

Estos es, $\frac{a-\bar{y}}{\sigma/\sqrt{n}}$ es el percentil $\alpha/2$ de la distribución normal estándar $z_{\alpha/2}$ o equivalentemente $-z_{1-\alpha/2}$. De esta forma, tenemos que $a=\bar{y}-z_{1-\alpha/2}\ \sigma/\sqrt{n}$. Análogamente tenemos que $b=\bar{y}+z_{1-\alpha/2}\ \sigma/\sqrt{n}$, y podemos concluir que un intervalo de credibilidad de $(1-\alpha)\times 100\%$ está dada por $\bar{y}\pm z_{1-\alpha/2}\ \sigma/\sqrt{n}$, el cual coincide con el intervalo de confianza para $\theta$ usando el enfoque de la inferencia clásica [@Zhang].

### Diferentes formas de hallar la distribución previa para $\theta$

En primer lugar, considere el caso para el cual la información previa se encuentra en un conjunto de datos $x_1,\cdots,x_{m}$ que corresponden a mediciones de la variable de estudio en otro punto del tiempo, en otro punto geográfico, o inclusive en otra población de estudio. En este caso, podemos tomar la media de la distribución previa $\mu$ como $\bar{x}$ y la varianza de la distribución previa $\tau^2$ como $S^2_x$.

En el caso en el que no se disponga de datos como información previa, sino que esta esté contenida en alguna estimación que se haya realizado anteriormente sobre $\theta$. Por ejemplo, si se dispone de algún modelamiento estadístico que se haya hecho previamente sobre $\theta$, es posible fácilmente obtener el valor estimado de $\theta$ y el error estándar de esta estimación y naturalmente estos dos valores serían los parámetros de la distribución previa: $\mu$ y $\tau^2$. 

Finalmente, si la estimación previa de $\theta$ se presenta en forma de un intervalo; por ejemplo, si se sabe que un intervalo de confianza para $\theta$ es $(15.3,\ 24.7)$, entonces es posible definir a $\mu$ como el punto medio de este intervalo, es decir, $\mu=20$ y para escoger el valor de $\tau^2$ se tiene en cuenta que en muchas ramas de la estadística, un intervalo de confianza se puede aproximar por $\hat{\theta}\pm 2\sqrt{var(\hat{\theta})}$. De esta forma, se puede definir $\tau^2=\left(\frac{24.7-20}{2}\right)^2\approx5.5$

### Distribuciones predictivas

Los siguientes resultados presentan las distribuciones predictivas previa y predictiva posterior para una observación o una nueva muestra.

\BeginKnitrBlock{proposition}<div class="proposition"><span class="proposition" id="prp:unnamed-chunk-69"><strong>(\#prp:unnamed-chunk-69) </strong></span>La distribución predictiva previa para una observación $y$ estáa dada por 

\begin{equation*}
y \sim Normal (\mu, \tau^2+\sigma^2)
\end{equation*}</div>\EndKnitrBlock{proposition}
<br>

\BeginKnitrBlock{proof}<div class="proof">\iffalse{} <span class="proof"><em>Prueba. </em></span>  \fi{}De la definición de función de distribución predictiva se tiene que

\begin{align*}
p(Y)&=\int p(Y \mid \theta)p(\theta \mid \mu,\tau^2)\ d\theta\\
&=\int_{-\infty}^{\infty} \frac{1}{\sqrt{2\pi\sigma^2}}\exp\left\{-\frac{1}{2\sigma^2}(y-\theta)^2\right\}
\frac{1}{\sqrt{2\pi\tau^2}}\exp\left\{-\frac{1}{2\tau^2}(\theta-\mu)^2\right\}d\theta
\end{align*}

@Berger desarrolló las siguientes igualdades

\begin{align*}
&\ \ \ \ \frac{1}{2}\left[\frac{(\theta-\mu)^2}{\tau^2}+\frac{(y-\theta)^2}{\sigma^2}\right]\\
&=\frac{1}{2}\left[\left(\frac{1}{\tau^2}+\frac{1}{\sigma^2}\right)\theta^2-2\left(\frac{\mu}{\tau^2}+\frac{y}{\sigma^2}\right)\theta+\left(\frac{\mu^2}{\tau^2}+\frac{y^2}{\sigma^2}\right)\right]\\
&=\frac{1}{2\tau_1^2}\left[\theta^2-2\tau_1^2\left(\frac{\mu}{\tau^2}+\frac{y}{\sigma^2}\right)\theta+\tau_1^4\left(\frac{\mu}{\tau^2}+\frac{y}{\sigma^2}\right)^2\right]+\frac{1}{2}\left(\frac{\mu^2}{\tau^2}+\frac{y^2}{\sigma^2}\right)-\frac{\tau_1^2}{2}\left(\frac{\mu}{\tau^2}+\frac{y}{\sigma^2}\right)^2\\
&=\frac{1}{2\tau_1^2}\left[\theta-\tau_1^2\left(\frac{\mu}{\tau^2}+\frac{y}{\sigma^2}\right)\right]^2+\frac{1}{2}\left[\left(\frac{1}{\sigma^2}-\frac{\tau_1^2}{\sigma^4}\right)y^2-2\frac{\mu\tau_1^2}{\tau^2\sigma^2}y+\left(\frac{\mu^2}{\tau^2}-\frac{\mu^2\tau_1^2}{\tau^4}\right)\right]\\
&=\frac{1}{2\tau_1^2}\left[\theta-\mu_1\right]^2+\frac{1}{2}\left[\frac{1}{\sigma^2+\tau^2}y^2-2\frac{\mu}{\sigma^2+\tau^2}y+\frac{\mu^2}{\sigma^2+\tau^2}\right]\\
&=\frac{1}{2\tau_1^2}\left[\theta-\mu_1\right]^2+\frac{1}{2(\sigma^2+\tau^2)}(y-\mu)^2
\end{align*}

Entonces

\begin{align*}
p(Y)&=\int_{-\infty}^{\infty} \frac{1}{2\pi\sigma\tau}\exp\left\{-\frac{1}{2\tau_1^2}(\theta-\mu_1)^2\right\}
\exp\left\{-\frac{1}{2(\tau^2+\sigma^2)}(y-\mu)^2\right\}d\theta\\
&= \frac{1}{\sqrt{2\pi\frac{\sigma^2\tau^2}{\tau_1^2}}}\exp\left\{-\frac{1}{2(\tau^2+\sigma^2)}(y-\mu)^2\right\}
\int_{-\infty}^{\infty} \frac{1}{\sqrt{2\pi\tau_1^2}}\exp\left\{-\frac{1}{2\tau_1^2}(\theta-\mu_1)^2\right\}d\theta\\
&= \frac{1}{\sqrt{2\pi(\tau^2+\sigma^2)}}\exp\left\{-\frac{1}{2(\tau^2+\sigma^2)}(y-\mu)^2\right\}
\end{align*}</div>\EndKnitrBlock{proof}
<br>

Una vez recolectados los datos $\mathbf{Y}=\{Y_1,\cdots,Y_n\}$, se obtiene la distribución predictiva posterior dada en el siguiente resultado. La demostración es similar al del resultado anterior.

\BeginKnitrBlock{proposition}<div class="proposition"><span class="proposition" id="prp:unnamed-chunk-71"><strong>(\#prp:unnamed-chunk-71) </strong></span>La distribución predictiva posterior para una nueva observación $\tilde{y}$ es
\begin{equation*}
\tilde{y} \mid \mathbf{Y} \sim Normal (\mu_n, \tau_n^2+\sigma^2)
\end{equation*}</div>\EndKnitrBlock{proposition}
<br>

En algunas situaciones, se quiere conocer el comportamiento probabilístico de más de una nueva observación, digamos $Y_1^*,\cdots,Y_{n^*}^*$, en este caso, lo ideal sería obtener la distribución conjunta predictiva posterior de la nueva muestra, $p(Y_1^*,\cdots,Y_{n^*}^*|\mathbf{Y})$. Sin embargo, esta distribución no es fácil de hallar, por lo que el énfasis se pondrá en la distribución predictiva posterior de la media de esta nueva muestra $\bar{Y}^*$, la cual está dada en el siguiente resultado. 

\BeginKnitrBlock{proposition}<div class="proposition"><span class="proposition" id="prp:PredNorm"><strong>(\#prp:PredNorm) </strong></span>La distribución predictiva posterior para la media muestral $\bar{Y}^*$ de una nueva muestra es 

\begin{equation*}
\bar{Y}^*|\mathbf{Y}\ \sim N\left(\mu_n, \frac{\sigma^2}{n^*}+\tau^2_n\right)
\end{equation*}

donde $\mu_n$ y $\tau^2_n$ fueron definidos en \@ref(eq:TauSigman).</div>\EndKnitrBlock{proposition}
<br>

\BeginKnitrBlock{proof}<div class="proof">\iffalse{} <span class="proof"><em>Prueba. </em></span>  \fi{}\begin{align*}
p(\bar{Y}^*|\mathbf{Y})&=\int_{-\infty}^\infty p(\bar{Y}^*|\theta)p(\theta|\mathbf{Y})\ d\theta\\
&=\int_{-\infty}^\infty (2\pi\frac{\sigma^2}{n^*})^{-1/2}\exp\left\{-\frac{n^*}{2\sigma^2}(\bar{y}^*-\theta)^2\right\}
(2\pi\tau_n^2)^{-1/2}\exp\left\{-\frac{1}{2\tau_n^2}(\theta-\mu_n)^2\right\}\ d\theta\\
&=\int_{-\infty}^\infty (2\pi)^{-1}(\frac{\sigma^2}{n^*}\tau_n^2)^{-1/2}\exp\left\{-\frac{1}{2}\left[\frac{(\bar{y}^*-\theta)^2}{\sigma^2/n^*}+\frac{(\theta-\mu_n)^2}{\tau^2_n}\right]\right\}\ d\theta\\
&=\underbrace{\int_{-\infty}^\infty(2\pi\frac{1}{n^*/\sigma^2+1/\tau^2_n})^{-1/2}\exp\left\{-\frac{1}{2}\left(\frac{n^*}{\sigma^2}+\frac{1}{\tau^2_n}\right)\left(\theta-\frac{\bar{y}^*/(\sigma^2/n^*)+\mu_n/\tau^2_n}{n^*/\sigma^2+1/\tau^2_n}\right)^2\right\}\ d\theta}_{\text{igual a 1}}\\
&\ \ \ \ \ \ \ (2\pi)^{-1/2}(\frac{\sigma^2}{n^*}\tau_n^2)^{-1/2}(\frac{n^*}{\sigma^2}+\frac{1}{\tau^2_n})^{-1/2}\exp\left\{-\frac{1}{2(\sigma^2/n^*+\tau^2_n)}(\bar{y}^*-\mu_n)^2\right\}\\
&=(2\pi)^{-1/2}(\frac{\sigma^2}{n^*}\tau_n^2)^{-1/2}(\frac{n^*}{\sigma^2}+\frac{1}{\tau^2_n})^{-1/2}\exp\left\{-\frac{1}{2(\sigma^2/n^*+\tau^2_n)}(\bar{y}^*-\mu_n)^2\right\}\\
&=(2\pi)^{-1/2}(\frac{\sigma^2}{n^*}+\tau^2_n)^{-1/2}\exp\left\{-\frac{1}{2(\sigma^2/n^*+\tau^2_n)}(\bar{y}^*-\mu_n)^2\right\}
\end{align*}</div>\EndKnitrBlock{proof}
<br>

Del anterior resultado, podemos ver que la esperanza de la distribución de $\bar{Y}^*|\mathbf{Y}$ es igual a la esperanza de $\theta|\mathbf{Y}$. A diferencia de la varianza de $\theta|\mathbf{Y}$, la varianza de $\bar{Y}^*|\mathbf{Y}$ tiene un componente adicional dado por $\sigma^2/n^*$. De esta forma, existirán tres fuentes de incertidumbre al momento de pronosticar $\bar{Y}^*$: la incertidumbre en la información previa, la incertidumbre en la muestra observada y la incertidumbre en la nueva muestra.

\BeginKnitrBlock{example}<div class="example"><span class="example" id="exm:unnamed-chunk-73"><strong>(\#exm:unnamed-chunk-73) </strong></span>En @Zhang[Ej. 2.3.6], se reportan datos sobre el grosor de láminas de vidrio templado de 3 cm. para controlar la calidad de una línea de producción. Estos datos son 3.56, 3.36, 2.99, 2.71, 3.31, 3.68, 2.78, 2.95, 2.82, 3.45, 3.42 y 3.15, con promedio de 3.18 cm. Suponga que, por especificaciones técnicas, se conoce que la varianza del grosor es de $0.1 cm^2$. Por otro lado, como información previa, se conoce que en la última inspección de calidad el grosor promedio fue de 2.8 cm. con una desviación estándar de 0.23 cm.

De la anterior información, se puede decir que el parámetro de interés $\theta$ sería el grosor promedio de las láminas. También podemos afirmar que $\sigma^2=0.1cm^2$, $\bar{y}=3.18cm$, $n=12$, y los parámetros de la distribución previa estarían dados por $\mu=2.8cm$ y $\tau=0.45cm$. De esta forma, podemos calcular los parámetros de la distribución posterior

\begin{align*}
\mu_n
&=\dfrac{\frac{12}{0.1}3.18+\frac{1}{0.23^2}2.8}{\frac{12}{0.1}+\frac{1}{0.23^2}}=3.13cm\\
\tau^2_n
&=\left(\frac{12}{0.1}+\frac{1}{0.23^2}\right)^{-1}=0.007cm^2
\end{align*}

Entonces, la distribución posterior del grosor promedio será $N(\mu_n=3.13cm,\ \tau^2_n=0.007cm^2)$. Es posible concluir que la estimación bayesiana del parámetro de interés corresponde a $3.13cm$, mientras que para calcular un intervalo de credibilidad de 95% para el parámetro de interés, se debe calcular los percentiles 2.5% y 97.5% de la distribución posterior de $\theta$, dados por $(2.966cm,\ 3.293cm)$.

A continuación se ilustra el uso de `STAN` para obtener la estimación bayesiana del parámetro $\theta$.</div>\EndKnitrBlock{example}



```r
NormalMedia <- '
data {
  int<lower=0> n;
  real y[n];
}
parameters {
  real theta;
}
model {
  y ~ normal(theta, 0.1);
  theta ~ normal(2.8, 0.23);
}
'

n <- 12
y <- c(3.56, 3.36, 2.99, 2.71, 3.31, 3.68, 
       2.78, 2.95, 2.82, 3.45, 3.42, 3.15)

sample_data <- list(y = y, n = n)
NormalMfit <- stan(model_code = NormalMedia,
                   data = sample_data, verbose = FALSE)
```


Después de la convergencia del proceso inferencial, la estimación bayesiana de $\theta$ es 3.1745 cm, mientras que un intervalo de credibilidad del 95\% es $(3.117 cm,\ 3.232 cm)$, resultados muy similares a lo obtenido calculando directamente $\mu_n$ y $\tau^2_n$. 


```r
print(NormalMfit, digits = 4, 
      pars = "theta", probs = c(0.025, 0.975))
```

```
## Inference for Stan model: 1a46eb2f3eff113b3e52c62bed60aed2.
## 4 chains, each with iter=2000; warmup=1000; thin=1; 
## post-warmup draws per chain=1000, total post-warmup draws=4000.
## 
##         mean se_mean     sd   2.5%  97.5% n_eff   Rhat
## theta 3.1765   8e-04 0.0289 3.1195 3.2323  1426 1.0016
## 
## Samples were drawn using NUTS(diag_e) at Sun Jun  6 23:46:02 2021.
## For each parameter, n_eff is a crude measure of effective sample size,
## and Rhat is the potential scale reduction factor on split chains (at 
## convergence, Rhat=1).
```

Las figuras \@ref(fig:posNormalMStan) muestra la distribución posterior para este ejemplo, junto con la estimación puntual, correspondiente a la media. 


```r
bayesplot::mcmc_areas(NormalMfit, pars = "theta", 
                      prob = 0.95)
```

<div class="figure" style="text-align: center">
<img src="3Uniparametricos_files/figure-html/posNormalMStan-1.svg" alt="Distribución posterior." width="576" />
<p class="caption">(\#fig:posNormalMStan)Distribución posterior.</p>
</div>

En el siguiente ejemplo se ilustra el uso de la distribución predictiva posterior. 

\BeginKnitrBlock{example}<div class="example"><span class="example" id="exm:unnamed-chunk-76"><strong>(\#exm:unnamed-chunk-76) </strong></span>Suponga que la fábrica debe hacer un despacho de 8 láminas, y se quiere conocer sobre el grosor promedio del despacho $\bar{y}^*$. Usando el resultado \@ref(prp:PredNorm), se tiene que la disribución de $\bar{Y}^*$ condicionado en los 12 datos observados está dada por 

\begin{equation*}
\bar{Y}^*|\mathbf{Y}\ \sim N\left(\mu_n,\ \frac{\sigma^2}{n^*}+\tau^2_n\right) = N\left(3.13cm,\ \frac{0.1}{8}+0.007\right) = N(3.13cm,\ 0.0195cm^2)
\end{equation*}

De esta forma, se puede afirmar que el grosor promedio del despacho es de 3.13cm con un intervalo del 95% dado por (2.85cm, 3.40cm). Nótese que el intervalo para $\bar{Y}^{*}$ es más ancho que el intervalo para $\theta$, pues este tiene una varianza mayor a la varianza de la distribución posterior de $\theta$.</div>\EndKnitrBlock{example}

## Modelo Normal con varianza desconocida

En esta sección se presentan los fundamentos necesario para realizar inferencia bayesiana en un modelo normal para el cual sí se conoce la media, pero no su varianza. En casi todas las aplicaciones prácticas, es común que ambos parámetros sean desconocidos. Sin embargo, estas secciones servirán de base teórica para desarrollar modelos más complejos. Supóngase además que se cuenta con información previa, la cual puede basarse sin pérdida de generalidad en observaciones anteriores de alguna muestra de tamaño $n_0$ cuya varianza estimada fue $\sigma^2_0$.

En este orden de ideas, y siguiendo la argumentación de las secciones anteriores, dado que la varianza de la distribución es un parámetro que toma valores positivos únicamente, es plausible plantear que su distribución previa sea
\begin{equation*}
\sigma^2 \sim Inversa-Gamma(n_0/2, n_0\sigma^2_0/2)
\end{equation*}

Siguiendo la regla de Bayes, y después de factorizar convenientemente, se encuentra el siguiente resultado que muestra que la distribución posterior es conjugada con respecto a la previa. 

\BeginKnitrBlock{proposition}<div class="proposition"><span class="proposition" id="prp:unnamed-chunk-77"><strong>(\#prp:unnamed-chunk-77) </strong></span>La distribución posterior condicional de $\sigma^2$ es 
\begin{equation}
\sigma^2  \mid  \theta, \mathbf{Y} \sim Inversa-Gamma\left(\dfrac{n_0+n}{2},\dfrac{v_0}{2}\right)
\end{equation}
con $v_0=n_0\sigma^2_0+(n-1)S^2+n(\bar{y}-\theta)^2$.</div>\EndKnitrBlock{proposition}

A continuación se presenta un ejemplo de modelación de la varianza en una distribución normal. Como se vió en los anteriores resultados, la distribución previa conjugada toma la forma de una Gamma-inversa, sin embargo el investigador debe sopesar las ventajas de considerar esta distribución, puesto que las características del modelo bayesiano no sólo recaen en la definición de la verosimilitud, sino también en la estructura de la distribución previa. Es decir, en algunas ocasiones es mejor ponderar las bondades predictivas y de ajuste de los modelos, antes que escoger una distribución conjugada. Por ejemplo, @GelmanVariance2006 presentan diferentes opciones de distribuciones previas para modelar la varianza; sin embargo, el desarrollo teórico de estas otras posibles escogencias necesitan de herramientas metodológicas que se considerarán en los posteriores capítulos.

\BeginKnitrBlock{example}<div class="example"><span class="example" id="exm:8escuelas"><strong>(\#exm:8escuelas) </strong></span>@Gelman95[ , sección 5.5.] presentan un estudio realizado para analizar los efectos de algunos programas de entrenamiento especial en preparación para una prueba estandariza de opción múltiple llamada *SAT*, la cual es utilizada por las universidades para tomar decisiones con respecto a la adminsión de sus estudiantes. Quienes presenten mejores puntuaciones tienen una mayor probabilidad de ser admitidos para cursas sus estudios de educación superior. 

Las puntuaciones del *SAT* pueden variar entre 200 y 800, con una media de 500 y una desviación estándar de 100. Los exámenes SAT están diseñados para reflejar los conocimientos adquiridos y las habilidades desarrolladas durante muchos años de educación. Sin embargo, cada una de las ocho escuelas en este estudio consideró que su programa de entrenamiento a corto plazo fue muy exitoso para aumentar los puntajes de la prueba. Además, no había ninguna razón previa para creer que alguno de los ocho programas fuera más efectivo que cualquier otro o que algunos fueran más similares en efecto entre sí que cualquier otro.

La variable de interés en este ejemplo es el efecto del programa de entrenamiento en ocho escuelas secuendarias. Este efecto se define como la diferencia entre el puntaje promedio de la escuela con el promedio de la escala de la prueba (500 puntos). Además,  en una administración especial de la prueba para en ocho escuelas secundarias. La figura \@ref(fig:hist8escuelas) muestra la distribución de los datos observados.</div>\EndKnitrBlock{example}


```r
Escuelas <- data.frame(
  row.names=c("A", "B", "C", "D", "E", "F", "G", "H"),
  efecto = c(28.39, 7.94, -2.75, 6.82,
             -0.64, 0.63, 18.01, 12.16))
```

<div class="figure" style="text-align: center">
<img src="3Uniparametricos_files/figure-html/hist8escuelas-1.svg" alt="Histograma de los efectos en las ocho escuelas." width="576" />
<p class="caption">(\#fig:hist8escuelas)Histograma de los efectos en las ocho escuelas.</p>
</div>

En el siguiente apartado se hará uso de `STAN` para realizar la inferencia bayesiana del parámetro de interés. Aunque, como ya se vio en las secciones anteriores, es posible utilizar simplemente `R`, acudiendo a los percentiles de la distribución posterior. Para este ejemplo, no se asume un conocimiento específico del fenómeno en la distribución previa. Por consiguiente se definirán parámetros no informativos sobre la distribución previa. 

Notando la forma funcional de la distribución Gamma-inversa, al hacer que sus parámetros tomen valores muy pequeños $(\alpha \rightarrow \infty, \ \ \beta  \rightarrow \infty)$ se tiene que 

\begin{align*}
p(\sigma^2|\alpha,\beta) 
&\propto (\sigma^2)^{-\alpha-1}
\exp(-\frac{\beta}{\sigma^2}) \\
&\propto \frac{1}{\sigma^2}
\end{align*}

La cual coincide con la distribución previa de Jeffreys. De esta forma, para que la distribución previa sea propia (integral definida de la función de densidad), se escogen valores cercanos a cero, pero no nulos; por ejemplo $\alpha = 0.001, \ \ \beta = 0.001$. La figura \@ref(fig:jefgammainv) muestra la densidad previa no informativa para el parámetro de interés.

<div class="figure" style="text-align: center">
<img src="3Uniparametricos_files/figure-html/jefgammainv-1.svg" alt="Distribución previa no informativa para la varianza de una distribución Normal" width="576" />
<p class="caption">(\#fig:jefgammainv)Distribución previa no informativa para la varianza de una distribución Normal</p>
</div>

Es posible notar que esta distribución provee un intervalo de confianza del 95% entre 176079526 e infinito. Por lo que se puede concluir que hay una muy baja probabilidad de que el parámetro tome valores muy bajos o muy altos. Por ejemplo, $Pr(\sigma^2 < 10) = 0.00859688$.


```r
library(pscl)
sapply(c(0.025, 0.975), 
       function(x) qigamma(x, 0.001, 0.001))
```

```
## [1] 176079526       Inf
```

```r
pigamma(10, 0.001, 0.001)
```

```
## [1] 0.00859688
```


```r
NormalVar <- '
data {
  int<lower=0> n;
  real mu;
  real y[n];
}
parameters {
  real sigma;
}
transformed parameters {
  real sigma2;
  sigma2 = pow(sigma, 2);
}
model {
  y ~ normal(mu, sigma);
  sigma2 ~ inv_gamma(0.001, 0.001);
}
generated quantities {
  real y_test[n];
  for(i in 1:n) {
    y_test[i] = normal_rng(mu, sigma);
  }
}
'

sample_data <- list(y = Escuelas$efecto, 
                    n = nrow(Escuelas),
                    mu = mean(Escuelas$efecto))
NormalVfit <- stan(model_code = NormalVar,
                   data = sample_data, verbose = FALSE)
```

Después de la convergencia del proceso inferencial, la estimación bayesiana de $\sigma^2$ es 10.8, y de $\sigma^2$ es 127.7 cm. De la misma forma, un intervalo de credibilidad del 95\% para la desviación estándar es $(6.659,\ 19.172)$. 


```r
print(NormalVfit, digits = 4, 
      pars = c("sigma", "sigma2"), probs = c(0.025, 0.975))
```

```
## Inference for Stan model: db1d737f5b63e0b7dcdc9498c360b64f.
## 4 chains, each with iter=2000; warmup=1000; thin=1; 
## post-warmup draws per chain=1000, total post-warmup draws=4000.
## 
##            mean se_mean      sd    2.5%    97.5% n_eff  Rhat
## sigma   10.9047   0.116  3.1332  6.6574  18.8534   730 1.014
## sigma2 128.7262   3.116 82.6527 44.3207 355.4494   704 1.016
## 
## Samples were drawn using NUTS(diag_e) at Sat Jun 12 11:27:49 2021.
## For each parameter, n_eff is a crude measure of effective sample size,
## and Rhat is the potential scale reduction factor on split chains (at 
## convergence, Rhat=1).
```

Las figuras \@ref(fig:posNormalVStan) muestra la distribución posterior para este ejemplo, junto con la estimación puntual, correspondiente a la desviación estándar. 


```r
bayesplot::mcmc_areas(NormalVfit, pars = "sigma", 
                      prob = 0.95)
```

<div class="figure" style="text-align: center">
<img src="3Uniparametricos_files/figure-html/posNormalVStan-1.svg" alt="Distribución posterior." width="576" />
<p class="caption">(\#fig:posNormalVStan)Distribución posterior.</p>
</div>


En cualquier caso, es posible sugerir mejores distribuciones no informativas eligiendo apropiadamente un límite superior $U$ y uno inferior $L$ en lugar de parámetros $\alpha$ y $\beta$ en la distribución Gamma-inversa. Por lo general, los límites se pueden establecer con bastante facilidad luego de un poco de reflexión sobre lo que $\sigma^2$ realmente significa en el mundo real. Por ejemplo, si fuese el error en algún tipo de cantidad física, $L$ no puede ser más pequeño que el tamaño de un átomo, o el tamaño más pequeño que pueda observar en su experimento. Es más, $U$ no podría ser más grande que la tierra (o el sol si se quisiera ser realmente conservador). De esta manera, una forma de mantener las propiedades de invariancia, es definir $\nu \sim \mathrm{Uniform}(\ln(L),\ln(U))$ para que $\sigma\sim \exp(U(\ln(L),\ln(U))$. El siguiente código en `STAN` permite establecer esta opción inferencial.



```r
NormalVar2 <- '
data {
  int<lower = 0> n;
  real mu;
  real y[n];
}
parameters {
  real <lower = 0> nu;
}
transformed parameters {
  real sigma;
  real sigma2;
  sigma = exp(nu);
  sigma2 = pow(sigma, 2);
}
model {
  y ~ normal(mu, sigma);
  nu ~ uniform(0.5, 20);
}
generated quantities {
  real y_test[n];
  for(i in 1:n) {
    y_test[i] = normal_rng(mu, sigma);
  }
}
'

sample_data <- list(y = Escuelas$efecto, 
                    n = nrow(Escuelas),
                    mu = mean(Escuelas$efecto))
NormalVfit2 <- stan(model_code = NormalVar2,
                   data = sample_data, verbose = FALSE)
```

Después de la convergencia del proceso inferencial, la estimación bayesiana de $\sigma^2$ es 10.8, y de $\sigma^2$ es 128.2 cm. De la misma forma, un intervalo de credibilidad del 95\% para la desviación estándar es $(6.596,\ 18.795)$. 


```r
print(NormalVfit2, digits = 4, 
      pars = c("sigma", "sigma2"), probs = c(0.025, 0.975))
```

```
## Inference for Stan model: 9fedbb4c644064e4fdb958182968aebc.
## 4 chains, each with iter=2000; warmup=1000; thin=1; 
## post-warmup draws per chain=1000, total post-warmup draws=4000.
## 
##            mean se_mean       sd    2.5%    97.5% n_eff   Rhat
## sigma   10.9963  0.0896   3.3723  6.6836  18.8918  1417 1.0000
## sigma2 132.2874  2.7314 103.1860 44.6711 356.9005  1427 1.0003
## 
## Samples were drawn using NUTS(diag_e) at Sat Jun 12 11:43:07 2021.
## For each parameter, n_eff is a crude measure of effective sample size,
## and Rhat is the potential scale reduction factor on split chains (at 
## convergence, Rhat=1).
```

Las figuras \@ref(fig:posNormalVStan2) muestra la distribución posterior para este ejemplo, junto con la estimación puntual, correspondiente a la desviación estándar. 


```r
bayesplot::mcmc_areas(NormalVfit2, pars = "sigma", 
                      prob = 0.95)
```

<div class="figure" style="text-align: center">
<img src="3Uniparametricos_files/figure-html/posNormalVStan2-1.svg" alt="Distribución posterior." width="576" />
<p class="caption">(\#fig:posNormalVStan2)Distribución posterior.</p>
</div>


