

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

\BeginKnitrBlock{proposition}
<span class="proposition" id="prp:unnamed-chunk-1"><strong>(\#prp:unnamed-chunk-1) </strong></span>La distribución posterior del parámetro $\theta$ sigue una distribución
\begin{equation*}
\theta \mid Y \sim Beta(y+\alpha,\beta-y+1)
\end{equation*}
\EndKnitrBlock{proposition}

\BeginKnitrBlock{proof}
\iffalse{} <span class="proof"><em>Prueba. </em></span>  \fi{}\begin{align*}
p(\theta \mid Y)&\propto 
p(Y \mid \theta)p(\theta \mid \alpha,\beta)\\
&=\frac{I_{\{0,1\}}(y)}{Beta(\alpha,\beta)}\theta^y\theta^{\alpha-1}(1-\theta)^{\beta-1}(1-\theta)^{1-y}I_{[0,1]}(\theta)\\
&\propto \theta^{y+\alpha-1}(1-\theta)^{\beta-y+1-1}I_{[0,1]}(\theta)
\end{align*}
    
Por lo tanto, factorizando convenientemente, se encuentra una expresión idéntica a la función de distribución de una variable aleatoria con distribución $Beta(y+\alpha,\beta-y+1)$.
\EndKnitrBlock{proof}
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

![(\#fig:jefber1)Distribución previa no informativa de Jeffreys para el parámetro de una distribución Bernoulli](3Uniparametricos_files/figure-latex/jefber1-1.pdf) 

\BeginKnitrBlock{proposition}
<span class="proposition" id="prp:unnamed-chunk-3"><strong>(\#prp:unnamed-chunk-3) </strong></span>La distribución predictiva previa para una observación $y$ está dada por
\begin{equation}
(\#eq:Predipreviabernou)
p(Y)=
\frac{Beta(y+\alpha,\beta-y+1)}{Beta(\alpha,\beta)}I_{\{0,1\}}(y)
\end{equation}
La cual define una auténtica función de densidad de probabilidad continua.
\EndKnitrBlock{proposition}
<br>

\BeginKnitrBlock{proof}
\iffalse{} <span class="proof"><em>Prueba. </em></span>  \fi{}De la definición de función de distribución predictiva se tiene que

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

Lo cual se verifica fácilmente teniendo en cuenta las propiedades de la función matemática Beta y de la función matemática Gamma.
\EndKnitrBlock{proof}
<br>

La distribución predictiva dada en \@ref(eq:Predipreviabernou) está
basada únicamente en la distribución previa del parámetro $\theta$. Una
vez observada la variable $Y$ se puede pensar en actualizar la
distribución predictiva basando la inferencia en la distribución
posterior del parámetro; esta distribución se da en el siguiente
resultado.

\BeginKnitrBlock{proposition}
<span class="proposition" id="prp:unnamed-chunk-5"><strong>(\#prp:unnamed-chunk-5) </strong></span>Después de la recolección de los datos, la distribución predictiva posterior para una nueva observación $\tilde{y}$ está dada por

\begin{equation}
p(\tilde{y} \mid Y)=
  \frac{Beta(\tilde{y}+y+\alpha,\beta-\tilde{y}-y+2)}{Beta(y+\alpha,\beta-y+1)}I_{\{0,1\}}(\tilde{y}),
\end{equation}
\EndKnitrBlock{proposition}
<br>

\BeginKnitrBlock{proof}
\iffalse{} <span class="proof"><em>Prueba. </em></span>  \fi{}De la definición de función de distribución predictiva se tiene que

\begin{align*}
p(\tilde{y} \mid Y)&=\int p(\tilde{y} \mid \theta)p(\theta \mid Y)\ d\theta\\
&=\int_0^1\theta^{\tilde{y}}(1-\theta)^{1-\tilde{y}}I_{\{0,1\}}(\tilde{y})\frac{\theta^{y+\alpha-1}(1-\theta)^{\beta-y+1-1}}{Beta(y+\alpha,\beta-y+1)}\ d\theta\\
&=\frac{Beta(\tilde{y}+y+\alpha,\beta-\tilde{y}-y+2)}{Beta(y+\alpha,\beta-y+1)}I_{\{0,1\}}(\tilde{y})\\
&\hspace{2cm}\times \int_0^1\frac{\theta^{\tilde{y}+y+\alpha-1}(1-\theta)^{\beta-\tilde{y}-y+2-1}}{Beta(\tilde{y}+y+\alpha,\beta-\tilde{y}-y+2)}\ d\theta\\
&=\frac{Beta(\tilde{y}+y+\alpha,\beta-\tilde{y}-y+2)}{Beta(y+\alpha,\beta-y+1)}I_{\{0,1\}}(\tilde{y})
\end{align*}
\EndKnitrBlock{proof}
<br>

En la práctica rara vez se observa la realización de una única variable
aleatoria Bernoulli $Y$, sino una muestra de variables aleatorias $Y_1$,
$\cdots$, $Y_n$. En este caso, la distribución posterior del parámetro
$\theta$ está dada en el siguiente resultado.

\BeginKnitrBlock{proposition}
<span class="proposition" id="prp:unnamed-chunk-7"><strong>(\#prp:unnamed-chunk-7) </strong></span>Cuando se tiene una muestra aleatoria $Y_1,\ldots,Y_n$ de variables con distribución Bernoulli de parámetro $\theta$, entonces la distribución posterior del parámetro de interés es

\begin{equation*}
\theta \mid Y_1,\ldots,Y_n \sim Beta\left(\sum_{i=1}^ny_i+\alpha,\beta-\sum_{i=1}^ny_i+n\right)
\end{equation*}
\EndKnitrBlock{proposition}

\BeginKnitrBlock{example}
<span class="example" id="exm:bernoelectoral"><strong>(\#exm:bernoelectoral) </strong></span>Es común en muchos países del mundo que se presenten encuestas de opinión electoral unas semanas antes de las elecciones presidenciales. Dentro de este tipo de encuestas se acostumbra a indagar acerca del favoritismo de los candidatos involucrados en la contienda electoral. Suponga que el candidato presidencial A está interesado en conocer su intención de voto previa a las elecciones. Para esto, él contrata a una firma encuestadora para la realización de una encuesta entre la población votante. El resultado de este estudio puede hacer cambiar o afirmar las estrategias publicitarias y la redefinición de la campaña electoral. La firma encuestadora decide implementar una estrategia de muestreo con un tamaño de muestra de doce mil personas. A cada respondiente se le realiza la siguiente pregunta: 
  
  > Si las elecciones presidenciales fueran mañana. ¿Usted votaría por el candidato A?
    
Las respuestas a esta pregunta son realizaciones de una muestra aleatoria de doce mil variables con densidad Bernoulli. Los resultados del estudio arrojan que 6360 personas de las personas entrevistadas, es decir un 53%, votarían por el suscrito candidato. Técnicamente se debe analizar esta cifra puesto que las implicaciones de ganar en una primera vuelta son grandes en el sentido económico, logístico y administrativo. Claramente, el dato 53% asegura una ventaja dentro de la muestra de doce mil personas. Sin embargo, es necesario realizar un estudio más profundo acerca de la caracterización estructural de la intención de voto del candidato en el electorado.
    
Con base en lo anteriormente expuesto, se decide utilizar la inferencia bayesiana puesto que existe información previa de un estudio anterior, contratado por el mismo candidato unos meses atrás en donde se entrevistaron a mil personas, con un favoritismo que estaba alrededor del 35 por ciento. Esta situación conlleva a la utilización de la metodología bayesiana que incorpora la información pasada acerca del mismo fenómeno.
    
El estadístico de la firma encuestadora decide utilizar una distribución previa Beta, definiendo los parámetros de la distribución previa como $\alpha$ igual al número de votantes a favor y $\beta$ igual al número de votantes en contra. Es decir, $Beta(\alpha=350, \beta=650)$. Por lo anterior, la distribución posterior del parámetro de interés, que representa la probabilidad de éxito en las elecciones presidenciales, es $Beta(6360+350, 650-6360+12000)=Beta(6710, 6290)$. Por lo tanto, utilizando la distribución posterior, se estima que la intención de voto por el candidato es de $\frac{6710}{6710+6290}=\frac{6710}{13000}=0.516$ y este valor equivale a la media de la distribución posterior. 
    
Sin embargo, si no se tuviese información previa como la suministrada por el estudio de meses anteriores, el análisis bayesiano sugeriría trabajar con una distribución previa no informativa, que en este caso, correspondería a una $Beta(\alpha=0.5, \beta=0.5)$. siguiendo el mismo análisis, se tiene que la distribución posterior es $Beta(6360.5, 5640.5)$. Finalmente, se estimaría que la intención de voto por el candidato es de $\frac{6350.5}{12001}=0.529$. 

La figuras \@ref(fig:BernoEj1) muestra el comportamiento de las distribuciones previas y posteriores en ambos escenarios. Nótese que la distribución no informativa influye muy poco en el comportamiento de la distribución posterior.
\EndKnitrBlock{example}

![(\#fig:BernoEj1)Distribuciones previas (línea punteada) y posteriores (línea sólida) para el ejemplo de las encuestas electorales.](3Uniparametricos_files/figure-latex/BernoEj1-1.pdf) 

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
##         mean se_mean     sd   2.5% 97.5% n_eff   Rhat
## theta 0.5155   1e-04 0.0042 0.5071 0.524  1591 1.0009
## 
## Samples were drawn using NUTS(diag_e) at Sat Jun  5 00:10:46 2021.
## For each parameter, n_eff is a crude measure of effective sample size,
## and Rhat is the potential scale reduction factor on split chains (at 
## convergence, Rhat=1).
```

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


\BeginKnitrBlock{proposition}
<span class="proposition" id="prp:unnamed-chunk-11"><strong>(\#prp:unnamed-chunk-11) </strong></span>La distribución posterior del parámetro $\theta$ sigue una distribución

\begin{equation*}
\theta \mid S \sim Beta(s+\alpha,\beta-s+n)
\end{equation*}
\EndKnitrBlock{proposition}
<br>

\BeginKnitrBlock{proof}
\iffalse{} <span class="proof"><em>Prueba. </em></span>  \fi{}\begin{align*}
p(\theta \mid S)&\propto p(S \mid \theta)p(\theta \mid \alpha,\beta)\\
&=\frac{\binom{n}{s}I_{\{0,1,\ldots,n\}}(s)}{Beta(\alpha,\beta)}
\theta^s\theta^{\alpha-1} (1-\theta)^{\beta-1}(1-\theta)^{n-s}I_{[0,1]}(\theta)\\
&\propto \theta^{s+\alpha-1} (1-\theta)^{\beta-s+n-1}I_{[0,1]}(\theta)
\end{align*}

Por lo tanto, factorizando convenientemente, se llega a una expresión idéntica a la función de distribución de una variable aleatoria con distribución $Beta(s+\alpha,\beta-s+n)$.
\EndKnitrBlock{proof}
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

![(\#fig:beta3)Funciones de verosimilitud, previa y posterior para $\alpha=2$, $\beta=5$, $s=8$ y $n=10$.](3Uniparametricos_files/figure-latex/beta3-1.pdf) 

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

![(\#fig:betan)Estimación posterior de $\theta$ para diferentes valores de $n$ y $s$ con $\alpha=\beta=5$.](3Uniparametricos_files/figure-latex/betan-1.pdf) 

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


\BeginKnitrBlock{proposition}
<span class="proposition" id="prp:unnamed-chunk-13"><strong>(\#prp:unnamed-chunk-13) </strong></span>La distribución predictiva previa para la observación particular de la suma de variables aleatorias Bernoulli, $s$, está dada por una distribución Beta-Binomial

\begin{equation}
p(S)=\binom{n}{s}\frac{Beta(s+\alpha,\beta-s+n)}{Beta(\alpha,\beta)}I_{\{0,1,\ldots,n\}}(s)
\end{equation}
\EndKnitrBlock{proposition}
<br>

\BeginKnitrBlock{proof}
\iffalse{} <span class="proof"><em>Prueba. </em></span>  \fi{}De la definición de función de distribución predictiva previa se tiene que

\begin{align*}
p(S)
&=\int p(S \mid \theta)p(\theta \mid \alpha,\beta)\ d\theta\\
&=\int_0^1\binom{n}{s}\theta^s(1-\theta)^{n-s}I_{\{0,1,\ldots,n\}}(s)
\frac{1}{Beta(\alpha,\beta)}\theta^{\alpha-1}(1-\theta)^{\beta-1}\ d\theta\\
&=\binom{n}{s}\frac{Beta(s+\alpha,\beta-s+n)}{Beta(\alpha,\beta)}I_{\{0,1,\ldots,n\}}(s)\\
&\hspace{2cm}\times\int_0^1\frac{\theta^{s+\alpha-1}(1-\theta)^{\beta-s+n-1}}{Beta(s+\alpha,\beta-s+n)}\ d\theta\\
&=\binom{n}{s}\frac{Beta(s+\alpha,\beta-s+n)}{Beta(\alpha,\beta)}I_{\{0,1,\ldots,n\}}(s)
\end{align*}
\EndKnitrBlock{proof}
<br>

Una vez observados los valores muestrales, podemos encontrar la
distribución predictiva posterior para una nueva variable binomial
$\tilde{S}$ en una muestra de tamaño $\tilde{n}$. Esta distribución se
encuentra en el siguiente resultado.

\BeginKnitrBlock{proposition}
<span class="proposition" id="prp:ResPredBinom"><strong>(\#prp:ResPredBinom) </strong></span>Después de la recolección de los datos $y_1$, $\cdots$, $y_n$, la distribución predictiva posterior para una nueva variable $\tilde{S}$ en una muestra del tamaño $\tilde{n}$ está dada por

\begin{equation}
(\#eq:Binompredict)
p(\tilde{s} \mid S)=\binom{\tilde{n}}{\tilde{s}}\frac{Beta(\tilde{s}+s+\alpha,\beta-\tilde{s}-s+n+\tilde{n})}{Beta(s+\alpha,\beta-s+n)}I_{\{0,1,\ldots,\tilde{n}\}}(\tilde{s}),
\end{equation}
\EndKnitrBlock{proposition}
<br>

\BeginKnitrBlock{proof}
\iffalse{} <span class="proof"><em>Prueba. </em></span>  \fi{}De la definición de función de distribución predictiva se tiene que

\begin{align*}
p(\tilde{s} \mid S)&=\int p(\tilde{s} \mid \theta)p(\theta \mid S)\ d\theta\\
&=\int_0^1 \binom{\tilde{n}}{\tilde{s}} \theta^{\tilde{s}}(1-\theta)^{\tilde{n}-\tilde{s}}I_{\{0,1,\ldots,\tilde{n}\}}(\tilde{s})
\frac{\theta^{s+\alpha-1}(1-\theta)^{\beta-s+n-1}}{Beta(s+\alpha,\beta-s+n)}\ d\theta\\
&=\binom{\tilde{n}}{\tilde{s}}\frac{Beta(\tilde{s}+s+\alpha,\beta-\tilde{s}-s+n+\tilde{n})}{Beta(s+\alpha,\beta-s+n)}I_{\{0,1,\ldots,\tilde{n}\}}(\tilde{s})\\
& \hspace{2cm}\times
\int_0^1\frac{\theta^{\tilde{s}+s+\alpha-1}(1-\theta)^{\beta-\tilde{s}-s+n+\tilde{n}-1}}
{Beta(\tilde{s}+s+\alpha,\beta-\tilde{s}-s+n+\tilde{n})}\ d\theta\\
&=\binom{\tilde{n}}{\tilde{s}}\frac{Beta(\tilde{s}+s+\alpha,\beta-\tilde{s}-s+n+\tilde{n})}{Beta(s+\alpha,\beta-s+n)}I_{\{0,1,\ldots,\tilde{n}\}}(\tilde{s})
\end{align*}
\EndKnitrBlock{proof}
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

Nótese que esta probabilidad es la misma contenida en la posición 501 del objeto `res` igual a \ensuremath{5.969157\times 10^{-4}}. Finalmente, se observa que la distribución predictiva \@ref(eq:Binompredict) corresponde a una distribución
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

Podemos observar que la posición 501 del objeto `res2` es igual a \ensuremath{5.969157\times 10^{-4}}, el cual es idéntico a lo obtenido en `res`. Adicionalmente, al escribir la distribución predictiva de
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

\BeginKnitrBlock{proposition}
<span class="proposition" id="prp:postbinom"><strong>(\#prp:postbinom) </strong></span>Cuando se tiene una sucesión de variables aleatorias $S_1,\ldots,S_i, \ldots,S_k$ independientes y con distribución $Binomial(n_i,\theta)$ para $i=1,\ldots,k$, entonces la distribución posterior del parámetro de interés $\theta$ es

\begin{equation*}
\theta \mid S_1,\ldots,S_k \sim Beta\left(\sum_{i=1}^ks_i+\alpha,\beta+\sum_{i=1}^k n_i-\sum_{i=1}^k s_i\right)
\end{equation*}
\EndKnitrBlock{proposition}
<br>

\BeginKnitrBlock{proof}
\iffalse{} <span class="proof"><em>Prueba. </em></span>  \fi{}\begin{align*}
p(\theta \mid S_1,\ldots,S_k)&\propto \prod_{i=1}^kp(S_i \mid \theta)p(\theta \mid \alpha,\beta)\\
&\propto \prod_{i=1}^k\theta^{\sum_{i=1}s_i}\theta^{\alpha-1}(1-\theta)^{\beta-1}
(1-\theta)^{\sum_{i=1}^kn_i-\sum_{i=1}^ks_i}I_{[0,1]}(\theta)\\
&= \theta^{\sum_{i=1}^ks_i+\alpha-1}(1-\theta)^{\sum_{i=1}^kn_i-\sum_{i=1}^ks_i+\beta}I_{[0,1]}(\theta)
\end{align*}

Por lo tanto, factorizando convenientemente, se encuentra una expresión idéntica a la función de densidad de la distribución $Beta\left(\sum_{i=1}^ks_i+\alpha,\beta+\sum_{i=1}^k n_i-\sum_{i=1}^n s_i\right)$.
\EndKnitrBlock{proof}
<br>

\BeginKnitrBlock{example}
<span class="example" id="exm:unnamed-chunk-23"><strong>(\#exm:unnamed-chunk-23) </strong></span>El siguiente conjunto de datos fue estudiado inicialmente por @Efron75 y se ha convertido en uno de los ejemplos prácticos más citados en la historia de la estadística moderna. Se trata de los porcentajes de bateo en una muestra de 18 jugadores profesionales en la temporada regular de béisbol en Estados Unidos en el año 1970. @wikiBat establece que, en términos generales, este valor representa la razón entre la cantidad de *hits* y el número de turnos al bate^[Un *hit* es la conexión efectuada por el bateador que coloca la pelota dentro del terreno de juego, permitiéndole alcanzar al menos una base, sin que se produzca un error de defensa del equipo contrario. Para lograr un hit, el bateador debe llegar a primera base antes de que ningún jugador defensivo lo toque con la bola en el trayecto del home a la inicial, o que el jugador de la defensa que tenga la bola pise la primera base antes que el bateador llegue a la misma.]. La fórmula para calcular esta estadística es $s/n$, donde $s$ es el número de *hits* y $n$ es el total de turnos. Este conjunto de datos está disponible en el paquete `pscl` de `R` y se puede cargar mediante el siguiente código computacional.
\EndKnitrBlock{example}



```r
library(pscl)
data(EfronMorris)
```


\begin{tabular}{l|l|l|r|r|r|r}
\hline
name & team & league & r & y & n & p\\
\hline
Roberto Clemente & Pitts & NL & 18 & 0.400 & 367 & 0.346\\
\hline
Frank Robinson & Balt & AL & 17 & 0.378 & 426 & 0.298\\
\hline
Frank Howard & Wash & AL & 16 & 0.356 & 521 & 0.276\\
\hline
Jay Johnstone & Cal & AL & 15 & 0.333 & 275 & 0.222\\
\hline
Ken Berry & Chi & AL & 14 & 0.311 & 418 & 0.273\\
\hline
Jim Spencer & Cal & AL & 14 & 0.311 & 466 & 0.270\\
\hline
Don Kessinger & Chi & NL & 13 & 0.289 & 586 & 0.263\\
\hline
Luis Alvarado & Bos & AL & 12 & 0.267 & 138 & 0.210\\
\hline
Ron Santo & Chi & NL & 11 & 0.244 & 510 & 0.269\\
\hline
Ron Swoboda & NY & NL & 11 & 0.244 & 200 & 0.230\\
\hline
Del Unser & Wash & AL & 10 & 0.222 & 277 & 0.264\\
\hline
Billy Williams & Chi & AL & 10 & 0.222 & 270 & 0.256\\
\hline
George Scott & Bos & AL & 10 & 0.222 & 435 & 0.303\\
\hline
Rico Petrocelli & Bos & AL & 10 & 0.222 & 538 & 0.264\\
\hline
Ellie Rodriguez & KC & AL & 10 & 0.222 & 186 & 0.226\\
\hline
Bert Campaneris & Oak & AL & 9 & 0.200 & 558 & 0.285\\
\hline
Thurman Munson & NY & AL & 8 & 0.178 & 408 & 0.316\\
\hline
Max Alvis & Mil & NL & 7 & 0.156 & 70 & 0.200\\
\hline
\end{tabular}

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
##         mean se_mean     sd  2.5%  97.5% n_eff  Rhat
## theta 0.2736   1e-04 0.0051 0.264 0.2838  1311 1.003
## 
## Samples were drawn using NUTS(diag_e) at Sat Jun  5 00:11:22 2021.
## For each parameter, n_eff is a crude measure of effective sample size,
## and Rhat is the potential scale reduction factor on split chains (at 
## convergence, Rhat=1).
```

Por otro lado, el mismo intervalo de credibilidad del 95% correspondiente se puede hallar mediante el siguiente código computacional de `R`.


```r
qbeta(c(0.025, 0.975), 2040, 5419)
```

```
## [1] 0.2634379 0.2836674
```

La figura \@ref(fig:BinomEj1) muestra el comportamiento de las distribuciones previa y posterior para este ejemplo. Nótese que, con un análisis frecuentista, se hubiese llegado a una estimación cercana de $\frac{1825}{6649}=0.274$. 



```r
p <- ggplot(data = data.frame(x = 0), 
            mapping = aes(x = x)) + 
  ylab("") + xlab(expression(theta))
previa <- function(x) dbeta(x, 215, 595)
posterior <- function(x) dbeta(x, 2040, 5419)
p + 
  stat_function(stat="function", fun = previa,
                mapping = aes(linetype="solid")) +
  stat_function(stat="function", fun = posterior, 
                mapping = aes(linetype="dashed")) +
	xlim(0.2,0.35) +
	scale_linetype_manual(name="", values=c("solid", "dashed"),
		label=c("Posterior","Previa")) 
```

![(\#fig:BinomEj1)Función de densidad previa y función de densidad posterior para el ejemplo de bateo.](3Uniparametricos_files/figure-latex/BinomEj1-1.pdf) 

Es posible analizar este conjunto de datos desde otra perspectiva al suponer que los jugadores no constituyen una muestra aleatoria y cada uno de ellos tiene un promedio de bateo diferente. Sin embargo, este análisis se deja como ejercicio en un capítulo posterior.


\BeginKnitrBlock{example}
<span class="example" id="exm:unnamed-chunk-29"><strong>(\#exm:unnamed-chunk-29) </strong></span>Continuando con el conjunto de datos de Efron y Morris, suponga que el entrenador de un equipo de las ligas inferiores está interesado en adquirir los servicios de Max Alvis. Este jugador no tuvo un buen promedio de bateo en la temporada y no tuvo muchos turnos al bate. El entrenador quiere conocer cuál será el número más probable de *hits* que anotará en la siguiente temporada. Teniendo en cuenta que es un jugador que viene de la liga profesional, lo más conveniente es que tenga muchos turnos al bate, digamos 400.

Para resolver este cuestionamiento, es conveniente recurrir a la función predictiva posterior, dada en el resultado \@ref(prp:ResPredBinom). Para este análisis, se define la caracterización estructural de la distribución previa del jugador que está dada por una $Beta(\alpha=7, \beta=38)$. La siguiente función en `R` permite obtener la distribución predictiva para este jugador, que se muestra en la figura \@ref(fig:PredBinom).
\EndKnitrBlock{example}


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


```r
predic <- data.frame(time = c(1:(n.ast + 1)), 
                     y = predictiva)
ggplot(predic, aes(time, y)) + 
  geom_line() + 
  xlab(expression(tilde(y))) + 
  ylab("Distribución posterior predictiva")
```

![(\#fig:PredBinom)Función de densidad predictiva posterior para el jugador Max Alvis.](3Uniparametricos_files/figure-latex/PredBinom-1.pdf) 

