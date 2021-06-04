

# Modelos uniparamétricos

Los modelos que están definidos en términos de un solo parámetro que pertenece al conjunto de los números reales se definen como modelos *uniparamétricos*. Este capítulo estudia modelos, discretos y continuos, que son comunes de implementar en la práctica. Dado que todos ellos son inducidos por familias de probabilidad conjugadas, entonces las estimaciones posteriores para los parámetros pueden hallarse sin necesidad de sofisticaciones computacionales. Es decir, con el uso de una simple calculadora de bolsillo, es posible realizar inferencia bayesiana propiamente dicha. Por lo tanto, en este capítulo, será menor el uso de software estadístico. Sin embargo, para cada modelo se incluye la sintaxis de programación en `R` y en `STAN` junto con ejemplos prácticos que permiten la familiarización e interiorización del ambiente computacional de este software que será indispensable en el desarrollo de capítulos posteriores.
    
## Modelo Bernoulli
    
Suponga que $Y$ es una variable aleatoria con distribución Bernoulli dada por:
\begin{equation}
p(Y \mid \theta)=\theta^y(1-\theta)^{1-y}I_{\{0,1\}}(y),
\end{equation}
    
Como el parámetro $\theta$ está restringido al espacio $\Theta=[0,1]$, entonces es posible formular varias opciones para la distribución previa del parámetro. En particular, la distribución uniforme restringida al intervalo $[0,1]$ o la distribución Beta parecen ser buenas opciones. Puesto que la distribución uniforme es un caso particular de la distribución Beta, entonces se iniciará con ésta. Por lo tanto la distribución previa del parámetro $\theta$ estará dada por

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
    
Del anterior resultado, podemos ver que la familia de distribuciones Beta es conjugada con respecto a la familia de distribuciones Bernoulli. Ahora consideremos cuál sería la distribución previa no informativa de Jeffreys para el parámetro $\theta$. De acuerdo a la definición \@ref(def:Jeffreys), se tiene que

\begin{equation*}
p(\theta)\propto I(\theta) ^{1/2}
\end{equation*}
    
En donde $I(\theta)$ es la información de Fisher del parámetro $\theta$, que en este caso está dada por

\begin{align*}
I(\theta)&=
-E\left\{\dfrac{\partial^2}{\partial\theta^2}\log{p(\mathbf{Y}\mid\theta)}\right\}\\
&=-E\left\{\dfrac{\partial^2}{\partial\theta^2}\{Y\log\theta+(1-Y)\log(1-\theta)\}\right\}\\
&=E\left\{\frac{Y}{\theta^2}+\frac{1-Y}{(1-\theta)^2}\right\}\\
&=\frac{1}{\theta(1-\theta)}
\end{align*}
    
De esta forma, la distribución previa no informativa de Jeffreys debe ser proporcional a $\theta^{-1/2}(1-\theta)^{-1/2}$, que asimismo corresponde a la distribución $Beta(1/2,1/2)$, cuya función de densidad se muestra en la figura \@ref(fig:jefber1) la cual asigna iguales pesos a los valores extremos del parámetro de interés y su característica de ser no informativa se representa en la simetría de la función alrededor del valor 0.5.

<div class="figure">
<img src="3Uniparametricos_files/figure-html/jefber1-1.png" alt="Distribución previa no informativa de Jeffreys para el parámetro de una distribución Bernoulli" width="672" />
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
    
La distribución predictiva dada en \@ref(eq:Predipreviabernou) está basada únicamente en la distribución previa del parámetro $\theta$. Una vez observada la variable $Y$ se puede pensar en actualizar la distribución predictiva basando la inferencia en la distribución posterior del parámetro; esta distribución se da en el siguiente resultado.
    
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
    
En la práctica rara vez se observa la realización de una única variable aleatoria Bernoulli $Y$, sino una muestra de variables aleatorias $Y_1$, $\cdots$, $Y_n$. En este caso, la distribución posterior del parámetro $\theta$ está dada en el siguiente resultado.
    
\BeginKnitrBlock{proposition}<div class="proposition"><span class="proposition" id="prp:unnamed-chunk-7"><strong>(\#prp:unnamed-chunk-7) </strong></span>Cuando se tiene una muestra aleatoria $Y_1,\ldots,Y_n$ de variables con distribución Bernoulli de parámetro $\theta$, entonces la distribución posterior del parámetro de interés es

\begin{equation*}
\theta \mid Y_1,\ldots,Y_n \sim Beta\left(\sum_{i=1}^ny_i+\alpha,\beta-\sum_{i=1}^ny_i+n\right)
\end{equation*}</div>\EndKnitrBlock{proposition}

\BeginKnitrBlock{example}<div class="example"><span class="example" id="exm:unnamed-chunk-8"><strong>(\#exm:unnamed-chunk-8) </strong></span>Es común en muchos países del mundo que se presenten encuestas de opinión electoral unas semanas antes de las elecciones presidenciales. Dentro de este tipo de encuestas se acostumbra a indagar acerca del favoritismo de los candidatos involucrados en la contienda electoral. Suponga que el candidato presidencial A está interesado en conocer su intención de voto previa a las elecciones. Para esto, él contrata a una firma encuestadora para la realización de una encuesta entre la población votante. El resultado de este estudio puede hacer cambiar o afirmar las estrategias publicitarias y la redefinición de la campaña electoral. La firma encuestadora decide implementar una estrategia de muestreo con un tamaño de muestra de doce mil personas. A cada respondiente se le realiza la siguiente pregunta: 
  
  > Si las elecciones presidenciales fueran mañana. ¿Usted votaría por el candidato A?
    
Las respuestas a esta pregunta son realizaciones de una muestra aleatoria de doce mil variables con densidad Bernoulli. Los resultados del estudio arrojan que 6360 personas de las personas entrevistadas, es decir un 53%, votarían por el suscrito candidato. Técnicamente se debe analizar esta cifra puesto que las implicaciones de ganar en una primera vuelta son grandes en el sentido económico, logístico y administrativo. Claramente, el dato 53% asegura una ventaja dentro de la muestra de doce mil personas. Sin embargo, es necesario realizar un estudio más profundo acerca de la caracterización estructural de la intención de voto del candidato en el electorado.
    
Con base en lo anteriormente expuesto, se decide utilizar la inferencia bayesiana puesto que existe información previa de un estudio anterior, contratado por el mismo candidato unos meses atrás en donde se entrevistaron a mil personas, con un favoritismo que estaba alrededor del 35 por ciento. Esta situación conlleva a la utilización de la metodología bayesiana que incorpora la información pasada acerca del mismo fenómeno.
    
El estadístico de la firma encuestadora decide utilizar una distribución previa Beta, definiendo los parámetros de la distribución previa como $\alpha$ igual al número de votantes a favor y $\beta$ igual al número de votantes en contra. Es decir, $Beta(\alpha=350, \beta=650)$. Por lo anterior, la distribución posterior del parámetro de interés, que representa la probabilidad de éxito en las elecciones presidenciales, es $Beta(6360+350, 650-6360+12000)=Beta(6710, 6290)$. Por lo tanto, utilizando la distribución posterior, se estima que la intención de voto por el candidato es de $\frac{6710}{6710+6290}=\frac{6710}{13000}=0.516$ y este valor equivale a la media de la distribución posterior. 
    
Sin embargo, si no se tuviese información previa como la suministrada por el estudio de meses anteriores, el análisis bayesiano sugeriría trabajar con una distribución previa no informativa, que en este caso, correspondería a una $Beta(\alpha=0.5, \beta=0.5)$. siguiendo el mismo análisis, se tiene que la distribución posterior es $Beta(6360.5, 5640.5)$. Finalmente, se estimaría que la intención de voto por el candidato es de $\frac{6350.5}{12001}=0.529$. 

La figuras \@ref(fig:BernoEj1) muestra el comportamiento de las distribuciones previas y posteriores en ambos escenarios. Nótese que la distribución no informativa influye muy poco en el comportamiento de la distribución posterior.</div>\EndKnitrBlock{example}


<div class="figure">
<img src="3Uniparametricos_files/figure-html/BernoEj1-1.png" alt="Distribuciones previas (línea punteada) y posteriores (línea sólida) para el ejemplo de las encuestas electorales." width="672" />
<p class="caption">(\#fig:BernoEj1)Distribuciones previas (línea punteada) y posteriores (línea sólida) para el ejemplo de las encuestas electorales.</p>
</div>

Utilizando el siguiente código en R, es posible conocer los intervalos de credibilidad para las dos distribuciones posteriores. Además, es posible concluir que, en ambos escenarios, el candidato aventaja significativamente a sus contrincantes y, salvo algún cambio drástico en el comportamiento del electorado, ganará las elecciones. Lo anterior se deduce puesto que el intervalo de credibilidad al 95 % no contiene ningún valor menor a 0.5


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

Por otro lado, el siguiente código en `STAN` permite obtener el mismo tipo de inferencia creando cuatro cadenas cuya distribución de probabilidad coincide con la distribución posterior del ejemplo.
    

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

```r
print(Berfit, digits = 2)
```

```
## Inference for Stan model: 62d7c91114fcc6227deb68f6059b2b09.
## 4 chains, each with iter=2000; warmup=1000; thin=1; 
## post-warmup draws per chain=1000, total post-warmup draws=4000.
## 
##           mean se_mean   sd     2.5%      25%      50%      75%    97.5% n_eff
## theta     0.52    0.00 0.00     0.51     0.51     0.52     0.52     0.52  1714
## lp__  -9005.24    0.01 0.68 -9007.14 -9005.39 -9004.97 -9004.81 -9004.76  2124
##       Rhat
## theta    1
## lp__     1
## 
## Samples were drawn using NUTS(diag_e) at Thu Jun  3 23:52:27 2021.
## For each parameter, n_eff is a crude measure of effective sample size,
## and Rhat is the potential scale reduction factor on split chains (at 
## convergence, Rhat=1).
```


    
    
    
    
    
    
