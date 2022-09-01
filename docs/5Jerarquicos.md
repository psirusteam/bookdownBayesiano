


# Modelos empíricos y jerárquicos

En las últimas décadas la formulación de modelos estadísticos ha evolucionado rápidamente, en parte, gracias a la capacidad de procesamiento de los equipos computacionales. En un principio, los modelos establecidos obedecían a reglas estándares que se suponían ciertas para toda la población. Sin embargo, el estado de la naturaleza de la mayoría de los problemas prácticos no sigue una regla común para todos y cada uno de los elementos de una población aleatoria. De hecho el sentido común establece que, para una misma población, pueden existir tendencias comunes entre diferentes miembros de la misma y la estructura de dispersión de los elementos puede obedecer comportamientos disímiles a través de éstos.

Lo anterior ha permitido que el investigador pueda proponer modelos que siguen comportamientos estructurales distintos y en algunos casos que se encuentran anidados en modelos más complejos. En el caso bayesiano, es claro que el momento de coyuntura en el cual el investigador no contempla un punto de retorno está dado en la formulación de la distribución previa para el vector de parámetros de interés $\btheta$. Más aún, la influencia de la distribución previa en la resultante distribución posterior está dada por la asignación del vector de hiperparámetros $\bEta$ que parametriza la distribución previa. Cuando los valores exactos de los hiperparámetros se desconocen o cuando no se tiene plena certeza del comportamiento estructural de la distribución previa, entonces es necesario estimarlos. En otras palabras, una mala asignación de los valores de los hiperparámetros conduce a una distribución previa que no es acorde con la realidad y esto puede conllevar a su vez a que la distribución posterior no sea apropiada, produciendo así resultados engañosos.

Siguiendo los fundamentos filosóficos de la estadística bayesiana, tener que estimar el vector de hiperparámetros envuelve al investigador en una paradoja cuya solución no siempre es sencilla. En primer lugar, nótese la forma de la distribución previa del vector de parámetros de interés: 

$$p(\btheta \mid \bEta)$$

A simple vista se puede concluir que $\bEta$ hace parte de la distribución previa, la cual, según la lógica de la filosofía bayesiana, involucra el conocimiento del investigador antes de la recolección de los datos. Por tanto la pregunta directa que surge es *¿por qué estimar algo que se supone conocido?*. En segundo lugar, y si se concibe tal estimación, la otra pregunta natural es *¿se deben utilizar los datos para estimar tales hiperparámetros?*. Las posibles respuestas a las anteriores preguntas definen caminos alternos en la estadística bayesiana.

1. Por un lado está la llamada *corriente bayesiana empírica* que utiliza los métodos de estimación puntual frecuentista para estimar los hiperparámetros y, por consiguiente, definir la distribución previa del vector de parámetros de interés. @Carlin96 menciona que en el análisis empírico se estima el vector de hiper-parámetros $\bEta$ con los datos observados, contradiciendo de alguna manera el espíritu y la filosofía de la corriente bayesiana radical y esta estimación se realiza con métodos frecuentistas ya sean paramétricos o no-paramétricos.
1. Por el otro lado se tiene la *corriente bayesiana jerárquica* que asume una posición totalmente bayesiana desde su concepción y establece un modelo posterior para los hiperparámetros.

En este capítulo se suponndrá que la variable de interés sigue un modelo común a toda la población aunque parametrizado por parámetros que toman distintos valores para cada individuo y que está regido por la siguiente expresión
\begin{equation*}
Y_i\sim p(Y_i \mid \theta_i)
\end{equation*}

## Análisis empírico

Este enfoque, criticado por algunos bayesianos radicales, se centra en la escogencia de una estimación $\hat{\bEta}$ para $\bEta$, obtenida como el valor que maximiza la verosimilitud marginal previa dada por

\begin{align}
(\#eq:ecua1)
p(\mathbf{Y} \mid \bEta)=\int p(\mathbf{Y} \mid \btheta)p(\btheta \mid \bEta) \ d\btheta
\end{align}

Por lo tanto todo el andamiaje inferencial está supeditado a la distribución posterior estimada, $p(\btheta \mid Y,\hat{\bEta})$. Una vez que esta está bien definida, el proceso de estimación puntual, estimación por intervalo y pruebas de hipótesis sigue su curso bayesiano idénticamente como en los capítulos anteriores.

En términos prácticos suponga que se tiene un modelo en dos etapas para cada una de las observaciones. Se asume que existen $n$ observaciones que, si bien no conforman una muestra aleatoria, conservan la característica de intercambiabilidad y están definidas en los siguientes términos
\begin{equation*}
Y_i \sim p(Y_i \mid \theta_i) \ \ \ \ \ \ \ i=1,\ldots,n
\end{equation*}

La segunda etapa comienza con la asignación de una distribución^[En esta etapa la distribución previa no está completamente especificada puesto que se desconocen los hiperparámetros que la indexan.] previa para los parámetros de interés $\theta_i$.
\begin{equation*}
\theta_i \sim p(\theta_i \mid \bEta) \ \ \ \ \ \ \ i=1,\ldots,n
\end{equation*}

Nótese que detrás de la asignación de la estructura probabilística para cada uno de los parámetros $\theta_i$, se supone que éstos últimos determinan una muestra aleatoria de la distribución $p(\btheta \mid \bEta)$. El objetivo de este enfoque es encontrar estimadores que maximicen la verosimilitud marginal previa la cual, para este caso particular y considerando independencia marginal entre las observaciones y el vector de hiperparámetros, es
\begin{align}
p(Y_i \mid \bEta)&=\int p(Y_i,\theta_i \mid \bEta) \ d\theta_i \notag \\
&=\int p(Y_i \mid \theta_i,\bEta)p(\theta_i \mid \bEta) \ d\theta_i \notag \\
&=\int p(Y_i \mid \theta_i)p(\theta_i \mid \bEta) \ d\theta_i
\end{align}

De lo anterior, la verosimilitud marginal previa del vector de observaciones dada por la expresión \@ref(eq:ecua1) queda convertida en
\begin{align}
p(Y \mid \bEta)&=\prod_{i=1}^np(Y_i \mid \bEta) \notag \\
&=\prod_{i=1}^n\int p(Y_i \mid \theta_i)p(\theta_i \mid \bEta) \ d\theta_i
\end{align}

A continuación, examinamos algunas casos que son de interés para los investigadores.  

### Modelo Binomial-Beta

Suponga el siguiente modelo binomial (intercambiable) en una primera etapa
\begin{equation*}
Y_i \mid \theta_i \sim Binomial(n_i,\theta_i)  \ \ \ \ \ \ \ i=1,\ldots,P.
\end{equation*}

Para la segunda etapa, se supone una muestra aleatoria proveniente de una misma distribución, tal que
\begin{equation*}
\theta_i \sim Beta(\alpha, \beta)  \ \ \ \ \ \ \ i=1,\ldots,p
\end{equation*}

Puesto que cada $\theta_i$ se encuentra en el intervalo $(0,1)$, sería apropiado asignarle una distribución Beta.

#### Análisis preliminar {-}

Es bien sabido que la distribución posterior para cada uno de los parámetros de interés involucrados en el anterior contexto está dada por
\begin{equation*}
\theta_i \mid Y_i \sim Beta(\alpha+Y_i, \beta+n_i-y_i)
\end{equation*}

para todo $i=1,\cdots,p$. Sin embargo, como se desconoce totalmente el valor de los hiperparámetros $\alpha$ y $\beta$, entonces se debe encontrar una estimación de estos para poder proseguir normalmente con la inferencia bayesiana, pero esta vez enfocados en la estimación de la distribución posterior dada por
\begin{equation*}
\theta_i \mid Y_i \sim Beta(\hat{\alpha}+Y_i, \hat{\beta}+n_i-y_i)
\end{equation*}

Para tal fin, nótese que la esperanza y la varianza previa de $\theta_i$ están dadas por las siguientes expresiones
\begin{align}
(\#eq:ecua3)
E(\theta_i)&=\frac{\alpha}{\alpha+\beta}\label{ecua2}\\
Var(\theta_i)&=\frac{\alpha\beta}{(\alpha+\beta)^2(\alpha+\beta+1)}
\end{align}

De donde se tiene que
\begin{align}
(\#eq:ecua5)
\alpha=E(\theta_i)(\alpha+\beta)
\end{align}

y también que
\begin{align}
(\#eq:ecua4)
1-E(\theta_i)=\frac{\beta}{\alpha+\beta}
\end{align}

Por lo tanto
\begin{align}
(\#eq:ecua6)
\beta=(1-E(\theta_i))(\alpha+\beta)
\end{align}

Reemplazando \@ref(eq:ecua5) y \@ref(eq:ecua4) en \@ref(eq:ecua3), se concluye que:
\begin{align*}
Var(\theta_i)&=\frac{E(\theta_i)(1-E(\theta_i))}{(\alpha+\beta+1)}
\end{align*}

Por consiguiente, 
\begin{align}
\alpha+\beta=\frac{E(\theta_i)(1-E(\theta_i))}{Var(\theta_i)}-1
\end{align}

Con el anterior razonamiento, es posible encontrar los estimadores de los hiperparámetros utilizando el método frecuentista de los momentos, los cuales corresponden a
\begin{align}
\widehat{\alpha+\beta}&=\frac{\bar{Y}(1-\bar{Y})}{S^2}-1
\end{align}

Donde $\bar{Y}$ y $S^2$ es el promedio y la varianza de las cantidades $Y_1/n_1, Y_2/n_2,\ldots, Y_P/n_P$, respectivamente. Ahora, teniendo en cuenta \@ref(eq:ecua5) y \@ref(eq:ecua6), se tiene que:
\begin{align}
\hat{\alpha}&=\left(\frac{\bar{Y}(1-\bar{Y})}{S^2}-1\right)\bar{Y}\\
\hat{\beta}&=\left(\frac{\bar{Y}(1-\bar{Y})}{S^2}-1\right)(1-\bar{Y})
\end{align}

Con las anteriores estimaciones es posible ahora conectarlas a la distribución posterior de $\theta_i$.

#### Análisis legítimo {-}

Según @Gelman03[, pág. 119], el anterior análisis no implica simplemente un punto de partida que da pie a la exploración de la idea de la estimación de los parámetros de la distribución posterior y, de ninguna manera, constituye un cálculo bayesiano, puesto que no está basado en ningún modelo de probabilidad. Sin embargo, el análisis empírico de esta situación, hace uso del esperanza y varianza condicional de la distribución beta de los parámetros $\theta_i$ ($i=1,\ldots, p$).

Para realizar este tipo de análisis, vamos a suponer que contamos con una variable $Y$, distribuida de forma binomial en $n$ ensayos y con probabilidad de éxito $\theta$. De esta manera, se tiene que el primer momento está dado por

\begin{align}
E_{binom}\left(\frac{Y}{n}\right)&=E_{beta}\left(E_{binom}\left(\frac{Y}{n} \mid \theta\right)\right) \notag \\
&=E_{beta}\left(\theta\right) \notag \\
(\#eq:ecua7)
&=\frac{\alpha}{\alpha+\beta}
\end{align}

Por otro lado, se tiene que la varianza, que es función del primer y segundo momento, está dada por

\begin{align*}
Var_{binom}\left(\frac{Y}{n}\right)
&=E_{beta}\left(Var_{binom}\left(\frac{Y}{n} \mid \theta\right)\right)
+ Var_{beta}\left(E_{binom}\left(\frac{Y}{n} \mid \theta\right)\right)  \\
&=E_{beta}\left( \frac{1}{n}\theta(1-\theta)\right)
+ Var_{beta}\left(\theta\right)  \\
&=\frac{1}{n}E_{beta}(\theta) - \frac{1}{n}E_{beta}(\theta^2)+ Var_{beta}\left(\theta\right)  \\
&=\frac{1}{n}E_{beta}(\theta) - \frac{1}{n}Var_{beta}(\theta)- \frac{1}{n}(E_{beta}\theta)^2+ Var_{beta}\left(\theta\right)  \\
&=\frac{n-1}{n}Var_{beta}(\theta) + \frac{1}{n}E_{beta}(\theta)(1-E_{beta}(\theta))  \\
&=\frac{n-1}{n}\frac{\alpha\beta}{(\alpha+\beta)^2(\alpha+\beta+1)} +
\frac{1}{n}\frac{\alpha\beta}{(\alpha+\beta)^2}   \\
&=\frac{1}{n}\frac{\alpha}{\alpha+\beta}\frac{\beta}{\alpha+\beta}
\left(\frac{n-1}{\alpha+\beta+1}+1\right)   \\
&=\frac{1}{n}E_{binom}\left(\frac{Y}{n}\right)\left(1-E_{binom}\left(\frac{Y}{n}\right)\right)
\left(\frac{n-1}{\alpha+\beta+1}+1\right)
\end{align*}

De esta última expresión, y despejando $\alpha +\beta$, se tiene que

\begin{align}
\alpha+\beta&=\frac{(n-1)E_{binom}\left(\frac{Y}{n}\right)\left(1-E_{binom}\left(\frac{Y}{n}\right)\right)}
{nVar_{binom}\left(\frac{Y}{n}\right)-
E_{binom}\left(\frac{Y}{n}\right)\left(1-E_{binom}\left(\frac{Y}{n}\right)\right)}-1\notag \\
&=\frac{E_{binom}\left(\frac{Y}{n}\right)\left(1-E_{binom}\left(\frac{Y}{n}\right)\right)-Var_{binom}\left(\frac{Y}{n}\right)}
{Var_{binom}\left(\frac{Y}{n}\right)-\frac{1}{n}E_{binom}\left(\frac{Y}{n}\right)\left(1-E_{binom}\left(\frac{Y}{n}\right)\right)}
\end{align}

Ahora, despejando $\alpha$ de la expresión \@ref(eq:ecua7) se tiene que
\begin{align}
\alpha=E_{binom}\left(\frac{Y}{n}\right)\frac{E_{binom}\left(\frac{Y}{n}\right)\left(1-E_{binom}\left(\frac{Y}{n}\right)\right)
-Var_{binom}\left(\frac{Y}{n}\right)}{Var_{binom}\left(\frac{Y}{n}\right)
-\frac{1}{n}E_{binom}\left(\frac{Y}{n}\right)\left(1-E_{binom}\left(\frac{Y}{n}\right)\right)}
\end{align}

Además, también despejando $\beta$ de \@ref(eq:ecua7) se tiene que
\begin{align}
\beta&=\frac{\alpha\left(1-E_{binom}\left(\frac{Y}{n}\right)\right)}{E_{binom}\left(\frac{Y}{n}\right)} \notag \\
&=\frac{E_{binom}\left(\frac{Y}{n}\right)\left(1-E_{binom}\left(\frac{Y}{n}\right)\right)-Var_{binom}\left(\frac{Y}{n}\right)}
{Var_{binom}\left(\frac{Y}{n}\right)-\frac{1}{n}E_{binom}\left(\frac{Y}{n}\right)\left(1-E_{binom}\left(\frac{Y}{n}\right)\right)}\left(1-E_{binom}\left(\frac{Y}{n}\right)\right)
\end{align}

El anterior enfoque nos ha llevado a poder expresar los parámetros de interés en términos de $E_{binom}\left(\frac{Y}{n}\right)$, $Var_{binom}\left(\frac{Y}{n}\right)$ y $n$. Una vez que podamos estimar las anteriores cantidades, es posible realizar la inferencia bayesiana empírica de la manera correcta. Para lo anterior, es necesario observar la naturaleza de las variables que, aunque no representan una muestra aleatoria, sí son una sucesión de variables aleatorias intercambiables. Por lo tanto, teniendo en cuenta que la inferencia se realiza con las cantidades $Y_1/n_1, Y_2/n_2, \ldots, Y_P/n_P$, es posible proponer los siguientes estimadores
\begin{align}
\widehat{E}_{binom}\left(\frac{Y}{n}\right)&=\bar{Y} \\
\widehat{Var}_{binom}\left(\frac{Y}{n}\right)&=S^2 \\
\widehat{n}&=\frac{1}{p}\sum_{i=1}^pn_i
\end{align}

Con base en lo anterior, las estimaciones empíricas de los parámetros $\alpha$ y $\beta$ son
\begin{align}
\hat{\alpha}=\bar{Y}\left(\frac{\bar{Y}\left(1-\bar{Y}\right)
-S^2}{S^2-\frac{1}{\hat{n}}\bar{Y}\left(1-\bar{Y}\right)}\right)
\end{align}

y
\begin{align}
\hat{\beta}&=(1-\bar{Y})\frac{\bar{Y}\left(1-\bar{Y}\right)-S^2}{S^2-\frac{1}{\hat{n}}\bar{Y}\left(1-\bar{Y}\right)}
\end{align}

respectivamente. Cuando la cantidad de ensayos $n_i$ es diferente en cada experimento, existen otras formas de obtener estimaciones para los parámetros $\alpha$ y $\beta$ [@Carlin96, pág. 81].

\BeginKnitrBlock{example}<div class="example"><span class="example" id="exm:BeisbolJerarquico"><strong>(\#exm:BeisbolJerarquico) </strong></span>En el ejemplo \@ref(exm:Beisbol) se estudiaron datos que correspondían al porcentaje de bateo de 18 jugadores profesionales de beisbol. A continuación se representa los cálculos necesarios para obtener una inferencia bayesiana empírica apropiada. </div>\EndKnitrBlock{example}


```r
library(pscl)
data(EfronMorris)
attach(EfronMorris)

y <- p # Porcentaje de bateo de los 18 jugadores
y
```

```
##  [1] 0.346 0.298 0.276 0.222 0.273 0.270 0.263 0.210 0.269 0.230 0.264 0.256
## [13] 0.303 0.264 0.226 0.285 0.316 0.200
```

```r
y.bar <- mean(y)
S2 <- var(y)
n.hat <- mean(n)
alfa <- y.bar * (y.bar * (1 - y.bar) - S2)/
  (S2 - y.bar * (1 - y.bar)/n.hat)
beta <- (1 - y.bar) * (y.bar * (1 - y.bar) - S2)/
  (S2 - y.bar * (1 - y.bar)/n.hat)
alfa
```

```
## [1] 56.99536
```

```r
beta
```

```
## [1] 158.0364
```

De donde podemos concluir que la distribución previa para cada $\theta_i$ es la distribución $Beta(57, 158)$, Nótese que la esperanza de esta distribución es 0.26, que coincide con el porcentaje de bateo promedio de los datos de los 18 jugadores. Ahora, podemos calcular los parámetros de la distribución para cada $\theta_i$ con $i=1,\cdots,18$ como sigue:


```r
alfa.new <- alfa + p * n
beta.new <- beta + (1 - p) * n
head(alfa.new)
```

```
## [1] 183.9774 183.9434 200.7914 118.0454 171.1094 182.8154
```

```r
head(beta.new)
```

```
## [1] 398.0544 457.0884 535.2404 371.9864 461.9224 498.2164
```
De esta manera, podemos realizar inferencias para cualquier $\theta_i$. Por ejemplo, para primer jugador, Roberto Clemente, la distribución posterior para su porcentaje de bateo es $Beta(184, 398)$; por consiguiente la estimación para el porcentaje de bateo de este jugador es $184/(184+398)=0.3162$ y su intervalo de credibilidad será $(0.279,0.354)$. 


























