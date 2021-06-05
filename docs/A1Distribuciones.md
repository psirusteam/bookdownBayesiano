# (APPENDIX) Apéndice {-} 

# Elementos de probabilidad

ss

## Distribuciones discretas

### Distribución uniforme discreta

\BeginKnitrBlock{definition}<div class="definition"><span class="definition" id="def:unnamed-chunk-1"><strong>(\#def:unnamed-chunk-1) </strong></span>Una variable aleatoria $Y$ tiene distribución uniforme discreta sobre el conjunto $\{1,2,\cdots,N\}$ si su función de densidad está dada por:
  
\begin{equation}
f_Y(y)=\frac{1}{N}I_{\{1,2,\cdots,N\}}(y)
\end{equation}</div>\EndKnitrBlock{definition}

Esta distribución describe situaciones donde los resultados de un experimento aleatorio tienen la misma probabilidad de ocurrencia. Entre los ejemplos de la distribución uniforme discreta en la vida práctica están el lanzamiento de una moneda corriente, el lanzamiento de un dado corriente, la extracción de una urna que contiene bolas enumeradas de 1 a $N$.

\BeginKnitrBlock{proposition}<div class="proposition"><span class="proposition" id="prp:unnamed-chunk-2"><strong>(\#prp:unnamed-chunk-2) </strong></span>Si $Y$ es una variable aleatoria con distribución uniforme discreta sobre el conjunto $\{1,2,\cdots,N\}$, entonces:

*  $E(Y)=\frac{N+1}{2}$.
*  $Var(Y)=\frac{N^2-1}{12}$.
*  $m_Y(t)=\sum_{i=1}^N\frac{e^{ti}}{N}$.
</div>\EndKnitrBlock{proposition}

### Distribución hipergeométrica

\BeginKnitrBlock{definition}<div class="definition"><span class="definition" id="def:unnamed-chunk-3"><strong>(\#def:unnamed-chunk-3) </strong></span>Una variable aleatoria $Y$ tiene distribución hipergeométrica con parámetros $n$, $R$ y $N$ si su función de densidad está dada por:

\begin{equation}
f_Y(y)=\frac{\binom{R}{y}\binom{N-R}{n-y}}{\binom{N}{n}}I_{\{0,1,\cdots,n\}}(y),
\end{equation}
y se nota como $Y\sim Hg(n,R,N)$.</div>\EndKnitrBlock{definition}
<br>

Suponga que en una urna hay $N$ bolas en total, donde $R$ de ellas son del color negro y los $N-R$ son del color blanco, se extrae aleatoriamente
$n$ bolas de la urna ($n<N$), entonces la variable "número de bolas negras extraídas" tiene distribución hipergeométrica con parámetros $n$, $R$ y
$N$. Otro uso de la distribución hipergeométrica es el problema de captura-recaptura.

\BeginKnitrBlock{proposition}<div class="proposition"><span class="proposition" id="prp:unnamed-chunk-4"><strong>(\#prp:unnamed-chunk-4) </strong></span>Si $Y$ es una variable aleatoria con distribución hipergeométrica con parámetros $n$, $R$ y $N$, entonces:

*  $E(Y)=\frac{nR}{N}$.
*  $Var(Y)=\frac{nR(N-R)(N-n)}{N^2(N-1)}$.
</div>\EndKnitrBlock{proposition}

El anterior resultado no incluye la función generadora de momentos, pues éste no ha resultado ser útil en la teoría relacionada con la distribución hipergeométrica.

### Distribución Bernoulli

La distribución Bernoulli debe su nombre al matemático suizo Jacob Bernoulli (1654-1705) que describe el éxito o fracaso de un evento.

\BeginKnitrBlock{definition}<div class="definition"><span class="definition" id="def:unnamed-chunk-5"><strong>(\#def:unnamed-chunk-5) </strong></span>Una variable aleatoria $Y$ tiene distribución Bernoulli con parámetro $p\in (0,1)$ si su función de densidad está dada por:
  
\begin{equation}
f_Y(y)=p^y(1-p)^{1-y}I_{\{0,1\}}(y),
\end{equation}

y se nota como $Y\sim Ber(p)$.</div>\EndKnitrBlock{definition}

\BeginKnitrBlock{proposition}<div class="proposition"><span class="proposition" id="prp:unnamed-chunk-6"><strong>(\#prp:unnamed-chunk-6) </strong></span>Si $Y$ es una variable aleatoria con distribución Bernoulli con parámetro $p$, entonces:
    
*  $E(Y)=p$.
*  $Var(Y)=p(1-p)$.
*  $m_Y(t)=pe^t+1-p$.
</div>\EndKnitrBlock{proposition}

### Distribución binomial

\BeginKnitrBlock{definition}<div class="definition"><span class="definition" id="def:unnamed-chunk-7"><strong>(\#def:unnamed-chunk-7) </strong></span>Una variable aleatoria $Y$ tiene distribución binomial con los parámetros $n\in \mathbb{N}$ y $p\in (0,1)$ si su función de densidad está dada por:
  
\begin{equation}
f_Y(y)=\binom{n}{y}p^y(1-p)^{n-y}I_{\{0,1,\cdots,n\}}(y),
\end{equation}

y se nota como $Y\sim Bin(n,p)$.</div>\EndKnitrBlock{definition}

\BeginKnitrBlock{proposition}<div class="proposition"><span class="proposition" id="prp:unnamed-chunk-8"><strong>(\#prp:unnamed-chunk-8) </strong></span>Sea $Y_1$, $\cdots$, $Y_n$ variables aleatorias independientes e idénticamente distribuidas con distribución Bernoulli con parámetro $p$,
entonces la variable $\sum_{i=1}^nY_i$ tiene distribución $Bin(n,p)$. Por ende, la distribución Bernoulli es un caso particular de la distribución binomial cuando $n=1$.</div>\EndKnitrBlock{proposition}
<br>

\BeginKnitrBlock{proof}<div class="proof">\iffalse{} <span class="proof"><em>Prueba. </em></span>  \fi{}La demostración radica en el hecho de que la función generadora de momentos caracteriza la distribución probabilística, entonces basta demostrar que la función generadora de momentos de $\sum_{i=1}^nX_i$ es la de una distribución $Bin(n,p)$. Tenemos lo siguiente:

\begin{align*}
m_{\sum Y_i}(t)=E(e^{\sum tY_i})&=E(\prod_{i=1}^ne^{tY_i})\\
               &=\prod_{i=1}^nE(e^{tY_i})\ \ \ \ (\text{por independencia})\\
               &=\prod_{i=1}^n(pe^t+1-p)\ \ \ \ (\text{definición de $m_{Y_i}(t)$})\\
               &=(pe^t+1-p)^n
\end{align*}</div>\EndKnitrBlock{proof}
<br>

Una aplicación de esta distribución es cuando tenemos un número $n$ de repeticiones independientes de un experimento donde cada uno tiene dos
posibles resultados que se podrían llamarse como éxito o fracaso y donde la probabilidad de éxito $p$ es constante en cada una de las repeticiones. Por tanto, la variable número de éxitos obtenidos en las $n$ repeticiones tiene distribución $Bin(n,p)$. La distribución binomial tiene dos parámetros, $n$ y $p$; sin embargo, cuando $n$ es conocido, la distribución dependerá sólo del valor $p$ que sería el único parámetro con espacio paramétrico $\Theta=(0,1)$.

\BeginKnitrBlock{proposition}<div class="proposition"><span class="proposition" id="prp:unnamed-chunk-10"><strong>(\#prp:unnamed-chunk-10) </strong></span>Si $Y$ es una variable aleatoria con distribución binomial con parámetros $n$ y $p$, entonces
    
*  $E(Y)=np$.
*  $Var(Y)=np(1-p)$.
*  $m_Y(t)=(pe^t+1-p)^n$.
</div>\EndKnitrBlock{proposition}

### Distribución Binomial negativa

\BeginKnitrBlock{definition}<div class="definition"><span class="definition" id="def:unnamed-chunk-11"><strong>(\#def:unnamed-chunk-11) </strong></span>Una variable aleatoria $Y$ tiene distribución Binomial negativa con parámetros $(\theta, r)$ si su función de densidad está dada por:

\begin{equation}
P(y\mid \theta, r)=\frac{\Gamma(r+y_i)}{y_i!\Gamma(r)}\theta^r(1-\theta)^{1-y_i}I_{(0,1,2,\ldots)}(y)
\end{equation}</div>\EndKnitrBlock{definition}

Esta distribución siempre ha tenido lugar al resolver el problema del número de ensayos necesarios para lograr un número específico de éxitos. Por supuesto, si $r$ es el número de éxitos necesarios y se conoce que la probabilidad de éxito es $\theta$, entonces la distribución binomial negativa corresponde a un modelo probabilístico, afianzado durante siglos, que permite la resolución de este tipo de situaciones. 

Por otro lado, es posible asignar al parámetro $r$ valores que sean reales; en este caso no hay ninguna interpretación práctica en el contexto del número de ensayos necesarios para determinados éxitos. Sin embargo, en términos de distribución, $r$ es un parámetro más. Esto nos lleva a uno de los verdaderos usos prácticos de esta distribución: la sobredispersión. Dado que la forma funcional de arriba corresponde a una generalización de la función de distribución Poisson, entonces es posible suponer que los datos de conteo vienen de una distribución binomial negativa. 

Lo anterior trae ventajas puesto que, si la media de los datos recolectados no corresponde con la varianza (característica esencial de la Poisson), entonces cualquier modelo que de allí surgiese sería altamente cuestionable. Si lo anterior se presenta es mejor acudir a la distribución binomial negativa dando valores reales al parámetro $r$.

\BeginKnitrBlock{proposition}<div class="proposition"><span class="proposition" id="prp:unnamed-chunk-12"><strong>(\#prp:unnamed-chunk-12) </strong></span>Si $Y$ es una variable aleatoria con distribución binomial-negativa con parámetros $(\theta, r)$, entonces
    
*  $E(Y)=\frac{r\theta}{1-\theta}$.
*  $Var(Y)=\frac{r\theta}{(1-\theta)^2}$.
*  $m_Y(t)=\left(\frac{1-\theta}{1-\theta e^t}\right)^r$.
</div>\EndKnitrBlock{proposition}

### Distribución de Poisson

La distribución de Poisson debe su nombre al francés Siméon-Denis Poisson (1781-1840) quien descubrió esta distribución en el año 1838, cuando
la usó para describir el número de ocurrencias de algún evento durante un intervalo de tiempo de longitud dada.

\BeginKnitrBlock{definition}<div class="definition"><span class="definition" id="def:unnamed-chunk-13"><strong>(\#def:unnamed-chunk-13) </strong></span>Una variable aleatoria $Y$ tiene distribución de Poisson con parámetros $\lambda>0$ si su función de densidad está dada por:
  
\begin{equation}
f_Y(y)=\frac{e^{-\lambda}\lambda^y}{y!}I_{\{0,1,\cdots\}}(y)
\end{equation}

y se nota como $Y\sim P(\lambda)$.</div>\EndKnitrBlock{definition}
<br>

Nótese que la distribución Poisson tiene solo un parámetro $\theta=\lambda$, y el espacio paramétrico es $\Theta=(0,\infty)$.

\BeginKnitrBlock{proposition}<div class="proposition"><span class="proposition" id="prp:unnamed-chunk-14"><strong>(\#prp:unnamed-chunk-14) </strong></span>Si $Y$ es una variable aleatoria con distribución Poisson con parámetro $\lambda$, entonces
    
*  $E(Y)=\lambda$.
*  $Var(Y)=\lambda$.
*  $m_Y(t)=\exp\{\lambda(e^t-1)\}$.
</div>\EndKnitrBlock{proposition}

\BeginKnitrBlock{proposition}<div class="proposition"><span class="proposition" id="prp:unnamed-chunk-15"><strong>(\#prp:unnamed-chunk-15) </strong></span>Sea $Y_1$, $\cdots$, $Y_n$ variables aleatorias independientes con distribución $P(\lambda_i)$ para $i=1,\cdots,n$, entonces la variable $\sum_{i=1}^nX_i$ tiene distribución $P(\sum_{i=1}^n\lambda_i)$.</div>\EndKnitrBlock{proposition}

## Distribuciones continuas

### Distribución Uniforme Continua

\BeginKnitrBlock{definition}<div class="definition"><span class="definition" id="def:unnamed-chunk-16"><strong>(\#def:unnamed-chunk-16) </strong></span>Una variable aleatoria $Y$ tiene distribución uniforme continua sobre el intervalo $[a,b]$ con $a<b$ si su función de densidad está dada por:
  
\begin{equation}
f_Y(y)=\frac{1}{b-a}I_{[a,b]}(y)
\end{equation}</div>\EndKnitrBlock{definition}

\BeginKnitrBlock{proposition}<div class="proposition"><span class="proposition" id="prp:unnamed-chunk-17"><strong>(\#prp:unnamed-chunk-17) </strong></span>Si $Y$ es una variable aleatoria con distribución uniforme continua sobre $[a,b]$, entonces
    
*  $E(Y)=\frac{a+b}{2}$.
*  $Var(Y)=\frac{(b-a)^2}{12}$.
*  $m_Y(t)=\frac{e^{bt}-e^{at}}{(b-a)t}$.
</div>\EndKnitrBlock{proposition}

### Distribución Weibull

\BeginKnitrBlock{definition}<div class="definition"><span class="definition" id="def:unnamed-chunk-18"><strong>(\#def:unnamed-chunk-18) </strong></span>Una variable aleatoria $Y$ tiene distribución uniforme continua sobre los reales positivos si su función de densidad está dada por:
  
\begin{equation}
p(Y\mid \theta, \gamma)=\frac{\theta}{\gamma^\theta}y^{\theta-1}exp\left\{-\frac{y^\theta}{\gamma^\theta}  \right\}I_{[0,\infty)}(y)
\end{equation}</div>\EndKnitrBlock{definition}

\BeginKnitrBlock{proposition}<div class="proposition"><span class="proposition" id="prp:unnamed-chunk-19"><strong>(\#prp:unnamed-chunk-19) </strong></span>Si $Y$ es una variable aleatoria con distribución Weibull, entonces
    
*  $E(Y)=\gamma \Gamma\left(1+\frac{1}{\theta}\right)$.
*  $Var(Y)=\gamma^2 \left[ \Gamma\left(1+\frac{2}{\theta}\right) + \Gamma^2\left(1+\frac{1}{\theta}\right)\right]$.
*  $m_Y(t)=\sum_{n=0}^\infty \frac{t^n\gamma^n}{n!}\Gamma\left(1+\frac{n}{\theta}\right), \ \theta\geq 1$.
</div>\EndKnitrBlock{proposition}

### Distribución valor-extremo

\BeginKnitrBlock{definition}<div class="definition"><span class="definition" id="def:unnamed-chunk-20"><strong>(\#def:unnamed-chunk-20) </strong></span>Una variable aleatoria $Y$ tiene distribución valor-extremo si su función de densidad está dada por:
  
\begin{equation}
p(y \mid \theta, \lambda )= \theta \exp(\theta y)\exp\left\{\lambda-\exp(\lambda+\theta y)\right\}
\end{equation}</div>\EndKnitrBlock{definition}

\BeginKnitrBlock{proposition}<div class="proposition"><span class="proposition" id="prp:unnamed-chunk-21"><strong>(\#prp:unnamed-chunk-21) </strong></span>Si $Y$ es una variable aleatoria con distribución valor-extremo, entonces
    
*  $E(Y)=-\frac{\lambda}{\theta}-\frac{\epsilon}{\theta}$.
*  $Var(Y)=\frac{\pi^2}{6\theta^2}$.
    
Donde $\pi\approx 3.1416$ es el número Pi y $\epsilon=0.5772$ es la constante de Euler.</div>\EndKnitrBlock{proposition}

### Distribución Gamma

\BeginKnitrBlock{definition}<div class="definition"><span class="definition" id="def:unnamed-chunk-22"><strong>(\#def:unnamed-chunk-22) </strong></span>Una variable aleatoria $Y$ tiene distribución Gamma con parámetro de forma $\alpha>0$ y parámetro de escala $\theta>0$ si su función de densidad está dada por:
  
\begin{equation}
p(\theta \mid \alpha,\beta)=\frac{\beta^\alpha}{\Gamma(\alpha)}\theta^{\alpha-1} e^{-\beta\theta}I_{(0,\infty)}(\theta).
\end{equation}

donde $\Gamma(k)=\int_0^{\infty}u^{k-1}\exp(-u)\ du$.</div>\EndKnitrBlock{definition}
<br>

La distribución Gamma tiene dos parámetros: $\alpha$ y $\beta$, en este caso, el vector de hiper-parámetros es $\btheta=(\alpha,\theta)'$ donde el espacio paramétrico está dado por $\Theta=(0,\infty)\times(0,\infty)$.

\BeginKnitrBlock{proposition}<div class="proposition"><span class="proposition" id="prp:unnamed-chunk-23"><strong>(\#prp:unnamed-chunk-23) </strong></span>Si $Y$ es una variable aleatoria con distribución Gamma con parámetro de forma $\alpha$ y parámetro de escala $\theta$, entonces
    
*  $E(Y)=\alpha/\beta$.
*  $Var(Y)=\alpha/\theta^2$.
</div>\EndKnitrBlock{proposition}

\BeginKnitrBlock{proposition}<div class="proposition"><span class="proposition" id="prp:unnamed-chunk-24"><strong>(\#prp:unnamed-chunk-24) </strong></span>Sea $Y_1$, $\cdots$, $Y_n$ variables aleatorias independientes con distribución Gamma con parámetro de forma $\alpha_i$ y parámetro de escala $\beta$ para $i=1,\cdots,n$, entonces la variable $\sum_{i=1}^nX_i$ tiene distribución Gamma con parámetro de forma $\sum_{i=1}^n\alpha_i$ y parámetro de escala $\theta$.</div>\EndKnitrBlock{proposition}

### Distribución Gamma-inversa

\BeginKnitrBlock{definition}<div class="definition"><span class="definition" id="def:unnamed-chunk-25"><strong>(\#def:unnamed-chunk-25) </strong></span>Una variable aleatoria $Y$ tiene distribución Gamma-inversa con parámetro de forma $\alpha>0$ y parámetro de escala $\beta>0$ si su función de densidad está dada por:
  
\begin{equation}
p(y \mid \alpha,\beta)=\frac{\beta^\alpha}{\Gamma(\alpha)}y^{-\alpha-1} e^{-\beta/y}I_{(0,\infty)}(y).
\end{equation}

donde $\Gamma(k)=\int_0^{\infty}u^{k-1}\exp(-u)\ du$.</div>\EndKnitrBlock{definition}
<br>

La distribución Gamma-inversa tiene dos parámetros: $\alpha$ y $\beta$; en este caso, el vector de hiper-parámetros es $\btheta=(\alpha,\beta)'$ donde el espacio paramétrico está dado por $\Theta=(0,\infty)\times(0,\infty)$.

\BeginKnitrBlock{proposition}<div class="proposition"><span class="proposition" id="prp:unnamed-chunk-26"><strong>(\#prp:unnamed-chunk-26) </strong></span>Si $Y$ es una variable aleatoria con distribución Gamma-inversa con parámetro de forma $\alpha$ y parámetro de escala $\beta$, entonces
    
*  $E(Y)=\beta/(\alpha-1)$.
*  $Var(Y)=\theta^2/(\alpha-1)^2(\alpha-2)$.
</div>\EndKnitrBlock{proposition}

\BeginKnitrBlock{proposition}<div class="proposition"><span class="proposition" id="prp:gammainver"><strong>(\#prp:gammainver) </strong></span>Si $X$ es una variable aleatoria con distribución $Gamma(\alpha,\beta)$, entonces $1/X$ tiene distribución $Gamma-inversa(\alpha,1/\beta)$.</div>\EndKnitrBlock{proposition}

### Distribución exponencial

\BeginKnitrBlock{definition}<div class="definition"><span class="definition" id="def:unnamed-chunk-27"><strong>(\#def:unnamed-chunk-27) </strong></span>Una variable aleatoria $Y$ tiene distribución exponencial con parámetro de escala $\theta>0$ si su función de densidad está dada por:
  
\begin{equation}
f_Y(y)=\frac{1}{\theta}e^{-y/\theta}I_{(0,\infty)}(y)
\end{equation}</div>\EndKnitrBlock{definition}

La distribución exponencial es un caso particular de la distribución Gamma cuando el parámetro de forma $k$ toma el valor 1, y usualmente se utiliza para describir la vida útil de un componente eléctrico o el tiempo necesario para la ocurrencia de algún evento.

\BeginKnitrBlock{proposition}<div class="proposition"><span class="proposition" id="prp:unnamed-chunk-28"><strong>(\#prp:unnamed-chunk-28) </strong></span>Si $Y$ es una variable aleatoria con distribución exponencial con parámetro $\theta$, entonces
    
*  $E(Y)=\theta$.
*  $Var(Y)=\theta^2$.
*  $m_Y(t)=\frac{1}{1-\theta t}$ para $t<1/\theta$, y no existe para otros valores de $t$.
</div>\EndKnitrBlock{proposition}

\BeginKnitrBlock{proposition}<div class="proposition"><span class="proposition" id="prp:unnamed-chunk-29"><strong>(\#prp:unnamed-chunk-29) </strong></span>Sea $Y_1$, $\cdots$, $Y_n$ variables aleatorias independientes e idénticamente distribuidas con distribución exponencial con parámetro de escala $\theta$, entonces la variable $\sum_{i=1}^nX_i$ tiene distribución Gamma con parámetro de forma $n$ y parámetro de escala $\theta$.</div>\EndKnitrBlock{proposition}

### Distribución Beta

\BeginKnitrBlock{definition}<div class="definition"><span class="definition" id="def:unnamed-chunk-30"><strong>(\#def:unnamed-chunk-30) </strong></span>Una variable aleatoria $Y$ tiene distribución Beta con parámetro de forma $\alpha>0$ y parámetro de escala $\beta>0$ si su función de densidad está dada por:
  
\begin{equation}
f_Y(y)=\frac{1}{Beta(\alpha,\beta)}y^{\alpha-1}(1-y)^{\beta-1}I_{[0,1]}(y).
\end{equation}

donde $Beta(\alpha,\beta)=\dfrac{\gamma(\alpha)\gamma(\beta)}{\gamma(\alpha+\beta)}$.</div>\EndKnitrBlock{definition}
<br>

La distribución Beta tiene dos parámetros: $\alpha$ y $\beta$; en este caso, el vector de parámetros es $\btheta=(\alpha,\beta)'$ donde el espacio paramétrico está dado por $\Theta=(0,\infty)\times(0,\infty)$. Pero cuando uno de los dos parámetros es fijo, por ejemplo $\theta$, entonces la distribución tendría un sólo parámetro: $k$.

\BeginKnitrBlock{proposition}<div class="proposition"><span class="proposition" id="prp:unnamed-chunk-31"><strong>(\#prp:unnamed-chunk-31) </strong></span>Si $Y$ es una variable aleatoria con distribución Gamma con parámetro de forma $k$ y parámetro de escala $\theta$, entonces
    
*  $E(Y)=\frac{\alpha}{\alpha + \beta}$.
*  $Var(Y)=\frac{\alpha\beta}{(\alpha+\beta)^2(\alpha+\beta+1)}$
  </div>\EndKnitrBlock{proposition}

### Distribución normal 

La distribución normal también es llamada la distribución gaussiana, rindiendo homenaje al matemático alemán Carl Friedrich Gauss (1777-1855). La distribución normal es, sin duda, una de las distribuciones más importantes, puesto que una gran parte de la teoría estadística fue desarrollada inicialmente para variables con esta distribución; por el otro lado, gracias al teorema central del límite, muchas distribuciones ajenas a la normal puede ser aproximadas por esta.

\BeginKnitrBlock{definition}<div class="definition"><span class="definition" id="def:unnamed-chunk-32"><strong>(\#def:unnamed-chunk-32) </strong></span>Una variable aleatoria $Y$ tiene distribución normal con parámetros $\mu$ y $\sigma^2$ si su función de densidad está dada por:
  
\begin{equation}
f_Y(y)=\frac{1}{\sqrt{2\pi\sigma^2}}\exp\left\{-\frac{1}{2\sigma^2}(y-\mu)^2\right\}I_\mathbb{R}(y),
\end{equation}

donde $\sigma>0$ y se nota como $Y\sim N(\mu,\sigma^2)$.</div>\EndKnitrBlock{definition}
<br>

La distribución normal tiene dos parámetros, representado como $\btheta=(\mu,\sigma^2)$, mientras que su espacio paramétrico es  $\Theta=\mathbb{R}\times(0,\infty)$.

\BeginKnitrBlock{proposition}<div class="proposition"><span class="proposition" id="prp:unnamed-chunk-33"><strong>(\#prp:unnamed-chunk-33) </strong></span>Si $Y$ es una variable aleatoria con distribución normal con parámetros $\mu$ y $\sigma^2$, entonces
    
*  $E(Y)=\mu$.
*  $Var(Y)=\sigma^2$.
*  $m_Y(t)=\exp\{\mu t+\frac{1}{2}\sigma^2t^2\}$.
</div>\EndKnitrBlock{proposition}

Cuando $\mu=0$ y $\sigma=1$, se dice que $Y$ tiene distribución normal estándar y usualmente se denota por $Z$.

\BeginKnitrBlock{proposition}<div class="proposition"><span class="proposition" id="prp:unnamed-chunk-34"><strong>(\#prp:unnamed-chunk-34) </strong></span>Si $Y\sim N(\mu,\sigma^2)$, y $\alpha$, $\beta$ son constantes, entonces la variable $\alpha Y+\beta$ tiene distribución $N(\alpha\mu+\beta,\alpha^2\sigma^2)$.</div>\EndKnitrBlock{proposition}
<br>

\BeginKnitrBlock{proof}<div class="proof">\iffalse{} <span class="proof"><em>Prueba. </em></span>  \fi{}Se usará el hecho de que la función generadora de momentos caracteriza la distribución probabilística. Se tiene que:

\begin{align*}
m_{\alpha Y+\beta}(t)&=E(e^{t(\alpha Y+\beta)})\\
                     &=E(e^{\alpha tY})e^{\beta t}\\
                     &=m_Y(\alpha t)e^{\beta t}\\
                     &=e^{\mu\alpha t+\sigma^2\alpha^2t/2}e^{\beta t}\\
                     &=e^{(\alpha\mu+\beta)t+\sigma^2\alpha^2t/2}
\end{align*}

la cual es la función generadora de momentos de una distribución $N(\alpha\mu+\beta,\alpha^2\sigma^2)$, y el resultado queda demostrado.</div>\EndKnitrBlock{proof}
<br>

Como consecuencia inmediata del anterior resultado, se define la estandarización, que es fundamental en la teoría relacionada con las distribuciones normales. Si $Y\sim N(\mu,\sigma^2)$, entonces la variable $Z=\frac{Y-\mu}{\sigma}$ tiene distribución normal estándar, y la anterior transformación se conoce como la normal estandarizada.

\BeginKnitrBlock{proposition}<div class="proposition"><span class="proposition" id="prp:unnamed-chunk-36"><strong>(\#prp:unnamed-chunk-36) </strong></span>Sea $Y_1$, $\cdots$, $Y_n$ variables aleatorias independientes, donde $Y_i\sim N(\mu_i,\sigma^2_i)$ con $i=1,\cdots,n$, entonces la variable $\sum_{i=1}^nY_i$ tiene distribución $N(\sum_{i=1}^n\mu_i,\sum_{i=1}^n\sigma_i^2)$.</div>\EndKnitrBlock{proposition}

### Distribución log-normal

\BeginKnitrBlock{definition}<div class="definition"><span class="definition" id="def:unnamed-chunk-37"><strong>(\#def:unnamed-chunk-37) </strong></span>Una variable aleatoria $Y$ tiene distribución log-normal si su función de densidad está dada por:
  
\begin{equation}
p(Y\mid \mu, \sigma^2)=\frac{1}{y\sqrt{2\pi\sigma^2}}\exp\{\frac{-1}{2\sigma^2(\ln(y)-\mu)^2}\}
\end{equation}

Nótese que si $\mu$ y $\sigma^2$ son la media y la varianza de $\ln(Y)$, entonces $\ln(Y)$ tiene distribución normal de media $\mu$ y varianza $\sigma^2$.</div>\EndKnitrBlock{definition}

\BeginKnitrBlock{proposition}<div class="proposition"><span class="proposition" id="prp:unnamed-chunk-38"><strong>(\#prp:unnamed-chunk-38) </strong></span>Si $Y$ es una variable aleatoria con distribución log-normal, entonces
    
*  $E(Y)=\exp(\mu+\sigma^2/2)$.
*  $Var(Y)=(\exp(\sigma^2-1)) \exp(2\mu+\sigma^2)$.
</div>\EndKnitrBlock{proposition}

### Distribución Ji-cuadrado

\BeginKnitrBlock{definition}<div class="definition"><span class="definition" id="def:unnamed-chunk-39"><strong>(\#def:unnamed-chunk-39) </strong></span>Una variable aleatoria $Y$ tiene distribución Ji-cuadrado con $n$ grados de libertad, con $n$ entero positivo, si su función de densidad está dada por:
  
\begin{equation}
f_Y(y)=\frac{y^{(n/2)-1}e^{-y/2}}{2^{n/2}\Gamma(n/2)}I_{(0,\infty)}(y),
\end{equation}

y se nota como $Y\sim\chi^2_n$.</div>\EndKnitrBlock{definition}
<br>

La distribución Ji-cuadrado con $n$ grados de libertad es un caso particular de la distribución Gamma cuando el parámetro de forma $k$ toma el valor $n/2$ y el parámetro de escala toma el valor 2. También, en la literatura estadística existe la siguiente definición para la distribución Ji-cuadrado.

\BeginKnitrBlock{definition}<div class="definition"><span class="definition" id="def:unnamed-chunk-40"><strong>(\#def:unnamed-chunk-40) </strong></span>Si $Z_1$, $\cdots$, $Z_n$ son variables aleatorias independientes e idénticamente distribuidas con distribución normal estándar, entonces la variable $\sum_{i=1}^nZ_i^2$ tiene distribución Ji-cuadrado con $n$ grados de libertad.</div>\EndKnitrBlock{definition}

\BeginKnitrBlock{proposition}<div class="proposition"><span class="proposition" id="prp:unnamed-chunk-41"><strong>(\#prp:unnamed-chunk-41) </strong></span>Si $Y$ es una variable aleatoria con distribución Ji-cuadro con $n$ grados de libertad, entonces
    
*  $E(Y)=n$.
*  $Var(Y)=2n$.
*  $m_Y(t)=\left(\frac{1}{1-2t}\right)^{n/2}$ para $t<1/2$, y no existe para otros valores de $t$.
</div>\EndKnitrBlock{proposition}

\BeginKnitrBlock{proposition}<div class="proposition"><span class="proposition" id="prp:unnamed-chunk-42"><strong>(\#prp:unnamed-chunk-42) </strong></span>Sea $Z_1$, $\cdots$, $Z_m$ variables aleatorias independientes con distribución $\chi^2_{n_i}$ para $i=1,\cdots,m$, entonces la variable $\sum_{i=1}^mZ_i$ tiene distribución Ji-cuadrado con $\sum_{i=1}^mn_i$ grados de libertad.</div>\EndKnitrBlock{proposition}

### Distribución t-student

El descubrimiento de la distribución t-student fue publicado por el estadístico inglés William Sealy Gosset (1876-1937) en el año 1908 cuando trabajaba en la famosa empresa cervecera *Guinness*. La publicación lo hizo de forma anónimo bajo el nombre de Student, pues Guinness le prohibía la publicación por ser el descubrimiento parte de resultados de investigación realizado por la empresa.

\BeginKnitrBlock{definition}<div class="definition"><span class="definition" id="def:unnamed-chunk-43"><strong>(\#def:unnamed-chunk-43) </strong></span>Una variable aleatoria $Y$ tiene distribución t-student con $n$ grados de libertad si su función de densidad está dada por:
  
\begin{equation}
f_Y(y)=\frac{\Gamma(\frac{n+1}{2})}{\sqrt{\pi n}\ \Gamma(\frac{n}{2})}\left(1+\frac{y^2}{n}\right)^{-(n+1)/2}I_\mathbb{R}(y),
\end{equation}

donde $n>0$ y se nota como $Y\sim t_n$.</div>\EndKnitrBlock{definition}
<br>

Otra definición que se encuentra frecuentemente en la literatura estadística es la siguiente.

\BeginKnitrBlock{definition}<div class="definition"><span class="definition" id="def:unnamed-chunk-44"><strong>(\#def:unnamed-chunk-44) </strong></span>Sea $Z$ una variable aleatoria con distribución normal estándar y $Y$ una variable aleatoria con distribución Ji-cuadrado con $n$ grados de libertad, si $Z$ y $Y$ son independientes, entonces la variable $\frac{Z}{\sqrt{Y/n}}$ tiene distribución t-student con $n$ grados de libertad.</div>\EndKnitrBlock{definition}
<br>

La función de densidad de la distribución t-student es muy parecida a la de distribución normal estándar, entre más grande sea el grado de libertad, más se parece a la distribución normal estándar.

\BeginKnitrBlock{proposition}<div class="proposition"><span class="proposition" id="prp:unnamed-chunk-45"><strong>(\#prp:unnamed-chunk-45) </strong></span>Si $Y$ es una variable aleatoria con distribución t-student con $n$ grados de libertad, entonces
    
*  $E(Y)=0$ para $n>1$.
*  $Var(Y)=\frac{n}{n-2}$ para $n>2$.

La distribución t-student no tiene función generadora de momentos.</div>\EndKnitrBlock{proposition}

### Distribución t-student generalizada

\BeginKnitrBlock{definition}<div class="definition"><span class="definition" id="def:unnamed-chunk-46"><strong>(\#def:unnamed-chunk-46) </strong></span>Una variable aleatoria $Y$ tiene distribución t-student con $n$ grados de libertad, parámetro de centralidad $\theta$ y parámetro de escala $\sigma^2$, si su función de densidad está dada por:

\begin{equation}
f_Y(y)=\frac{\Gamma((n+1)/2)}{\Gamma(n/2)\sqrt{n\pi}\sigma}\left[1+\frac{1}{n}\left(\frac{y-\theta}{\sigma}\right)^2\right]^{-(n+1)/2}
I_\mathbb{R}(y),
\end{equation}

donde $n>0$ y se nota como $Y\sim t_n(\theta,\sigma^2)$.</div>\EndKnitrBlock{definition}

\BeginKnitrBlock{proposition}<div class="proposition"><span class="proposition" id="prp:unnamed-chunk-47"><strong>(\#prp:unnamed-chunk-47) </strong></span>Si $Y$ es una variable aleatoria con distribución t-student generalizada, entonces
    
*  $E(Y)=\theta$ para $n>1$.
*  $Var(Y)=\frac{n}{n-2}\sigma^2$ para $n>2$.
</div>\EndKnitrBlock{proposition}

### Distribución F

La distribución F también se conoce como la distribución F de Fisher o distribución de Fisher-Snedecor, refiriendo al gran estadístico Ronald Aylmer Fisher (1890-1962) y el fundador del primer departamento de estadística en los Estados Unidos, George Waddel Snedecor (1881-1974).

\BeginKnitrBlock{definition}<div class="definition"><span class="definition" id="def:unnamed-chunk-48"><strong>(\#def:unnamed-chunk-48) </strong></span>Una variable aleatoria $Y$ tiene distribución F con $m$ grados de libertad en el numerador y $n$ grados de libertado en el denominador si su función de densidad está dada por:
  
\begin{equation}
f_Y(y)=\frac{\Gamma(\frac{m+n}{2})}{\Gamma(\frac{m}{2})\Gamma(\frac{n}{2})}\left(\frac{m}{n}\right)^{m/2}\frac{z^{\frac{m}{2}-1}}{\left(1+\frac{m}{n}z\right)^{\frac{m+n}{2}}},
\end{equation}

y se nota como $Y\sim F^m_n$.</div>\EndKnitrBlock{definition}
<br>

Otra definición para la distribución F es como sigue:

\BeginKnitrBlock{definition}<div class="definition"><span class="definition" id="def:unnamed-chunk-49"><strong>(\#def:unnamed-chunk-49) </strong></span>Sea $Y$ y $Y$ variables aleatorias independientes con distribuciones Ji-cuadrado con $m$ y $n$ grados de libertad, respectivamente, entonces la variable $\frac{Y/m}{Y/n}$ tiene distribución F con $m$ grados de libertad en el numerador y $n$ grados de libertado en el denominador.</div>\EndKnitrBlock{definition}

\BeginKnitrBlock{proposition}<div class="proposition"><span class="proposition" id="prp:unnamed-chunk-50"><strong>(\#prp:unnamed-chunk-50) </strong></span>Si $Y$ es una variable aleatoria con distribución F con $m$ grados de libertad en el numerador y $n$ grados de libertado en el denominador, entonces
    
*  $E(Y)=\frac{n}{n-2}$ para $n>2$.
*  $Var(Y)=\frac{2n^2(m+n-2)}{m(n-2)^2(n-4)}$ para $n>4$.

La distribución F no tiene función generadora de momentos.</div>\EndKnitrBlock{proposition}

## Distribuciones multivariadas

### Distribución Multinomial

\BeginKnitrBlock{definition}<div class="definition"><span class="definition" id="def:unnamed-chunk-51"><strong>(\#def:unnamed-chunk-51) </strong></span>Un vector aleatorio $\mathbf{Y}=(Y_1,\ldots,Y_p')$ tiene distribución multinomial si su función de densidad está dada por:

\begin{equation}
p(\mathbf{Y} \mid \btheta)=\binom{n}{y_1,\ldots,y_p}\theta_1^{y_1}\cdots\theta_p^{y_p} \ \ \ \ \ \theta_i>0 \texttt{, } \sum_{i=1}^p\theta_i=1  \texttt{ y } \sum_{i=1}^py_i=p
\end{equation}

donde

\begin{equation}
\binom{p}{y_1,\ldots,y_p}=\frac{p!}{y_1!\cdots y_p!}.
\end{equation}</div>\EndKnitrBlock{definition}

Como @Gelman03, afirma esta distribución es una generalización de la distribución binomial. La distribución marginal de una sola variable $Y_i$ es $Binomial(p,\theta_i)$

\BeginKnitrBlock{proposition}<div class="proposition"><span class="proposition" id="prp:unnamed-chunk-52"><strong>(\#prp:unnamed-chunk-52) </strong></span>Si $\mathbf{Y}$ es una vector aleatorio con distribución multinomial, entonces
    
*  $E(\mathbf{Y})=p(\theta_1,\ldots,\theta_p)'$.
*  $Var(\mathbf{Y})_{ij}=
    \begin{cases}
    p\theta_i(1-\theta_i) & \text{si $i=j$}\\
    -p\theta_i\theta_j & \text{si $i\neq j$}
    \end{cases}$
</div>\EndKnitrBlock{proposition}

### Distribución Dirichelt

\BeginKnitrBlock{definition}<div class="definition"><span class="definition" id="def:unnamed-chunk-53"><strong>(\#def:unnamed-chunk-53) </strong></span>Un vector aleatorio $\mathbf{Y}=(Y_1,\ldots,Y_p')$ tiene distribución Dirichelt si su función de densidad está dada por:

\begin{equation}
p(\mathbf{Y} \mid \btheta)=\frac{\Gamma(\theta_1+\cdots+\theta_p)}{\Gamma(\theta_1)\cdots\Gamma(\theta_p)}
y^{\theta_1-1}\cdots y^{\theta_p-1} \ \ \ \ \ \theta_i>0 \texttt{ y } \sum_{i=1}^p\theta_i=1.
\end{equation}</div>\EndKnitrBlock{definition}

Esta distribución es una generalización de la distribución beta. La distribución marginal de una sola variable $Y_i$ es $Beta(\theta_i,(\sum_{i=1}^p\theta_i)-\theta_i)$

\BeginKnitrBlock{proposition}<div class="proposition"><span class="proposition" id="prp:unnamed-chunk-54"><strong>(\#prp:unnamed-chunk-54) </strong></span>Si $\mathbf{Y}$ es una vector aleatorio con distribución Dirichlet, entonces
    
*  $E(\mathbf{Y})=(\sum_{i=1}^p\theta_i)^{-1}(\theta_1,\ldots,\theta_p)'$.
*  $Var(\mathbf{Y})_{ij}=
   \begin{cases}
   \frac{\theta_i(\sum_{i=1}^p\theta_i-\theta_i)}{(\sum_{i=1}^p\theta_i)^2(\sum_{i=1}^p\theta_i+1)} & \text{si $i=j$}\\
   -\frac{\theta_i\theta_j)}{(\sum_{i=1}^p\theta_i)^2(\sum_{i=1}^p\theta_i+1)} & \text{si $i\neq j$}
    \end{cases}$
    </div>\EndKnitrBlock{proposition}

### Distribución Normal Multivariante

\BeginKnitrBlock{definition}<div class="definition"><span class="definition" id="def:unnamed-chunk-55"><strong>(\#def:unnamed-chunk-55) </strong></span>Un vector aleatorio $\mathbf{Y}=(Y_1,\ldots,Y_p')$ tiene distribución normal multivariante de orden $p$, denotada como $\mathbf{Y}\sim N_p(\btheta,\bSigma)$, si su función de densidad está dada por:

\begin{equation}
p(\mathbf{Y} \mid \btheta,\bSigma)=(2\pi)^{-p/2} \mid \bSigma \mid ^{-1/2}
\exp\left\{-\frac{1}{2}(\mathbf{y}-\btheta)'\bSigma(\mathbf{y}-\btheta)\right\}
\end{equation}

donde $\mid \bSigma \mid$ se refiere al determinante de la matriz $\bSigma$, la cual es simétrica y definida positiva de orden $p\times p$.</div>\EndKnitrBlock{definition}
<br>

La distribución Normal Multivariante es el baluarte de una gran cantidad de técnicas y métodos estadísticos como son los modelos lineales, los modelos lineales generalizados, el análisis factorial, etc. Algunas de sus propiedades se citan a continuación.

\BeginKnitrBlock{proposition}<div class="proposition"><span class="proposition" id="prp:normalmulti"><strong>(\#prp:normalmulti) </strong></span>Si $\mathbf{Y}=(Y_1,\ldots,Y_p')$ es una vector aleatorio con distribución normal multivariante, entonces

*  La distribución marginal de cualquier subconjunto de componentes de $\mathbf{Y}$ es también normal multivariante. Por ejemplo si $\mathbf{Y}$ es particionado en $\mathbf{Y}=(\mathbf{Y}_1',\mathbf{Y}_2')$, entonces $p(\mathbf{Y}_1)$ seguiría una distribución normal multivariante, al igual que $p(\mathbf{Y}_2)$.
*  Cualquier transformación lineal de $\mathbf{Y}$ es normal multivariante y su dimensión equivale al rango de la transformación. en particular, la suma de las componentes del vector, dada por $\sum_{i=1}^pY_i$ sigue una distribución normal univariada.
*  La distribución condicional de $\mathbf{Y}$, restringida a un subespacio lineal es normal.
*  La distribución condicional de cualquier sub-vector de elementos de $\mathbf{Y}$ dados los restantes elementos es normal multivariante. Más aún, si $\mathbf{Y}$ es particionado en $\mathbf{Y}=(\mathbf{Y}_1',\mathbf{Y}_2')$, entonces $p(\mathbf{Y}_1 \mid \mathbf{Y}_2)$ es normal multivariada con
\begin{align*}
E(\mathbf{Y}_1 \mid \mathbf{Y}_2)&=E(\mathbf{Y}_1)+Cov(\mathbf{Y}_1,\mathbf{Y}_2)(Var(\mathbf{Y}_2))^{-1}(\mathbf{Y}_2-E(\mathbf{Y}_2))\\
Var(\mathbf{Y}_1 \mid \mathbf{Y}_2)&=Var(\mathbf{Y}_1)-Cov(\mathbf{Y}_1,\mathbf{Y}_2)(Var(\mathbf{Y}_2))^{-1}Cov(\mathbf{Y}_2,\mathbf{Y}_1)
\end{align*}
*  Si $\mathbf{X}$ es un vector con distribución normal multivariante, entonces $\mathbf{X}+\mathbf{Y}$ tiene una distribución normal multivariante. En particular si $\mathbf{X}$ es independiente de $\mathbf{Y}$, comparten el mismo orden $p$ y $\mathbf{X}\sim N_p(\bmu,\bGamma)$, entonces $\mathbf{X}+\mathbf{Y}\sim N_p(\bmu+\btheta,\bGamma+\bSigma)$.
</div>\EndKnitrBlock{proposition}

\BeginKnitrBlock{proposition}<div class="proposition"><span class="proposition" id="prp:unnamed-chunk-56"><strong>(\#prp:unnamed-chunk-56) </strong></span>Si $\mathbf{Y}$ es una vector aleatorio con distribución Normal Multivariante, entonces
    
*  $E(\mathbf{Y})=\btheta=(\theta_1,\ldots,\theta_n)'$.
*  $Var(\mathbf{Y})=\bSigma$
</div>\EndKnitrBlock{proposition}

\BeginKnitrBlock{proposition}<div class="proposition"><span class="proposition" id="prp:unnamed-chunk-57"><strong>(\#prp:unnamed-chunk-57) </strong></span>Dado $\mathbf{Y}$ un vector aleatorio particionado como $\mathbf{Y}=(\mathbf{Y}_1',\mathbf{Y}_2')$ con esperanza $\btheta=(\btheta_1',\btheta_2')$ y matrix de varianzas y covarianzas

\begin{equation*}
\bSigma=\begin{pmatrix}
\bSigma_{11}&\bSigma_{12}\\
\bSigma_{21}&\bSigma_{22}
\end{pmatrix}.
\end{equation*}

Si $\mathbf{Y}_1 \mid \mathbf{Y}_2\sim N(\btheta_1+\bSigma_{12}\bSigma_{22}^{-1}(\mathbf{Y}_2-\btheta_2),\bSigma_{11}-\bSigma_{12}\bSigma_{22}^{-1}\bSigma_{21})$ y $\mathbf{Y}_2\sim N(\btheta_2,\bSigma_{22})$, entonces se tiene que

\begin{equation*}
\mathbf{Y}\sim N(\btheta,\bSigma).
\end{equation*}</div>\EndKnitrBlock{proposition}

\BeginKnitrBlock{proposition}<div class="proposition"><span class="proposition" id="prp:unnamed-chunk-58"><strong>(\#prp:unnamed-chunk-58) </strong></span>Si $\mathbf{Y}_1,\ldots,\mathbf{Y}_n$ es una muestra aleatoria de vectores con distribución Normal Multivariante, entonces la verosimilitud de la muestra se puede escribir como

\begin{equation}
\prod_{i=1}^np(\mathbf{Y}_i \mid \btheta,\bSigma)\propto \mid \bSigma \mid ^{-n/2}\exp\left\{-\frac{1}{2}traza(\bSigma^{-1}\mathbf{S}_{\btheta})\right\}
\end{equation}

Donde $\mathbf{S}_{\btheta}=\sum_{i=1}^n(\mathbf{Y}_i-\btheta)(\mathbf{Y}_i-\btheta)'$.</div>\EndKnitrBlock{proposition}
<br>

\BeginKnitrBlock{proof}<div class="proof">\iffalse{} <span class="proof"><em>Prueba. </em></span>  \fi{}La verosimilitud de la muestra aleatoria está dada por

\begin{align*}
\prod_{i=1}^np(\mathbf{Y}_i \mid \btheta,\bSigma)
&\propto \mid \bSigma \mid ^{-n/2}\exp\left\{-\frac{1}{2}\sum_{i=1}^n
(\mathbf{Y}_i-\btheta)'\bSigma^{-1}(\mathbf{Y}_i-\btheta)\right\}\\
&= \mid \bSigma \mid ^{-n/2}\exp\left\{-\frac{1}{2}traza\left(\bSigma^{-1}\mathbf{S}_{\btheta}\right)\right\}
\end{align*}

Puesto que, por las propiedades del operador $traza$, se tiene que

*  Si $c$ es un escalar, entonces $c=traza(c)$.
*  Si $\mathbf{A}$ y $\mathbf{B}$ son dos matrices, entonces $traza(\mathbf{AB})=traza(\mathbf{BA})$
*  Si $\mathbf{A}_i$ (i=1,\ldots,n) son matrices del mismo tamaño, entonces $\sum_{i=1}^ntraza(\mathbf{A}_i)=traza\left(\sum_{i=1}^n\mathbf{A}_i\right)$

Por lo anterior,

\begin{align*}
\sum_{i=1}^n(\mathbf{Y}_i-\btheta)'\bSigma^{-1}(\mathbf{Y}_i-\btheta)
&=traza\left[\sum_{i=1}^n(\mathbf{Y}_i-\btheta)'\bSigma^{-1}(\mathbf{Y}_i-\btheta)\right]\\
&=\sum_{i=1}^ntraza[\bSigma^{-1}(\mathbf{Y}_i-\btheta)(\mathbf{Y}_i-\btheta)')]\\
&=traza\left[\bSigma^{-1}\sum_{i=1}^n(\mathbf{Y}_i-\btheta)(\mathbf{Y}_i-\btheta)')\right]\\
&=traza(\bSigma^{-1}\mathbf{S}_{\btheta})
\end{align*}</div>\EndKnitrBlock{proof}

### Distribución Wishart

\BeginKnitrBlock{definition}<div class="definition"><span class="definition" id="def:unnamed-chunk-60"><strong>(\#def:unnamed-chunk-60) </strong></span>Sea $\bSigma$ una matriz aleatoria simétrica y definida positiva de tamaño $p\times p$. Se dice que $\bSigma$ tiene distribución Wishart con $v$ grados de libertad, denotada como $\mathbf{Y}\sim Wishart_v(\bLambda)$, si su función de densidad está dada por:
  
\begin{align}
p(\bSigma)&=\left( 2^{vp/2}\pi^{p(p-1)/4}\prod_{i=1}^p \Gamma\left(\frac{v+1-i}{2}\right)    \right)^{-1} \notag \\
&\hspace{2cm}\times
 \mid \bLambda \mid ^{-v/2} \mid \bSigma \mid ^{(v-p-1)/2}
\exp\left\{ -\frac{1}{2}traza(\bLambda^{-1}\bSigma)\right\}
\end{align}

donde $\mid \bLambda \mid$ se refiere al determinante de la matriz $\bLambda$, la cual es simétrica y definida positiva de orden $p\times p$.</div>\EndKnitrBlock{definition}

\BeginKnitrBlock{proposition}<div class="proposition"><span class="proposition" id="prp:unnamed-chunk-61"><strong>(\#prp:unnamed-chunk-61) </strong></span>Si $\bSigma$ es una matriz aleatoria con distribución Wishart con $v$ grados de libertad, entonces $E(\bSigma)=v\bLambda$</div>\EndKnitrBlock{proposition}

### Distribución inversa-Wishart

\BeginKnitrBlock{definition}<div class="definition"><span class="definition" id="def:unnamed-chunk-62"><strong>(\#def:unnamed-chunk-62) </strong></span>Sea $\bSigma$ una matriz aleatoria simétrica y definida positiva de tamaño $p\times p$. Se dice que $\bSigma$ tiene distribución Wishart con $v$ grados de libertad, denotada como $\mathbf{Y}\sim Wishart_v(\bLambda)$, si su función de densidad está dada por:

\begin{align}
p(\bSigma)&=\left( 2^{vp/2}\pi^{p(p-1)/4}\prod_{i=1}^p \Gamma\left(\frac{v+1-i}{2}\right)    \right)^{-1} \notag \\
&\hspace{2cm}\times
 \mid \bLambda \mid ^{v/2} \mid \bSigma \mid ^{-(v+p+1)/2}
\exp\left\{ -\frac{1}{2}traza(\bLambda\bSigma^{-1})\right\}
\end{align}

donde $\mid \bLambda \mid$ se refiere al determinante de la matriz $\bLambda$, la cual es simétrica y definida positiva de orden $p\times p$.</div>\EndKnitrBlock{definition}

\BeginKnitrBlock{proposition}<div class="proposition"><span class="proposition" id="prp:unnamed-chunk-63"><strong>(\#prp:unnamed-chunk-63) </strong></span>Si $\bSigma$ es una matriz aleatoria con distribución inversa-Wishart con $v$ grados de libertad, entonces $E(\bSigma)=\dfrac{1}{v-p-1}\bLambda$</div>\EndKnitrBlock{proposition}

\BeginKnitrBlock{proposition}<div class="proposition"><span class="proposition" id="prp:unnamed-chunk-64"><strong>(\#prp:unnamed-chunk-64) </strong></span>Si $\bSigma^{-1}$ es una matriz aleatoria con distribución inversa-Wishart, entonces con $\bSigma$ tiene distribución Wishart.</div>\EndKnitrBlock{proposition}

