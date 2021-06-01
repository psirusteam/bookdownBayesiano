# Simulación de distribuciones de probabilidad

Como lo afirma @Gelman03 la simulación numérica es parte central del análisis bayesiano puesto que la generación de datos provenientes de una distribución de probabilidad se puede realizar fácilmente, incluso cuando la forma estructural de ésta no es conocida o es muy complicada computacionalmente. A lo largo de la historia del desarrollo de la teoría estadística, la simulación de distribuciones de probabilidad ha jugado un papel importante. Aunque son innumerables los métodos de generación de datos, en este apartado, se da cuenta de unos pocos, quizás lo más usados en este auge computacional.

## Método de la transformación uniforme

Al momento de la simulación estocástica de observaciones provenientes de alguna distribución de interés, la distribución uniforme es quizás la más usada y la más importante. El siguiente resultado adaptado de @Casella así lo confirma.

\BeginKnitrBlock{proposition}
<span class="proposition" id="prp:unnamed-chunk-1"><strong>(\#prp:unnamed-chunk-1) </strong></span>Si $U$ es una variable aleatoria con distribución uniforme en el intervalo $(0,1)$, entonces la la variable aleatoria $F^{-1}(U)$ tiene distribución F.
\EndKnitrBlock{proposition}
<br>

Aunque la función $F$ no necesariamente es una función uno a uno (por lo menos no lo es en el caso discreto) sí se puede verificar que $F^{-1}(U)$ es única con probabilidad uno. Una definición general, que encaja en el caso continuo o discreto, de la función $F$ inversa es la siguiente

\BeginKnitrBlock{definition}
<span class="definition" id="def:unnamed-chunk-2"><strong>(\#def:unnamed-chunk-2) </strong></span>Para cualquier función $F$ definida sobre $\mathbb{R}$, se define la función inversa generalizada de $F$ como

\begin{equation}
F^{-1}(u)=\inf \{x \mid F(x)\geq u\}
\end{equation}
\EndKnitrBlock{definition}


\BeginKnitrBlock{example}
<span class="example" id="exm:unnamed-chunk-3"><strong>(\#exm:unnamed-chunk-3) </strong></span>Suponga que $X$ es una variable aleatoria con distribución exponencial. De esta forma, su función de densidad acumulativa viene dada por

\begin{equation*}
F(x)=1-\exp\{ -\theta x\}
\end{equation*}

Del anterior resultado se tiene que si $u$ es una realización de una variable $U \sim Uniforme(0,1)$, entonces $F^{-1}(u)$ es una realización de una variable con distribución exponencial. Como $x=F^{-1}(u)$, entonces $F(x)=u$ y despejando $x$, se llega a que la siguiente expresión

\begin{equation*}
F^{-1}(u)=-\frac{\ln(1-u)}{\theta}
\end{equation*}

entrega una forma diáfana para la simulación de una observación con distribución exponencial. Para simular una muestra de $n$ observaciones, simplemente se repite el anterior procedimiento $n$ veces. En `R`, el código necesario para la simulación de una muestra de tamaño 1000 proveniente de una distribución exponencial con parámetro $\theta=5$ es
\EndKnitrBlock{example}


```r
theta <- 5
u <- runif(1000)
rexpo <- log(1 - u)/(-theta)
1/mean(rexpo)
```

```
## [1] 5.145846
```

```r
hist(rexpo, freq = FALSE)
lines(density(rexpo), col = 2)
```

![(\#fig:unnamed-chunk-4)Histograma de observaciones con distribución exponencial](A3Simulacion_files/figure-latex/unnamed-chunk-4-1.pdf) 

## Método de la grilla

Existen distribuciones de probabilidad cuya forma estructural es muy compleja. Mas aún, existen distribuciones de probabilidad conocidas para las cuales la inversa de la función de de densidad acumulativa es difícil de solucionar analíticamente. En los anteriores casos, el método analítico dado por el teorema de la transformación integral de probabilidad no siempre resulta efectivo. Sin embargo, es posible realizar una variante, manteniendo el espíritu de la anterior técnica.

El presente método utiliza una distribución discreta para aproximar cualquier tipo de distribución (discreta o continua) sin importar su nivel de complejidad. El algoritmo que enmarca este método se da a continuación:

1. Escribir la densidad de interés como $f(\cdot)$ y establecer el rango de la variable aleatoria de interés.
2. Fijar un conjunto de $n$ valores $x_1<\cdots<x_n$ equiespaciados que cubran una gran parte del rango de la variable aleatoria.
3. Para $x_k$ $(k=1,\ldots,n)$ calcular $f(x_k)$ que equivale al valor de la densidad en el punto $x_k$. Nótese que si $f(\cdot)$ es una función de densidad continua, entonces $f(x_k)$ no corresponde a una probabilidad;
4. Calcular la probabilidad asociada al punto $x_k$ definida por la aproximación discreta a $f(\cdot)$ y dada por
\begin{equation*}
p(x_k)=\frac{f(x_k)}{\sum_{k=1}^{n}f(x_k)}
\end{equation*}
5. Calcular la función de densidad acumulativa aproximada definida como
\begin{equation*}
F(x)=
\begin{cases}
0, \ \ \ \text{si $x<x_1$}\\
\sum_{l=1}^kp(x_l), \ \ \ \text{si $x_k\leq x<x_{k+1}$}\\
1, \ \ \ \text{si $x>x_n$}
\end{cases}
\end{equation*}
6. Simular una observación $u$ proveniente de una distribución uniforme continua en el intervalo $(0,1)$.
7. Si $F(x_k)<u\leq F(x_{k+1})$, entonces $F^{-1}(u)=x_{k+1}$ y por consiguiente el valor $x_{k+1}$ es una pseudo-observación proveniente de la densidad de interés.

Nótese que en el anterior proceso, la unidad $x_{k+1}$ es seleccionada con probabilidad $p_{k+1}$; puesto que

\begin{align*}
P(F(x_k)< U \leq F(x_{k+1}))&=F(x_{k+1})-F(x_{k})\\
&=\sum_{l=1}^{k+1}p(x_l)-\sum_{l=1}^kp(x_l)=p_{k+1}
\end{align*}

Si se quiere extraer una muestra aleatoria de $N$ observaciones provenientes de la distribución de interés, entonces basta con repetir el anterior proceso $N$ veces. Por supuesto, como se trata de una muestra aleatoria cada selección se debe realizar con repetición; de esta manera no importa si $N>n$. Suponiendo que el conjunto $x_1,\ldots,x_n$ conforma una grilla de puntos lo suficientemente cercanos y que no sucede nada importante entre cada uno de ellos, entonces esta técnica debe tener un buen funcionamiento.

\BeginKnitrBlock{example}
<span class="example" id="exm:unnamed-chunk-5"><strong>(\#exm:unnamed-chunk-5) </strong></span>El siguiente código computacional permite utilizar el método de la grilla para simular mil valores provenientes de una distribución exponencial con parámetro $\theta = 5$.
\EndKnitrBlock{example}


```r
theta <- 5
x.grid <- seq(0, 100, by=0.01)
p.exp <- theta * exp(-theta * x.grid)
r.exp <- sample(x.grid, 1000, prob = p.exp, replace = T)
1/mean(r.exp)
```

```
## [1] 5.05459
```

```r
hist(r.exp, freq = FALSE)
lines(density(rexpo), col = 2)
```

![(\#fig:unnamed-chunk-6)Histograma de observaciones con distribución exponencial](A3Simulacion_files/figure-latex/unnamed-chunk-6-1.pdf) 
\BeginKnitrBlock{example}
<span class="example" id="exm:unnamed-chunk-7"><strong>(\#exm:unnamed-chunk-7) </strong></span>De la misma manera, el método de la grilla permite simular valores de un distribución discreta. El siguiente código computacional permite utilizar el método de la grilla para simular mil valores provenientes de una distribución Poisson con parámetro $\theta = 2$.
\EndKnitrBlock{example}


```r
p.poisson <-function(theta, x.grid){
  N <- length(x.grid)
  res <- rep(NA, N)
  for(k in 1:N){
    P1 <- exp(-theta) * theta^(x.grid[k])
    P2 <- factorial(x.grid[k])
    res[k] <- P1/P2
  }
  return(res)
}

theta <- 2
x.grid <- seq(0, 100, by = 1)
f.x <- p.poisson(theta, x.grid)
p.x <- f.x/sum(f.x)
sum(p.x)
```

```
## [1] 1
```

```r
rpois <- sample(x.grid, 1000, prob=p.x, replace = T)
mean(rpois)
```

```
## [1] 2.059
```

```r
var(rpois)
```

```
## [1] 2.095615
```

```r
hist(rpois, freq = FALSE)
```

![(\#fig:unnamed-chunk-8)Histograma de observaciones con distribución Poisson](A3Simulacion_files/figure-latex/unnamed-chunk-8-1.pdf) 

\BeginKnitrBlock{example}
<span class="example" id="exm:unnamed-chunk-9"><strong>(\#exm:unnamed-chunk-9) </strong></span>El método de la grilla también puede utilizarse para simular observaciones de una distribución biparamétrica, univariada y continua. El siguiente código computacional permite utilizar el método de la grilla para simular mil valores provenientes de una distribución Gamma con parámetros $\alpha = 4$, $\beta = 2$.
\EndKnitrBlock{example}


```r
p.gamma <- function(a, b, x.grid){
  N <-length(x.grid)
  res <- rep(NA, N)
  for(k in 1:N){
    P1 <- (b^a)/gamma(a)
    P2 <- x.grid[k]^(a - 1)
    P3 <- exp(-b * x.grid[k])
    res[k] <- P1*P2*P3
  }
  return(res)
}

alpha <- 4
beta <- 2
x.grid <- seq(0, 100, by = 0.1)
f.x <- p.gamma(alpha, beta, x.grid)
p.x <- f.x / sum(f.x)

rgamma <- sample(x.grid, 1000, prob = p.x, replace = T)
mean(rgamma)
```

```
## [1] 2.0458
```

```r
var(rgamma)
```

```
## [1] 0.9759383
```

```r
hist(rgamma, freq = F)
lines(density(rgamma), col = 2)
```

![(\#fig:unnamed-chunk-10)Histograma de observaciones con distribución Gamma.](A3Simulacion_files/figure-latex/unnamed-chunk-10-1.pdf) 

\BeginKnitrBlock{example}
<span class="example" id="exm:unnamed-chunk-11"><strong>(\#exm:unnamed-chunk-11) </strong></span>Para comprobar el poder de este método de simulación, se presenta el siguiente código que permite simular valores de una distribución multiparamétrica, bivariada y continua. En particular, se simulan valores de la distribución Normal multivariante con vector de medias $\bmu = (2, 4)'$ y matriz de varianzas covarianzas 
$\bSigma = 
\begin{bmatrix} 
25 & 30 \\
30 & 16 
\end{bmatrix}$ 
\EndKnitrBlock{example}


```r
p.normal2 <- function(mu, Sigma, x, y){
  P1  <- 1/(2 * pi)
  P2  <- 1/sqrt(det(Sigma))
  P3a <- t((c(x, y) - mu)) %*% solve(Sigma) %*% (c(x, y) - mu)
  P3  <- exp((-1/2) * P3a)
  res <- P1 * P2 * P3
  return(res)
}

grilla <- function(a, b){
  A <- seq(1:length(a))
  unoA <- rep(1, length(A))
  B <- seq(1:length(b))
  unoB <- rep(1, length(B))
  P1 <- kronecker(A, unoB)
  P2 <- kronecker(unoA, B)
  grid <- cbind(a[P1], b[P2])
  return(grid)
}

mu1 <- c(2, 4)
Sigma1 <- matrix(c(25, 10, 10, 16), nrow=2)

x.grid <- seq(mu1[1] - 3 * sqrt(Sigma1[1, 1]),
              mu1[1] + 3 * sqrt(Sigma1[1, 1]),
              by = 0.5)
y.grid <- seq(mu1[2] - 3 * sqrt(Sigma1[2, 2]),
              mu1[2] + 3 * sqrt(Sigma1[2, 2]),
              by = 0.5)
xy.grid <- grilla(x.grid, y.grid)
N.grid <- dim(xy.grid)[1]

f.xy <- rep(NA, N.grid)
for(j in 1:N.grid){
  f.xy[j] <- p.normal2(mu1, Sigma1, 
                       xy.grid[j, 1], 
                       xy.grid[j, 2])
}

p.xy <- as.vector(f.xy/sum(f.xy))
sum(p.xy)
```

```
## [1] 1
```

```r
rnormal2 <- sample(N.grid, 1000, prob = p.xy, replace = T)
rxy.normal2 <- xy.grid[rnormal2, ]
rx.normal <- rxy.normal2[, 1]
ry.normal <- rxy.normal2[, 2]

colMeans(rxy.normal2)
```

```
## [1] 1.6305 3.7500
```

```r
var(rxy.normal2)
```

```
##          [,1]     [,2]
## [1,] 25.67640 10.63701
## [2,] 10.63701 15.69820
```

```r
hist(rx.normal, freq = F)
lines(density(rx.normal), col = 2)
```

![(\#fig:unnamed-chunk-12-1)Histogramas de observaciones con distribución Normal bivariada.](A3Simulacion_files/figure-latex/unnamed-chunk-12-1.pdf) 

```r
hist(ry.normal, freq = F)
lines(density(ry.normal), col = 2)
```

![(\#fig:unnamed-chunk-12-2)Histogramas de observaciones con distribución Normal bivariada.](A3Simulacion_files/figure-latex/unnamed-chunk-12-2.pdf) 