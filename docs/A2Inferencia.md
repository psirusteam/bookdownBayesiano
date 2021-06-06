# Matriz de información

\BeginKnitrBlock{definition}<div class="definition"><span class="definition" id="def:unnamed-chunk-1"><strong>(\#def:unnamed-chunk-1) </strong></span>Dada $X$ una variable aleatoria con función de densidad $f(x,\theta)$, donde $\theta$ es el parámetro de la distribución, y además existe $\dfrac{\partial}{\partial\theta}\ln{f(x,\theta)}$, entonces se define la información contenida en $X$ acerca de $\theta$ como

\begin{equation}
I_X(\theta)=E\left\{\left[\frac{\partial}{\partial\theta}\ln{f(X,\theta)}\right]^2\right\}.
\end{equation}</div>\EndKnitrBlock{definition}


\BeginKnitrBlock{proposition}<div class="proposition"><span class="proposition" id="prp:unnamed-chunk-2"><strong>(\#prp:unnamed-chunk-2) </strong></span>En la anterior definición, si además existe $\dfrac{\partial^2}{\partial\theta^2}\ln{f(x,\theta)}$, entonces se tiene que

\begin{equation}
I_X(\theta)=-E\left\{\dfrac{\partial^2}{\partial\theta^2}\ln{f(X,\theta)}\right\}.
\end{equation}</div>\EndKnitrBlock{proposition}
<br>

Las anteriores definiciones introducen la información contenida en una variable; sin embargo, cuando tenemos disponible una muestra aleatoria, es necesario definir la información contenida en una muestra aleatoria acerca de algún parámetro.

\BeginKnitrBlock{definition}<div class="definition"><span class="definition" id="def:unnamed-chunk-3"><strong>(\#def:unnamed-chunk-3) </strong></span>Dada $X_1$, $\cdots$, $X_n$ variables aleatorias con función de densidad $f(x_i,\theta)$, donde $\theta$ es el parámetro de la distribución, y además existe $\dfrac{\partial}{\partial\theta}\ln{\prod_{i=1}^nf(x_i,\theta)}$, entonces se define la información contenida en la muestra aleatoria acerca de $\theta$ como

\begin{equation}
I_{X_1,\cdots,X_n}(\theta)=E\left\{\left[\frac{\partial}{\partial\theta}\ln{\prod_{i=1}^nf(X_i,\theta)}\right]^2\right\}.
\end{equation}</div>\EndKnitrBlock{definition}


\BeginKnitrBlock{proposition}<div class="proposition"><span class="proposition" id="prp:unnamed-chunk-4"><strong>(\#prp:unnamed-chunk-4) </strong></span>Dada $X_1$, $\cdots$, $X_n$ una muestra aleatoria, entonces

\begin{equation*}
I_{X_1,\cdots,X_n}(\theta)=nI_X(\theta),
\end{equation*}

donde $I_X(\theta)=I_{X_i}(\theta)$, con $i=1,\cdots,n$. Es decir, en una muestra aleatoria, cada variable aporta la misma cantidad de información, y la cantidad total de información en la muestra es la suma de la información en cada variable.</div>\EndKnitrBlock{proposition}
<br>

\BeginKnitrBlock{proof}<div class="proof">\iffalse{} <span class="proof"><em>Prueba. </em></span>  \fi{}\begin{align*}
I_{X_1,\cdots,X_n}(\theta)&=E\left\{\left[\frac{\partial}{\partial\theta}\ln{\prod_{i=1}^nf(X_i,\theta)}\right]^2\right\}\\                   &=E\left\{\left[\sum_{i=1}^n\frac{\partial}{\partial\theta}\ln{f(X_i,\theta)}\right]^2\right\}\\                       &=E\left\{\sum_{i=1}^n\left[\frac{\partial}{\partial\theta}\ln{f(X_i,\theta)}\right]^2\right\}+\\
                          &\ \ \ \ \ \ \ \ \ \ \ \ \underbrace{E\left\{\sum_{\substack{i,j=1\\i\neq j}}^n\left[\frac{\partial}{\partial\theta}\ln{f(X_i,\theta)}\frac{\partial}{\partial\theta}\ln{f(X_j,\theta)}\right]\right\}}_{=0,\ \text{por la independencia entre}\ X_i\ \text{y}\ X_j}\\                      &=\sum_{i=1}^nE\left\{\left[\frac{\partial}{\partial\theta}\ln{f(X_i,\theta)}\right]^2\right\}\\
                          &=\sum_{i=1}^nI_X(\theta)=nI_X(\theta).
\end{align*}</div>\EndKnitrBlock{proof}
<br>

\BeginKnitrBlock{example}<div class="example"><span class="example" id="exm:unnamed-chunk-6"><strong>(\#exm:unnamed-chunk-6) </strong></span>Sea $X_1$, $\cdots$, $X_n$ una muestra aleatoria proveniente de la distribución $N(\mu,\sigma^2)$, la información contenida en la muestra acerca de $\mu$ es $n/\sigma^2$. Para verificar esta afirmación, calculamos la información acerca de $\mu$ en una variable $X$ con distribución $N(\mu,\sigma^2)$. Tenemos:
  
\begin{align*}
I_X(\mu)&=-E\left\{\dfrac{\partial^2}{\partial\mu^2}\ln{f(X,\theta)}\right\}\\
        &=-E\left\{\dfrac{\partial^2}{\partial\mu^2}\left[-\frac{1}{2}\ln2\pi\sigma^2-\frac{1}{2\sigma^2}(X-\mu)^2\right]\right\}\\
        &=-E\left\{\frac{\partial}{\partial\mu}\left[\frac{X-\mu}{\sigma^2}\right]\right\}\\
        &=-E\left\{-\frac{1}{\sigma^2}\right\}\\
        &=\frac{1}{\sigma^2}.
\end{align*}

Ahora, usando el Resultado 2.3.4, se tiene que $I_{X_1,\cdots,X_n}(\mu)=n/\sigma^2$.</div>\EndKnitrBlock{example}
<br>

Nótese que esta información, en primer lugar, depende del tamaño $n$ de manera que entre más grande sea la muestra, hay mayor información acerca de $\mu$; en segundo lugar, entre más pequeña sea la varianza $\sigma^2$, la cantidad de información acerca de $\mu$ también incrementa, esto es natural, puesto que si $\sigma^2$ es pequeña, los datos de la muestra están muy concentrados alrededor de $\mu$, entonces estos datos aportan más información que otros datos con más dispersión.

\BeginKnitrBlock{definition}<div class="definition"><span class="definition" id="def:unnamed-chunk-7"><strong>(\#def:unnamed-chunk-7) </strong></span>Dada una variable aleatoria $X$ con función de densidad $f(x,\btheta)$, la matriz de información contenida en $X$ acerca de $\btheta$ se define como

\begin{equation}
I_X(\btheta)=E\left\{\frac{\partial\ln f(X,\btheta)}{\partial\btheta}\left(\frac{\partial\ln f(X,\btheta)}{\partial\btheta}\right)'\right\}
\end{equation}</div>\EndKnitrBlock{definition}


\BeginKnitrBlock{definition}<div class="definition"><span class="definition" id="def:unnamed-chunk-8"><strong>(\#def:unnamed-chunk-8) </strong></span>Dada una muestra aleatoria $X_1$, $\cdots$, $X_n$ con función de densidad $f(x_i,\btheta)$, la matriz de información contenida en la muestra acerca de $\btheta$ se define como

\begin{equation*}
I_{X_1,\cdots,X_n}(\btheta)=E\left\{\frac{\partial\ln \prod_{i=1}^nf(X_i,\btheta)}{\partial\btheta}\left(\frac{\partial\ln \prod_{i=1}^nf(X_i,\btheta)}{\partial\btheta}\right)'\right\}
\end{equation*}</div>\EndKnitrBlock{definition}


\BeginKnitrBlock{example}<div class="example"><span class="example" id="exm:unnamed-chunk-9"><strong>(\#exm:unnamed-chunk-9) </strong></span>Dada una muestra aleatoria $X_1$, $\cdots$, $X_n$ con distribución común $N(\mu,\sigma^2)$, vamos a hallar la matriz de información contenida en la muestra acerca del vector de parámetros $(\mu,\sigma^2)$. Tenemos que

\begin{align*}
&\ \ \ \ \ \ \ I_{X_1,\cdots,X_n}(\mu,\sigma^2)\\
&=E\left\{
\begin{pmatrix}
\dfrac{\partial\ln \prod_{i=1}^nf(X_i,\mu,\sigma^2)}{\partial\mu}\\
\dfrac{\partial\ln \prod_{i=1}^nf(X_i,\mu,\sigma^2)}{\partial\sigma^2}
\end{pmatrix}
\begin{pmatrix}
\dfrac{\partial\ln \prod_{i=1}^nf(X_i,\mu,\sigma^2)}{\partial\mu}&
\dfrac{\partial\ln \prod_{i=1}^nf(X_i,\mu,\sigma^2)}{\partial\sigma^2}
\end{pmatrix}
\right\}\\
&=E\left\{
\begin{pmatrix}
\dfrac{\sum_{i=1}^nX_i-n\mu}{\sigma^2}\\
\dfrac{\sum_{i=1}^n(X_i-\mu)^2-n\sigma^2}{2\sigma^4}
\end{pmatrix}
\begin{pmatrix}
\dfrac{\sum_{i=1}^nX_i-n\mu}{\sigma^2}&
\dfrac{\sum_{i=1}^n(X_i-\mu)^2-n\sigma^2}{2\sigma^4}
\end{pmatrix}
\right\}\\
&=E\left\{\begin{pmatrix}
\dfrac{(\sum_{i=1}^nX_i-n\mu)^2}{\sigma^4}&\dfrac{(\sum_{i=1}^nX_i-n\mu)(\sum_{i=1}^n(X_i-\mu)^2-n\sigma^2)}{2\sigma^6}\\
\dfrac{(\sum_{i=1}^nX_i-n\mu)(\sum_{i=1}^n(X_i-\mu)^2-n\sigma^2)}{2\sigma^6}&\dfrac{(\sum_{i=1}^n(X_i-\mu)^2-n\sigma^2)^2}{4\sigma^8}
\end{pmatrix}\right\}
\end{align*}

Donde el primer elemento diagonal de la anterior matriz está dada por
\begin{align*}
E\left\{\dfrac{(\sum_{i=1}^nX_i-n\mu)^2}{\sigma^4}\right\}&=\left[Var\left(\sum_{i=1}^nX_i-n\mu\right)+(E\left(\sum_{i=1}^nX_i-n\mu\right))^2\right]/\sigma^4\\
&=n\sigma^2/\sigma^4=n/\sigma^2.
\end{align*}

El segundo elemento diagonal está dada por
\begin{align}\label{feo}
&\ \ \ \ \ E\left\{\dfrac{(\sum_{i=1}^n(X_i-\mu)^2-n\sigma^2)^2}{4\sigma^8}\right\}\\
&=\frac{1}{4\sigma^8}E\left\{\left[\sum_{i=1}^n(X_i-\mu)^2\right]^2+n^2\sigma^4-2n\sigma^2\sum_{i=1}^n(X_i-\mu)^2\right\}\\
&=\frac{1}{4\sigma^8}\left\{Var(\sum_{i=1}^n(X_i-\mu)^2)+\left[E(\sum_{i=1}^n(X_i-\mu)^2)\right]^2+n^2\sigma^4-2n\sigma^2E\left[\sum_{i=1}^n(X_i-\mu)^2\right]\right\}
\end{align}

Usando el hecho de que
\begin{equation*}
\frac{\sum_{i=1}^n(X_i-\mu)^2}{\sigma^2}\sim\chi^2_n
\end{equation*}

y la esperanza y varianza de la distribución $\chi^2_n$, tenemos que la expresión (\ref{feo}) está dada por
\begin{equation*}
\frac{1}{4\sigma^8}\left\{2n\sigma^4+\left[n\sigma^2\right]^2+n^2\sigma^4-2n\sigma^2n\sigma^2\right\}=\frac{n}{2\sigma^4}.
\end{equation*}

Finalmente, el elemento fuera de la diagonal de la matriz $I_{X_1,\cdots,X_n}(\mu,\sigma^2)$ está dado por
\begin{align*}
&\ \ \ \ \ \ E\left\{\left(\sum_{i=1}^nX_i-n\mu\right)\left(\sum_{i=1}^n(X_i-\mu)^2-n\sigma^2\right)\right\}\\
&=E\left\{\sum_{i=1}^nX_i\left(\sum_{i=1}^n(X_i-\mu)^2-n\sigma^2\right)-n\mu\left(\sum_{i=1}^n(X_i-\mu)^2-n\sigma^2\right)\right\}\\
&=E\left\{\sum_{i=1}^nX_i\sum_{i=1}^n(X_i-\mu)^2\right\}-n\sigma^2E\left(\sum_{i=1}^nX_i\right)-n\mu E\left(\sum_{i=1}^n(X_i-\mu)^2\right)+n^2\mu\sigma^2\\
&=E\left(\sum_{i=1}^nX_i\sum_{i=1}^nX_i^2\right)-2\mu E\left[(\sum_{i=1}^nX_i)^2\right]+n^2\mu^3-n^2\mu\sigma^2-n^2\mu\sigma^2+n^2\mu\sigma^2\\
&=E\left(\sum_{i=1}^nX_i^3+\sum_{i\neq j}X_iX_j^2\right)-2\mu(n\sigma^2+n^2\mu^2)+n^2\mu^3-n^2\mu\sigma^2\\
&=\sum_{i=1}^n\left[3\mu E(X_i^2)-2\mu^3\right]+\sum_{i\neq j}E(X_i)E(X_j^2)-2n\mu\sigma^2-2n^2\mu^3+n^2\mu^3-n^2\mu\sigma^2\\
&=3n\mu(\sigma^2+\mu^2)-2n\mu^3+\mu(\sigma^2+\mu^2)(n^2-n)-2n\mu\sigma^2-2n^2\mu^3+n^2\mu^3-n^2\mu\sigma^2\\
&=0
\end{align*}

De donde obtenemos finalmente la matriz de información $I_{X_1,\cdots,X_n}(\mu,\sigma^2)$ dada por
\begin{equation*}
I_{X_1,\cdots,X_n}(\mu,\sigma^2)=\begin{pmatrix}
\dfrac{n}{\sigma^2}&0\\
0&\dfrac{n}{2\sigma^4}
\end{pmatrix}
\end{equation*}</div>\EndKnitrBlock{example}
