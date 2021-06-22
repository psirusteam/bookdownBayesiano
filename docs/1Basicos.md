
# Tópicos básicos

Para poder entender la racionalidad del paradigma bayesiano es importante reconocer que las innovaciones de las metodologías bayesianas descansan en principios clásicos de la teoría de probabilidad. 

## Teoría de la decisión

El problema estadístico de estimar un parámetro se puede ver dentro del contexto de la teoría de decisión: la estimación que proveemos, sea en el ámbito de la estadística clásica o la estadística bayesiana, depende de los datos muestrales, $\mathbf{X}$, de tal forma que si éstos cambian, la estimación también cambia. De esta manera, el proceso de estimación puede ser representado como una función que toma un conjunto de datos muestrales y los convierte en una estimación ($A(\mathbf{X})$ o simplemente $A$) del parámetro de interés. En la teoría de decisión, la anterior función se conoce como una regla de decisión.

Así como en la vida cotidiana, por la incertidumbre del futuro (en el ámbito estadístico, por la incertidumbre acerca del parámetro), toda acción que se tome (toda estimación que se provea) puede traer consigo un grado de falla o riesgo. Y es necesario escoger la acción óptima que de alguna forma minimice ese riesgo. Formalizando esta idea intuitiva, se define la función de pérdida $L$ que asocia a cada dupla conformada por la acción tomada y el parámetro de interés $\theta$, $(A, \ \theta)$ con un número no negativo que cuantifica la pérdida que ocasiona la acción (o la estimación) $A$ con respecto al parámetro $\theta$.

Es claro que se desea escoger aquella acción que minimice de alguna forma la pérdida que ésta ocasiona, pero la función $L$ no se puede minimizar directamente, puesto que:

  * En el ámbito de la estadística clásica, el parámetro $\theta$ se considera fijo, y los datos muestrales $\mathbf{X}$ aleatorios. Como la función de pérdida $L$ depende de $\mathbf{X}$, entonces ésta también será una variable aleatoria, y no se puede minimizar directamente. Por lo tanto se define el riesgo o la pérdida promedio como la esperanza matemática de $L$; denotando el riesgo como $R$, éste está definido como $R=E(L)$ (la esperanza se toma con respecto a la distribución probabilística de $\mathbf{X}$).
  
  * En el ámbito de la estadística bayesiana, $\theta$ sigue siendo una cantidad fija, pero la incertidumbre que tiene el investigador sobre la localización del parámetro se puede modelar mediante funciones de probabilidad. La herramienta fundamental para conocer características de $\theta$ es su función de densidad posterior $p(\theta|\mathbf{X})$. En este caso, el riesgo $R$ se define como

\begin{equation*}
R=E(L)=\int L(A, \theta)p(\theta|\mathbf{X})d\theta
\end{equation*}

En cualquiera de los dos casos anteriores, se busca la estimación que minimice el riesgo $R$. Ilustramos los anteriores conceptos en los siguientes ejemplos tanto en la estadística clásica como en la estadística bayesiana.

\BeginKnitrBlock{example}<div class="example"><span class="example" id="exm:unnamed-chunk-1"><strong>(\#exm:unnamed-chunk-1) </strong></span>Sea $X_i$ con $i=1,\cdots, n$ una muestra aleatoria con media $\theta$ y varianza $\sigma^2$, ambas fijas, y suponga que se desea encontrar el mejor estimador de $\theta$ bajo la función de pérdida cuadrática dada por

\begin{equation*}
L(A,\theta)=(A-\theta)^2
\end{equation*}

cuyo riesgo asociado está dado por $R=E(A-\theta)^2$. En primer lugar, buscaremos dicho estimador dentro de todas las formas lineales de $X_i$, es decir, los estimadores de la forma $A=\sum_{i=1}^nc_iX_i$. Por tanto, el riesgo se puede expresar como
\begin{align*}
R=E(A-\theta)^2&=Var(A)+(E(A)-\theta)^2\\
&=\sum_{i=1}^nc_i^2\sigma^2+\theta^2(\sum_{i=1}^nc_i-1)^2
\end{align*}

Y al buscar los coeficientes $c_i$ que minimizan la anterior expresión, encontramos que $c_i=\theta^2/(\sigma^2+n\theta^2)$ para todo $i$. Como estos coeficientes conducen a un estimador que depende del parámetro desconocido, concluimos que no hay ningún estimador que minimiza el riesgo.

Para encontrar una solución, es necesario restringir aún más el rango de estimadores; para eso, se impone la restricción de que $\sum_{i=1}^n c_i=1$. De esta forma, el riesgo está dado por $R=\sum c_i^2\sigma^2$. Dado que $\sigma^2$ es fijo, al minimizar $\sum c_i^2$ sujeto a la restricción, se tiene que la solución es $c_i=1/n$ para todo $i$, y así encontramos que el mejor estimador (en el sentido de minimizar el riesgo de la función de pérdida cuadrática) dentro de todas las formas lineales con $\sum c_i=1$ es la media muestral $\bar{X}$.</div>\EndKnitrBlock{example}

\BeginKnitrBlock{example}<div class="example"><span class="example" id="exm:unnamed-chunk-2"><strong>(\#exm:unnamed-chunk-2) </strong></span>Suponga que se desea estimar un parámetro de interés $\theta$ en el contexto de la estadística bayesiana y denotamos la función de densidad posterior de $\theta$ como $p(\theta|\mathbf{X})$, entonces si utilizamos la función de pérdida cuadrática, el riesgo asociado será

\begin{align*}
R&=E(L(A,\theta))=E (A-\theta)^2=Var(\theta)+(E(\theta)-A)^2
\end{align*}

que es minimizado si $A=E(\theta)$. Es decir, la mejor acción para estimar $\theta$ es utilizar su tomada con respecto a la distribución posterior $p(\theta|\mathbf{X})$.</div>\EndKnitrBlock{example}

\BeginKnitrBlock{example}<div class="example"><span class="example" id="exm:unnamed-chunk-3"><strong>(\#exm:unnamed-chunk-3) </strong></span>En el mismo contexto del ejemplo anterior, si cambiamos la función de pérdida a la siguiente
\begin{equation*}
L(A,\theta)=|A-\theta|=(A-\theta)I_{(A\geq\theta)}+(\theta-A)I_{(\theta>A)}
\end{equation*}

El riesgo estará dado por
\begin{align*}
R&=E(L(A,\theta))\\
&=\int L(A,\theta)p(\theta|\mathbf{X})d\theta\\
&=\int_{(A\geq\theta)}(A-\theta)p(\theta|\mathbf{X})d\theta+\int_{(\theta>A)}(\theta-A)p(\theta|\mathbf{X})d\theta
\end{align*}

Derivando el riesgo con respecto a la acción $A$, se tiene que
\begin{equation*}
\frac{\partial R}{\partial A}=\int_{(A\geq\theta)}p(\theta|\mathbf{X})d\theta-\int_{(\theta>A)}p(\theta|\mathbf{X})d\theta
\end{equation*}

Igualando a cero, tenemos que
\begin{equation*}
\int_{(A\geq\theta)}p(\theta|\mathbf{X})d\theta=\int_{(\theta>A)}p(\theta|\mathbf{X})d\theta=0.5
\end{equation*}

Y concluimos que la acción $A$ que induce menor riesgo corresponde al percentil 50% o la mediana de la distribución posterior de $\theta$.</div>\EndKnitrBlock{example}
<br>

De los anteriores ejemplos se observa que, bajo un mismo contexto, cuando se utilizan diferentes funciones de pérdida, también se obtienen distintas estimaciones, y distintas acciones que optimizan el riesgo.

## Algunos resultados de probabilidad

Antes de entrar en el repaso de estos conceptos fundamentales, se definen los conceptos de **parámetro** y **espacio paramétrico** asociados a una distribución de probabilidad.

1. Un parámetro es aquella cantidad que define la forma funcional de una distribución de probabilidad; es decir, cuando el parámetro cambia de valor, la función de densidad y la función de distribución cambian. Las distribuciones de probabilidad pueden tener más de un parámetro. Cuando una distribución tiene solo un parámetro, éste se denota usualmente por $\theta$, cuando se presenta más de un parámetro, la notación se cambia a $\btheta$, representando el vector de parámetros. 
2. El espacio paramétrico, $\Theta$, es el conjunto que contiene todos los posibles valores que puede tomar el parámetro o el vector de parámetros. Para distribuciones con un solo parámetro, $\Theta$ será un subconjunto de $\mathbb{R}$, mientras que para distribuciones con dos o más parámetros, $\Theta$ será un subconjunto de $\mathbb{R}\times\mathbb{R}$.

Para entender los fundamentos de la modelación bayesiana, es necesario
recordar algunas definiciones y resultados de la teoría de probabilidad
que ayudarán a hacer más expedito este periplo por la estadística
bayesiana. En términos de notación, se utilizará indistintamente la
expresión de integral, $\int$, indicando la sumatoria, en el caso de las
variables aleatorias discretas o la integral de Riemann-Stieltjes en el
caso de las variables aleatorias continuas.

\BeginKnitrBlock{definition}<div class="definition"><span class="definition" id="def:unnamed-chunk-4"><strong>(\#def:unnamed-chunk-4) </strong></span>Sean $\mathbf{X}=(X_1,\ldots,X_p)'$, $\mathbf{Y}=(Y_1,\ldots,Y_q)'$ dos vectores aleatorios definidos sobre los espacios de  muestreo $\mathcal{X}$, $\mathcal{Y}$, respectivamente. Suponga que la distribución conjunta de estos vectores aleatorios está dada por $p(\mathbf{X},\mathbf{Y})$. La distribución marginal de $\mathbf{X}$ está dada por
\begin{equation}
p(\mathbf{X})=\int p(\mathbf{X},\mathbf{Y})\ d\mathbf{Y}
\end{equation}
y la distribución condicional de $\mathbf{X}$ dado $\mathbf{Y}$ como
\begin{equation}
p(\mathbf{X} \mid \mathbf{Y})
=\frac{p(\mathbf{X},\mathbf{Y})}{p(\mathbf{Y})}
\end{equation}</div>\EndKnitrBlock{definition}


\BeginKnitrBlock{proposition}<div class="proposition"><span class="proposition" id="prp:Res121"><strong>(\#prp:Res121) </strong></span>Suponga los vectores $\mathbf{X}$, $\mathbf{Y}$ y un tercer vector $\mathbf{Z}=(Z_1,\ldots,Z_r)'$ definido sobre el espacio de muestreo  $\mathcal{Z}$. Entonces se tiene que
\begin{equation}
p(\mathbf{X} \mid \mathbf{Z})=\int p(\mathbf{X},\mathbf{Y} \mid \mathbf{Z})\ d\mathbf{Y}
\end{equation}
y
\begin{equation}
p(\mathbf{X} \mid \mathbf{Y},\mathbf{Z})=\frac{p(\mathbf{X},\mathbf{Y} \mid \mathbf{Z})}{p(\mathbf{Y} \mid \mathbf{Z})}
\end{equation}</div>\EndKnitrBlock{proposition}
<br>

\BeginKnitrBlock{proof}<div class="proof">\iffalse{} <span class="proof"><em>Prueba. </em></span>  \fi{}En primer lugar, nótese que
\begin{align*}
\int p(\mathbf{X},\mathbf{Y} \mid \mathbf{Z})\ d\mathbf{Y}&=
\int \frac{p(\mathbf{X},\mathbf{Y},\mathbf{Z})}{p(\mathbf{Z})}\ d\mathbf{Y}\\
&=\frac{1}{p(\mathbf{Z})} \int p(\mathbf{X},\mathbf{Y},\mathbf{Z}) \ d\mathbf{Y}\\
&=\frac{1}{p(\mathbf{Z})} p(\mathbf{X},\mathbf{Z})=p(\mathbf{X} \mid \mathbf{Z})
\end{align*}

Por otro lado,

\begin{align*}
\frac{p(\mathbf{X},\mathbf{Y} \mid \mathbf{Z})}{p(\mathbf{Y} \mid \mathbf{Z})}=
\frac{p(\mathbf{X},\mathbf{Y},\mathbf{Z})}{p(\mathbf{Z})} \diagup
\frac{p(\mathbf{Y},\mathbf{Z})}{p(\mathbf{Z})}
=\frac{p(\mathbf{X},\mathbf{Y},\mathbf{Z})}{p(\mathbf{Y},\mathbf{Z})}=p(\mathbf{X} \mid \mathbf{Y},\mathbf{Z})
\end{align*}</div>\EndKnitrBlock{proof}
<br>

\BeginKnitrBlock{definition}<div class="definition"><span class="definition" id="def:unnamed-chunk-6"><strong>(\#def:unnamed-chunk-6) </strong></span>Sean $\mathbf{X}$, $\mathbf{Y}$, $\mathbf{Z}$ vectores aleatorios, se dice que $\mathbf{X}$ es condicionalmente independiente de $\mathbf{Y}$ con respecto a $\mathbf{Z}$ si satisfacen la siguiente expresión
\begin{equation}
p(\mathbf{X},\mathbf{Y} \mid \mathbf{Z})=p(\mathbf{X} \mid \mathbf{Z})p(\mathbf{Y} \mid \mathbf{Z})
\end{equation}</div>\EndKnitrBlock{definition}


\BeginKnitrBlock{proposition}<div class="proposition"><span class="proposition" id="prp:Res122"><strong>(\#prp:Res122) </strong></span>Si $\mathbf{X}$ es condicionalmente independiente de $\mathbf{Y}$ con respecto a $\mathbf{Z}$, entonces se tiene que
\begin{equation}
p(\mathbf{X} \mid \mathbf{Y},\mathbf{Z})=p(\mathbf{X} \mid \mathbf{Z})
\end{equation}</div>\EndKnitrBlock{proposition}
<br>

\BeginKnitrBlock{proof}<div class="proof">\iffalse{} <span class="proof"><em>Prueba. </em></span>  \fi{}Como $p(\mathbf{X},\mathbf{Y} \mid \mathbf{Z})=\dfrac{p(\mathbf{X},\mathbf{Y},\mathbf{Z})}{p(\mathbf{Z})}$, entonces

\begin{align*}
p(\mathbf{X} \mid \mathbf{Y},\mathbf{Z})=\frac{p(\mathbf{X},\mathbf{Y},\mathbf{Z})}{p(\mathbf{Y},\mathbf{Z})}
=\frac{p(\mathbf{X},\mathbf{Y} \mid \mathbf{Z})p(\mathbf{Z})}{p(\mathbf{Y},\mathbf{Z})}
=\frac{p(\mathbf{X} \mid \mathbf{Z})p(\mathbf{Y} \mid \mathbf{Z})}{p(\mathbf{Y} \mid \mathbf{Z})}=p(\mathbf{X} \mid \mathbf{Z})
\end{align*}</div>\EndKnitrBlock{proof}
<br>

\BeginKnitrBlock{proposition}<div class="proposition"><span class="proposition" id="prp:unnamed-chunk-8"><strong>(\#prp:unnamed-chunk-8) </strong></span>Si $\mathbf{X}$ es independiente de $\mathbf{Y}$, entonces $\mathbf{X}$ es condicionalmente independiente de $\mathbf{Y}$ dado cualquier otro vector $\mathbf{Z}$.</div>\EndKnitrBlock{proposition}
<br>

\BeginKnitrBlock{proof}<div class="proof">\iffalse{} <span class="proof"><em>Prueba. </em></span>  \fi{}Nótese que
\begin{equation*}
p(\mathbf{X},\mathbf{Y}\mid \mathbf{Z})=p(\mathbf{X} \mid \mathbf{Y},\mathbf{Z})p(\mathbf{Y} \mid \mathbf{Z})=p(\mathbf{X} \mid \mathbf{Z})p(\mathbf{Y} \mid \mathbf{Z})
\end{equation*}

puesto que, utilizando la hipótesis de independencia, se tiene que
\begin{equation*}
p(\mathbf{X} \mid \mathbf{Y})=p(\mathbf{X})
\end{equation*}</div>\EndKnitrBlock{proof}
<br>

## Teorema de Bayes

Desde la revolución estadística de Pearson y Fisher, la inferencia
estadística busca encontrar los valores que parametrizan a la
distribución desconocida de los datos. El primer enfoque, propuesto por
Pearson, afirmaba que si era posible observar a la variable de interés
en todos y cada uno de los individuos de una población, entonces era
posible calcular los parámetros de la distribución de la variable de
interés; por otro lado, si solo se tenía acceso a una muestra
representativa, entonces era posible calcular una estimación de tales
parámetros. Sin embargo, Fisher discrepó de tales argumentos, asumiendo
que las observaciones están sujetas a un error de medición y por lo
tanto, así se tuviese acceso a toda la población, sería imposible calcular
los parámetros de la distribución de la variable de interés.

Del planteamiento de Fisher resultaron una multitud de métodos
estadísticos para la estimación de los parámetros poblacionales. Es
decir, si la distribución de $\mathbf{Y}$ está parametrizada por
$\btheta=(\theta_1,\ldots,\theta_K)$, $\btheta \in \Theta$ con $\Theta$
el espacio paramétrico inducido por el comportamiento de la variable de
interés, el objetivo de la teoría estadística inferencial es calcular
una estimación $\hat{\btheta}$ del parámetro $\btheta$, por medio de los
datos observados. En este enfoque, los parámetros se consideran
cantidades fijas y constantes. Sin embargo, en la última mitad del siglo
XX, algunos investigadores estadísticos comenzaron a reflexionar acerca
de la naturaleza de $\btheta$ y enfocaron la inferencia estadística de
una manera distinta: *asumiendo que la distribución de la variable de
interés está condicionada a valores específicos de los parámetros*. Es
decir, en términos de notación, si la variable de interés es
$\mathbf{Y}$, su distribución condicionada a los parámetros toma la
siguiente forma $p(\mathbf{Y} \mid \btheta)$. Esto implica claramente
que en este nuevo enfoque la naturaleza de los parámetros no es
constante.

En términos de inferencia para $\btheta$, es necesario encontrar la
distribución de los parámetros condicionada a la observación de los
datos. Para este fin, es necesario definir la distribución conjunta de
la variable de interés con el vector de parámetros. 
\begin{equation*}
p(\btheta,\mathbf{Y})=p(\btheta)p(\mathbf{Y} \mid \btheta)
\end{equation*}

A la distribución $p(\btheta)$ se le conoce con el nombre de
distribución *previa* y en ella se enmarcan todas y cada una de las
creencias que se tienen acerca del comportamiento estocástico del vector
de parámetros antes de que ocurra la recolección de los datos; $p(\mathbf{Y} \mid \btheta)$ es la distribución de muestreo,
verosimilitud o distribución de los datos. Por otro lado, la
distribución del vector de parámetros condicionada a los datos
observados está dada por 

\begin{equation}
(\#eq:Bayes)
p(\btheta \mid \mathbf{Y})=\frac{p(\btheta,\mathbf{Y})}{p(\mathbf{Y})}=\frac{p(\btheta)p(\mathbf{Y} \mid \btheta)}{p(\mathbf{Y})}
\end{equation}

A la distribución $p(\btheta \mid \mathbf{Y})$ se le conoce con el
nombre de distribución *posterior* y en ella se enmarcan las
creencias actualizadas acerca del comportamiento estocástico del vector
de parámetros teniendo en cuenta los datos observados $\mathbf{Y}$.
Nótese que la expresión \@ref(eq:Bayes) se compone de una fracción cuyo
denominador no depende del vector de parámetros y considerando a los
datos observados como fijos, corresponde a una constante y puede ser
obviada. Por lo tanto, otra representación de la regla de Bayes está
dada por 

\begin{align}
(\#eq:Bayes1)
p(\btheta \mid \mathbf{Y})\propto p(\mathbf{Y} \mid \btheta)p(\btheta)
\end{align}

@Gelman03 menciona que esta expresión se conoce como la
distribución *posterior no-normalizada* y encierra el núcleo
técnico de la inferencia bayesiana. La constante $p(\mathbf{Y})$
faltante en la expresión \@ref(eq:Bayes1) se da a continuación.

\BeginKnitrBlock{proposition}<div class="proposition"><span class="proposition" id="prp:Res131"><strong>(\#prp:Res131) </strong></span>La expresión $p(\mathbf{Y})$ corresponde a una constante $k$ tal que
\begin{equation*}
k=p(\mathbf{Y})=E_{\btheta}[p(Y \mid \btheta)]
\end{equation*}</div>\EndKnitrBlock{proposition}
<br>

\BeginKnitrBlock{proof}<div class="proof">\iffalse{} <span class="proof"><em>Prueba. </em></span>  \fi{}Nótese que
\begin{equation*}
k=p(\mathbf{Y})=\int p(\mathbf{Y},\btheta)\ d\btheta=\int p(\btheta)p(\mathbf{Y} \mid \btheta)\ d\btheta.
\end{equation*}
entonces
\begin{align*}
k&=\int p(\mathbf{Y} \mid \btheta)p(\btheta)\ d\btheta\\
&=E_{\btheta}[p(Y \mid \btheta)]
\end{align*}</div>\EndKnitrBlock{proof}
<br>

Curiosamente, el reverendo Thomas Bayes nunca publicó este resultado,
sino que después de su fallecimiento, su amigo el filósofo Richard
Price, encontró los escritos dentro de sus pertenencias, y éstos fueron
publicados en el 1764 en
*Philosophical Transactions of the Royal Society of London*. Aunque
el teorema de Bayes fue nombrado en honor de Thomas Bayes, es casi
seguro que él mismo no sospechaba del gran impacto de su resultado. De hecho, aproximadamente una década más tarde, Pierre-Simon Laplace también descrubrió el mismo principio, y dedicó gran parte de su vida extendiéndolo y formalizándolo. Más aún, él analizó grandes volumenes de datos relacionados a los nacimientos en diferentes paises para confirmar esta teoría, y sentó las bases de la estadística bayesiana.

A continuación se presenta un ejemplo simple de este sencillo pero
poderoso teorema.

\BeginKnitrBlock{example}<div class="example"><span class="example" id="exm:unnamed-chunk-11"><strong>(\#exm:unnamed-chunk-11) </strong></span>Suponga que una fábrica del sector industrial produce bolígrafos y que la producción está a cargo de tres máquinas. La primera máquina produce el 50% del total de bolígrafos en el año, la segunda máquina produce el 30% y la última maquina produce el restante 20%. Por supuesto, esta producción esta sujeta al error y por tanto, basados en la experiencia, es posible reconocer que, de los artículos producidos por la primera máquina, el 5% resultan defectuosos; de los artículos producidos por la segunda máquina, el 2% resultan defectuosos y, de los artículos producidos por la última máquina, el 6% resultan defectuosos.

Una pregunta natural que surge es acerca de la probabilidad de selección de un artículo defectuoso y para responder a esta pregunta con rigurosidad de probabilística es necesario enfocar la atención en los tópicos básicos que dejamos atrás. En primer lugar, el experimento en cuestión es la selección de un bolígrafo. Para este experimento, una terna $(\Omega, \mathfrak{F}, P)$ ^[$\Omega$ denota el conjunto de todos lo posibles resultados del experimento, $\mathfrak{F}$ denota una $\sigma$-álgebra y $P$ hace referencia ana medida de probabilidad propiamente definida.], llamada comúnmente espacio de medida o espacio de probabilidad, está dada por

1. El espacio muestral: $\Omega=\{\text{defectuoso}, \text{No defectouso}\}$
1. La $\sigma$-álgebra: $\mathfrak{F}=\{\Omega, \phi, \{\text{Defectuoso}\}, \{\text{No Defectuoso}\}\}$
1. La función de probabilidad:
  \begin{align*}
  p: \mathfrak{F} &\longrightarrow [0,1]\\
     \Omega &\longrightarrow 1\\
     \phi &\longrightarrow 0\\
     \{Defectuoso\}&\longrightarrow P(D)\\
     \{No Defectuoso\}&\longrightarrow 1-P(D)
  \end{align*}
  en donde, acudiendo al teorema de probabilidad total, se define
  \begin{equation*}
  p(D)=p(D \mid M1)P(M1)+p(D \mid M2)P(M2)+p(D \mid M3)P(M3)
  \end{equation*}

Sin embargo, también es posible plantearse otro tipo de preguntas que sirven para calibrar el proceso de producción de artículos defectuosos. Por ejemplo, cabe preguntarse acerca de la probabilidad de que, habiendo seleccionado un artículo defectuoso, éste provenga de la primera máquina^[Por supuesto que la pregunta también es válida al indagar por la probabilidad de que habiendo seleccionado un artículo defectuoso, éste provenga de la segunda o tercera máquina.]. En esta ocasión, el experimento ha cambiado y ahora se trata de seleccionar un artículo defectuoso y para responder a tal cuestionamiento, se debe establecer rigurosamente el espacio de probabilidad que puede estar dado por

1. El espacio muestral: $\Omega=\{M1, M2, M3 \}$
1. La $\sigma$-álgebra: $\mathfrak{F}^+=\{\Omega, \phi, \{M1\}, \{M2,M3\}\}$
1. La función de probabilidad:
  \begin{align*}
  p: \mathfrak{F}^+ &\longrightarrow [0,1]\\
     \Omega &\longrightarrow 1\\
     \phi &\longrightarrow 0\\
     \{M1\}&\longrightarrow p(M1 \mid D)\\
     \{M2,M3\}&\longrightarrow 1-p(M1 \mid D)
  \end{align*}
  en donde, acudiendo a la probabilidad condicional, se define
  \begin{equation*}
  p(M1 \mid D)=\frac{p(D \mid M1)P(M1)}{p(D \mid M1)P(M1)+p(D \mid M2)P(M2)+p(D \mid M3)P(M3)}
  \end{equation*}

La anterior función de probabilidad se conoce con el nombre de regla de probabilidad de Bayes y, aparte de ser el baluarte de la mayoría de investigaciones estadísticas que se plantean hoy en día, ha sido la piedra de tropiezo de muchos investigadores radicales que trataron de estigmatizar este enfoque tildando a sus seguidores de mediocres matemáticos y pobres probabilistas afirmando que la regla de probabilidad de Bayes es sólo un artilugio diseñado para divertirse en el tablero.

Pues bien, la interpretación de la regla de bayes se puede realizar en el sentido de actualización de la estructura probabilística que gobierna el experimento. Y esta actualización tiene mucho sentido práctico cuando se cae en la cuenta de que la vida real está llena de calibradores y que las situaciones generadas son consecuencia de algún cambio estructural. De esta forma, el conocimiento de la probabilidad de que el artículo sea producido por la primera máquina se actualiza al conocer que este artículo particular es defectuoso y de esta manera calibra la estructura aleatoria que existe detrás del contexto de la fábrica de bolígrafos. Aparte de servir para resolver problemas como el anteriormente mencionado, la regla de bayes ha marcado el comienzo de un nuevo enfoque de análisis de datos, no solamente porque hace explícitas las relaciones causales entre los procesos aleatorios, sino también porque facilita la inferencia estadística y la interpretación de los resultados.</div>\EndKnitrBlock{example}
<br>

En el campo de la medicina, también se ha visto un gran número de la
aplicación del teorema de Bayes. A continuación se enuncia uno de ellos:

\BeginKnitrBlock{example}<div class="example"><span class="example" id="exm:unnamed-chunk-12"><strong>(\#exm:unnamed-chunk-12) </strong></span>El Grupo de Trabajo de Servicios Preventivos de los Estados Unidos (USPSTF) hizo unas nuevas y controversiales recomendaciones [recomendaciones](https://www.uspreventiveservicestaskforce.org/uspstf/recommendation/breast-cancer-screening) sobre la detección del cáncer de mama dentro de los cuales no recomienda el examen de la mamografía en mujeres entre 40 y 49 años de edad, afirmando que la práctica bienal de este examen debe ser una decisión individual según el contexto particular de la paciente. Por otro lado, la USPSTF sí recomienda tal práctica de forma bienal en grupos de mujeres de entre 50 y 74 años de edad, puesto que no encontró suficiente evidencia de beneficio o daño adicional en realizar este examen en mujeres mayores a los 74 años. Además, también recomendó *no* realizar auto exámanes de senos, contrario a las recomendaciones y consejos que da la mayoría de los profesionales y organizaciones de la salud, incluyendo la *Amerian Cancer Society*. Como información adicional, se sabe que:
  
* Los expertos estiman que un 12.3% de las mujeres desarrollan formas invasivas del cáncer de mama durante la vida.
* La probabilidad de que una mujer desarrolle el cáncer de mama entre los 40 y los 49 años de edad es 1 en 69, y esta probabilidad aumenta a medida que envejezca, de tal forma que llega a ser de 1 en 38 en mujeres de entre 50 y 59 años.
* El cáncer de mama es más difícil de detectar en mujeres jóvenes puesto que el tejido mamario es más denso y fibroso. Los expertos estiman que la tasa de un falso positivo es de 97.8 por cada 1000 mujeres de 40 y 49 años, y esta tasa disminuye a 86.6 por cada 1000 mujeres entre 50 y 59 años. 
* La tasa de un falso negativo es de 1 por cada 1000 mujeres de 40 y 49 años, y es de 1.1 por cada 1000 mujeres entre 50 y 59 años.

Resumiendo las anteriores afirmaciones, tenemos las siguientes probabilidades

| Probabilidad         | 40 - 49      | 50 - 59 años  |
|----------------------|--------------|---------------|
| Cáncer               | 1/69=0.01449 | 1/38=0.02632  |
| No cáncer            | 68/69=0.9855 | 37/38=0.97368 |
| Positivo $\mid$ No cáncer | 0.0978       | 0.0866        |
| Negativo $\mid$ No cáncer | 0.9022       | 0.9134        |
| Positivo $\mid$ Cáncer    | 0.999        | 0.9989        |
| Negativo $\mid$ Cáncer    | 0.001        | 0.0011        |

Utilizando la regla de Bayes, se puede calcular las siguientes probabilidades para mujeres de 40 y 49 años: 
\begin{align*}
P(\text{Cáncer}|\text{Positivo})&=\frac{P(\text{Positivo}|\text{Cáncer})P(\text{Cáncer})}{P(\text{Positivo}|\text{Cáncer})P(\text{Cáncer})+P(\text{Positivo}|\text{No cáncer})P(\text{No cáncer})}\\
&=\frac{0.999*0.01449}{0.999*0.01449+0.0978*0.9855}\\
&=0.1305
\end{align*}

\begin{align*}
P(\text{Cáncer}|\text{Negativo})&=\frac{P(\text{Negativo}|\text{Cáncer})P(\text{Cáncer})}{P(\text{Negativo}|\text{Cáncer})P(\text{Cáncer})+P(\text{Negativo}|\text{No cáncer})P(\text{No cáncer})}\\
&=\frac{0.001*0.01449}{0.001*0.01449+0.9022*0.9855}\\
&=0.0000163
\end{align*}

Similarmente, se puede calcular estas dos probabilidades para las mujeres de 50 y 59 años.

| Probabilidad         | 40 - 49 años | 50 - 59 años |
|----------------------|--------------|--------------|
| Cáncer $\mid$ Positivo    | 0.1305985    | 0.23769      |
| No cáncer $\mid$ Positivo | 0.8694223    | 0.7623123    |
| Cáncer $\mid$ Negativo    | 0.0000163    | 0.0000326    |
| No cáncer $\mid$ Negativo | 0.9999837    | 0.9999674    |

Los anteriores resultados muestran cómo cambia la probabilidad de tener cáncer al condicionar en los resultados de la pruebe. Entre estos valores se puede ver que, con un resultado positivo en el examen, la probabilidad de tener efectivamente el cáncer es aproximadamente diez puntos porcentuales más bajo en mujeres de edad de 40 y 49 años, de donde se puede sustentar la recomendación de no efectuar este examen en mujeres de este rango de edad.</div>\EndKnitrBlock{example}
<br>
