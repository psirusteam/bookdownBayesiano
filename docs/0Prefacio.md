

# Antes de comenzar {-}

La estadística bayesiana ha tenido un auge bastante importante en las últimas décadas. La mayoría de publicaciones actuales tienen un componente bayesiano importante. Todos las áreas de la ciencia de datos han explotado este paradigma y la evolución de la estadística ha ido en aumento a medida que la revolución bayesiana va tomando forma. 

## Cuestionamientos sobre el enfoque bayesiano {-}

@GelmanObjections presenta algunos de los cuestionamientos que algunos estadísticos anti-bayesianos han argumentado en contra de este paradigma que, sin lugar a dudas, ha proporcionado una valiosa herramienta de modelación en la ciencia contemporanea. Revisemos algunos de estos argumentos: 

> La inferencia bayesiana es una teoría matemática coherente pero no brinda la suficiente confianza en usos científicos. Las distribuciones *previas* subjetivas no inspiran confianza porque ni siquiera existe algún principio objetivo para elegir una distribución previa no informativa. ¿De dónde vienen las distribuciones previas? No confío en ellas y no veo ninguna razón para recomendarlas a otra gente, apenas me siento cómodo acerca de su coherencia filosófica.

Este argumento es débil puesto que la teoría bayesiana es una teoría cinetífica apoyada en los axiomas matemáticos de la teoría de la medida y de probabilidad. De la mismaforma, nótese que tampoco existe un principio objetivo para escoger una verosimilitud. ¿De dónde vienen las regresiones logísticas? ¿quién dijo que los datos eran normales? Como toda ciencia, la estadística se basa en procedimientos subjetivos que inducen resultados que se pueden probar de una manera objetiva. Al decidir usar una determinada distribución previa, el investigador está haciendo uso de su conocimiento objetivo sobre el fenómeno de interés. Esto no dista mucho de la planificación de un estudio por muestreo o de un experimento, en donde se hace uso de la información auxiliar disponible para definir la mejor versión del estudio. Además, como se verá más adelante, sí existen principios objetivos que permiten decidir acerca de la elección de una distribución previa; por ejemplo, la invarianza de la distribución previa frente a transformaciones de los parámetros.

> La teoría bayesiana requiere un pensamiento mucho más profundo sobre la situación y recomendarle a los investigadores comunes el uso del teorema de Bayes es como darle al hijo del vecino la llave de un *F-16*. De veras que, yo comenzaría con algo de métodos probados y confiables, y entonces generalizaría la situación utilizando los principios estadísticos y la teoría del minimax, que no dependen de ninguna creencia subjetiva. Especialmente cuando las distribuciones previas que veo en la práctica toman formas conjugadas. ¡Qué coincidencia!

Como científicos e investigadores debemos tratar con el conocimiento objetivo y dejar a un lado las creencias subjetivas. Es por eso que las distribuciones previas que se manejan en la inferencia bayesiana son objetivas de la misma forma que lo son los métodos frecuentistas al asignar un modelo probabilístico a la verosimilitud de los datos. El resultado final sólo depende del modelo asumido y de los datos recolectados. A pesar de que algunos resultados de la inferencia bayesiana coinciden con el acercamiento frecuentista, esto no sucede en todos los casos. Si la distribución es conjugada, simplemente quiere decir que es posible utilizar un generador de números aleatorios conocido; sin embargo, en pleno siglo XXI, esto ya no constituye un problema.

> Dejando de lado las preocupaciones matemáticas, me gustan las estimaciones insesgadas, los intervalos de confianza con un nivel real de cobertura. Pienso que la manera correcta de inferir es acercarse al parámetro tanto como sea posible y desarrollar métodos robustos que trabajen con supuestos mínimos. El acercamiento bayesiano intenta aproximar el insesgamiento, mientras asume supuestos más y más fuertes. En los viejos tiempos, los métodos Bayesianos por lo menos tenían la virtud de estar matemáticamente limpios. Hoy en día, cualquier inferencia se realiza mediante el uso de las cadenas de Markov con métodos de Monte Carlo (MCMC). Lo anterior significa que, no sólo no se pueden evaluar las características estadísticas del método, sino que tampoco se puede asegurar su convergencia.

Los métodos bayesianos parecen moverse rápidamente hacia la computación elaborada. Para bien o para mal, la computación se está convirtiendo en una plataforma central para el desarrollo científico y estadístico. Por otro lado, estos mismos adelantos de computación científica permiten evaluar las características de los modelos bayesianos y la convergencia de las cadenas de la distribución posterior. Haciendo uso de la rigurosidad científica, el investigador debe conocer a profundidad el espíritu de los métodos MCMC y verificar que la distribución posterior conjunta sobre un vector de parámetros no sea impropia, y por supuesto verificar que las cadenas tienen propiedades estacionarias.

> La gente tiende a creer los resultados que apoyan sus preconceptos y descreen los resultados que los sorprenden, ésta es una forma errada y sesgada de pensar. Pues bien, los métodos bayesianos animan este modo indisciplinado de pensamiento. Estoy seguro que muchos estadísticos bayesianos están actuando de buena fe; sin embargo, al mismo tiempo, también están proporcionando estímulo a investigadores descuidados y poco éticos por todas partes, porque el investigador queda estancado al momento de escoger una distribución previa.

Si hay una seria diferenciación entre las creencias subjetivas y los resultados posteriores, debería ser un indicador de revaluar el modelo usado. Además, ante el desconocimiento del fenómeno, el investigador bayesiano puede utilizar una distribución previa débil y añadir más información si se necesita. Las verificaciones predictivas (previas y posteriores) son una parte esencial del método bayesiano que obliga a repensar las creencias del investigador con respecto al parámetro de interés. Este ejercicio redunda en el replanteamineto de la distribución previa mediante el estudio de las distribuciones predictivas, decantándose al final por el mejor modelo.

> Los cálculos de la teoría de la decisión guían a la idea de que el muestreo probabilístico y la asignación aleatoria de tratamientos son ineficaces, de que los mejores diseños y muestras son los determinísticos. No tengo ningún conflicto con estos cálculos matemáticos; el conflicto es más profundo, en los fundamentos filosóficos, en la idea de que el objetivo de la estadística consiste en tomar una decisión óptima. Un estimador bayesiano es un estimador estadístico que reduce al mínimo el riesgo promedio. Sin embargo, cuando hacemos estadística, no estamos intentando *reducir al mínimo el riesgo promedio*, estamos intentando hacer estimación y juzgamiento de hipótesis.

Un estimador bayesiano es un estimador estadístico que minimiza el riesgo promedio. Uno de los primeros tópicos que se presentan en este libro es el de la teoría de la decisión y funciones de perdida, como herramientas fundamentales del aprendizaje estadístico [@hastie]. Además, como se verá más adelante, la asignación de las unidades experimentales al tratamiento o la inclusión de las unidades muestrales en un estudio probabilístico debe y puede ser tenido en cuenta en los modelos bayesianos, mediante la inclusión en el modelo de las variables que intervinieron en la selección de las unidades. De la misma forma, el juzgamiento de hipótesis es una práctica que se extiende en la modelación bayesiana.

> No puedo estar al tanto de lo que están haciendo todos esos Bayesianos hoy en día. Desafortunadamente, toda clase de personas están siendo seducidas por las promesas de la inferencia automática con la *magia del MCMC*. Desearía que todos paráramos de una vez y por todas y empezáramos, de nuevo, a hacer estadística de la forma en que debe ser hecha: volviendo a los viejos tiempos en que un $p$-valor era utilizado para algo, cuando un intervalo de confianza tenía significado, y el sesgo estadístico era algo que se quería eliminar y no algo que se debiera abrazar.

Los métodos Bayesianos algunas veces son presentados como un motor de inferencia automática. Sin embargo, la inferencia bayesiana tiene tres etapas: formulación del modelo, ajuste del modelo a los datos, evaluación del ajuste. Así que el procedimiento no es mágico ni automático. Además, una de las ventajas de la estadística bayesiana es que deja de lado las sofisticaciones de la inferencia clásica en donde, por ejemplo, la simple interpretación de un intervalo de confianza se hace muy complicada a la luz del razonamiento lógico. De la misma forma los valores $p$ constituyen un paradigma cada vez más revalorado en la investigación social. 

## Acerca de la notación {-}

Antes de empezar las próximas secciones, es necesario revisar la
notación que se seguirá de ahora en adelante. Del teorema de Bayes
resultan tres grandes definiciones que constituyen la base de la
estadística Bayesiana y que a lo largo de este texto se mencionarán
diferenciándolas por medio de la notación. El símbolo más importante de
la estadística matemática es $p$, el cual indica que existe una
distribución de probabilidad para los datos, para el vector de
parámetros, condicional o no. De hecho todos las definiciones y
resultados anteriores han estado supeditadas al uso de esta monótona
notación. En el ámbito de la notación de investigación internacional es
común diferenciar las distribuciones con el fin de hacer más ameno el
estudio del enfoque Bayesiano. En este texto se seguirá esta distinción.
Un ejemplo claro en donde $p$ representa cuatro funciones distintas en
una sola ecuación es el siguiente:

$$p(\theta \mid y)=p(y \mid \theta)\frac{p(\theta)}{p(y)}$$

@Gelman95 explica por qué la notación simple, con el uso (a
veces abuso) de la letra $p$ es más rigurosa de lo que, a simple vista,
pueda parecer y comenta que,

> En realidad no me gusta la notación que la mayoría de los estadísticos usen $f$ para las distribuciones de muestreo; $\pi$, para las distribuciones previas y $L$, para las verosimilitudes. Este estilo de notación se desvía de lo que realmente es importante. La notación no debería depender del orden en que las distribuciones son especificadas. Todas ellas son distribuciones de probabilidad, eso es lo realmente importante.

Esto tiene sentido, aún más cuando se estudian las propiedades
estadísticas de los estimadores desde el punto de vista de la teoría de
la medida. Siendo así, el símbolo $p$ se refiere a una notación para una
medida de probabilidad, quizás inducida por un elemento aleatorio. De
hecho, en la ecuación que determina la regla de Bayes, cada una de las
$p$ son medidas de probabilidad que no comparten el mismo espacio de
medida (ni la misma $\sigma$-álgebra, ni el mismo espacio muestral).

De hecho, todo queda claro al realizar un diagrama que permita ver el
espacio de salida y el espacio de llegada de los elementos aleatorios
que inducen (si es el caso), cada una de las distribuciones de
probabilidad. Por otra parte, Bob Carpenter concluye que:

> Una vez resuelto el problema de identificación de los espacios, la notación estadística depende en gran manera del contexto y aunque la regla de Bayes no necesite de mucha explicación, es necesario conocerlo todo acerca del contexto para poder interpretar las funciones que la conforman... El problema se hace mucho más agudo para los estadísticos novatos, pero eso se resuelve con la práctica. Una vez que uno sabe lo que está haciendo, se vuelve obvia la referencia de la distribución $p$.

Por lo anterior, es natural que algunos de los textos clásicos de
estadística matemática, los autores asumen que el lector sigue la idea
de la referencia de la distribución $p$ en cuestión.
