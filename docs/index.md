--- 
title: "Modelos Bayesianos con R y STAN"
author: "Andrés Gutiérrez - Hanwen Zhang"
date: "2021-07-03"
documentclass: book
bibliography: [book.bib]
biblio-style: apalike
link-citations: yes
colorlinks: yes
lot: yes
lof: yes
fontsize: 10pt
github-repo: psirusteam/bookdownBayesiano
description: "Este es el repositorio del libro *Modelos Bayesianos con R y STAN*."
knit: "bookdown::render_book"
header-includes:
- \usepackage{setspace}
- \def\bPi{\boldsymbol \Pi }
- \def\bGamma{\boldsymbol \Gamma}
- \def\bgamma{\boldsymbol \gamma}
- \def\beps{\boldsymbol \varepsilon}
- \def\bbeta{\boldsymbol \beta}
- \def\bLambda{\boldsymbol \Lambda}
- \def\blambda{\boldsymbol \lambda}
- \def\bBeta{\boldsymbol \Beta}
- \def\bEta{\boldsymbol \eta}
- \def\balpha{\boldsymbol \alpha}
- \def\btheta{\boldsymbol \theta}
- \def\bmu{\boldsymbol \mu}
- \def\bSigma{\boldsymbol \Sigma}
- \def\bphi{\boldsymbol \phi}
- \def\bpi{\boldsymbol \pi}
- \def\bxi{\boldsymbol \xi}
# - \usepackage[spanish]{babel}
lang: es
linkcolor: blue
# bookdown::render_book("index.Rmd", "bookdown::pdf_book")
# bookdown::render_book("index.Rmd", "bookdown::tufte_html_book")
# bookdown::render_book("index.Rmd", "bookdown::gitbook")
#output: bookdown::html_book
output:
  pdf_document:
    toc: true
    toc_depth: 2
    keep_tex: true
  gitbook:
    df_print: kable
    css: "style.css"
---



# Prefacio {-}
<div class="figure">
<img src="Pics/CClicence.png" alt="Licencia de Creative Commons" width="100px" />
<p class="caption">(\#fig:unnamed-chunk-1)Licencia de Creative Commons</p>
</div>

La versión online de este libro estáa licenciaada bajo una [Licencia Internacinal de Creative Commons para compartir con atribución no comercial 4.0](http://creativecommons.org/licenses/by-nc-sa/4.0/). 
