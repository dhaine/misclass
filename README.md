<!-- README.md is generated from README.Rmd. Please edit that file -->
Biais de classement dans les études longitudinales de la santé mammaire
=======================================================================

(Misclassification Bias in Longitudinal Udder Health Studies)
=============================================================

[![Build Status](https://travis-ci.org/dhaine/misclass.svg?branch=master)](https://travis-ci.org/dhaine/misclass)

R package to assess misclassification bias in hierarchical longitudinal studies and the effect of various sampling strategies to control for it. This package is part of the scientific paper "Controlling Misclassification Bias in Longitudinal Udder Health Studies" (submitted to [Preventive Veterinary Medicine](https://www.journals.elsevier.com/preventive-veterinary-medicine)) and the [2017 SVEPM](http://www.svepm2017.org/) proceedings "Sampling Strategies to Control Misclassification Bias in Longitudinal Udder Health Studies" by Denis Haine, Ian Dohoo, Daniel Scholl, Henryk Stryhn, and Simon Dufour.

License
-------

This package is free and open source software, licensed under GPL-2.

Abstract
--------

Using imperfect tests may lead to biased estimates of disease frequency and of association with exposure. For instance in longitudinal udder health studies, both quarters at risk and incident intramammary infections (IMI) can be wrongly identified, resulting in selection and misclassification bias, respectively. Diagnostic accuracy can possibly be improved by using duplicate or triplicate samples for identifying quarters at risk and incident IMI.

The objectives of this study were to evaluate the relative impact of selection and misclassification biases resulting from IMI misclassification on measures of disease frequency (incidence) and of association with hypothetical exposures. The effect of improving the sampling strategy by collecting duplicate or triplicate samples at first or second sampling was also assessed.

Data sets from a hypothetical cohort study were simulated and analyzed based on a separate scenario for two common mastitis pathogens representing two distinct prevailing patterns. *Staphylococcus aureus*, a relatively uncommon pathogen with a low incidence, is identified with excellent sensitivity and almost perfect specificity. Coagulase negative staphylococci (CNS) are more prevalent, with a high incidence, and with milk bacteriological culture having fair Se but excellent Sp. The generated data sets for each scenario were emulating a longitudinal cohort study with two milk samples collected one month apart from each quarter of a random sample of 30 cows/herd, from 100 herds, with a herd-level exposure having a known strength of association. Incidence of IMI and measure of association with exposure (odds ratio; OR) were estimated using Markov Chain Monte Carlo (MCMC).

For *S. aureus* biases were small with an observed incidence of 0.29 versus a true incidence of 0.25 IMI/100 quarter-month. In the CNS scenario, diagnostic errors in the two samples led to important selection (40 IMI/100 quarter-month) and misclassification (23 IMI/100 quarter-month) biases for estimation of IMI incidence, respectively. These biases were in opposite direction and therefore the incidence measure obtained using single sampling on both the first and second test (29 IMI/100 quarter-month) was exactly the true value.

In the *S. aureus* scenario the OR for exposure showed little bias (observed OR of 3.1 versus true OR of 3.2). The CNS scenario revealed the presence of a large misclassification bias moving the association towards the null value (OR of 1.7 versus true OR of 2.6). Little improvement could be brought using different sampling strategies aiming at improving Se and/or Sp on first and/or second sampling and using a two out of three interpretation for IMI definition.

Increasing number of samples or tests can prevent bias in some situations but efforts can be spared by holding to a single sampling approach in others. When designing longitudinal studies, evaluating potential biases and best sampling strategy is as critical as the choice of test.

Installation
------------

You can get the most current version from **GitHub** with **devtools** package:

``` r
devtools::install_github('dhaine/misclass')
```
