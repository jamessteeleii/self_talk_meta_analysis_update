#   Cumulative evidence synthesis and consideration of "research waste" using Bayesian methods: An example updating a previous meta-analysis of self-talk interventions for sport/motor performance

## Abstract
  In the present paper we demonstrate the application of methods for cumulative evidence synthesis including Bayesian meta-analysis, and exploration of questionable research practices such as publication bias or *p*-hacking, in the sport and exercise sciences for the evaluation of experimental interventions. The use of such methods can aid in study planning and avoid "research waste". In demonstrating and discussing these methods we use the example of self-talk interventions and their effects upon sport/motor performance given a quantitative evidence synthesis has not been conducted on this topic, to the best of our knowledge, since 2011 when [Hatzigeorgiadis et al., (2011)](https://journals.sagepub.com/doi/abs/10.1177/1745691611413136) conducted their systematic review and meta-analysis. As such, this topic is ripe to use in demonstrating cumulative methods such as Bayesian updating. Therefore, our aim was to conduct an updated systematic review and Bayesian meta-analysis replicating the search, inclusion, and models of Hatzigeorgiadis et al. (2011) and demonstrate the application of cumulative evidence synthesis methods including; consideration of the initial probability that a new study of the effects of self-talk interventions would shift our prior belief in their effectiveness, the application of priors taken from the previous meta-analysis to be updated by new studies identified to a new posterior estimate of effect, and consideration of other possible sources of research waste from questionable research practices such as publication bias and *p*-hacking. Such methods as those demonstrated here, when used prospectively, can aid researchers in determining whether further research of a particular experimental intervention is in fact warranted. Considering the limited resources and time for conducting research we hope that highlighting the application of these methods might help researchers in the field to avoid research waste and more productively direct their research efforts.


## Reproducibility
This repository contains the necessary files and code to reproduce the analyses, figures, and the manuscript. 

## Usage
To reproduce the analyses, you will need to have R (https://cran.r-project.org/) and RStudio (https://www.rstudio.com/products/rstudio/download/#download) installed on your computer.

To help with reproducibility, this project uses the `renv` R package (see https://rstudio.github.io/renv/articles/renv.html). With `renv`, the state of this R project can be easily loaded as `renv` keeps track of the required R packages (including version), and (if known) the external source from which packages were retrieved (e.g., CRAN, Github). With `renv`, packages are installed to a project specific library rather than your user or system library. The `renv` package must be installed on your machine before being able to benefit from its features. The package can be installed using the following command:

``` r
install.packages("renv")
```

Once you have `renv` installed, you can get a copy of this repository on your machine by clicking the green Code button then choose Download zip. Save to your machine and extract. After extraction, double click the `self_talk_meta_analysis_update.Rproj` file in the root directory. This will automatically open RStudio. This will ensure all paths work on your system as the working directory will be set to the location of the `.Rproj` file. Upon opening, RStudio will recognize the `renv` files and you will be informed that the project library is out of sync with the lockfile. At shown in the console pane of RStudio, running `renv::restore()` will install the packages recorded in the lockfile. This could take some time depending on your machine and internet connection.

## Targets analysis pipeline

This project also uses a function based analysis pipeline using
[`targets`](https://books.ropensci.org/targets/). Instead of script based pipelines the `targets` package makes use of functions applied to targets specified within the pipeline. The targets can be viewed in the `_targets.R` file, and any user defined functions are available in `R/functions.r`.

You can view the existing targets pipeline by downloading the `targets_pipeline.html` file and opening it in your browser.

Useful console functions:

- `tar_edit()` opens the make file
- `tar_make()` to run targets
- `tar_visnetwork()` to view pipeline

Shield: [![CC BY-NC-SA 4.0][cc-by-nc-sa-shield]][cc-by-nc-sa]

This work is licensed under a
[Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License][cc-by-nc-sa].

[![CC BY-NC-SA 4.0][cc-by-nc-sa-image]][cc-by-nc-sa]

[cc-by-nc-sa]: http://creativecommons.org/licenses/by-nc-sa/4.0/
  [cc-by-nc-sa-image]: https://licensebuttons.net/l/by-nc-sa/4.0/88x31.png
[cc-by-nc-sa-shield]: https://img.shields.io/badge/License-CC%20BY--NC--SA%204.0-lightgrey.svg

