#   Self-talk interventions and sport/motor performance: An updated systematic review and Bayesian meta-analysis

This is a project where the earlier meta-analysis by [Hatzigeorgiadis et al., (2011)](https://journals.sagepub.com/doi/abs/10.1177/1745691611413136) of self-talk interventions and sports performance was updated using a Bayesian approach.

## Abstract

Self-talk has been researched as an aspect of mental preparation for performance in sports and other motor tasks since the late 1980s. In 2011 Hatzigeorgiadis and colleagues systematically reviewed and meta-analysed three decades of self-talk intervention research including 32 studies. Yet, despite the general proliferation of meta-analyses, this topic has not been meta-analysed in the decade since their review which, at least for the main effect, provided a reasonably precise standardised mean difference estimate (0.48 [95% confidence interval: 0.38, 0.58]). Several further studies on self-talk interventions have been conducted in that time and it is of interest to explore what additional evidence they offer regarding the effects of self-talk interventions on sport/motor performance. Bayesian approaches are well positioned to explore how additional evidence changes our understanding of an effect; to see whether it has changed in sign, magnitude, or precision, or whether further research has largely been a ‘waste’. Therefore, our aim was to conduct an updated systematic review and Bayesian meta-analysis replicating the search, inclusion, and models of Hatzigeorgiadis et al. Informative priors were taken directly from Hatzigeorgiadis et al. A total of 34 studies providing 128 effects nested in 64 groups across experiments 42 were included in the final updated meta-analysis representing data from 18761 participants. The overall posterior pooled estimate for the standardised mean difference was almost exactly the same as the prior: 0.47 [95% quantile interval: 0.39, 0.56]. Bayes factors were calculated for a range of effect sizes and indicated that the included studies largely reflected ‘Weak’ evidence against effects ranging from 0.30 to 0.59, and only provided ‘Decisive’ evidence or greater against more extreme effects: either very small (i.e., <0.02) or large (i.e., >0.81). Results were largely similar for all moderator analyses which were updated too; either providing relatively weak evidence against effects found in the previous meta-analysis or evidence suggesting smaller effects for certain moderators. The findings of our updated Bayesian meta-analyses reiterate the positive effect of self-talk interventions on sport/motor performance. However, they also suggest that cumulatively the past decade and more of research has done little to further our understanding of these effects. Considering the limited resources and time for conducting research, it may be worth moving onto to other more pertinent questions regarding psychological constructs impact upon sport/motor performance.

## Reproducibility

This project uses
[`renv`](https://rstudio.github.io/renv/articles/renv.html#reproducibility):

- `renv::snapshot()` save state
- `renv::restore()` load state

where state refers to package versions used by the project.

## Targets analysis pipeline

This project also uses a function based analysis pipeline using
[`targets`](https://books.ropensci.org/targets/)

Useful console functions:

- `tar_edit()` opens a the make file
- `tar_make()` to run targets
- `tar_visnetwork()` to view pipeline

You can view the current interactive pipeline by clicking [here](https://raw.githack.com/jamessteeleii/self_talk_meta_analysis/main/visnetwork.html).
