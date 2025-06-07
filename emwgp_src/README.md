# Code for "Identifying predictors of clinical outcomes using the projection-predictive feature selection - a proof of concept on the example of Crohn's disease"

Steps to run this code:

1.    Install R, RStudio, and all required R packages (see `session_info.txt` for a list of these, their versions, as well as the R version and the RStudio version used in the original computing environment). For the **cmdstanr** package and its CmdStan system dependency (the CmdStan version used here is also shown in `session_info.txt`), see the installation instructions in the **cmdstanr** ["Getting started"](https://mc-stan.org/cmdstanr/articles/cmdstanr.html) vignette.
2.    Open file `open_me.Rproj` in RStudio.
3.    Run (`source()`) file `descriptive_nums_genchar_outc.R`.
4.    Restart the R session.
5.    Run (`source()`) file `descriptive_nums_cat.R`.
6.    Restart the R session.
7.    Run (`source()`) file `descriptive_figures.R`.
8.    Restart the R session.
9.    Move folder `output` somewhere else.
10.    Run (`source()`) file `ppfs.R`.
11.   Restart the R session.
12.   Run (`source()`) file `comparison.R`.
13.   Restart the R session.
14.   Run (`source()`) file `comparison_smooth.R`.
