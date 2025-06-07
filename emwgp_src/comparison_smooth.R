#_______________________________________________________________________________
# Code for "Identifying predictors of clinical outcomes using the
# projection-predictive feature selection - a proof of concept on the example of
# Crohn's disease"
#
# Comparison of the predictive performance of the new (disease) activity index
# vs. that of existing (disease) activity indices, but allowing for nonlinear
# effects in the existing activity indices' ordinal regression models
#_______________________________________________________________________________

# User options ------------------------------------------------------------

console_oneFile <- TRUE
seed_overall <- 463684667

voutc <- "colo_score"
outc_single_noun <- "endoscopy"
vactidx <- c(
  "PCDAI",
  "abbrPCDAI",
  "modPCDAI",
  "shPCDAI",
  "wPCDAI",
  "MINI"
)
newactidx <- "SESM"

# Setup -------------------------------------------------------------------

# If necessary, create output folder:
if (!dir.exists("output")) dir.create("output")
out_folder <- "comparison_smooth"
if (!dir.exists(file.path("output", out_folder))) {
  dir.create(file.path("output", out_folder))
} else if (identical(Sys.getenv("CD_rm_output"), "TRUE")) {
  unlink(file.path("output", out_folder), recursive = TRUE)
  dir.create(file.path("output", out_folder))
} else {
  stop("Folder \"", file.path("output", out_folder), "\" already exists. ",
       "Stopping here to avoid overwriting files.")
}

# Start writing stdout (console output) and possibly stderr (errors, warnings,
# messages) to file(s):
sink_stderr <- !identical(Sys.getenv("CD_sink_stderr"), "FALSE")
if (isTRUE(console_oneFile)) {
  if (file.exists(file.path("output", out_folder, "console_stdout_stderr.txt"))) {
    file.remove(file.path("output", out_folder, "console_stdout_stderr.txt"))
  }
  console_con <- file(file.path("output", out_folder, "console_stdout_stderr.txt"),
                      open = "a")
  sink(console_con, type = "output", append = TRUE, split = TRUE)
  if (isTRUE(sink_stderr)) {
    sink(console_con, type = "message")
  }
} else {
  if (file.exists(file.path("output", out_folder, "console_stdout.txt"))) {
    file.remove(file.path("output", out_folder, "console_stdout.txt"))
  }
  if (file.exists(file.path("output", out_folder, "console_stderr.txt"))) {
    file.remove(file.path("output", out_folder, "console_stderr.txt"))
  }
  sink(file.path("output", out_folder, "console_stdout.txt"), type = "output",
       split = TRUE)
  if (isTRUE(sink_stderr)) {
    console_con_stderr <- file(file.path("output", out_folder, "console_stderr.txt"),
                               open = "wt")
    sink(console_con_stderr, type = "message")
  }
}

# Timestamp:
cat("\n-----\n")
cat("Timestamp at the start of the script:\n")
print(Sys.time())
cat("-----\n")

# Packages ----------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
  library(brms)
  library(ggplot2)
})

# Package options ---------------------------------------------------------

options(mc.cores = parallel::detectCores(logical = FALSE))

source("helpers/ggplot2_settings_global.R")

# Internal options --------------------------------------------------------

warn_orig <- options(warn = 1)
plot_prefix <- "PUB_"

# Post-process the user input ---------------------------------------------

stopifnot(is.character(vactidx), is.vector(vactidx), length(vactidx) > 0L)
vactidx <- setNames(nm = vactidx)

# Data --------------------------------------------------------------------

C_dat_orig <- readRDS("CD_data.rds")

cat("\n-----\n")
cat("Number of rows in the original dataset:\n")
print(nrow(C_dat_orig))
cat("-----\n")

#_______________________________________________________________________________
# Activity indices
#_______________________________________________________________________________

cat("\n")
cat("_______________________________________________________\n",
    "Separate analysis of each activity index:\n", sep = "")
cat("\n")

# New activity index ------------------------------------------------------

cat("\n")
cat("..............\nActivity index: ", newactidx, sep = "")
cat("\n")

## Missing values (NAs) ---------------------------------------------------

# Indicator for observations which were excluded for the new activity index:
N_hasNA <- readRDS(file.path("output", "ppfs", "hasNA.rds"))

cat("\n")
cat("Number of rows in the dataset used for ", newactidx, ":\n", sep = "")
print(sum(!N_hasNA))
cat("\n")

## MCMC diagnostics -------------------------------------------------------

N_bfit <- readRDS(file.path("output", "ppfs", "C_bfit.rds"))
source("helpers/check_MCMC_diagn.R")
N_MCMC_diagn <- check_MCMC_diagn(C_stanfit = N_bfit$fit, pars = "disc",
                                 include = FALSE)
stopifnot(N_MCMC_diagn$all_OK)
N_pfit <- readRDS(file.path("output", "ppfs", "proj", "C_pfit.rds"))
N_MCMC_diagn <- check_MCMC_diagn(C_stanfit = N_pfit$fit, pars = "disc",
                                 include = FALSE)
stopifnot(N_MCMC_diagn$all_OK)

## Predictive performance of the new activity index -----------------------

# Read-in CV-PPs:
N_cvpps <- readRDS(file.path("output", "ppfs", "proj", "PPs_from_vsel.rds"))

# Expand to the size of the original dataset:
N_cvpps_orig <- rep(NA_integer_, nrow(C_dat_orig))
stopifnot(identical(length(N_cvpps), sum(!N_hasNA)))
N_cvpps_orig[!N_hasNA] <- N_cvpps
rm(N_cvpps)

cat("..............\n")

#_______________________________________________________________________________
# Loop over each existing activity index and compare with the new activity index
#_______________________________________________________________________________

# The first of the existing activity indices ------------------------------
# (this one is handled separately to avoid model recompilations later)

vactidx_i <- vactidx[1]

cat("\n")
cat("..............\nActivity index: ", vactidx_i, sep = "")
cat("\n")

## Exclude rows with NAs from the dataset ---------------------------------

vfit <- c(voutc, vactidx_i)
vfit <- sub(".*\\| ", "", vfit)
vfit <- grep(":", vfit, value = TRUE, invert = TRUE)
vfit <- sub("^s\\(", "", vfit)
vfit <- sub("\\)$", "", vfit)
C_hasNA <- apply(C_dat_orig[, ..vfit], 1, anyNA)
C_dat <- C_dat_orig[!C_hasNA, ]
cat("\n")
cat("Number of rows in the dataset after excluding rows with at least one NA:\n")
print(nrow(C_dat))
cat("\n")

## Formula and family specification ---------------------------------------

# In our case, setting `k = 4` and `bs = 'cr'` in the mgcv::s() smooth term is
# not overly restrictive (checked by running the frequentist mgcv::gam()
# function for each existing activity index in turn, thereby varying arguments
# `k` and `bs`, and then checking the summary() of the resulting fits), but in a
# future (prospective) study, the mgcv::s() settings need to be checked again:
C_formula <- as.formula(paste(voutc, "~", paste0("s(", vactidx_i, ", k = 4, bs = 'cr')")))
C_family <- cumulative(link = "logit", link_disc = "log", threshold = "flexible")

## Prior specification ----------------------------------------------------

# In our case, the prior for the smooth-term SD needs to be made tighter (i.e.,
# more mass towards zero and a lighter tail compared to the default
# `student_t(3, 0, 2.5)` prior for `sds`) to avoid computational problems, but
# in a future (prospective) study, this might not be necessary:
C_prior <- set_prior("normal(0, 1)", class = "sds")

## Model fit --------------------------------------------------------------

cat("\n")
cat("Prior:\n")
print(C_prior)
cat("\n")
cat("\n")
cat("Seed: ", seed_overall, sep = "")
cat("\n")

if (packageVersion("cmdstanr") >= "0.5.3") {
  options(cmdstanr_write_stan_file_dir = ".")
} else {
  options(cmdstanr_write_stan_file_dir = getwd())
}
C_bfit <- brm(formula = C_formula,
              data = C_dat,
              family = C_family,
              prior = C_prior,
              iter = 4000,
              backend = "cmdstanr",
              adapt_delta = 0.99,
              max_treedepth = 15L,
              seed = seed_overall,
              refresh = 0)

## Check possibly data-dependent priors -----------------------------------
## (not really needed here (only when calling brms:::update.brmsfit() later),
## but do this nevertheless just to be warned if there are data-dependent
## priors)

C_prior_smmry <- prior_summary(C_bfit)
stopifnot(identical(
  C_prior_smmry$prior[C_prior_smmry$class == "Intercept" & C_prior_smmry$coef == ""],
  "student_t(3, 0, 2.5)"
))
stopifnot(identical(
  C_prior_smmry$prior[C_prior_smmry$class == "b" & C_prior_smmry$coef == ""],
  ""
))

## MCMC diagnostics -------------------------------------------------------

C_MCMC_diagn <- check_MCMC_diagn(C_stanfit = C_bfit$fit, pars = "disc",
                                 include = FALSE)

## Cross-validation (CV) --------------------------------------------------

seed_overall <<- c(seed_overall[1] - 1L, seed_overall)
set.seed(seed_overall[1])
seed_overall <<- c(seed_overall[1] - 1L, seed_overall)
C_kfold <- kfold(C_bfit, K = 25, seed = seed_overall[1])

# Extract CV-PPs:
C_cvpps <- exp(C_kfold$pointwise[, "elpd_kfold"])

# Expand to the size of the original dataset:
C_cvpps_orig <- rep(NA_integer_, nrow(C_dat_orig))
stopifnot(identical(length(C_cvpps), sum(!C_hasNA)))
C_cvpps_orig[!C_hasNA] <- C_cvpps
rm(C_cvpps)

## In case of problematic MCMC diagnostics: Save certain information ------

probl_diagn <- list()
if (!C_MCMC_diagn$all_OK) {
  probl_diagn <- c(
    probl_diagn,
    setNames(list(
      list("model_formula" = C_formula, "diagn" = C_MCMC_diagn,
           "kfold_res" = C_kfold)
    ), vactidx_i)
  )
}

## Prepare output ---------------------------------------------------------

C_res <- setNames(list(list("U_vactidx" = vactidx_i,
                            "U_bfit" = C_bfit,
                            "U_cvpps_orig" = C_cvpps_orig,
                            "U_prior" = C_prior)),
                  vactidx_i)

cat("..............\n")

# All remaining existing activity indices ---------------------------------

## Function for fitting a single updated model ----------------------------

U_i <- function(vactidx_i, P_formula, P_prior, P_bfit) {
  # Updated formula:
  U_formula <- update(P_formula,
                      paste(".", "~", paste0("s(", vactidx_i, ", k = 4, bs = 'cr')")))
  cat("\n")
  cat("..............\nActivity index: ", vactidx_i, sep = "")
  cat("\n")

  # Updated dataset:
  # Exclude rows with NAs from the dataset:
  vfit <- c(voutc, vactidx_i)
  vfit <- sub(".*\\| ", "", vfit)
  vfit <- grep(":", vfit, value = TRUE, invert = TRUE)
  vfit <- sub("^s\\(", "", vfit)
  vfit <- sub("\\)$", "", vfit)
  U_hasNA <- apply(C_dat_orig[, ..vfit], 1, anyNA)
  U_dat <- C_dat_orig[!U_hasNA, ]
  cat("\n")
  cat("Number of rows in the dataset after excluding rows with at least one NA:\n")
  print(nrow(U_dat))
  cat("\n")

  # Updated prior:
  U_prior <- P_prior
  cat("\n")
  cat("Prior:\n")
  print(U_prior)
  cat("\n")

  # Updated seed:
  seed_overall <<- c(seed_overall, seed_overall[length(seed_overall)] + 1L)
  cat("\n")
  cat("Seed: ", seed_overall[length(seed_overall)], sep = "")
  cat("\n")

  # Fit updated model:
  U_bfit <- update(P_bfit,
                   formula. = U_formula,
                   newdata = U_dat,
                   prior = U_prior,
                   seed = seed_overall[length(seed_overall)])

  # Check possibly data-dependent priors (because the problem is that
  # brms:::update.brmsfit() doesn't update data-dependent default priors):
  U_prior_smmry <- prior_summary(U_bfit)
  stopifnot(identical(
    U_prior_smmry$prior[U_prior_smmry$class == "Intercept" & U_prior_smmry$coef == ""],
    "student_t(3, 0, 2.5)"
  ))
  stopifnot(identical(
    U_prior_smmry$prior[U_prior_smmry$class == "b" & U_prior_smmry$coef == ""],
    ""
  ))

  # Updated MCMC diagnostics:
  U_MCMC_diagn <- check_MCMC_diagn(C_stanfit = U_bfit$fit, pars = "disc",
                                   include = FALSE)

  # Updated CV:
  seed_overall <<- c(seed_overall[1] - 1L, seed_overall)
  set.seed(seed_overall[1])
  seed_overall <<- c(seed_overall[1] - 1L, seed_overall)
  U_kfold <- kfold(U_bfit, K = 25, seed = seed_overall[1])

  U_cvpps <- exp(U_kfold$pointwise[, "elpd_kfold"])

  U_cvpps_orig <- rep(NA_integer_, nrow(C_dat_orig))
  stopifnot(identical(length(U_cvpps), sum(!U_hasNA)))
  U_cvpps_orig[!U_hasNA] <- U_cvpps
  rm(U_cvpps)

  # In case of problematic MCMC diagnostics: Save certain information:
  if (!U_MCMC_diagn$all_OK) {
    probl_diagn <<- c(
      probl_diagn,
      setNames(list(
        list("model_formula" = U_formula, "diagn" = U_MCMC_diagn,
             "kfold_res" = U_kfold)
      ), vactidx_i)
    )
  }

  cat("..............\n")

  return(list("U_vactidx" = vactidx_i,
              "U_bfit" = U_bfit,
              "U_cvpps_orig" = U_cvpps_orig,
              "U_prior" = U_prior))
}

## Loop -------------------------------------------------------------------

# Fit updated models:
U_res <- lapply(vactidx[-1], U_i, P_formula = C_formula, P_prior = C_prior,
                P_bfit = C_bfit)
U_res <- c(C_res, U_res)
rm(C_res)

# Comparison --------------------------------------------------------------

cat("\n")
cat("_______________________________________________________\n",
    "Comparison of the new activity index to the existing activity indices.\n",
    "NOTE:\n",
    "For the pairwise comparisons, only those observations are used for which\n",
    "it was possible to calculate *both* of the compared activity indices.\n",
    sep = "")
cat("\n")

## Separate evaluation of each activity index -----------------------------

U_cvpps_sep <- c(
  lapply(U_res, "[[", "U_cvpps_orig"),
  setNames(list(N_cvpps_orig), newactidx)
)
U_Nobs_sep <- lapply(
  U_cvpps_sep,
  function(lx) {
    sum(!is.na(lx))
  }
)
setDT(U_cvpps_sep)
setDT(U_Nobs_sep)
U_cvpps_sep_l <- melt(
  U_cvpps_sep,
  measure.vars = names(U_cvpps_sep),
  variable.name = "act_idx",
  value.name = "cvpp"
)
U_Nobs_sep_l <- melt(U_Nobs_sep,
                     measure.vars = names(U_Nobs_sep),
                     variable.name = "act_idx",
                     value.name = "N_obs")
U_Nobs_sep_l[, N_obs_char := paste("N =", N_obs)]
seed_overall <- c(seed_overall, seed_overall[length(seed_overall)] + 1L)
set.seed(seed_overall[length(seed_overall)])

ggprefixed1 <- ggplot(data = U_cvpps_sep_l,
                      mapping = aes(x = act_idx, y = cvpp)) +
  geom_boxplot(outlier.shape = NA, color = "gray50") +
  geom_jitter(width = 0.2, height = 0, alpha = 0.4) +
  xlab("Disease activity index") +
  ylab(paste("Predictive probability for observed", outc_single_noun)) +
  scale_y_continuous(labels = scales::label_percent()) +
  coord_cartesian(ylim = c(0, 1)) +
  scale_color_identity() +
  geom_label(aes(y = 1, label = N_obs_char), data = U_Nobs_sep_l)
print(ggprefixed1)

## Pairwise comparisons ---------------------------------------------------
## (always the new activity index vs. an existing activity index)

U_cvpps_paired <- U_cvpps_sep[
  ,
  lapply(.SD, function(x) {
    get(newactidx) - x
  }),
  .SDcols = unname(vactidx)
]
U_Nobs_paired <- U_cvpps_paired[
  ,
  lapply(.SD, function(lx) {
    sum(!is.na(lx))
  }),
  .SDcols = unname(vactidx)
]
U_cvpps_paired_l <- melt(U_cvpps_paired,
                         measure.vars = unname(vactidx),
                         variable.name = "act_idx",
                         value.name = "cvpp_diff")
U_Nobs_paired_l <- melt(U_Nobs_paired,
                        measure.vars = names(U_Nobs_paired),
                        variable.name = "act_idx",
                        value.name = "N_obs")
U_Nobs_paired_l[, N_obs_char := paste("N =", N_obs)]
seed_overall <- c(seed_overall, seed_overall[length(seed_overall)] + 1L)
set.seed(seed_overall[length(seed_overall)])
pretty_index <- "indices"
if (length(vactidx) == 1) {
  pretty_index <- "index"
}
saveRDS(
  U_cvpps_paired_l,
  file = file.path("output", out_folder, "U_cvpps_paired_l.rds")
)
ggprefixed2 <- ggplot(data = U_cvpps_paired_l,
                      mapping = aes(x = act_idx, y = cvpp_diff)) +
  geom_boxplot(outlier.shape = NA, color = "gray50") +
  geom_jitter(width = 0.2, height = 0, alpha = 0.4) +
  xlab("Disease activity index") +
  ylab(paste0("Difference (of ", newactidx, " vs. shown ", pretty_index, ") ",
              "of\npredictive probability for observed ", outc_single_noun)) +
  scale_y_continuous(labels = scales::label_percent()) +
  scale_color_identity() +
  geom_label(aes(y = U_cvpps_paired_l[, max(cvpp_diff, na.rm = TRUE) + 0.05],
                 label = N_obs_char),
             data = U_Nobs_paired_l)
print(ggprefixed2)

## Patchworked plot -------------------------------------------------------

library(patchwork)
print((ggprefixed1 / ggprefixed2) + plot_annotation(tag_levels = "A"))
ggsave(file.path("output", out_folder, paste0(plot_prefix,
                                              "all_CV-PPs_patchwork.jpeg")),
       width = 7, height = 2 * 7 * 0.618)
saveRDS(
  last_plot(),
  file = file.path("output", out_folder, paste0(plot_prefix,
                                                "all_CV-PPs_patchwork.rds"))
)

# Output ------------------------------------------------------------------

stopifnot(!any(duplicated(seed_overall)))
saveRDS(seed_overall, file.path("output", out_folder, "seeds_used.rds"))
cat("\n")
cat("_______________________________________________________\n",
    "Summary of problematic MCMC diagnostics:\n", sep = "")
cat("\n")
print(names(probl_diagn))
cat("\n")
saveRDS(probl_diagn, file.path("output", out_folder, "probl_MCMC_diagn.rds"))

saveRDS(U_res, file.path("output", out_folder, "U_res.rds"))

# Teardown ----------------------------------------------------------------

# Reset global options:
options(warn_orig)

# Timestamp:
cat("\n-----\n")
cat("Timestamp at the end of the script:\n")
print(Sys.time())
cat("-----\n")

# Stop writing stdout (console output) and possibly stderr (errors, warnings,
# messages) to file(s):
if (isTRUE(sink_stderr)) {
  sink(type = "message")
}
sink(type = "output")
