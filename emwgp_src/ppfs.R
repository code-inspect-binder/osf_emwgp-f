#_______________________________________________________________________________
# Code for "Identifying predictors of clinical outcomes using the
# projection-predictive feature selection - a proof of concept on the example of
# Crohn's disease"
#
# Projection-predictive feature selection (PPFS), with the reference model being
# the full model with a regularized horseshoe (RH) prior
#_______________________________________________________________________________

# User options ------------------------------------------------------------

console_oneFile <- TRUE

voutc <- "colo_score"
outc_nm_pretty <- "Endoscopic score"
outc_adjective <- "endoscopic"
outc_single_noun <- "endoscopy"
outc_nms <- c("remission", "mild", "moderate", "severe")
vpreds <- c(
  "age",
  "weight_gain",
  "condition",
  "activity_limitation",
  "stoolquantcategory",
  "stool_consistency",
  "stool_blood",
  "abdominal_pain",
  "abdominal_pain_night",
  "abdfinding",
  "abdfinding_pressurepain",
  "abdfinding_resistance",
  "analfinding",
  "ext.manifestation_sum",
  "lab_Hb",
  "lab_hk",
  "lab_mcv",
  "thr_log",
  "leuko_log",
  "crp_log",
  "lab_albumin_si",
  "calp_log",
  "cohort"
)
vpreds_PP <- c("patID")

# The prior assumption for the number of relevant coefficients:
p0 <- 5

# Passed to argument `nterms_max` of projpred::cv_varsel():
projpred_nterms_max <- 3
# The number of predictor terms to use when calling projpred::project():
proj_nterms <- 2
# The predictive performance statistic to use for projpred:
proj_evalstats <- c("mlpd")

# Setup -------------------------------------------------------------------

# If necessary, create output folder:
if (!dir.exists("output")) {
  dir.create("output")
} else if (identical(Sys.getenv("CD_rm_output"), "TRUE")) {
  unlink("output", recursive = TRUE)
  dir.create("output")
} else {
  stop("Folder \"output\" already exists. Stopping here to avoid overwriting ",
       "files.")
}
out_folder <- "ppfs"
if (!dir.exists(file.path("output", out_folder))) {
  dir.create(file.path("output", out_folder))
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

# Set seed:
set.seed(284867796)

# Packages ----------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
  library(brms)
  library(bayesplot)
  library(loo)
  library(projpred)
  library(ggplot2)
})

# Package options ---------------------------------------------------------

options(mc.cores = parallel::detectCores(logical = FALSE))

source("helpers/ggplot2_settings_global.R")

bayesplot_theme_set(theme_default(base_size = 12, base_family = "sans") +
                      grid_lines())
color_scheme_set("blue")

options(projpred.mlvl_pred_new = TRUE)

# Internal options --------------------------------------------------------

warn_orig <- options(warn = 1)
loo_forProj <- FALSE
plot_prefix <- ""

# Post-process the user input ---------------------------------------------

if (length(vpreds_PP) > 0L) {
  vpreds_PP <- paste0("(1 | ", vpreds_PP, ")")
}

# Data --------------------------------------------------------------------

C_dat <- readRDS("CD_data.rds")

# Modifications of the dataset which depend on the formula ----------------
# (the formula is implied by the user options above)

## Exclude rows with NAs from the dataset ---------------------------------

cat("\n-----\n")
cat("Number of rows in the dataset before excluding rows with at least one",
    "NA (in any of the columns relevant for the model):\n")
print(nrow(C_dat))
cat("-----\n")
vfit <- c(voutc, vpreds, vpreds_PP)
vfit <- sub(".*\\| ", "", vfit)
vfit <- grep(":", vfit, value = TRUE, invert = TRUE)
vfit <- sub("^s\\(", "", vfit)
vfit <- sub("\\)$", "", vfit)
hasNA <- apply(C_dat[, ..vfit], 1, anyNA)
saveRDS(hasNA, file.path("output", out_folder, "hasNA.rds"))
C_dat <- C_dat[!hasNA, ]
cat("\n-----\n")
cat("Number of rows in the dataset after excluding rows with at least one",
    "NA (in any of the columns relevant for the model):\n")
print(nrow(C_dat))
cat("-----\n")

## For the horseshoe prior: -----------------------------------------------
## Standardize the input variables with population-level effects (with the term
## "input variables" following the definition by gelman_weakly_2008, i.e. the
## input variables are the predictors before adding the intercept, building
## interactions etc.)

### Preparations ----------------------------------------------------------

vpreds_noInt <- grep(":", vpreds, value = TRUE, invert = TRUE)
stopifnot(length(vpreds_noInt) > 0L)

# Create matrix of input variables with population-level effects:
C_datMM <- model.matrix(
  as.formula(paste("~", paste(vpreds_noInt, collapse = " + "))),
  data = C_dat
)
# Handle attributes:
rownames(C_datMM) <- C_dat$obsID
attr(C_datMM, "contrasts") <- NULL

# Remove column "(Intercept)":
C_datMM <- C_datMM[, colnames(C_datMM) != "(Intercept)"]
# The removal of column "(Intercept)" also removes all attributes (including
# `assign`):
stopifnot(identical(names(attributes(C_datMM)), c("dim", "dimnames")))

# Modify the column names:
for (vpreds_i in sort(grep("_c$|cohort", vpreds_noInt, value = TRUE, invert = TRUE),
                      decreasing = TRUE)) {
  colnames(C_datMM) <- sub(paste0("^", vpreds_i, "(.+)"), "\\1", colnames(C_datMM))
}
colnames(C_datMM) <- gsub("\\+", "P", colnames(C_datMM))
colnames(C_datMM) <- gsub(",", ".", colnames(C_datMM))
colnames(C_datMM) <- gsub("\\(", ".", colnames(C_datMM))
colnames(C_datMM) <- gsub("\\[", ".", colnames(C_datMM))
colnames(C_datMM) <- gsub("\\]", ".", colnames(C_datMM))

### Standardization -------------------------------------------------------

C_datMM <- scale(C_datMM, center = TRUE)
colnames(C_datMM) <- paste0("STD_", colnames(C_datMM))

C_centers <- attr(C_datMM, "scaled:center")
saveRDS(C_centers, file.path("output", out_folder, "C_centers.rds"))
C_scales <- attr(C_datMM, "scaled:scale")
saveRDS(C_scales, file.path("output", out_folder, "C_scales.rds"))

# Post-processing:
# The number of population-level coefficients:
K_coefs <- ncol(C_datMM)
# Create a `data.table` and merge with `C_dat`:
C_dat <- merge(C_dat, as.data.table(C_datMM, keep.rownames = "obsID"),
               by = "obsID",
               all = TRUE,
               sort = FALSE)

# Model specification -----------------------------------------------------

## Formula (and distributional family) ------------------------------------

vpreds_STD <- grep("^STD_", names(C_dat), value = TRUE)

C_formula_tmp <- as.formula(paste(
  voutc, "~",
  paste(c(vpreds_STD, vpreds_PP), collapse = " + ")
))

C_bformula <- bf(
  formula = C_formula_tmp,
  family = cumulative(link = "logit", link_disc = "log", threshold = "flexible")
)

C_formula <- formula(C_bformula)
stopifnot(identical(voutc, as.character(C_formula)[2]))
stopifnot(identical(voutc, as.character(C_formula_tmp)[2]))
C_termlabs <- labels(terms(C_formula))
stopifnot(identical(C_termlabs, attr(terms(C_formula), "term.labels")))
stopifnot(identical(C_termlabs, attr(terms(C_formula_tmp), "term.labels")))

## Prior ------------------------------------------------------------------

### Default priors --------------------------------------------------------

cat("\n-----\n")
cat("Default priors:\n")
C_prior_def <- get_prior(C_bformula, data = C_dat)
print(C_prior_def)
cat("-----\n")

### Customized priors -----------------------------------------------------
### Here: Regularized horseshoe prior

# The number of coefficients:
cat("\n-----\n")
cat("Number of coefficients:\n")
print(K_coefs)
cat("-----\n")

# The pre-specified value of p0:
cat("\n-----\n")
cat("p0:\n")
print(p0)
cat("-----\n")

C_prior <- prior(horseshoe(df = 1, par_ratio = p0 / (K_coefs - p0)))

# Model fit ---------------------------------------------------------------

brm_file_arg <- Sys.getenv("CD_brms_file")
if (identical(brm_file_arg, "")) {
  brm_file_arg <- NULL
}
if (packageVersion("cmdstanr") >= "0.5.3") {
  options(cmdstanr_write_stan_file_dir = ".")
} else {
  options(cmdstanr_write_stan_file_dir = getwd())
}
C_bfit <- brm(formula = C_bformula,
              data = C_dat,
              prior = C_prior,
              iter = 4000,
              backend = "cmdstanr",
              adapt_delta = 0.99,
              max_treedepth = 15L,
              save_pars = save_pars(all = TRUE),
              file = brm_file_arg,
              seed = 244762505,
              refresh = 0)

C_bformula_ff <- formula(C_bfit)
stopifnot(identical(formula(C_bformula_ff), C_formula))

cat("\n-----\n")
cat("Priors that were used:\n")
print(prior_summary(C_bfit))
cat("-----\n")

# Save model fit object:
saveRDS(C_bfit, file = file.path("output", out_folder, "C_bfit.rds"))

# Further analyses based on the model -------------------------------------

run_loo <- !identical(Sys.getenv("CD_run_loo"), "FALSE")
source("helpers/analyze_model.R")

out_folder <- file.path(out_folder, "proj")
if (!dir.exists(file.path("output", out_folder))) {
  dir.create(file.path("output", out_folder))
}
plot_prefix <- "PUB_"
source("helpers/projpred.R")
cat("\n_________________________\n")
cat("The following output refers to the selected submodel (onto which the",
    "reference model was projected).\n\n")
loo_forProj <- TRUE
stopifnot(packageVersion("loo") >= package_version("2.5.1"))
source("helpers/analyze_model.R")

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
