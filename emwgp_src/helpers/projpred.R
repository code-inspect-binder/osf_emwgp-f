#_______________________________________________________________________________
# Code for "Identifying predictors of clinical outcomes using the
# projection-predictive feature selection - a proof of concept on the example of
# Crohn's disease"
#
# Helper file for the projection-predictive feature selection
#_______________________________________________________________________________

# Reference model object --------------------------------------------------

C_refModel <- get_refmodel(C_bfit, brms_seed = 739811025, latent = TRUE)

# Parallelization setup ---------------------------------------------------

ncores <- min(parallel::detectCores(logical = FALSE), 5L)
doParallel::registerDoParallel(ncores)
trigger_default <- options(projpred.prll_prj_trigger = 200)

# Run cv_varsel() ---------------------------------------------------------

C_cvvs <- cv_varsel(C_refModel,
                    method = "forward",
                    cv_method = "kfold",
                    K = 25,
                    ndraws_pred = n_draws %/% 3,
                    nterms_max = projpred_nterms_max,
                    parallel = TRUE,
                    seed = 1825482833)
saveRDS(C_cvvs, file.path("output", out_folder, "C_cvvs.rds"))

# Predictive performance plot ---------------------------------------------
# (= model size selection plot)

stat_pretty <- setNames(nm = proj_evalstats)
stat_pretty <- toupper(stat_pretty)
stopifnot(identical(proj_evalstats, "mlpd"))
ggeval <- plot(C_cvvs, stats = proj_evalstats, deltas = TRUE,
               ranking_nterms_max = NA)
ggeval <- ggeval + facet_null()
ggeval <- ggeval + scale_y_continuous(
  sec.axis = sec_axis(~ exp(.), name = bquote(Lambda*" "*"GMPD"))
)
ggeval <- ggeval + labs(y = bquote(Delta*" "*.(stat_pretty)))
print(ggeval)
ggsave(file.path("output", out_folder,
                 paste0(plot_prefix, "projpred_search_deltas.jpeg")),
       width = 7, height = 7 * 0.618)
saveRDS(last_plot(),
        file = file.path("output", out_folder,
                         paste0(plot_prefix, "projpred_search_deltas.rds")))

# Summary -----------------------------------------------------------------

for (deltas_i in c(TRUE, FALSE)) {
  if (deltas_i) {
    deltas_pretty <- " "
  } else {
    deltas_pretty <- " not "
  }
  cat("\n-----\n")
  cat("Summary of projpred results (using 'deltas = ", deltas_i,
      "', so", deltas_pretty, "based on the difference vs. ",
      "the baseline model):\n", sep = "")
  smmry_out <- summary(C_cvvs,
                       stats = c("elpd", proj_evalstats),
                       deltas = deltas_i,
                       type = c("mean", "se", "lower", "upper", "diff", "diff.se"))
  print(smmry_out, digits = 3)
  cat("-----\n")
  write.csv2(
    smmry_out$selection,
    file = file.path("output", out_folder,
                     paste0(plot_prefix, "projpred_search_deltas",
                            substr(deltas_i, 1, 1), "_smmry.csv")),
    row.names = FALSE
  )
}

# Predictor ranking(s) ----------------------------------------------------

rk <- ranking(C_cvvs)
cat("\n-----\n")
cat("Full-data predictor ranking (from projpred):\n", sep = "")
print(rk[["fulldata"]])
cat("-----\n")
cat("\n-----\n")
cat("CV ranking proportions:\n", sep = "")
pr_rk <- cv_proportions(rk)
print(pr_rk)
cat("-----\n")
write.csv2(
  pr_rk,
  file = file.path("output", out_folder,
                   paste0(plot_prefix, "projpred_pr_rk.csv")),
  row.names = FALSE
)

# Predictive probabilities ------------------------------------------------

# Needed for the manual extraction of the predictive probabilities:
stopifnot(length(unique(C_cvvs$summaries$sub[[1 + proj_nterms]]$oscale$wcv)) == 1)

# Manual extraction of the predictive probabilities:
logPPs_from_vsel <- C_cvvs$summaries$sub[[1 + proj_nterms]]$oscale$lppd
PPs_from_vsel <- exp(logPPs_from_vsel)
cat("\n-----\n")
cat("The MLPD at size `proj_nterms = ", proj_nterms,
    "` (should coincide with that from the summary() output):\n", sep = "")
print(mean(logPPs_from_vsel))
cat("-----\n")
cat("\n-----\n")
cat("The GMPD at size `proj_nterms = ", proj_nterms, "`:\n", sep = "")
print(exp(mean(logPPs_from_vsel)))
cat("-----\n")
saveRDS(PPs_from_vsel,
        file = file.path("output", out_folder, "PPs_from_vsel.rds"))

# Projection onto selected (i.e., final) submodel -------------------------

C_predictors_final <- rk[["fulldata"]][seq_len(proj_nterms)]
C_proj <- project(C_refModel, solution_terms = C_predictors_final, ndraws = n_draws)
C_proj_draws_mat <- as.matrix(C_proj)
saveRDS(C_proj_draws_mat, file.path("output", out_folder, "C_proj_draws_mat.rds"))
C_proj_pars <- colnames(C_proj_draws_mat)

# Update `C_centers` and `C_scales`:
C_centers <- C_centers[sub("^b_STD_", "", grep("^b_STD_", C_proj_pars, value = TRUE))]
saveRDS(C_centers, file.path("output", out_folder, "C_centers.rds"))
C_scales <- C_scales[sub("^b_STD_", "", grep("^b_STD_", C_proj_pars, value = TRUE))]
saveRDS(C_scales, file.path("output", out_folder, "C_scales.rds"))

# Create a "stanfit" object from the projected draws:
stopifnot("b_Intercept" %in% C_proj_pars)
C_proj_draws_mat <- cbind(
  C_draws_mat[, grep("^b_Intercept", C_pars), drop = FALSE] -
    C_proj_draws_mat[, "b_Intercept"],
  C_proj_draws_mat[, grep("^b_STD_", C_proj_pars), drop = FALSE],
  C_proj_draws_mat[, grep("^sd_patID__Intercept$", C_proj_pars), drop = FALSE],
  C_draws_mat[, "disc", drop = FALSE],
  C_proj_draws_mat[, grep("^r_patID", C_proj_pars), drop = FALSE]
)
C_proj_pars <- colnames(C_proj_draws_mat)
colnames(C_proj_draws_mat) <- gsub(",", "COMMA", C_proj_pars)
C_proj_pars <- colnames(C_proj_draws_mat)
colnames(C_proj_draws_mat) <- gsub("\\[", "SQBL", C_proj_pars)
C_proj_pars <- colnames(C_proj_draws_mat)
colnames(C_proj_draws_mat) <- gsub("\\]", "SQBR", C_proj_pars)
C_proj_pars <- colnames(C_proj_draws_mat)
colnames(C_proj_draws_mat) <- gsub("\\.", "DOT", C_proj_pars)
C_proj_pars <- colnames(C_proj_draws_mat)
colnames(C_proj_draws_mat) <- gsub("-", "HYPHEN", C_proj_pars)
C_proj_pars <- colnames(C_proj_draws_mat)
colnames(C_proj_draws_mat) <- gsub("/", "SLASH", C_proj_pars)
C_proj_pars <- colnames(C_proj_draws_mat)
C_proj_draws_arr <- array(
  C_proj_draws_mat, dim = c(n_drawsPerChain, n_chains, ncol(C_proj_draws_mat))
)
dimnames(C_proj_draws_arr) <- list("iterations" = NULL,
                                   "chains" = dimnames(C_draws_arr)[[2]],
                                   "parameters" = C_proj_pars)
dummy_stan_version <- strsplit(cmdstanr::cmdstan_version(), "\\.")[[1]]
dummy_sampler_params <- rstan::get_sampler_params(C_bfit$fit,
                                                  inc_warmup = FALSE)
for (j_idx in seq_len(n_chains)) {
  stan_args_j <- C_bfit$fit@stan_args[[j_idx]]
  out_j <- cbind("lp__" = C_draws_arr[, j_idx, "lp__"],
                 dummy_sampler_params[[j_idx]],
                 C_proj_draws_arr[, j_idx, ])
  colnames(out_j) <- gsub("^[[:digit:]]+", "",
                          gsub("\\|| |\\(|\\)", "", colnames(out_j)))
  # Add dummy warmup draws:
  out_j <- rbind(matrix(0, nrow = nrow(out_j), ncol = ncol(out_j)),
                 out_j)
  ### From ?rstan::read_stan_csv (alternatively, analogous code may be derived
  ### from the CmdStanR Bernoulli example):
  writeLines(
    c("# Samples Generated by Stan",
      paste0("# stan_version_major=", dummy_stan_version[1]),
      paste0("# stan_version_minor=", dummy_stan_version[2]),
      paste0("# stan_version_patch=", dummy_stan_version[3]),
      paste0("# init=", stan_args_j$init),
      paste0("# seed=", stan_args_j$seed),
      paste0("# chain_id=", stan_args_j$chain_id),
      paste0("# iter=", stan_args_j$iter),
      paste0("# warmup=", stan_args_j$warmup),
      "# save_warmup=1",
      paste0("# thin=", stan_args_j$thin),
      paste0("# refresh=", ifelse(is.null(stan_args_j$refresh),
                                  max(stan_args_j$iter / 10, 1),
                                  stan_args_j$refresh)),
      paste0("# stepsize=", ifelse(is.null(stan_args_j$stepsize),
                                   1,
                                   stan_args_j$stepsize)),
      paste0("# stepsize_jitter=", ifelse(is.null(stan_args_j$stepsize_jitter),
                                          0,
                                          stan_args_j$stepsize_jitter)),
      "# adapt_engaged=1",
      "# adapt_gamma=0.05",
      paste0("# adapt_delta=", stan_args_j$control$adapt_delta),
      "# adapt_kappa=0.75",
      "# adapt_t0=10",
      paste0("# max_treedepth=", stan_args_j$control$max_treedepth),
      "# sampler_t=NUTS(diag_e)",
      "# sample_file=./DUMMY_FILE.csv",
      "# append_samples=0",
      "#",
      paste(colnames(out_j), collapse = ",")),
    con = file.path("output", out_folder, paste0("proj_draws_chain", j_idx, ".csv"))
  )
  write.table(
    out_j,
    file = file.path("output", out_folder, paste0("proj_draws_chain", j_idx, ".csv")),
    sep = ",",
    append = TRUE,
    col.names = FALSE,
    row.names = FALSE
  )
  ###
  stan_csv_con_j <- file(
    file.path("output", out_folder, paste0("proj_draws_chain", j_idx, ".csv")),
    open = "a"
  )
  writeLines(
    c("# ",
      "#  Elapsed Time: 1.0 seconds (Warm-up)",
      "#                1.0 seconds (Sampling)",
      "#                1.0 seconds (Total)",
      "# "),
    con = stan_csv_con_j
  )
  close(stan_csv_con_j)
}
C_proj_stanfit <- rstan::read_stan_csv(
  file.path("output", out_folder,
            paste0("proj_draws_chain", seq_len(dim(C_proj_draws_arr)[2]), ".csv"))
)
# Replace the placeholders in the parameter names by their original symbols:
stopifnot(identical(C_proj_stanfit@model_pars, names(C_proj_stanfit)))
names(C_proj_stanfit) <- gsub("COMMA", ",", names(C_proj_stanfit))
names(C_proj_stanfit) <- gsub("SQBL", "[", names(C_proj_stanfit))
names(C_proj_stanfit) <- gsub("SQBR", "]", names(C_proj_stanfit))
names(C_proj_stanfit) <- gsub("DOT", ".", names(C_proj_stanfit))
names(C_proj_stanfit) <- gsub("HYPHEN", "-", names(C_proj_stanfit))
names(C_proj_stanfit) <- gsub("SLASH", "/", names(C_proj_stanfit))
C_proj_stanfit@model_pars <- names(C_proj_stanfit)

# Create a `brmsfit` object from the projected draws:
if (exists("vpreds")) rm(vpreds)
if (exists("vpreds_noInt")) rm(vpreds_noInt)
if (exists("C_datMM")) rm(C_datMM)
if (exists("K_coefs")) rm(K_coefs)
if (exists("p0")) rm(p0)
if (exists("C_formula_tmp")) rm(C_formula_tmp)
if (any(grepl("^sd_.*__", C_proj_pars))) {
  vpreds_PP <- sub("^sd_", "",
                   sub("__.*$", "", grep("^sd_.*__", C_proj_pars, value = TRUE)))
  stopifnot(identical(length(vpreds_PP), 1L))
  vpreds_PP <- paste0("(1 | ", vpreds_PP, ")")
} else {
  vpreds_PP <- character()
}
C_bformula <- update(
  C_bformula,
  paste(". ~",
        paste(setdiff(c(sub("^b_", "", grep("^b_STD_", C_proj_pars, value = TRUE)),
                        vpreds_PP),
                      c("Intercept", "sigma")),
              collapse = " + ")
  )
)
C_formula <- formula(C_bformula)
C_termlabs <- labels(terms(C_formula))
if (exists("vfit")) rm(vfit)
if (exists("hasNA")) rm(hasNA)
if (exists("C_family")) {
  stop("Assuming that the outcome family is included in the `brmsformula` ",
       "object.")
}
rm(C_prior_def)
C_prior <- NULL
if (exists("hasCoef")) {
  rm(hasCoef)
  rm(hasCoef_cat)
  rm(hasCoef_cont)
}
C_pfit <- brm(
  formula = C_bformula,
  data = C_dat,
  prior = C_prior,
  iter = 2 * n_drawsPerChain,
  backend = "cmdstanr",
  adapt_delta = 0.99,
  max_treedepth = 15L,
  save_pars = save_pars(all = TRUE),
  seed = 244762505,
  refresh = 0,
  empty = TRUE
)
C_pfit$fit <- C_proj_stanfit
C_pfit <- rename_pars(C_pfit)
saveRDS(C_pfit, file = file.path("output", out_folder, "C_pfit.rds"))
# Replace the existing object `C_bfit` and the remaining objects referring to
# it (apart from those which are created in `helpers/analyze_model.R`):
C_bfit <- C_pfit
rm(C_pfit)
C_bformula_ff <- formula(C_bfit)
stopifnot(identical(formula(C_bformula_ff), C_formula))

# Parallelization teardown ------------------------------------------------

options(trigger_default)
doParallel::stopImplicitCluster()
foreach::registerDoSEQ()
