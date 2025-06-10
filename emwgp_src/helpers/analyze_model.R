#_______________________________________________________________________________
# Code for "Identifying predictors of clinical outcomes using the
# projection-predictive feature selection - a proof of concept on the example of
# Crohn's disease"
#
# Helper file for analyzing a given Bayesian model
#_______________________________________________________________________________

# Preparations ------------------------------------------------------------

# Extract the posterior draws:
C_draws_arr <- as.array(C_bfit)
C_draws_mat <- as.matrix(C_bfit)

# Auxiliary objects:
n_drawsPerChain <- dim(C_draws_arr)[1]
n_chains <- dim(C_draws_arr)[2]
n_draws <- nrow(C_draws_mat)

C_pars <- colnames(C_draws_mat)
stopifnot(identical(C_pars, variables(C_bfit)))

# Posterior draws:
C_prep <- prepare_predictions(C_bfit)

# Thresholds:
C_thres <- C_draws_mat[,
                       grep("^b_Intercept", colnames(C_draws_mat), value = TRUE),
                       drop = FALSE]
stopifnot(identical(structure(C_thres, nchains = n_chains), C_prep$thres$thres))
saveRDS(C_thres, file.path("output", out_folder, "C_thres.rds"))

# Population-level coefficients:
C_plcoef <- C_draws_mat[,
                        setdiff(grep("^b_", colnames(C_draws_mat), value = TRUE),
                                colnames(C_thres)),
                        drop = FALSE]
stopifnot(identical(structure(C_plcoef, nchains = n_chains), C_prep$dpars$mu$fe$b))
saveRDS(C_plcoef, file.path("output", out_folder, "C_plcoef.rds"))

# Standard deviation of partially pooled effects (here: "sd_patID__Intercept"):
if (identical(vpreds_PP, "(1 | patID)")) {
  C_glsd <- C_draws_mat[,
                        grep("^sd_", colnames(C_draws_mat), value = TRUE),
                        drop = FALSE]
  stopifnot(identical(ncol(C_glsd), 1L))
  C_glsd <- C_glsd[, 1]
  saveRDS(C_glsd, file.path("output", out_folder, "C_glsd.rds"))
} else if (length(vpreds_PP) > 0) {
  warning("NOTE: Unknown partially pooled effects. Therefore, 'C_glsd' was not ",
          "created. In case of other partially pooled effects (not '(1 | patID)'), ",
          "the code from \"helpers/analyze_model.R\" needs to be adapted.",
          call. = FALSE)
}

# MCMC diagnostics --------------------------------------------------------

source("helpers/check_MCMC_diagn.R")
C_MCMC_diagn <- check_MCMC_diagn(C_stanfit = C_bfit$fit,
                                 pars = "disc",
                                 include = FALSE)
cat("\n-----\n")
cat("Are all MCMC diagnostics OK?:\n")
print(C_MCMC_diagn$all_OK)
cat("-----\n")

# Model diagnostics -------------------------------------------------------

if (run_loo && !isTRUE(loo_forProj)) {
  ## Pareto k-values from PSIS-LOO CV ---------------------------------------

  ### Run LOO CV (part 1) ---------------------------------------------------

  C_loo <- loo(
    C_bfit, save_psis = TRUE,
    cores = if (.Platform$OS.type == "windows") 1 else getOption("mc.cores", 1)
  )
  saveRDS(C_loo, file.path("output", out_folder, "C_loo_orig.rds"))

  ### Pareto k-values -------------------------------------------------------

  cat("\n-----\n")
  cat("Pareto k-values:\n")
  print(pareto_k_table(C_loo))
  cat("-----\n")

  ## Effective number of parameters from PSIS-LOO CV ------------------------

  ### Run LOO CV (part 2) ---------------------------------------------------

  # If necessary: LOO CV with refits for problematic observations:
  if (any(C_loo$diagnostics$pareto_k > 0.7)) {
    message("reloo() will need to re-fit.")
  }
  C_loo <- reloo(C_bfit, loo = C_loo, seed = 244762505 + 2L)
  saveRDS(C_loo, file.path("output", out_folder, "C_loo.rds"))

  ### Effective number of parameters and comparisons ------------------------

  cat("\n-----\n")
  cat("Effective number of parameters:\n")
  print(C_loo$estimates["p_loo", , drop = FALSE])
  cat("-----\n")

  # For comparison: The real (total) number of parameters:
  real_p <- length(
    grep("^Intercept|^lp__$|^lprior$|^disc$|^z_|^zb|^L_|^zs_",
         C_pars, value = TRUE, invert = TRUE)
  )
  cat("\n-----\n")
  cat("Real (total) number of parameters:\n")
  print(real_p)
  cat("-----\n")
  # Check:
  C_standata <- standata(C_bfit)
  real_p_ch <- C_standata$K +
    sum(grepl("^sigma$", C_pars))
  # Take the thresholds (of an ordinal regression) into account:
  if (!is.null(C_standata$nthres)) {
    real_p_ch <- real_p_ch +
      C_standata$nthres
  }
  # Take the regularized horseshoe prior into account:
  if (all(c("hs_global", "hs_slab") %in% C_pars)) {
    real_p_ch <- real_p_ch +
      sum(grepl("^hs_global$", C_pars)) +
      sum(grepl("^hs_slab$", C_pars)) +
      sum(grepl("^hs_local", C_pars))
  }
  # Take partially pooled effects into account:
  if (identical(sum(grepl("^M_", names(C_standata))), 1L)) {
    stopifnot(identical(C_standata$M_1, sum(grepl("^sd_", C_pars))))
    real_p_ch <- real_p_ch +
      C_standata$M_1 * C_standata$N_1 +
      C_standata$M_1 # For the standard deviation(s) of partially pooled effects.
    if ("NC_1" %in% names(C_standata)) {
      stopifnot(identical(C_standata$NC_1, sum(grepl("^cor_", C_pars))))
      real_p_ch <- real_p_ch +
        C_standata$NC_1
    }
  } else if (!identical(sum(grepl("^M_", names(C_standata))), 0L)) {
    stop("Unknown value for 'sum(grepl(\"^M_\", names(C_standata)))'.")
  }
  stopifnot(identical(real_p, real_p_ch))

  # Just in case that one wants to set `real_p` in relation to the number of
  # observations:
  real_N <- nrow(C_dat)
  cat("\n-----\n")
  cat("Number of observations (i.e. number of rows in the dataset)",
      "(for comparison with the effective and the real (total) number of",
      "parameters):\n")
  print(real_N)
  cat("-----\n")
}

# Conditional-effects plots -----------------------------------------------

if (!identical(plot_prefix, "")) {
  ggs_condEff <- lapply(
    setNames(
      nm = sub("^s\\(", "",
               sub("\\)$", "",
                   grep("\\|", C_termlabs, value = TRUE, invert = TRUE)))
    ),
    function(vpreds_i) {
      C_ceff <- conditional_effects(
        C_bfit,
        effects = vpreds_i,
        categorical = TRUE
      )
      if (!identical(plot_prefix, "")) {
        if (identical(vpreds_i, "STD_crp_log")) {
          pretty_xlab <- "standardized log C-reactive protein (standardized log mg/L)"
        } else if (identical(vpreds_i, "STD_calp_log")) {
          pretty_xlab <- "standardized log fecal calprotectin (standardized log mg/kg)"
        } else if (identical(vpreds_i, "STD_thr_log")) {
          pretty_xlab <- "standardized log platelets (standardized log Gpt/L)"
        } else if (identical(vpreds_i, "STD_leuko_log")) {
          pretty_xlab <- "standardized log white blood cells (standardized log Gpt/L)"
        }
      }
      if (exists("pretty_xlab")) {
        stopifnot(identical(length(C_ceff), 1L))
        levels(C_ceff[[1]]$effect2__) <- outc_nms
        attr(C_ceff[[1]], "response") <- outc_nm_pretty
        attr(C_ceff[[1]], "effects")[1] <- pretty_xlab
      }
      C_ceff_plot <- plot(C_ceff)
      if (exists("pretty_xlab") &&
          grepl("standardized ", pretty_xlab) &&
          grepl("log ", pretty_xlab)) {
        do_std_log <- function(x_breaks_orig) {
          (log(x_breaks_orig) - C_centers[sub("^STD_", "", vpreds_i)]) /
            C_scales[sub("^STD_", "", vpreds_i)]
        }
        doback_std_log <- function(x_breaks) {
          exp(C_scales[sub("^STD_", "", vpreds_i)] * x_breaks +
                C_centers[sub("^STD_", "", vpreds_i)])
        }
        cust_ggbreaks <- function(x_limits) {
          tmp_breaks_log <- scales::breaks_log(n = 8, base = 10)
          do_std_log(tmp_breaks_log(doback_std_log(x_limits)))
        }
        cust_gglabels <- doback_std_log
        print(
          C_ceff_plot[[1]] +
            ylab(paste("Probability for", outc_adjective, "score")) +
            scale_x_continuous(
              name = sub("^white", "White",
                         sub("^platelets", "Platelets",
                             sub("^fecal", "Fecal",
                                 gsub("standardized ", "",
                                      gsub("log ", "", pretty_xlab))))),
              breaks = cust_ggbreaks,
              labels = cust_gglabels
            ) +
            scale_color_manual(breaks = outc_nms, values = outc_colors_brewed) +
            scale_fill_manual(breaks = outc_nms, values = outc_colors_brewed)
        )
      }
      return(last_plot())
    }
  )
  library(patchwork)
  print(wrap_plots(ggs_condEff, ncol = 1) + plot_annotation(tag_levels = "A"))
  ggsave(file.path("output", out_folder, paste0(plot_prefix, "condEff_patchwork.jpeg")),
         width = 7, height = length(ggs_condEff) * 7 * 0.618)
  saveRDS(
    last_plot(),
    file = file.path("output", out_folder, paste0(plot_prefix, "condEff_patchwork.rds"))
  )
}
