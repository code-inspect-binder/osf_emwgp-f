#_______________________________________________________________________________
# Code for "Identifying predictors of clinical outcomes using the
# projection-predictive feature selection - a proof of concept on the example of
# Crohn's disease"
#
# Helper file specifying ggplot2 settings
#_______________________________________________________________________________

theme_set(theme_bw(base_size = 12))

label_fun_thous <- scales::label_number()

scale_x_continuous <- function(...) {
  dot_args <- list(...)
  if ("labels" %in% names(dot_args)) {
    label_fun_orig <- dot_args$labels
    stopifnot(is.function(label_fun_orig))
    dot_args$labels <- NULL
    if (!identical(label_fun_orig, scales::label_percent(), ignore.environment = TRUE)) {
      label_fun <- function(x_breaks) {
        label_fun_thous(label_fun_orig(x_breaks))
      }
    } else {
      label_fun <- label_fun_orig
    }
  } else {
    label_fun <- label_fun_thous
  }
  do.call(ggplot2::scale_x_continuous, args = c(dot_args, list(labels = label_fun)))
}
scale_y_continuous <- function(...) {
  dot_args <- list(...)
  if ("labels" %in% names(dot_args)) {
    label_fun_orig <- dot_args$labels
    stopifnot(is.function(label_fun_orig))
    dot_args$labels <- NULL
    if (!identical(label_fun_orig, scales::label_percent(), ignore.environment = TRUE)) {
      label_fun <- function(x_breaks) {
        label_fun_thous(label_fun_orig(x_breaks))
      }
    } else {
      label_fun <- label_fun_orig
    }
  } else {
    label_fun <- label_fun_thous
  }
  do.call(ggplot2::scale_y_continuous, args = c(dot_args, list(labels = label_fun)))
}

# Color palette for the outcome categories:
outc_colors_brewed <- rev(RColorBrewer::brewer.pal(4, "RdYlBu"))
