#_______________________________________________________________________________
# Code for "Identifying predictors of clinical outcomes using the
# projection-predictive feature selection - a proof of concept on the example of
# Crohn's disease"
#
# Descriptive statistics: Quantities (numbers, metrics) for general
# characteristics variables and the outcome
#_______________________________________________________________________________

# Packages ----------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
})

# Data --------------------------------------------------------------------

C_dat <- readRDS("CD_data.rds")

# Number of visits:
n_vis <- nrow(C_dat)
print(n_vis)

# Rename colo categories:
C_dat$colonoscopy[C_dat$colo_score == "colo0"] <- "remission"
C_dat$colonoscopy[C_dat$colo_score == "colo1"] <- "mild"
C_dat$colonoscopy[C_dat$colo_score == "colo2"] <- "moderate"
C_dat$colonoscopy[C_dat$colo_score == "colo3"] <- "severe"
C_dat$colonoscopy <- factor(
  C_dat$colonoscopy,
  levels = c("remission", "mild", "moderate", "severe")
)

# Dataset containing only the first visit for every patient:
C_first <- C_dat[visID == 1, ]

# Number of patients:
n_pat <- nrow(C_first)
print(n_pat)

# General characteristics -------------------------------------------------
# Only taking the first visit here

## Sex --------------------------------------------------------------------

print(C_first[, .(.N, Percentage = round(.N / n_pat, digits = 3)), sex])

## Age, weight, height ----------------------------------------------------

print(cbind(C_first[,
                    lapply(.SD, function(x) {
                      c(mean(x, na.rm = TRUE),
                        sd(x, na.rm = TRUE),
                        quantile(x, na.rm = TRUE))
                    }),
                    .SDcols = c("age", "weight", "height")],
            "stat" = c("mean",
                       "sd",
                       paste0("quantile_", names(quantile(0))))))

# Outcome -----------------------------------------------------------------
# Taking all visits here

colo_tab <- table(C_dat$colonoscopy)
print(colo_tab)
print(proportions(colo_tab))
