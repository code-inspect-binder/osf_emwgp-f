#_______________________________________________________________________________
# Code for "Identifying predictors of clinical outcomes using the
# projection-predictive feature selection - a proof of concept on the example of
# Crohn's disease"
#
# Descriptive statistics: Frequencies of categorical predictors (except for
# those from the "general characteristics" table)
#_______________________________________________________________________________

# Packages ----------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
})

# Data --------------------------------------------------------------------

C_dat <- readRDS("CD_data.rds")

# Number of visits --------------------------------------------------------

C_Nvis <- nrow(C_dat)

# Frequencies -------------------------------------------------------------

print(C_dat[, .(.N, Percentage = round(.N / C_Nvis, digits = 2)), weight_gain])

print(C_dat[, .(.N, Percentage = round(.N / C_Nvis, digits = 2)), condition])

print(C_dat[, .(.N, Percentage = round(.N / C_Nvis, digits = 2)), activity_limitation])

print(C_dat[, .(.N, Percentage = round(.N / C_Nvis, digits = 2)), stoolquantcategory])

print(C_dat[, .(.N, Percentage = round(.N / C_Nvis, digits = 2)), stool_consistency])

print(C_dat[, .(.N, Percentage = round(.N / C_Nvis, digits = 2)), stool_blood])

print(C_dat[, .(.N, Percentage = round(.N / C_Nvis, digits = 2)), abdominal_pain])

print(C_dat[, .(.N, Percentage = round(.N / C_Nvis, digits = 2)), abdominal_pain_night])

print(C_dat[, .(.N, Percentage = round(.N / C_Nvis, digits = 2)), abdfinding])

print(C_dat[, .(.N, Percentage = round(.N / C_Nvis, digits = 2)), abdfinding_pressurepain])

print(C_dat[, .(.N, Percentage = round(.N / C_Nvis, digits = 2)), abdfinding_resistance])

print(C_dat[, .(.N, Percentage = round(.N / C_Nvis, digits = 2)), analfinding])

print(C_dat[, .(.N, Percentage = round(.N / C_Nvis, digits = 2)), appetite])

print(C_dat[, .(.N, Percentage = round(.N / C_Nvis, digits = 2)), cohort])

print(C_dat[, .(.N, Percentage = round(.N / C_Nvis, digits = 2)), ext.manifestation_sum])

print(C_dat[, .(.N, Percentage = round(.N / C_Nvis, digits = 2)), height_gain])

print(C_dat[, .(.N, Percentage = round(.N / C_Nvis, digits = 2)), peri_eczema])

print(C_dat[, .(.N, Percentage = round(.N / C_Nvis, digits = 2)), weight_gain])
