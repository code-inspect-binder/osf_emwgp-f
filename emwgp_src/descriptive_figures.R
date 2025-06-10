#_______________________________________________________________________________
# Code for "Identifying predictors of clinical outcomes using the
# projection-predictive feature selection - a proof of concept on the example of
# Crohn's disease"
#
# Descriptive statistics: Figures
#_______________________________________________________________________________

# Packages ----------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

# Package options ---------------------------------------------------------

source("helpers/ggplot2_settings_global.R")

# Setup -------------------------------------------------------------------

# Set seed:
set.seed(1407698374)

# If necessary, create output folder:
if (!dir.exists("output")) dir.create("output")
if (!dir.exists(file.path("output", "figures"))) {
  dir.create(file.path("output", "figures"))
}

# Data --------------------------------------------------------------------

C_dat <- readRDS("CD_data.rds")

# Rename colo categories:
C_dat$colonoscopy[C_dat$colo_score == "colo0"] <- "remission"
C_dat$colonoscopy[C_dat$colo_score == "colo1"] <- "mild"
C_dat$colonoscopy[C_dat$colo_score == "colo2"] <- "moderate"
C_dat$colonoscopy[C_dat$colo_score == "colo3"] <- "severe"
C_dat$colonoscopy <- factor(
  C_dat$colonoscopy,
  levels = c("remission", "mild", "moderate", "severe")
)

# Long dataset for the lab variables:
C_dat_l <- melt(C_dat,
                id.vars = c("patID", "visID", "colonoscopy"),
                measure.vars = c("lab_calprotectin_si", "lab_Hb", "lab_hk",
                                 "lab_mcv", "lab_thr","lab_leuko",
                                 "lab_crp_si", "lab_albumin_si", "lab_bsg_si"),
                variable.name= "laboratory_parameter",
                variable.factor = FALSE,
                value.name= "lab_value")
C_dat_l$laboratory_parameter[
  C_dat_l$laboratory_parameter == "lab_calprotectin_si"
] <- "Fecal calprotectin (mg/kg)"
C_dat_l$laboratory_parameter[
  C_dat_l$laboratory_parameter == "lab_Hb"
] <- "Hemoglobin (mmol/L)"
C_dat_l$laboratory_parameter[
  C_dat_l$laboratory_parameter == "lab_hk"
] <- "Hematocrit (L/L)"
C_dat_l$laboratory_parameter[
  C_dat_l$laboratory_parameter == "lab_mcv"
] <- "Mean corpuscular volume (fL)"
C_dat_l$laboratory_parameter[
  C_dat_l$laboratory_parameter == "lab_thr"
] <- "Platelets (Gpt/L)"
C_dat_l$laboratory_parameter[
  C_dat_l$laboratory_parameter == "lab_leuko"
] <- "White blood cells (Gpt/L)"
C_dat_l$laboratory_parameter[
  C_dat_l$laboratory_parameter == "lab_crp_si"
] <- "C-reactive protein (mg/L)"
C_dat_l$laboratory_parameter[
  C_dat_l$laboratory_parameter == "lab_albumin_si"
] <- "Albumin (g/L)"
C_dat_l$laboratory_parameter[
  C_dat_l$laboratory_parameter == "lab_bsg_si"
] <- "Erythrocyte sedimentation rate (mm/hour)"

# Create a long dataset, containing the number of observations per colo
# category:
C_dat_obs_colo <- C_dat_l[,
                          .(N_obs_colo = sum(!is.na(lab_value))),
                          by = .(colonoscopy, laboratory_parameter)]
C_dat_l_colo <- merge(C_dat_l,
                      C_dat_obs_colo,
                      by = c("colonoscopy", "laboratory_parameter"),
                      all.x = TRUE, all.y = FALSE,
                      sort = FALSE)
colo_lev_bu <- levels(C_dat_l_colo$colonoscopy)
C_dat_l_colo$colonoscopy <- paste0(C_dat_l_colo$colonoscopy, "\n(N = ",
                                   C_dat_l_colo$N_obs_colo, ")")
colo_lev_new <- unique(unlist(lapply(colo_lev_bu, function(lev_x) {
  grep(paste0("^", lev_x), C_dat_l_colo$colonoscopy, value = TRUE)
})))
C_dat_l_colo$colonoscopy <- factor(C_dat_l_colo$colonoscopy,
                                   levels = colo_lev_new)

# Plots -------------------------------------------------------------------

# Lab variables vs. colo:
ggplot(C_dat_l_colo, aes(x = colonoscopy, y = lab_value)) +
  geom_boxplot() +
  xlab("") +
  ylab("") +
  scale_x_discrete() +
  theme_bw(base_size = 15) +
  facet_wrap("laboratory_parameter", ncol = 3, scales = "free") +
  theme(strip.text.x = element_text(size = 10))
ggsave(file.path("output", "figures", "colo_laboratory.jpeg"),
       width = 12, height = 12 * 0.618)

# CRP vs. colo (jittered boxplot):
gg_colo_crp <- ggplot(data = C_dat,
                      mapping = aes(x = colonoscopy, y = lab_crp_si)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, height = 0, alpha = 0.4) +
  xlab("Endoscopic score") +
  ylab("C-reactive protein (mg/L)") +
  scale_x_discrete(drop = FALSE) +
  scale_y_continuous(trans = "log10")

# FC vs. colo (jittered boxplot):
gg_colo_calp <- ggplot(data = C_dat,
                       mapping = aes(x = colonoscopy, y = lab_calprotectin_si)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, height = 0, alpha = 0.4) +
  xlab("Endoscopic score") +
  ylab("Fecal calprotectin (mg/kg)") +
  scale_x_discrete(drop = FALSE) +
  scale_y_continuous(trans = "log10")

print((gg_colo_crp / gg_colo_calp) + plot_annotation(tag_levels = "A"))
ggsave(file.path("output", "figures", "colo_CRP_FC_jittered_boxplots.jpeg"),
       width = 7, height = 2 * 7 * 0.618)

# Scatterplot of FC vs. CRP, colored by colo:
# Basic plot:
scatter_crp_calp_colo <- ggplot(
  C_dat, aes(x = lab_crp_si , y = lab_calprotectin_si, color = colonoscopy)
) +
  geom_point(size = 2, alpha = 0.6) +
  scale_color_manual(breaks = levels(C_dat$colonoscopy),
                     values = outc_colors_brewed)
scatter_crp_calp_colo <- scatter_crp_calp_colo +
  scale_x_continuous(trans = "log10") +
  scale_y_continuous(trans = "log10")+
  xlab("C-reactive protein (mg/L)") +
  ylab("Fecal calprotectin (mg/kg)") +
  labs(color = "Endoscopic score:") +
  theme(legend.position = "bottom")
# Add contour lines:
scatter_crp_calp_colo <- scatter_crp_calp_colo +
  geom_density_2d(bins = 3)
# Calculate colonoscopy-specific medians and add them to the plot:
C_dat_med_by_colo <- C_dat[
  !is.na(lab_crp_si) & !is.na(lab_calprotectin_si) & !is.na(colonoscopy),
  .("crp_med" = median(lab_crp_si),
    "calp_med" = median(lab_calprotectin_si)),
  by = colonoscopy
]
scatter_crp_calp_colo <- scatter_crp_calp_colo +
  geom_point(mapping = aes(x = crp_med, y = calp_med),
             data = C_dat_med_by_colo, size = 5, shape = 7, stroke = 1.2)
print(scatter_crp_calp_colo)
ggsave(file.path("output", "figures", "colo_CRP_FC_scatter.jpeg"),
       width = 7, height = 7 * 0.618)
