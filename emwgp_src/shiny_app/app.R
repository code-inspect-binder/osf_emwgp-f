#_______________________________________________________________________________
# Code for "Identifying predictors of clinical outcomes using the
# projection-predictive feature selection - a proof of concept on the example of
# Crohn's disease"
#
# Shiny application for calculating predictive probabilities for each outcome
# category
#
# Copyright (C) 2023  Frank Weber
#
# Licensed under CC BY-SA 3.0 (<https://creativecommons.org/licenses/by-sa/3.0/>):
#
# > This work is licensed under the
# > Creative Commons Attribution-ShareAlike 3.0 Unported License. To view a copy
# > of this license, visit http://creativecommons.org/licenses/by-sa/3.0/ or
# > send a letter to Creative Commons, PO Box 1866, Mountain View, CA 94042, USA.
#_______________________________________________________________________________

# Preparations and global definitions -------------------------------------

library(shiny)

outc_nms <- c("remission", "mild", "moderate", "severe")

empty_res <- as.data.frame(
  setNames(replicate(length(outc_nms), character()), outc_nms)
)

## Read-in the posterior draws --------------------------------------------

C_thres <- readRDS(file.path("data", "C_thres.rds"))
C_plcoef <- readRDS(file.path("data", "C_plcoef.rds"))

n_draws <- nrow(C_plcoef)
stopifnot(identical(n_draws, nrow(C_thres)))

stopifnot(all(grepl("^b_STD_", colnames(C_plcoef))))

## Read-in the objects needed for the standardization of the predictors ----

C_centers <- readRDS(file.path("data", "C_centers.rds"))
C_scales <- readRDS(file.path("data", "C_scales.rds"))

# UI ----------------------------------------------------------------------

ui <- fluidPage(
  titlePanel("Selected endoscopic submodel (SESM)"),
  h4("Disclaimer:"),
  p(span("The predictive probabilities are based on a proof-of-concept dataset",
         "and should not be used for clinical practice.", style = "color:red")),
  p("The predictive probabilities should not be the only means for judging the",
    "disease activity. In particular, other causes for elevated CRP or FC",
    "values need to be excluded."),
  h4("Reference:"),
  p("Wirthgen E, Weber F, Kubickova-Weber L, Schiller B, Schiller S, Radke M,",
    "D", HTML("&auml;", .noWS = "outside"), "britz J. Identifying predictors",
    "of clinical outcomes using the projection-predictive feature selection",
    HTML("&mdash;", .noWS = "outside"), "a proof of concept on the example of",
    "Crohn's disease.", em("Frontiers in Pediatrics", .noWS = "after"),
    ". 2023;11:1170563.", a("doi:10.3389/fped.2023.1170563",
                            href = "https://doi.org/10.3389/fped.2023.1170563",
                            target = "_blank", .noWS = "after")),
  sidebarLayout(
    sidebarPanel(
      helpText("Please specify a value for all predictors listed below."),
      numericInput("crp_sel",
                   label = "C-reactive protein (mg/L):",
                   value = NA, step = 0.5),
      numericInput("calp_sel",
                   label = "Fecal calprotectin (mg/kg):",
                   value = NA, step = 0.5)
    ),
    mainPanel(
      h4("Predictive probabilities for the endoscopic score categories"),
      tableOutput("prediction")
    )
  )
)

# Server ------------------------------------------------------------------

server <- function(input, output, session) {

  dat_new <- reactive({
    data.frame("crp_log" = log(input$crp_sel),
               "calp_log" = log(input$calp_sel))
  })

  predict_res <- reactive({
    isInvalid_dat_new <- any(sapply(dat_new(), function(x) {
      # In fact, `x` should always be of length 1 and much of this code is
      # probably redundant (we could probably simply return `!is.finite(x)`).
      anyNA(x) ||
        ### Should not be possible:
        any(!is.na(x) & x == "") ||
        ###
        any(!is.na(x) & is.infinite(x))
    }))
    if (isTRUE(isInvalid_dat_new)) {
      return(empty_res)
    } else {
      ## Standardize input data -------------------------------------------------

      dat_new_STD <- dat_new()
      dat_new_STD <- as.data.frame(
        scale(dat_new_STD,
              center = C_centers[names(dat_new_STD)],
              scale = C_scales[names(dat_new_STD)])
      )
      names(dat_new_STD) <- paste0("STD_", names(dat_new_STD))

      ## Perform prediction -----------------------------------------------------

      # Linear predictor:
      eta_man <- C_plcoef %*% t(as.matrix(
        dat_new_STD[, sub("^b_", "", colnames(C_plcoef)), drop = FALSE]
      ))

      # Transform to the original scale of the outcome:
      plogis_thres <- apply(C_thres, 2, plogis, location = eta_man)
      stopifnot(identical(
        colnames(plogis_thres),
        paste0("b_Intercept[", seq_len(ncol(plogis_thres)), "]")
      ))
      p_post <- t(apply(plogis_thres, 1, function(x) {
        c(x[1], diff(x), 1 - x[length(x)])
      }))
      names(dimnames(p_post))[2] <- "outc"
      dimnames(p_post)[[2]] <- outc_nms

      ## Aggregate over posterior iterations ------------------------------------

      p_mean <- apply(p_post, 2, mean)

      ## Output -----------------------------------------------------------------

      p_out_num <- 100 * p_mean
      # Function for rounding to integer values under a sum constraint, taken
      # from josliber's (<https://stackoverflow.com/users/3093387/josliber>)
      # answer
      # <https://stackoverflow.com/questions/32544646/round-vector-of-numerics-to-integer-while-preserving-their-sum/32544987#32544987>
      # to topic "Round vector of numerics to integer while preserving their sum"
      # (<https://stackoverflow.com/questions/32544646/round-vector-of-numerics-to-integer-while-preserving-their-sum>),
      # licensed under CC BY-SA 3.0
      # (<https://creativecommons.org/licenses/by-sa/3.0/>):
      round_sum_constr <- function(x) {
        y <- floor(x)
        indices <- tail(order(x - y), round(sum(x)) - sum(y))
        y[indices] <- y[indices] + 1
        y
      }
      p_out_int <- round_sum_constr(p_out_num)
      p_out <- setNames(paste(round(p_out_int), "%"), names(p_mean))
      return(as.data.frame(as.list(p_out)))
    }
  })

  output$prediction <- renderTable({
    predict_res()
  }, align = "c")

  session$onSessionEnded(function() {
    stopApp()
  })

}

# Call to shinyApp() ------------------------------------------------------

shinyApp(ui = ui, server = server)
