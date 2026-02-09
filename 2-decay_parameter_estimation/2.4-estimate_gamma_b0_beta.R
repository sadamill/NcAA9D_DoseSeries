# calculate γ, β, and B0 for RADDOSE input

library(tidyverse)
library(cowplot)

source("scripts/functions.R")

scaling_params <- dplyr::bind_cols(
  readr::read_csv(file = "2-decay_parameter_estimation/input/r_input/fwds.csv",
                  col_select = dplyr::select(dataset_num, dose),
                  col_names = c("dataset_num", "fwd")),
  readr::read_csv(file = "2-decay_parameter_estimation/input/r_input/ccp4_wilson_b.csv",
                  col_names = "wilson_b"),
  readr::read_csv(file = "2-decay_parameter_estimation/input/r_input/ccp4_wilson_scale.csv",
                  col_names = "scale")
) |>
  dplyr::mutate(k = 1/scale) |> 
  dplyr::mutate(
    b_used = c(rep(TRUE, 9), rep(FALSE, 16)),
    k_used = c(rep(TRUE, 17), rep(FALSE, 8))
  )

wilson_b_fit <- lm(wilson_b ~ fwd, scaling_params, subset = 1:9)
scale_fit <- lm(log(k) ~ {fwd^2}, scaling_params, subset = 1:17)$coefficients |> 
  setNames(c("c", "g"))
scale_fit[1] <- exp(scale_fit[1])
scale_fit[2] <- sqrt(-scale_fit[2])

parameters <- tibble::tibble(
  gamma = scale_fit[2],
  b_naught = wilson_b_fit$coefficients[1],
  beta = wilson_b_fit$coefficients[2]
)

wilson_b_plot <- ggplot2::ggplot(scaling_params, ggplot2::aes(x = fwd, y = wilson_b, color = b_used)) +
  ggplot2::geom_point() +
  ggplot2::geom_abline(
    slope = wilson_b_fit$coefficients[2],
    intercept = wilson_b_fit$coefficients[1],
    inherit.aes = FALSE,
    linetype = 2
  ) +
  ggplot2::scale_color_manual(
    "Used in Regression?",
    labels = c("Yes", "No"),
    breaks = c(TRUE, FALSE),
    values = c("black", "gray50")
  ) +
  ggplot2::theme_bw() +
  ggplot2::theme(
    legend.position = "inside",
    legend.position.inside = c(0.05, 0.95),
    legend.justification = c(0, 1),
    legend.background = ggplot2::element_rect(color = "black")
  ) +
  ggplot2::labs(
    x = "FWD (MGy)",
    y = bquote("Wilson B-Factor ("*Å^2*")")
  )

scales_plot <- ggplot2::ggplot(scaling_params, ggplot2::aes(x = fwd, y = k, color = k_used)) +
  ggplot2::geom_point() +
  ggplot2::geom_function(
    fun = \(fwd) {scale_fit[1] * exp(-(scale_fit[2]^2) * (fwd^2))},
    inherit.aes = FALSE,
    linetype = 2
  ) +
  ggplot2::scale_color_manual(
    "Used in Regression?",
    labels = c("Yes", "No"),
    breaks = c(TRUE, FALSE),
    values = c("black", "gray50")
  ) +
  ggplot2::theme_bw() +
  ggplot2::theme(
    legend.position = "inside",
    legend.position.inside = c(0.05, 0.05),
    legend.justification = c(0, 0),
    legend.background = ggplot2::element_rect(color = "black")
  ) +
  ggplot2::labs(
    x = "FWD (MGy)",
    y = "Scale Factor K"
  )

ggplot2::ggsave("2-decay_parameter_estimation/output/r_output/wilson_b_plot.svg",
                plot = wilson_b_plot,
                height = 5, width = 8,
                unit = "cm",
                create.dir = TRUE)
ggplot2::ggsave("2-decay_parameter_estimation/output/r_output/scales_plot.svg",
                plot = scales_plot,
                height = 5, width = 8,
                unit = "cm",
                create.dir = TRUE)
readr::write_csv(parameters, "2-decay_parameter_estimation/output/r_output/decay_params.csv")
