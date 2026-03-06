# calculate γ, β, and B0 for RADDOSE input

library(propagate)
library(tidyverse)
library(cowplot)

source("global_functions.R")

scaling_params <- dplyr::bind_cols(
  readr::read_csv(file = "2-decay_parameter_estimation/input/r_input/fwds.csv",
                  col_select = dataset_number | dose)[1:25,],
  readr::read_csv(file = "2-decay_parameter_estimation/input/r_input/ccp4_wilson_b.csv",
                  col_names = "wilson_b"),
  readr::read_csv(file = "2-decay_parameter_estimation/input/r_input/ccp4_wilson_scale.csv",
                  col_names = "scale")
) |>
  dplyr::rename(fwd = dose) |> 
  dplyr::mutate(k = 1/scale) |> 
  dplyr::mutate(
    b_used = c(rep(TRUE, 9), rep(FALSE, 16)),
    k_used = c(rep(TRUE, 17), rep(FALSE, 8))
  )

wilson_b_fit <- lm(wilson_b ~ fwd, scaling_params, subset = 1:9)
scale_fit <- nls(k ~ c * exp(-(fwd^2) * (gamma^2)), scaling_params, start = c(c = 0.03, gamma = 0.02), subset = k_used == TRUE)

parameters <- tibble::tibble(
  gamma = coef(scale_fit)[2],
  b_naught = wilson_b_fit$coefficients[1],
  beta = wilson_b_fit$coefficients[2]
)

scales_points <- predictNLS(scale_fit, newdata = select(filter(scaling_params, k_used == TRUE), fwd))$summary[, c(2, 5, 6)] |> 
  as_tibble() |> 
  rename(pred = mean.2, lower = `2.5%`, upper = `97.5%`) |>
  bind_cols(select(filter(scaling_params, k_used == TRUE), fwd),
            current = _)

wilson_b_plot <- ggplot2::ggplot(scaling_params, ggplot2::aes(x = fwd, y = wilson_b, color = b_used)) +
  ggplot2::geom_point() +
  ggplot2::geom_smooth(data = filter(scaling_params, b_used == TRUE),
                       method = "lm",
                       show.legend = FALSE,
                       linetype = 5, linewidth = 0.5) +
  ggplot2::scale_color_manual("Point Used",
                              labels = c("Yes", "No"),
                              breaks = c(TRUE, FALSE),
                              values = c("black", "gray80")) +
  ggtheme_light() +
  ggplot2::theme(legend.position = "inside",
                 legend.position.inside = c(0.95, 0.05),
                 legend.justification = c(1, 0),
                 legend.background = ggplot2::element_rect(color = "black")) +
  ggplot2::labs(
    x = "Average FWD (MGy)",
    y = bquote("Wilson B-Factor ("*Å^2*")")
  )

scales_plot <- ggplot2::ggplot(scaling_params, ggplot2::aes(x = fwd, y = k, color = k_used)) +
  ggplot2::geom_ribbon(data = scales_points, 
                       aes(x = fwd, ymin = lower, ymax = upper), 
                       inherit.aes = FALSE, 
                       fill = "gray", alpha = 0.5) +
  ggplot2::geom_line(data = scales_points, 
                     aes(x = fwd, y = pred), 
                     inherit.aes = FALSE,
                     color = "black", linetype = 5, linewidth = 0.5) +
  ggplot2::geom_point() +
  ggplot2::scale_color_manual("Point Used", 
                              labels = c("Yes", "No"), 
                              breaks = c(TRUE, FALSE),
                              values = c("black", "gray80")) +
  ggtheme_light() +
  ggplot2::theme(
    legend.position = "inside",
    legend.position.inside = c(0.05, 0.05),
    legend.justification = c(0, 0),
    legend.background = ggplot2::element_rect(color = "black")
  ) +
  ggplot2::labs(x = "Average FWD (MGy)",
                y = "Scale Factor K")

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
