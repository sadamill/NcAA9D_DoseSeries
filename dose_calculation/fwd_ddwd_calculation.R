source("scripts/functions.R")

doses <- dplyr::bind_rows(
  calculate_dose("fwd"),
  calculate_dose("ddwd")
)

ggplot2::ggplot(doses, ggplot2::aes(x = start_angle, y = dose, color = dose_type)) +
  ggplot2::geom_point() +
  ggplot2::scale_color_manual(
    "Dose Type",
    labels = c("Fluence-Weighted Dose", "Diffraction Decay-Weighted Dose"),
    breaks = c("fwd", "ddwd"),
    values = c("#fa8a15", "#c688ff")
  ) +
  ggplot2::facet_wrap(~ dataset_type) +
  ggplot2::labs(x = "Start Angle (φ, °)", y = "Dose (MGy)") +
  ggplot2::theme_bw() +
  ggplot2::theme(legend.position = "top")
