library(tidyverse)

source("global_functions.R")

doses <- dplyr::bind_rows(
  readr::read_csv("3-ddwd_calculation/input/r_input/fwds.csv"),
  calculate_dose("ddwd")
)

doses_plot <- ggplot2::ggplot(doses, ggplot2::aes(x = start_angle, y = dose, color = dose_type)) +
  ggplot2::geom_point() +
  ggplot2::scale_color_manual(
    "Dose Type",
    labels = c("FWD", "DDWD"),
    breaks = c("fwd", "ddwd"),
    values = c("#fa8a15", "#c688ff")
  ) +
  ggplot2::facet_wrap(~ dataset_type, ncol = 1) +
  ggplot2::labs(x = "Start Angle (φ, °)", y = "Average Dose (MGy)") +
  ggplot2::theme_bw() +
  ggplot2::theme(legend.position = "top")

filter(doses, dose_type == "ddwd") |> 
  readr::write_csv("3-ddwd_calculation/output/r_output/ddwds.csv")
ggplot2::ggsave("3-ddwd_calculation/output/r_output/fwd_ddwd.svg",
                                plot = doses_plot,
                                height = 10, width = 8,
                                unit = "cm",
                                create.dir = TRUE)
