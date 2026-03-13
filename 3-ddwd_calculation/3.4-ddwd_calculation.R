library(tidyverse)

source("global_functions.R")

doses <- dplyr::bind_rows(
  readr::read_csv("3-ddwd_calculation/input/r_input/fwds.csv"),
  calculate_dose("ddwd")
)

doses_plot <- ggplot2::ggplot(doses, ggplot2::aes(x = start_angle, y = dose, color = dose_type, linetype = dataset_type)) +
  ggplot2::geom_line() +
  ggplot2::scale_color_manual(
    "Dose Type",
    labels = c("DDWD", "FWD"),
    breaks = c("ddwd", "fwd"),
    values = c("#001e7a", "#e2b8ff")
  ) +
  ggplot2::scale_linetype_manual(
    "Dataset Type",
    labels = c("Pseudohelix", "Wedge"),
    breaks = c("pseudohelix", "wedge"),
    values = c(1, 2)
  ) +
  ggplot2::scale_x_continuous(breaks = seq(5, 185, 45)) +
  ggplot2::scale_y_continuous(breaks = seq(0, 18, 3)) +
  ggplot2::labs(x = "Start φ Angle (°)", y = "Average Dose (MGy)") +
  ggtheme_light() +
  theme(legend.position = "right") +
  coord_cartesian(xlim = c(5, 185), ylim = c(0, 17), expand = FALSE)

filter(doses, dose_type == "ddwd") |> 
  readr::write_csv("3-ddwd_calculation/output/r_output/ddwds.csv")
ggplot2::ggsave("3-ddwd_calculation/output/r_output/fwd_ddwd.svg",
                                plot = doses_plot,
                                height = 5, width = 8,
                                unit = "cm",
                                create.dir = TRUE)
