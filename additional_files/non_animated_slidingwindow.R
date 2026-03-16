library(tidyverse)

source("global_functions.R")

samples_doses <- readr::read_csv("4-sampling/output/all_datasets.csv") |> 
  arrange(sampled)
samples_doses_anim <- samples_doses |> 
  mutate(pseudohelix = c(rep(0, 140), 1:36))

calculate_dose <- function(type) {
  if (!type %in% c("fwd", "ddwd")) {stop("invalid dose type provided")}
  
  directory <- switch(type,
                      "fwd" = "1-fwd_calculation/input/r_input",
                      "ddwd" = "3-ddwd_calculation/input/r_input")
  
  raddose_output <- lapply(1:38, function(i) {
    read.csv(stringr::str_glue("{directory}/wedge{i}.csv"), header = TRUE) |> 
      dplyr::rename(angle = DWD.Angle)
  })
  
  # calculate the respective subwedge average doses and average these across wedges
  # 2-37 to calculate the average DWD for a pseudhelix
  subwedge_matrix <- sapply(5:184, function(start_angle) {
    subwedge_averages <- lapply(2:37, function(subwedge_number) {
      wedge <- raddose_output[[subwedge_number]]
      subwedge_start <- start_angle + ((subwedge_number - 2) * 5) # Start angle depends on the subwedge you are calculating
      subwedge_average <- dplyr::filter(wedge, angle >= subwedge_start & angle <= subwedge_start + 1)$DWD |>  # Extract the DWD column for angles within a target subwedge
        mean() # Average the extracted DWDs
      return(c(subwedge_average, subwedge_start))
    }) # Creates a vector of subwedge average DWDs
    return(subwedge_averages)
  })
  
  subwedges <- subwedge_matrix |> t() |> 
    as_tibble() |> 
    setNames(1:36) |> 
    pivot_longer(1:36, names_to = "wedge") |> 
    unnest_wider(value, names_sep = "_") |> 
    setNames(c("wedge_number", "ddwd", "angle")) |> 
    mutate(
      wedge_number = as.integer(wedge_number),
      angle = as.integer(angle)
    ) |> 
    arrange(wedge_number) |> 
    bind_cols(rep(0:179, 36)) |> 
    rename(delta_angle = ...4)
  
  return(subwedges)
}

ddwds <- calculate_dose("ddwd")

starts <- sapply(filter(samples_doses, sampled == TRUE)$start_angle, function(angle) {
  seq(angle, angle + 175, 5)
}) |> as.vector()
pseudohelix <- rep(1:36, 36) |> sort()

dose_rects <- tibble::tibble(
  pseudohelix = pseudohelix,
  wedge = rep(1:36, 36),
  start = starts
) |> group_by(pseudohelix) |> 
  mutate(min_start = min(start)) |> 
  ungroup() |> 
  filter(pseudohelix %in% c(20))

angles_plot <- ggplot2::ggplot() +
  geom_rect(
    data = ddwds, 
    aes(xmin = wedge_number-0.5, 
        xmax = wedge_number+0.5, 
        ymin = angle, 
        ymax = angle+1, 
        fill = ddwd)
  ) +
  scale_fill_viridis_c(name = "DDWD",
                       option = "inferno",
                       breaks = c(1, 4, 7, 10)) +
  scale_x_continuous(
    breaks       = seq(1, 36, 3),
    minor_breaks = seq(0.5, 36.5, 1)
  ) +
  scale_y_continuous(
    breaks = c(5, 45, 90, 135, 180, 225, 270, 315, 360),
    minor_breaks = NULL,
    position = "right"
  ) +
  coord_cartesian(expand = FALSE,
                  ylim = c(0, 450)) + 
  ggtheme_light() +
  theme(
    panel.grid.minor.x = element_line(color = "gray80", linewidth = 0.2),
    panel.grid.major.x = element_blank(),
    legend.position    = "none",
    legend.justification = c(1, 0),
    legend.position.inside = c(0.98, 0.02),
    legend.key.width = unit(0.3, "cm"),
    legend.key.height = unit(0.2, "cm"),
    panel.background = element_blank(),
    panel.border = element_blank(),
    axis.ticks = element_blank()
  ) +
  labs(x = "Wedge Number", y = "φ Angle (°)") +
  geom_rect(
    data = dose_rects, 
    aes(xmin = wedge - 0.5, 
        xmax = wedge + 0.5, 
        ymin = start, 
        ymax = start + 5), 
    alpha = 0.2, 
    fill = "white", 
    color = "white", 
    linewidth = 0.1
  )

angles_plot

ggsave("additional_scripts/output/graph.tiff", height = 7, width = 8, units = "cm", dpi = 600)
  