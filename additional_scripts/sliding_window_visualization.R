library(tidyverse)
library(gganimate)
library(magick)

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
  ungroup()

ggplot2::ggplot() +
  geom_tile(data = ddwds, aes(x = wedge_number, y = delta_angle, fill = ddwd)) +
  scale_fill_viridis_c(option = "inferno", name = "DDWD") +
  coord_cartesian(expand = FALSE) + 
  theme_classic() +
  labs(x = "Wedge Number", y = "Δφ Angle (°)")

samples_plot <- ggplot2::ggplot(
  samples_doses, aes(
    x = start_angle,
    y = ddwd, 
    color = sample_type, 
    shape = sample_type
  )
) +
  geom_line(
    aes(x = start_angle, y = ddwd), 
    inherit.aes = FALSE, 
    color = "gray"
  ) +
  geom_point(data = filter(samples_doses, sampled == TRUE)) +
  geom_point(data = filter(samples_doses_anim, pseudohelix >= 1),
             size = 3) +
  scale_color_manual(
    "Sample Type",
    breaks = c("deterministic", "random"),
    labels = c("Deterministic", "Random"),
    values = c("#933032", "#97d775")
  ) +
  scale_shape_discrete(
    "Sample Type",
    breaks = c("deterministic", "random"),
    labels = c("Deterministic", "Random")
  ) +
  ggtheme_light() +
  labs(
    x = "Start Angle (φ, °)",
    y = "DDWD (MGy)"
  ) +
  theme(
    legend.position = "inside",
    legend.justification = c(1, 0),
    legend.position.inside = c(0.95, 0.1),
    legend.background = element_blank()
  ) +
  gganimate::transition_states(pseudohelix, transition_length = 0)

angles_plot <- ggplot2::ggplot() +
  geom_rect(
    data = ddwds, 
    aes(xmin = wedge_number-0.5, 
        xmax = wedge_number+0.5, 
        ymin = angle, 
        ymax = angle+1, 
        fill = ddwd)
  ) +
  scale_fill_viridis_c(option = "inferno", name = "DDWD") +
  scale_x_continuous(
    breaks       = seq(1, 36, 3),
    minor_breaks = seq(0.5, 36.5, 1)
  ) +
  coord_cartesian(expand = FALSE) + 
  ggtheme_light() +
  theme(
    panel.grid.minor.x     = element_line(color = "gray80", linewidth = 0.2),
    panel.grid.major.x     = element_blank(),
    legend.position        = "inside",
    legend.justification   = c(1, 0),
    legend.position.inside = c(1, 0),
    legend.background      = element_blank()
  ) +
  labs(x = "Wedge Number", y = "φ Angle (°)")

angles_anim <- angles_plot +
  geom_rect(
    data = dose_rects, 
    aes(xmin = wedge - 0.5, 
        xmax = wedge + 0.5, 
        ymin = start, 
        ymax = start + 5), 
    alpha = 0.2, 
    fill = "white", 
    color = "white", 
    linewidth = 0.2
  ) +
  geom_text(
    data = dose_rects,
    aes(label = stringr::str_glue("Pseudohelix {pseudohelix}\nφ = {min_start}-{min_start+180}")),
    stat = "unique",
    x = 4,
    y = 300,
    color = "black",
    hjust = 0,
    size = 7
  ) +
  gganimate::transition_states(pseudohelix, transition_length = 0)

samples_gif <- gganimate::animate(plot = samples_plot,
                                  height = 8, width = 12, units = "cm", res = 300,
                                  fps = 100, nframes = 500,
                                  renderer = magick_renderer())
angles_gif <- gganimate::animate(plot = angles_anim, 
                                 height = 8, width = 12, units = "cm", res = 300,
                                 fps = 100, nframes = 500,
                                 renderer = magick_renderer())

new_gif <- magick::image_append(c(samples_gif[1], angles_gif[1]), stack = TRUE)
for(i in 2:500) {
  combined <- magick::image_append(c(samples_gif[i], angles_gif[i]), stack = TRUE)
  new_gif <- c(new_gif, combined)
}

ggsave("additional_scripts/output/angles.tiff", angles_plot, height = 10, width = 16, unit = "cm")
image_write(new_gif, "additional_scripts/output/samples_anim.gif")
