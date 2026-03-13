library(tidyverse)
source("global_functions.R")

n_samples <- 36
max_distance <- 25 
set.seed(12)

pseudohelix_doses <- readr::read_csv("4-sampling/input/ddwds.csv") |> 
  dplyr::filter(dataset_type == "pseudohelix") |> 
  dplyr::rename(ddwd = dose) |> 
  dplyr::select(dataset_number, start_angle, ddwd) |> 
  dplyr::mutate(weight = {diff(c(0, ddwd)) |> abs()})

deterministic_points <- pseudohelix_doses |>
  filter(dataset_number %in% seq(1, 176, max_distance))

candidate_points <- pseudohelix_doses |>
  filter(!dataset_number %in% seq(1, 176, max_distance))

random_points <- candidate_points %>% slice_sample(
  n = n_samples-length(seq(1, 176, max_distance)),
  weight_by = weight
) |> 
  arrange(dataset_number)

samples <- bind_rows(random_points, deterministic_points) |> 
  arrange(dataset_number)

pseudohelix_doses <- mutate(
  pseudohelix_doses, 
  sampled = dataset_number %in% c(random_points$dataset_number, deterministic_points$dataset_number),
  sample_type = case_when(dataset_number %in% random_points$dataset_number ~ "random",
                          dataset_number %in% deterministic_points$dataset_number ~ "deterministic",
                          .default = "none")
)

p <- ggplot(pseudohelix_doses, aes(x = start_angle, y = ddwd)) +
  geom_line(color = "gray95") +
  geom_point(
    data = filter(pseudohelix_doses, sampled == TRUE), aes(color = sample_type, shape = sample_type)
  ) +
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
  scale_x_continuous(breaks = seq(5, 185, 90)) +
  scale_y_continuous(breaks = seq(0, 9, 3)) +
  ggtheme_light() +
  theme(legend.justification = c(1, 0),
        legend.position = c(0.95, 0.05)) +
  labs(
    x = "Start Angle (φ, °)",
    y = "DDWD (MGy)"
  ) +
  coord_cartesian(
    xlim = c(5, 185),
    ylim = c(0, 9),
    expand = FALSE
  )

samples_plot <- ggExtra::ggMarginal(
  p, 
  type = "histogram",
  xparams = list(binwidth = 10),
  yparams = list(binwidth = 0.5)
)

ggplot2::ggsave("4-sampling/output/samples.svg", samples_plot, 
                height = 5, width = 8, unit = "cm")
write_csv(samples, "4-sampling/output/samples.csv")
write_csv(pseudohelix_doses, "4-sampling/output/all_datasets.csv")
