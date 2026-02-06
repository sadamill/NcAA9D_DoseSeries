library(tidyverse)
source("functions.R")

n_samples <- 36
max_distance <- 25 
set.seed(12)

calculate_dose("ddwd")

pseudohelix_doses <- calculate_dose("ddwd") |> 
  filter(dataset_type == "pseudohelix") |> 
  rename(ddwd = dose) |> 
  select(dataset_number, start_angle, ddwd) |> 
  mutate(weight = {diff(c(0, ddwd)) |> abs()})

assured_datasets <- pseudohelix_doses |>
  filter(dataset_number %in% seq(1, 176, max_distance))

candidate_points <- pseudohelix_doses |>
  filter(!dataset_number %in% seq(1, 176, max_distance))
samples <- candidate_points %>% slice_sample(
  n = n_samples-length(seq(1, 176, max_distance)),
  weight_by = weights
) |> 
  bind_rows(assured_datasets) |> 
  arrange(dataset_number)

pseudohelix_doses <- mutate(pseudohelix_doses, sampled = dataset_number %in% samples$dataset_number)

p <- ggplot(samples, aes(x = start_angle, y = ddwd)) +
  geom_point(data = filter(pseudohelix_doses, sampled == FALSE), color = "gray95") +
  geom_point(color = "red3") +
  theme_bw() +
  labs(
    x = "Start Angle (φ, °)",
    y = "DDWD (MGy)"
  )

ggExtra::ggMarginal(
  p, 
  type = "histogram",
  xparams = list(binwidth = 10),
  yparams = list(binwidth = 0.5)
)

write_csv(samples, "sampling/output/samples.csv")
