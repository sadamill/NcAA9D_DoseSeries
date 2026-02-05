n_samples <- 36
max_distance <- 25 
set.seed(12)

pseudohelix_doses <- filter(doses, dataset_type == "pseudohelix") %>% 
  select(dataset_number, start_angle, ddwd) %>% 
  mutate(weight = {diff(c(0, ddwd)) %>% abs()})

assured_datasets <- pseudohelix_doses %>% filter(dataset_number %in% seq(1, 176, max_distance))

candidate_points <- pseudohelix_doses %>% filter(!dataset_number %in% seq(1, 176, max_distance))
samples <- candidate_points %>% slice_sample(
  n = n_samples-length(seq(1, 176, max_distance)),
  weight_by = weight
) %>% 
  bind_rows(assured_datasets) %>% 
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

write_csv(samples, "samples.csv")
