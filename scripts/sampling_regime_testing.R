n_samples <- 36
max_distance <- 25 
doses <- filter(dwds, dataset_type == "pseudohelix") %>% 
  mutate(weight = {diff(c(0, ddwd)) %>% abs()})

assured_datasets <- doses %>% filter(dataset_number %in% seq(1, 176, max_distance))

candidate_points <- doses %>% filter(!dataset_number %in% seq(1, 176, max_distance))
samples <- candidate_points %>% slice_sample(
  n = n_samples-length(seq(1, 176, max_distance)),
  weight_by = weight
) %>% 
  bind_rows(assured_datasets) %>% 
  arrange(dataset_number)

p <- ggplot(samples, aes(x = start_angle, y = ddwd)) +
  geom_point() +
  theme_classic()

ggExtra::ggMarginal(
  p, 
  type = "histogram",
  xparams = list(binwidth = 10),
  yparams = list(binwidth = 0.5)
)

write_csv(samples, "samples.csv")
