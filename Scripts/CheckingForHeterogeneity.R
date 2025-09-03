# Occupancy Pairwise Comparisons ------------------------------------------
combinations <- combn(seq(36), 2) %>% t()
pseudohelix_occ <- sapply(pseudohelixAtoms, function(atoms) {atoms$o})
pseudohelix_occ_rmsd <- apply(combinations, 1, function(pair) {
  {(pseudohelix_occ[,pair[1]] - pseudohelix_occ[,pair[2]]) ** 2} %>% mean() %>% sqrt()
}) %>% as_tibble_col(column_name = "rmsd") %>%
  mutate(type = "pseudohelices") %>% cbind(combinations)
wedge_occ <- sapply(wedgeAtoms, function(atoms) {atoms$o})
wedge_occ_rmsd <- apply(combinations, 1, function(pair) {
  {(wedge_occ[,pair[1]] - wedge_occ[,pair[2]]) ** 2} %>% mean() %>% sqrt()
}) %>% as_tibble_col(column_name = "rmsd") %>% 
  mutate(type = "wedges") %>% cbind(combinations)
rmsd_occ <- rbind(pseudohelix_occ_rmsd, wedge_occ_rmsd) %>% 
  rename(ref_dataset = "1", comp_dataset = "2") %>% mutate(parameter = "occupancies")

# B-factor Pairwise Comparisons -------------------------------------------
pseudohelix_bf <- sapply(pseudohelixAtoms, function(atoms) {atoms$b})
pseudohelix_bf_rmsd <- apply(combinations, 1, function(pair) {
  {(pseudohelix_bf[,pair[1]] - pseudohelix_bf[,pair[2]]) ** 2} %>% mean() %>% sqrt()
}) %>% as_tibble_col(column_name = "rmsd") %>%
  mutate(type = "pseudohelices") %>% cbind(combinations)
wedge_b <- sapply(wedgeAtoms, function(atoms) {atoms$b})
wedge_bf_rmsd <- apply(combinations, 1, function(pair) {
  {(wedge_b[,pair[1]] - wedge_b[,pair[2]]) ** 2} %>% mean() %>% sqrt()
}) %>% as_tibble_col(column_name = "rmsd") %>% 
  mutate(type = "wedges") %>% cbind(combinations)
rmsd_b <- rbind(pseudohelix_bf_rmsd, wedge_bf_rmsd) %>% 
  rename(ref_dataset = "1", comp_dataset = "2") %>% mutate(parameter = "b_factors")

# XYZ Pairwise Comparisons ------------------------------------------------
pseudohelix_xyz <- sapply(pseudohelixList, function(pdb) {
  pdb$xyz %>% as.numeric()
}) %>% t()
wedge_xyz <- sapply(wedgeList, function(pdb) {
  pdb$xyz %>% as.numeric()
}) %>% t()
pseudohelix_rmsd <- rmsd(pseudohelix_xyz)
wedge_rmsd <- rmsd(wedge_xyz)

pseudohelix_xyz_rmsd <- apply(combinations, 1, function(pair) {
  row <- pair[1]
  col <- pair[2]
  return(pseudohelix_rmsd[row,col])
}) %>% as_tibble_col(column_name = "rmsd") %>% 
  mutate(type = "pseudohelices") %>% cbind(combinations)
wedge_xyz_rmsd <- apply(combinations, 1, function(pair) {
  row <- pair[1]
  col <- pair[2]
  return(wedge_rmsd[row,col])
}) %>% as_tibble_col(column_name = "rmsd") %>% 
  mutate(type = "wedges") %>% cbind(combinations)
rmsd_xyz <- rbind(pseudohelix_xyz_rmsd, wedge_xyz_rmsd) %>% 
  rename(ref_dataset = "1", comp_dataset = "2") %>% mutate(parameter = "xyz")

# RMSD Plotting -----------------------------------------------------------
all_rmsds <- rbind(rmsd_occ, rmsd_b, rmsd_xyz)
rmsd_plot <- ggplot() +
  geom_tile(
    data = all_rmsds %>% filter(parameter == "occupancies"),
    mapping = aes(x = ref_dataset, y = comp_dataset, fill = rmsd)
  ) +
  scale_fill_viridis_c(name = "Occupancy RMSD") +
  new_scale_fill() +
  geom_tile(
    data = all_rmsds %>% filter(parameter == "b_factors"),
    mapping = aes(x = ref_dataset, y = comp_dataset, fill = rmsd)
  ) +
  scale_fill_viridis_c(name = "B-Factor RMSD") +
  new_scale_fill() +
  geom_tile(
    data = all_rmsds %>% filter(parameter == "xyz"),
    mapping = aes(x = ref_dataset, y = comp_dataset, fill = rmsd)
  ) +
  scale_fill_viridis_c(name = "Coordinate RMSD") +
  facet_grid(
    type ~ parameter,
    labeller = labeller(
      .default = str_to_title,
      parameter = c(b_factors = "B-Factors", xyz = "Coordinates")
    )
  ) +
  theme_bw() +
  theme(panel.spacing = unit(0, "lines"),
    text = element_text(color = "black"), 
    strip.background = element_rect(fill = "white"), 
    plot.background = element_blank()
  ) +
  scale_x_continuous(expand = c(0, 0), breaks = seq(0, 36, 3), minor_breaks = seq(0, 36, 1)) +
  scale_y_reverse(expand = c(0, 0), breaks = seq(0, 36, 3), minor_breaks = seq(0, 36, 1)) +
  labs(x = "Reference Dataset", y = "Comparison Dataset")

# general shit ------------------------------------------------------------
avg_occ <- apply(wedge_occ, 2, mean)
