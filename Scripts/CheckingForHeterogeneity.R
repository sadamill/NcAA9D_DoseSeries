# Occupancy Pairwise Comparisons ------------------------------------------
combinations <- expand_grid(reference = 1:36, comparison = 1:36)
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
  rename(ref_dataset = "reference", comp_dataset = "comparison") %>% mutate(parameter = "occupancies")

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
  rename(ref_dataset = "reference", comp_dataset = "comparison") %>% mutate(parameter = "b_factors")

# XYZ Pairwise Comparisons ------------------------------------------------
pseudohelix_coord <- sapply(pseudohelixList, function(pdb) {
  pdb$xyz %>% as.numeric()
}) %>% t()
wedge_coord <- sapply(wedgeList, function(pdb) {
  pdb$xyz %>% as.numeric()
}) %>% t()
pseudohelix_rmsd <- rmsd(pseudohelix_coord)
wedge_rmsd <- rmsd(wedge_coord)

pseudohelix_coord_rmsd <- apply(combinations, 1, function(pair) {
  row <- pair[1]
  col <- pair[2]
  return(pseudohelix_rmsd[row,col])
}) %>% as_tibble_col(column_name = "rmsd") %>% 
  mutate(type = "pseudohelices") %>% cbind(combinations)
wedge_coord_rmsd <- apply(combinations, 1, function(pair) {
  row <- pair[1]
  col <- pair[2]
  return(wedge_rmsd[row,col])
}) %>% as_tibble_col(column_name = "rmsd") %>% 
  mutate(type = "wedges") %>% cbind(combinations)
rmsd_coord <- rbind(pseudohelix_coord_rmsd, wedge_coord_rmsd) %>% 
  rename(ref_dataset = "reference", comp_dataset = "comparison") %>% mutate(parameter = "coordinates")

all_rmsds <- rbind(rmsd_occ, rmsd_b, rmsd_coord) %>% 
  mutate(parameter = factor(parameter, levels = c("b_factors", "occupancies", "coordinates")))

# Local RMSDs ---------------------------------------------------------
make_rmsd_df <- function(parameter, window_size = 5) {
  wedges_data <- filter(all_rmsds, parameter == !!parameter & type == "wedges")
  pseudohelices_data <- filter(all_rmsds, parameter == !!parameter & type == "pseudohelices")
  
  wedges_df <- tibble(
    parameter = rep(parameter, 36 - window_size + 1),
    type = rep("wedges", 36 - window_size + 1),
    window_middle = seq((window_size+1)/2, 36+1-(window_size+1)/2),
    local_rmsd = sapply(1:(36 - window_size + 1), function(i) {
      start_idx <- i
      end_idx <- i + window_size - 1
      combinations <- combn(start_idx:end_idx, 2)
      window_rmsd <- sapply(1:ncol(combinations), function(col) {
        ref <- combinations[1, col]
        comp <- combinations[2, col]
        rmsd <- filter(wedges_data, ref_dataset == ref & comp_dataset == comp)$rmsd
        return(rmsd)
      })
      window_average <- mean(window_rmsd)
      return(window_average)
    })
  )
  pseudohelices_df <- tibble(
    parameter = rep(parameter, 36 - window_size + 1),
    type = rep("pseudohelices", 36 - window_size + 1),
    window_middle = seq((window_size+1)/2, 36+1-(window_size+1)/2),
    local_rmsd = sapply(1:(36 - window_size + 1), function(i) {
      start_idx <- i
      end_idx <- i + window_size - 1
      combinations <- combn(start_idx:end_idx, 2)
      window_rmsd <- sapply(1:ncol(combinations), function(col) {
        ref <- combinations[1, col]
        comp <- combinations[2, col]
        rmsd <- filter(pseudohelices_data, ref_dataset == ref & comp_dataset == comp)$rmsd
        return(rmsd)
      })
      window_average <- mean(window_rmsd)
      return(window_average)
    })
  )
  
  combined_df <- rbind(wedges_df, pseudohelices_df)
  return(combined_df)
}

local_rmsds <- list()
local_rmsds$occ <- make_rmsd_df("occupancies")
local_rmsds$bf <- make_rmsd_df("b_factors")
local_rmsds$coord <- make_rmsd_df("coordinates")
local_rmsds$combined <- rbind(local_rmsds$occ, local_rmsds$bf, local_rmsds$coord)

local_rmsd_plot <- local_rmsds$combined %>% 
  ggplot(aes(x = window_middle, y = local_rmsd, color = type)) +
  geom_point() +
  geom_line() +
  facet_wrap(
    . ~ parameter,
    scales = "free",
    labeller = as_labeller(
      c(
        b_factors = "B-Factor",
        coordinates = "Coordinates",
        occupancies = "Occupancy"
      )
    )
  ) +
  labs(
    x = "Window Middle",
    y = "Local Average RMSD"
  ) +
  ggtheme_light() +
  scale_color_manual(
    name = "Dataset Type", 
    labels = c("Wedges", "Pseudohelices"), 
    breaks = c("wedges", "pseudohelices"),
    values = c("#0096c5", "#b8008c")
  ) +
  theme(legend.position = "bottom")
