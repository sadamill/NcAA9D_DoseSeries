# Occupancy Pairwise Comparisons ------------------------------------------
combinations <- tidyr::expand_grid(reference = 1:36, comparison = 1:36)
pseudohelix_occ <- sapply(pseudohelixAtoms, function(atoms) {atoms$o})
pseudohelix_occ_rmsd <- apply(combinations, 1, function(pair) {
  {(pseudohelix_occ[,pair[1]] - pseudohelix_occ[,pair[2]]) ** 2} %>% mean() %>% sqrt()
}) %>% tibble::as_tibble_col(column_name = "rmsd") %>%
  dplyr::mutate(type = "pseudohelices") %>% cbind(combinations)
wedge_occ <- sapply(wedgeAtoms, function(atoms) {atoms$o})
wedge_occ_rmsd <- apply(combinations, 1, function(pair) {
  {(wedge_occ[,pair[1]] - wedge_occ[,pair[2]]) ** 2} %>% mean() %>% sqrt()
}) %>% tibble::as_tibble_col(column_name = "rmsd") %>% 
  dplyr::mutate(type = "wedges") %>% cbind(combinations)
rmsd_occ <- rbind(pseudohelix_occ_rmsd, wedge_occ_rmsd) %>% 
  dplyr::rename(ref_dataset = "reference", comp_dataset = "comparison") %>% dplyr::mutate(parameter = "occupancies")

# B-factor Pairwise Comparisons -------------------------------------------
pseudohelix_bf <- sapply(pseudohelixAtoms, function(atoms) {atoms$b})
pseudohelix_bf_rmsd <- apply(combinations, 1, function(pair) {
  {(pseudohelix_bf[,pair[1]] - pseudohelix_bf[,pair[2]]) ** 2} %>% mean() %>% sqrt()
}) %>% tibble::as_tibble_col(column_name = "rmsd") %>%
  dplyr::mutate(type = "pseudohelices") %>% cbind(combinations)
wedge_b <- sapply(wedgeAtoms, function(atoms) {atoms$b})
wedge_bf_rmsd <- apply(combinations, 1, function(pair) {
  {(wedge_b[,pair[1]] - wedge_b[,pair[2]]) ** 2} %>% mean() %>% sqrt()
}) %>% tibble::as_tibble_col(column_name = "rmsd") %>% 
  dplyr::mutate(type = "wedges") %>% cbind(combinations)
rmsd_b <- rbind(pseudohelix_bf_rmsd, wedge_bf_rmsd) %>% 
  dplyr::rename(ref_dataset = "reference", comp_dataset = "comparison") %>% dplyr::mutate(parameter = "b_factors")

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
}) %>% tibble::as_tibble_col(column_name = "rmsd") %>% 
  dplyr::mutate(type = "pseudohelices") %>% cbind(combinations)
wedge_coord_rmsd <- apply(combinations, 1, function(pair) {
  row <- pair[1]
  col <- pair[2]
  return(wedge_rmsd[row,col])
}) %>% tibble::as_tibble_col(column_name = "rmsd") %>% 
  dplyr::mutate(type = "wedges") %>% cbind(combinations)
rmsd_coord <- rbind(pseudohelix_coord_rmsd, wedge_coord_rmsd) %>% 
  dplyr::rename(ref_dataset = "reference", comp_dataset = "comparison") %>% dplyr::mutate(parameter = "coordinates")

all_rmsds <- rbind(rmsd_occ, rmsd_b, rmsd_coord) %>% 
  dplyr::mutate(parameter = factor(parameter, levels = c("b_factors", "occupancies", "coordinates")))
