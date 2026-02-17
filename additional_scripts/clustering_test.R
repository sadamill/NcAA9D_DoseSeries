library(factoextra)
library(cowplot)

wide_wedges <- wedgeAtoms |> purrr::map(\(pdb) {
  pdb |> as_tibble() |> 
    select(c(eleno, x:b)) |> 
    pivot_wider(
      names_from = c(eleno),
      values_from = c(x, y, z, b, o)
    )
}) |> purrr::list_rbind() |> 
  dplyr::select(where(~ var(.) > 0))

pca <- prcomp(wide_wedges, scale = TRUE)
factoextra::fviz_eig(pca, addlabels = TRUE, ncp = 36)
factoextra::fviz_cos2(pca, choice = "var", axes = 1:4)

pca_filtered <- pca$x[,1:2] |> dplyr::as_tibble()


cowplot::plot_grid(
  factoextra::fviz_nbclust(pca_filtered, kmeans, "gap_stat"),
  factoextra::fviz_nbclust(pca_filtered, kmeans, "silhouette"),
  factoextra::fviz_nbclust(pca_filtered, kmeans, "wss"),
  ncol = 2
)
  
factoextra::eclust(pca_filtered, FUNcluster = "kmeans", k = 3)

pca_filtered |> 
  bind_cols(number = 1:36) |> 
  ggplot(aes(x = PC1, y = PC2, color = number, label = number)) +
  geom_path() +
  geom_text()
