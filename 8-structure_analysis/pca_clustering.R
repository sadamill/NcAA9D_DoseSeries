
# shape data to wide format
wide_wedges <- wedgeAtoms |> purrr::map(\(pdb) {
  pdb |> as_tibble() |> 
    filter(!elesy == "H") |> # filter out hydrogens
    select(c(eleno, x:b)) |> 
    pivot_wider(
      names_from = c(eleno),
      values_from = c(x, y, z, o, b)
    )
}) |> purrr::list_rbind() |> 
  dplyr::select(where(~ var(.) > 0))
wide_pseudohelices <- pseudohelixAtoms |> purrr::map(\(pdb) {
  pdb |> as_tibble() |> 
    filter(!elesy == "H") |> # filter out hydrogens
    select(c(eleno, x:b)) |> 
    pivot_wider(
      names_from = c(eleno),
      values_from = c(x, y, z, o, b)
    )
}) |> purrr::list_rbind() |> 
  dplyr::select(where(~ var(.) > 0))

# perform pca
pca_wedge <- prcomp(wide_wedges, scale. = TRUE)
pca_pseudohelix <- prcomp(wide_pseudohelices, scale. = TRUE)

# determine the number of "significant" dimensions by broken stick metho
broken_stick <- function(pca) {
  
  n <- length(pca$sdev)
  
  # calculate the length distribution of a stick broken into n pieces
  sticks <- PCDimension::brokenStick(1:n, n)

  # eigenvalues are the square of standard deviations
  eigenvalues <- factoextra::get_eig(pca) |> as_tibble() |> 
    bind_cols(dimension = 1:36,
              middle = _,
              broken_stick_percent = PCDimension::brokenStick(1:n, n)*100) |>
    clean_names()
  
  # plot the percentages against each other
  plot <- eigenvalues[1:10,] |> 
    pivot_longer(cols = c(variance_percent, broken_stick_percent),
                 names_to = "percent_type",
                 names_ptypes = factor()) |> 
    ggplot(aes(x = dimension, group = percent_type)) +
    geom_col(data = ~ filter(.x, percent_type == "variance_percent"), 
             aes(y = value, fill = percent_type),
             color = "black") +
    geom_line(data = ~ filter(.x, percent_type == "broken_stick_percent"),
              aes(y = value, color = percent_type)) +
    geom_point(data = ~ filter(.x, percent_type == "broken_stick_percent"),
               aes(y = value, fill = percent_type),
               size = 3, shape = 23, color = "black") +
    scale_fill_manual("Percent Type",
                      values = c("broken_stick_percent" = "plum4", "variance_percent" = "thistle2"),
                      aesthetics = c("color", "fill"),
                      labels = c("Variance\nExplained", "Broken Stick")) +
    ggtheme_light() +
    theme(legend.justification = c(1, 1),
          legend.position.inside = c(0.9, 0.9),
          legend.key.spacing.y = unit(2, "points")) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
    scale_x_continuous(breaks = 1L:10L, minor_breaks = NULL) +
    labs(x = "Principal Component",
         y = "Percentage")
  
  # broken stick analysis to determine minimum number of dimensions
  ndim <- PCDimension::bsDimension(eigenvalues$eigenvalue)
  # plotting is odd in one dimension, and it doesn't hurt to include one more pc
  if(ndim < 2) {ndim <- 2}
  
  return(list(plot = plot, ndim = ndim))
}
bs_wedge <- broken_stick(pca_wedge)
bs_pseudohelix <- broken_stick(pca_pseudohelix)
ggplots$Light$pca$broken_stick <- cowplot::plot_grid(bs_wedge$plot, bs_pseudohelix$plot,
                                        ncol = 1,
                                        labels = c("Wedges", "Pseudohelices"),
                                        label_x = 0.5, hjust = 0.5, vjust = 2)

pca_wedge_filtered <- pca_wedge$x[,1:bs_wedge$ndim] |> dplyr::as_tibble()
pca_pseudohelix_filtered <- pca_pseudohelix$x[,1:bs_pseudohelix$ndim] |> dplyr::as_tibble()

wedge_cluster_choice <- NbClust::NbClust(pca_wedge_filtered, method = "kmeans", index = "alllong")
pseudohelix_cluster_choice <- NbClust::NbClust(pca_pseudohelix_filtered, method = "kmeans", index = "alllong")

wedge_clusters <- bind_cols(dataset_number = 1:36,
                            cluster = factor(wedge_cluster_choice$Best.partition),
                            pca_wedge_filtered)
pseudohelix_clusters <- bind_cols(dataset_number = 1:36,
                            cluster = factor(pseudohelix_cluster_choice$Best.partition),
                            pca_pseudohelix_filtered)

plot_clusters <- function(data) {
  outline_pts <- tibble(
    x = c(0,250,250,0,0,
          0,-50,-50,0,
          -50,200,250)/33,
    y = c(0,0,1188,1188,0,
          0,50,1238,1188,
          1238,1238,1188)/33+0.5,
    id = c(1,1,1,1,1,
           2,2,2,2,
           3,3,3)/33
  )
  
  deltas <- tibble(
    x = c(250, rep(c(0,-250,-50,0,50,250), 35), 
          0,-50,-250,0,50),
    y = c(0, rep(c(33,0,50,-33,-50,33), 35), 
          33,50,0,-33,-50)
  )
  
  cluster_vec <- data$cluster
  
  polygon_pts <- purrr::accumulate(
    seq(1, nrow(deltas)), 
    \(acc_deltas, idx) {
      
      delta_next <- deltas[idx,]
      x <- acc_deltas$x + delta_next$x
      y <- acc_deltas$y + delta_next$y
      
      point <- tibble(x = x, y = y)
      return(point)
    },
    .init = tibble(x = 0, y = 0)
  )[-1,] |> 
    bind_cols(id = sort(rep(1:36, 6)), 
              cluster = {purrr::map(cluster_vec, \(x) rep(x, 6)) |> list_c()}) |> 
    mutate(x = x/33,
           y = y/33+0.5)
  
  palette <- "PNWColors::Starfish"
  long_p <- ggplot(polygon_pts, aes(x = x, y = y, fill = cluster, color = cluster, group = id)) +
    geom_polygon() +
    geom_path(data = outline_pts, aes(x = x, y = y), inherit.aes = FALSE) +
    scale_fill_manual(labels = 1:3, breaks = 1:3, values = c("#01608c", "#9462ff", "#ee8cab")) +  
    scale_color_manual(labels = 1:3, breaks = 1:3, values = c("#01608c", "#9462ff", "#ee8cab")) +  
    coord_cartesian(expand = FALSE) +
    scale_x_continuous(breaks = NULL) +
    scale_y_continuous(breaks = seq(1, 36, 2), position = "right") +
    labs(x = NULL, y = "Dataset Number") +
    ggtheme_light() +
    theme(legend.position = "none", panel.grid = element_blank(), panel.border = element_blank())
  envelope_p <- data |> 
    ggplot(aes(x = PC1, y = PC2, color = cluster, fill = cluster)) +
    geom_point(show.legend = FALSE) +
    ggforce::geom_mark_hull(expand = 0.01, radius = 0.01) +
    scale_fill_manual(labels = 1:3, breaks = 1:3, values = c("#01608c", "#9462ff", "#ee8cab")) +  
    scale_color_manual(labels = 1:3, breaks = 1:3, values = c("#01608c", "#9462ff", "#ee8cab")) +   
    ggtheme_light() +
    theme(legend.position = "none")
  legend <- get_legend(
    data |> 
      ggplot(aes(x = PC1, y = PC2, color = cluster, fill = cluster)) +
      geom_point(show.legend = FALSE) +
      ggforce::geom_mark_hull(expand = 0, radius = 0) +
      scale_fill_manual(labels = 1:3, breaks = 1:3, values = c("#01608c", "#9462ff", "#ee8cab")) +  
      scale_color_manual(labels = 1:3, breaks = 1:3, values = c("#01608c", "#9462ff", "#ee8cab")) +  
      labs(fill = "Cluster", color = "Cluster") +
      theme(legend.position = "right", legend.key.height = unit(4, "mm"), rect = element_blank())
  )
  
  right_col <- plot_grid(legend, long_p, ncol = 1, rel_heights = c(0.25,1))
  combined_plot <- cowplot::plot_grid(envelope_p, right_col, ncol = 2, rel_widths = c(1, 0.25), align = "h")
  
  return(combined_plot)
}

ggplots$Light$pca$wedge_clusters <- plot_clusters(wedge_clusters)
ggplots$Light$pca$pseudohelix_clusters <- plot_clusters(pseudohelix_clusters)
