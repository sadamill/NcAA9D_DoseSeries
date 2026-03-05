
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

plot_clusters <- function(data, dose = pseudohelixDose) {
  
  arg <- deparse(substitute(data))
  
  if(arg == "wedge_clusters") {
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
    shadow_pts <- tibble(
      x = c(0, 0, -50, -50,
            0, 250, 200, -50)/33,
      y = c(0, 1188, 1238, 50,
            1188, 1188, 1238, 1238)/33+0.5,
      id = c(1, 1, 1, 1,
             2, 2, 2, 2)
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
    
    long_p <- ggplot(polygon_pts, aes(x = x, y = y, fill = cluster, color = cluster, group = id)) +
      geom_polygon() +
      geom_polygon(data = shadow_pts, aes(x = x, y = y, group = id), inherit.aes = FALSE, color = "black", alpha = 0.4) +
      geom_path(data = outline_pts, aes(x = x, y = y), inherit.aes = FALSE) +
      scale_fill_manual(labels = 1:3, breaks = 1:3, values = c("#01608c", "#9462ff", "#ee8cab")) +  
      scale_color_manual(labels = 1:3, breaks = 1:3, values = c("#01608c", "#9462ff", "#ee8cab")) +  
      coord_cartesian(expand = FALSE) +
      scale_x_continuous(breaks = NULL) +
      scale_y_continuous(breaks = seq(1, 36, 2), position = "right") +
      labs(x = NULL, y = "Dataset Number") +
      ggtheme_light() +
      theme(legend.position = "none", panel.grid = element_blank(), panel.border = element_blank())
    
    xlab <- "PC1 (42.7%)"
    ylab <- "PC2 (5.5%)"
  }
  if(arg == "pseudohelix_clusters") {
    
    line <- tibble(
      x = c(0.75, 1.25, 1, 1, 0.75, 1.25,
            1.75, 2.25, 2, 2, 1.75, 2.25,
            2.75, 3.25, 3, 3, 2.75, 3.25),
      y = c(1, 1, 1, 9, 9, 9,
            1, 1, 1, 9, 9, 9,
            1, 1, 1, 9, 9, 9),
      id = c(1, 1, 2, 2, 3, 3,
             4, 4, 5, 5, 6, 6,
             7, 7, 8, 8, 9, 9)
    )
    
    points <- tibble(
      x = as.numeric(data$cluster),
      y = dose,
      cluster = data$cluster
    )
    
    long_p <- ggplot(points, aes(x = x, y = y, fill = cluster, color = cluster)) +
      geom_path(data = line, aes(x = x, y = y, group = id), inherit.aes = FALSE, color = "black") +
      geom_point(size = 2) +
      coord_cartesian(xlim = c(0, 4)) +
      scale_fill_manual(labels = 1:3, breaks = 1:3, values = c("#01608c75", "#9462ff75", "#ee8cab75")) +  
      scale_color_manual(labels = 1:3, breaks = 1:3, values = c("#01608c75", "#9462ff75", "#ee8cab75")) +  
      scale_x_continuous(breaks = NULL) +
      scale_y_continuous(breaks = 0:10, position = "right") +
      labs(x = NULL, y = "Average DDWD") +
      ggtheme_light() +
      theme(legend.position = "none", panel.border = element_blank(), axis.ticks = element_blank())
    
    xlab <- "PC1 (16.6%)"
    ylab <- "PC2 (12.1%)"
  }
  
  envelope_p <- data |> 
    ggplot(aes(x = PC1, y = PC2, color = cluster, fill = cluster)) +
    geom_point() +
    ggforce::geom_mark_hull(expand = 0.015, radius = 0.015) +
    scale_fill_manual(labels = 1:3, breaks = 1:3, values = c("#01608c", "#9462ff", "#ee8cab")) +  
    scale_color_manual(labels = 1:3, breaks = 1:3, values = c("#01608c", "#9462ff", "#ee8cab")) +   
    ggtheme_light() +
    theme(legend.position = "none") +
    labs(x = xlab, y = ylab)
  legend <- get_legend(
    data |> 
      ggplot(aes(x = PC1, y = PC2, color = cluster, fill = cluster)) +
      ggforce::geom_mark_hull() +
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

pseudohelix_pca_3d <- plot_ly(pseudohelix_clusters, x = ~PC1, y = ~PC2, z = ~PC3, color = ~cluster, colors = c("#01608c", "#9462ff", "#ee8cab")) |> 
  add_markers() |> 
  layout(scene = list(camera = list(projection = list(type = "orthographic")),
                      xaxis = list(title = "PC1 (16.6%)"),
                      yaxis = list(title = "PC2 (12.1%)"),
                      zaxis = list(title = "PC3 (8.8%)")))

htmlwidgets::saveWidget(as_widget(pseudohelix_pca_3d), "8-structure_analysis/output/plots/Light/pca_pseudohelix_clusters_3d.html")
system("rm -r 8-structure_analysis/output/plots/Light/pca_pseudohelix_clusters_3d_files/")
