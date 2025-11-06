# Function setup ----------------------------------------------------------

ggtheme_dark <- function() {
  list(
    ggplot2::theme_dark(base_size = 8, base_family = "ArialMT"), 
    ggplot2::theme(
      #Overall elements
      rect = ggplot2::element_blank(), 
      text = ggplot2::element_text(color = "white"), 
      line = ggplot2::element_line(color = "black"), 
      
      #Legend positioning
      legend.position = "inside", 
      legend.position.inside = c(0.85, 0.25), 
      
      #Manual override of desired ggplot2::theme elements
      panel.background = ggplot2::element_rect(fill = "black", color = "white"), 
      legend.key = ggplot2::element_blank(), #panel.background automatically maps to legend.key; I want to override this
      legend.background = ggplot2::element_rect(fill = "black", color = "white"), 
      strip.background = ggplot2::element_rect(fill = "black", color = "white"),
      axis.text = ggplot2::element_text(color = 'gray'),
      legend.key.height = unit(3, "mm")
    )
  )
} # Define a custom ggplot2::ggplot dark ggplot2::theme
ggtheme_light <- function() {
  list(
    ggplot2::theme_bw(base_size = 8, base_family = "ArialMT"), 
    ggplot2::theme(
      #Overall elements
      text = ggplot2::element_text(color = "black"), 
      
      #Legend positioning
      legend.position = "inside", 
      legend.position.inside = c(0.85, 0.25), 
      
      #Manual override of desired ggplot2::theme elements
      legend.background = ggplot2::element_rect(fill = "white", color = "black"), 
      strip.background = ggplot2::element_rect(fill = "white"), 
      plot.background = ggplot2::element_blank(),
      legend.key.height = unit(3, "mm")
    )
  )
} # Define a custom ggplot2::ggplot light ggplot2::theme
faceting <- function(facetVar, datasetType) {
  # Create a list containing all the necessary facet variables
  mapping_list <- list(
    Residue = c(
      "CO2" = "CO[2]", 
      "Ax"  = "H[2]*O[ax]", 
      "Eq"  = "H[2]*O[eq]", 
      "Glu" = "'Intact Glu30'", 
      "Oxy" = "'Dioxygen'"
    ),
    AtomPair = c(
      "Cu-Tyr" = "'Cu-Tyr168Oη'",
      "Cu-NTerm" = "Cu-N[term]",
      "Cu-His1ND1" = "'Cu-His1Nδ'",
      "Cu-His84NE2" = "'Cu-His84Nε'",
      "Cu-Eq" = "Cu-H[2]*O[eq]",
      "Cu-Ax" = "Cu-H[2]*O[ax]"
    ),
    AngleID = c(
      "T1" = "θ[1]", 
      "T2" = "θ[2]", 
      "T3" = "θ[3]", 
      "TT" = "θ[T]", 
      "THH" = "θ[H-H]", 
      "TH1" = "θ[H1]", 
      "THN" = "θ[HN]"
    )
  )
  
  list(
    # Dynamically create a faceting layer based on the input facetVar string
    ggplot2::facet_wrap(
      as.formula(paste("~", facetVar)), # Create a formula based on the input facetVar string
      scales = 'free',
      labeller = ggplot2::labeller(
        # Use a custom ggplot2::labeller from the mapping list; use tidy eval (!! to inject the user-supplied string and := to allow for LHS evaluation; R coerces the input string into a symbol)
        !!facetVar := ggplot2::as_labeller(mapping_list[[facetVar]], label_parsed), # 
        .default = label_parsed
      )),
    # Dynamically assign axis labels depending on the data type and facet used
    ggplot2::labs(
      x = if(datasetType == 'Pseudohelix') {
        "Diffraction-Weighted Dose (MGy)"
      } else if(datasetType == "Wedge") {
        "Wedge Number"
      },
      y = if(facetVar == "Residue") {
        "Occupancy"
      } else if(facetVar == 'AtomPair') {
        "Distance (Å)"
      } else if(facetVar == "AngleID") {
        "Angle (°)"
      }
    )
  )
}
scales <- function(facetVar) {
  
  # Manually alter the scales of individual facets
  angle_window <- 7
  occupancy_window <- 0.4
  bfactor_window <- 8
  distance_window <- 0.25
  
  ggh4x::facetted_pos_scales(
    y = if(facetVar == "AngleID") {
      list(
        ggplot2::scale_y_continuous(limits = c(94 - angle_window / 2, 94 + angle_window / 2)),
        ggplot2::scale_y_continuous(limits = c(94 - angle_window / 2, 94 + angle_window / 2)),
        ggplot2::scale_y_continuous(limits = c(171 - angle_window / 2, 171 + angle_window / 2)),
        ggplot2::scale_y_continuous(limits = c(0 - angle_window / 2, 0 + angle_window / 2)),
        ggplot2::scale_y_continuous(limits = c(68 - angle_window / 2, 68 + angle_window / 2)),
        ggplot2::scale_y_continuous(limits = c(-4 - angle_window / 2, -4 + angle_window / 2)),
        ggplot2::scale_y_continuous(limits = c(5 - angle_window / 2, 5 + angle_window / 2))
    )} else if(facetVar == "Residue") {
      list(
        ggplot2::scale_y_continuous(limits = c(0.6 - occupancy_window / 2, 0.6 + occupancy_window / 2)),
        ggplot2::scale_y_continuous(limits = c(0.55 - occupancy_window / 2, 0.55 + occupancy_window / 2)),
        ggplot2::scale_y_continuous(limits = c(0.35 - occupancy_window / 2, 0.35 + occupancy_window / 2)),
        ggplot2::scale_y_continuous(limits = c(0.55 - occupancy_window / 2, 0.55 + occupancy_window / 2)),
        ggplot2::scale_y_continuous(limits = c(0.65 - occupancy_window / 2, 0.65 + occupancy_window / 2))
    )} else if(facetVar == "AtomPair") {
      list(
        ggplot2::scale_y_continuous(limits = c(2.375 - distance_window / 2, 2.375 + distance_window / 2)),
        ggplot2::scale_y_continuous(limits = c(1.99 - distance_window / 2, 1.99 + distance_window / 2)),
        ggplot2::scale_y_continuous(limits = c(1.95 - distance_window / 2, 1.95 + distance_window / 2)),
        ggplot2::scale_y_continuous(limits = c(1.98 - distance_window / 2, 1.98 + distance_window / 2)),
        ggplot2::scale_y_continuous(limits = c(2.18 - distance_window / 2, 2.18 + distance_window / 2)),
        ggplot2::scale_y_continuous(limits = c(2.66 - distance_window / 2, 2.66 + distance_window / 2))
    )} else {NULL}) # If input facetVar isn't known, just automatically scale the axes
}

scatter_dark <- function(data, mapping, facetVar, datasetType) {
  ggplot2::ggplot(data, mapping) +
    scale_shape_manual(
      "Chain", 
      labels = c("A", "B"), 
      breaks = c("A", "B"),
      values = c(16, 17), 
    ) +
    ggplot2::scale_color_manual(
      "Chain", 
      labels = c("A", "B"), 
      breaks = c("A", "B"),
      values = c("#8100b6", "#5cb344"), 
    ) +
    geom_smooth(
      method = "lm", 
      linewidth = 0, 
      fill = "gray", 
      show.legend = FALSE,
      na.rm = TRUE,
      fullrange = TRUE
    ) + #Standard error plotting
    geom_smooth(
      method = "lm", 
      linetype = 5, 
      linewidth = 0.5,
      show.legend = FALSE,
      se = FALSE,
      na.rm = TRUE,
      fullrange = TRUE
    ) + #Linear regression line
    ggplot2::geom_point(
      size = 1, 
    ) + #Point for each occupancy value
    faceting(facetVar = facetVar, datasetType = datasetType) +
    ggtheme_dark() +
    coord_cartesian(expand = FALSE) +
    scales(facetVar = facetVar) +
    if(facetVar == "AngleID") {
      ggplot2::theme(legend.position.inside = c(0.85, 0.15))
    } else if(facetVar == "AtomPair") {
      ggplot2::theme(legend.position = "right")
    }
} #Make scatter plot with fitted linear regression
scatter_light <- function(data, mapping, facetVar, datasetType) {
  ggplot2::ggplot(data, mapping) +
    scale_shape_manual(
      "Chain", 
      labels = c("A", "B"), 
      breaks = c("A", "B"),
      values = c(16, 17), 
    ) +
    ggplot2::scale_color_manual(
      "Chain", 
      labels = c("A", "B"), 
      breaks = c("A", "B"),
      values = c("#8100b6", "#5cb344"), 
    ) +
    geom_smooth(
      method = "lm", 
      linewidth = 0, 
      fill = "gray", 
      show.legend = FALSE,
      na.rm = TRUE,
      fullrange = TRUE
    ) + #Standard error plotting
    geom_smooth(
      method = "lm", 
      linetype = 5, 
      linewidth = 0.5,
      show.legend = FALSE,
      se = FALSE,
      na.rm = TRUE,
      fullrange = TRUE
    ) + #Linear regression line
    ggplot2::geom_point(
      size = 1, 
    ) + #Point for each occupancy value
    faceting(facetVar = facetVar, datasetType = datasetType) +
    ggtheme_light() +
    coord_cartesian(expand = FALSE) +
    scales(facetVar = facetVar) +
    if(facetVar == "AngleID") {
      ggplot2::theme(legend.position.inside = c(0.85, 0.15))
    } else if(facetVar == "AtomPair") {
      ggplot2::theme(legend.position = "right")
    }
} #Make scatter plot with fitted linear regression

ggdarklight <- function(data, key) {
  yvar <- if(deparse(substitute(data)) == "stackedOccupancies") {
    'Occupancy'
  } else if(deparse(substitute(data)) == "stackedBFactors") {
    'bFactor'
  } else if(deparse(substitute(data)) == "stackedDistances") {
    'Distance'
  } else if(deparse(substitute(data)) == "stackedAngles") {
    'Angle'
  } else if(deparse(substitute(data)) == "dwds") {
    'Dataset Number'
  } else {stop("Invalid data input")}
  
  facetStr <- if(deparse(substitute(data)) == "stackedOccupancies") {
    "Residue"
  } else if(deparse(substitute(data)) == "stackedDistances") {
    "AtomPair"
  } else if(deparse(substitute(data)) == "stackedAngles") {
    "AngleID"
  } else if(deparse(substitute(data)) == 'dwds') {
    NA
  } else{stop("Invalid data input")}
  
  lightplots <- list(
    Pseudohelices = scatter_light(
      data = dplyr::filter(data$Pseudohelices, !is.na(!!yvar)),
      mapping = ggplot2::aes(
        x = Dose,
        y = .data[[yvar]],
        color = Molecule,
        shape = Molecule
      ),
      facetVar = facetStr,
      datasetType = "Pseudohelix"
    ),
    Wedges = scatter_light(
      data = dplyr::filter(data$Wedges, !is.na(!!yvar)),
      mapping = ggplot2::aes(
        x = WedgeNumber,
        y = .data[[yvar]],
        color = Molecule,
        shape = Molecule
      ),
      facetVar = facetStr,
      datasetType = "Wedge"
    )
  )
  
  darkplots <- list(
    Pseudohelices = scatter_dark(
      data = dplyr::filter(data$Pseudohelices, !is.na(!!yvar)),
      mapping = ggplot2::aes(
        x = Dose,
        y = .data[[yvar]],
        color = Molecule,
        shape = Molecule
      ),
      facetVar = facetStr,
      datasetType = "Pseudohelix"
    ),
    Wedges = scatter_dark(
      data = dplyr::filter(data$Wedges, !is.na(!!yvar)),
      mapping = ggplot2::aes(
        x = WedgeNumber,
        y = .data[[yvar]],
        color = Molecule,
        shape = Molecule
      ),
      facetVar = facetStr,
      datasetType = "Wedge"
    )
  )
  
  ggplots$Light[[key]] <<- lightplots
  ggplots$Dark[[key]] <<- darkplots
} #Combination of scatter_dark and scatter_light to ease the plotting of multiple datasets with multiple themes

# Scatter plots -----------------------------------------------------------

ggdarklight(stackedOccupancies, "Occupancies")
ggdarklight(stackedDistances, "Distances")
ggdarklight(stackedAngles, "Angles")

ggplots$Dark$Dose$DWDs <- ggplot2::ggplot(dwds, ggplot2::aes(x = datasetNumber, y = dwd_MGy, color = datasetType, shape = datasetType)) +
  geom_line(data = dplyr::filter(dwds, datasetType == "Wedges")) +
  ggplot2::geom_point(data = dplyr::filter(dwds, datasetType == "Wedges")) +
  geom_line(data = dplyr::filter(dwds, datasetType == "Pseudohelices")) +
  ggplot2::geom_point(data = dplyr::filter(dwds, datasetType == "Pseudohelices")) +
  ggtheme_dark() +
  ggplot2::scale_color_manual(
    "Dataset Type", 
    labels = c("Wedges", "Pseudohelices"), 
    breaks = c("Wedges", "Pseudohelices"), 
    values = c("#016ad6", "#c2db4d")
  ) +
  scale_shape_manual(
    "Dataset Type", 
    labels = c("Wedges", "Pseudohelices"), 
    breaks = c("Wedges", "Pseudohelices"), 
    values = c(16, 15)
  ) +
  coord_cartesian(expand = FALSE) +
  ggplot2::theme(
    legend.position.inside = c(0.75, 0.22)
  ) +
  ggplot2::labs(
    x = "Dataset Number",
    y = "Diffraction-Weighted Dose (MGy)"
  )

ggplots$Light$Dose$DWDs <- ggplot2::ggplot(dwds, ggplot2::aes(x = datasetNumber, y = dwd_MGy, color = datasetType, shape = datasetType)) +
  geom_line(data = dplyr::filter(dwds, datasetType == "Wedges")) +
  ggplot2::geom_point(data = dplyr::filter(dwds, datasetType == "Wedges")) +
  geom_line(data = dplyr::filter(dwds, datasetType == "Pseudohelices")) +
  ggplot2::geom_point(data = dplyr::filter(dwds, datasetType == "Pseudohelices")) +
  ggtheme_light() +
  ggplot2::scale_color_manual(
    "Dataset Type", 
    labels = c("Wedges", "Pseudohelices"), 
    breaks = c("Wedges", "Pseudohelices"), 
    values = c("#016ad6", "#c2db4d")
  ) +
  scale_shape_manual(
    "Dataset Type", 
    labels = c("Wedges", "Pseudohelices"), 
    breaks = c("Wedges", "Pseudohelices"), 
    values = c(16, 15)
  ) +
  coord_cartesian(expand = FALSE) +
  ggplot2::theme(
    legend.position.inside = c(0.75, 0.22)
  ) +
  ggplot2::labs(
    x = "Dataset Number",
    y = "Diffraction-Weighted Dose (MGy)"
  )

ggplots$Dark$Stats$CrystalStats <- ggplot2::ggplot(crystal_stats$combined, ggplot2::aes(x = dataset_number, y = value, color = dataset_type)) +
  ggplot2::geom_point() +
  ggtheme_dark() +
  ggplot2::scale_color_manual(
    "Dataset Type",
    labels = c("Wedges", "Pseudohelices"), 
    breaks = c("Wedges", "Pseudohelices"), 
    values = c("#016ad6", "#c2db4d")
  ) +
  ggplot2::facet_wrap(
    . ~ statistic,
    scales = "free",
    ncol = 3,
    labeller = ggplot2::as_labeller(
      c(
        cc1_2 = "CC[1/2]",
        multiplicity = "Multiplicity",
        r_free = "R[free]",
        r_work = "R[work]",
        completeness_percent = "'Completeness (%)'",
        mean_i_sigma_i = "'Mean '*I/σ[I]",
        wilson_b_factor = "'Wilson B-factor'",
        rms_angles = "RMS[angles]",
        rms_bonds = "RMS[bonds]",
        average_b_factor = "'Average B-factor'"
      ),
      label_parsed
    )
  ) +
  ggplot2::theme(axis.title.y = ggplot2::element_blank(), legend.position.inside = c(0.85, 0.1)) +
  ggplot2::labs(x = "Dataset Number")
ggplots$Light$Stats$CrystalStats <- ggplot2::ggplot(crystal_stats$combined, ggplot2::aes(x = dataset_number, y = value, color = dataset_type)) +
  ggplot2::geom_point() +
  ggtheme_light() +
  ggplot2::scale_color_manual(
    "Dataset Type",
    labels = c("Wedges", "Pseudohelices"), 
    breaks = c("Wedges", "Pseudohelices"), 
    values = c("#016ad6", "#c2db4d")
  ) +
  ggplot2::facet_wrap(
    . ~ statistic,
    scales = "free",
    ncol = 3,
    labeller = ggplot2::as_labeller(
      c(
        cc1_2 = "CC[1/2]",
        multiplicity = "Multiplicity",
        r_free = "R[free]",
        r_work = "R[work]",
        completeness_percent = "'Completeness (%)'",
        mean_i_sigma_i = "'Mean '*I/σ[I]",
        wilson_b_factor = "'Wilson B-factor'",
        rms_angles = "RMS[angles]",
        rms_bonds = "RMS[bonds]",
        average_b_factor = "'Average B-factor'"
      ),
      label_parsed
    )
  ) +
  ggplot2::theme(axis.title.y = ggplot2::element_blank(), legend.position.inside = c(0.85, 0.1)) +
  ggplot2::labs(x = "Dataset Number")

# RMSD Plotting -----------------------------------------------------------
rmsd_plots <- list()

base_plot <- function() {
  ggplot2::ggplot(all_rmsds, ggplot2::aes(x = ref_dataset, y = comp_dataset, fill = rmsd)) +
    ggplot2::geom_tile(
      data = all_rmsds %>% dplyr::filter(parameter == "occupancies"),
      mapping = ggplot2::aes(x = ref_dataset, y = comp_dataset, fill = rmsd)
    ) +
    ggplot2::scale_fill_viridis_c(name = "Occupancy RMSD", limits = c(0.01, 0.025), oob = scales::squish, n.breaks = 4) +
    ggnewscale::new_scale_fill() +
    ggplot2::geom_tile(
      data = all_rmsds %>% dplyr::filter(parameter == "b_factors"),
      mapping = ggplot2::aes(x = ref_dataset, y = comp_dataset, fill = rmsd)
    ) +
    ggplot2::scale_fill_viridis_c(name = "B-Factor RMSD", limits = c(0, 1.2), oob = scales::squish, breaks = c(0, 0.4, 0.8, 1.2)) +
    ggnewscale::new_scale_fill() +
    ggplot2::geom_tile(
      data = all_rmsds %>% dplyr::filter(parameter == "coordinates"),
      mapping = ggplot2::aes(x = ref_dataset, y = comp_dataset, fill = rmsd)
    ) +
    ggplot2::scale_fill_viridis_c(name = "Coordinate RMSD", limits = c(0.04, 0.16), oob = scales::squish, n.breaks = 4) +
    ggplot2::facet_grid(
      parameter ~ type,
      labeller = ggplot2::labeller(
        .default = str_to_title,
        parameter = c(b_factors = "B-Factors", coordinates = "Coordinates", occupancies = "Occupancy")
      )
    )
}

rmsd_plots$base_plots$light <- base_plot() +
  ggplot2::theme_bw(base_size = 8, base_family = "ArialMT") +
  ggplot2::theme(
    panel.spacing = unit(0, "mm"),
    text = ggplot2::element_text(color = "black"), 
    panel.border = ggplot2::element_rect(fill = NA, color = "black", linewidth = 1),
    strip.background = ggplot2::element_rect(fill = "white"), 
    plot.background = ggplot2::element_blank(),
    legend.position = "none"
  ) +
  ggplot2::scale_x_continuous(expand = c(0, 0), breaks = seq(0, 36, 3)) +
  ggplot2::scale_y_reverse(expand = c(0, 0), breaks = seq(0, 36, 3)) +
  ggplot2::labs(x = "Reference Dataset", y = "Comparison Dataset")

rmsd_plots$base_plots$dark <- base_plot() +
  ggtheme_dark() + 
  ggplot2::theme(
    panel.spacing = unit(0, "mm"),
    rect = ggplot2::element_blank(),
    text = ggplot2::element_text(color = "white"), 
    line = ggplot2::element_line(color = "black"), 
    panel.border = ggplot2::element_rect(fill = NA, color = "white", linewidth = 1),
    strip.background = ggplot2::element_rect(fill = "black", color = "white"),
    axis.text = ggplot2::element_text(color = 'gray'),
    axis.ticks = ggplot2::element_line(color = "gray"),
    legend.position = "none"
  ) +
  ggplot2::scale_x_continuous(expand = c(0, 0), breaks = seq(0, 36, 3)) +
  ggplot2::scale_y_reverse(expand = c(0, 0), breaks = seq(0, 36, 3)) +
  ggplot2::labs(x = "Reference Dataset", y = "Comparison Dataset")

make_legend <- function(parameter, theme) {
  p <- all_rmsds %>% dplyr::filter(parameter == !!parameter) %>% 
    ggplot2::ggplot(ggplot2::aes(x = ref_dataset, y = comp_dataset, fill = rmsd)) + 
    ggplot2::geom_tile()
  
  p <- switch(
    parameter,
    "occupancies" = p + ggplot2::scale_fill_viridis_c(name = "Occupancy", limits = c(0.01, 0.025), oob = scales::squish, n.breaks = 4),
    "b_factors" = p + ggplot2::scale_fill_viridis_c(name = "B-Factor", limits = c(0, 1.2), oob = scales::squish, breaks = c(0, 0.4, 0.8, 1.2)),
    "coordinates" = p + ggplot2::scale_fill_viridis_c(name = "Coordinate", limits = c(0.04, 0.16), oob = scales::squish, n.breaks = 4)
  )
  
  if(theme == "dark") {
    p <- p + ggplot2::theme(
      text = ggplot2::element_text(color = "white"),
      rect = ggplot2::element_blank()
    )
  }
  
  legend <- cowplot::get_legend(p)
  
  return(legend)
}

rmsd_plots$legends$light$occ <- make_legend("occupancies", "light")
rmsd_plots$legends$dark$occ <- make_legend("occupancies", "dark")
rmsd_plots$legends$light$bfact <- make_legend("b_factors", "light")
rmsd_plots$legends$dark$bfact <- make_legend("b_factors", "dark")
rmsd_plots$legends$light$coord <- make_legend("coordinates", "light")
rmsd_plots$legends$dark$coord <- make_legend("coordinates", "dark")

rmsd_plots$legends$light$combined <- cowplot::plot_grid(rmsd_plots$legends$light$bfact, rmsd_plots$legends$light$occ, rmsd_plots$legends$light$coord, ncol = 1) +
  cowplot::draw_label("RMSD", fontface = "bold", y = 0.99, vjust = 1)
rmsd_plots$legends$dark$combined <- cowplot::plot_grid(rmsd_plots$legends$dark$bfact, rmsd_plots$legends$dark$occ, rmsd_plots$legends$dark$coord, ncol = 1) +
  cowplot::draw_label("RMSD", fontface = "bold", y = 0.99, vjust = 1, color = "white")

ggplots$Light$Comparisons$RMSDs <- suppressWarnings(cowplot::plot_grid(rmsd_plots$base_plots$light, rmsd_plots$legends$light$combined, ncol = 2, rel_widths = c(1, 0.15)))
ggplots$Dark$Comparisons$RMSDs <- suppressWarnings(cowplot::plot_grid(rmsd_plots$base_plots$dark, rmsd_plots$legends$dark$combined, ncol = 2, rel_widths = c(1, 0.15)))

# Trend plotting ----------------------------------------------------------

dark.trend <- function(trend, datasetType) {
  regressor <- ifelse(datasetType == "Pseudohelices", "MGy", "Wedge Number")
  
  ggplot2::ggplot(
    dplyr::filter(longData[[datasetType]], Estimate != "Contrast", Measurement == trend), # Plot only trends A and B (exclude contrast coefficients)
    ggplot2::aes(x = Residue, y = Coefficient) # Start off with inverted axes to allow for asterisk offset
  ) + 
    ggplot2::geom_hline(yintercept = 0, color = 'gray') + # Vertical line to show zero mark
    ggplot2::geom_segment( # Line portion of barbell
      data = dplyr::filter(wideData[[datasetType]], Measurement == trend),
      ggplot2::aes(
        x = Residue, 
        xend = Residue, 
        y = Coefficient_TrendA, 
        yend = Coefficient_TrendB, 
        color = PValue_Contrast <= 0.05
      ), 
      show.legend = FALSE, 
      inherit.aes = FALSE, 
      linewidth = 3,
      alpha = 0.5
    ) +
    ggplot2::scale_color_manual( # Significant contrasts get colored light blue
      values = c("gray", "skyblue2")
    ) +
    ggnewscale::new_scale_color() + # Need new color scale for dots
    ggplot2::scale_color_manual(
      "Chain", 
      labels = c("A", "B"), 
      breaks = c("TrendA", "TrendB"),
      values = c("#8100b6", "#5cb344"), 
    ) +
    ggtheme_dark() +
    ggplot2::geom_point( # Blocks out ggplot2::geom_segment to allow for transparent dots
      color = "black", 
      size = 3,
    ) +
    ggplot2::geom_point( # Individual trend plotting
      ggplot2::aes(color = Estimate), 
      size = 3,
      alpha = 0.6
    ) +
    ggplot2::geom_text( # Asterisks for significance
      data = dplyr::filter(longData[[datasetType]], Estimate == "TrendA", Measurement == trend),
      ggplot2::aes(
        label = Significance,
        x = Residue, # Offset vertically for only trend A (proportional to number of data points)
        y = Coefficient,
      ),
      size = 3, 
      position = ggplot2::position_nudge(x = 0.25),
      color = '#bd9adb',
      inherit.aes = FALSE
    ) +
    ggplot2::geom_text( # Asterisks for significance
      data = dplyr::filter(longData[[datasetType]], Estimate == "TrendB", Measurement == trend),
      ggplot2::aes(
        label = Significance,
        x = Residue, # Offset vertically for only trend A (proportional to number of data points)
        y = Coefficient,
      ),
      size = 3, 
      position = ggplot2::position_nudge(x = 0.3),
      color = '#b2deab',
      inherit.aes = FALSE
    ) +
    ggplot2::theme(
      legend.position.inside = c(0.87, 0.2), 
      panel.grid.major.y = ggplot2::element_blank()
    ) +
    ggplot2::scale_x_discrete(labels = c(
      "T1" = bquote(θ[1]), 
      "T2" = bquote(θ[2]), 
      "T3" = bquote(θ[3]), 
      "TH1" = bquote(θ[H1]), 
      "THN" = bquote(θ[HN]), 
      "THH" = bquote(θ[H-H]), 
      "TT" = bquote(θ[T]), 
      "CuNterm" = bquote(Cu-N[term]),
      "CuAx" = bquote(Cu-H[2]*O[ax]), 
      "CuEq" = bquote(Cu-H[2]*O[eq]), 
      "CuHis1ND1" = bquote(Cu-His1Nδ), 
      "CuHis84NE2" = bquote(Cu-His84Nε),
      "CuTyr" = bquote(Cu-Tyr168Oη),
      "Ax" = bquote(H[2]*O[ax]), 
      "Eq" = bquote(H[2]*O[eq]), 
      "CO2" = bquote(CO[2]), 
      "Glu" = bquote(Intact*" "*Glu30)
    )) +
    ggplot2::labs(
      x = if(trend %in% c("Occupancies", "BFactors")) {
        "Residue"
      } else if(trend == "Distances") {
        "Atom Pair"
      } else if(trend == "Angles") {
        "Angle ID"
      } else {stop("Invalid trend input for trend visualization")},
      y = if(trend == "Occupancies") {
        stringr::str_glue("Occupancy Trend (Δ/{regressor})")
      } else if(trend == "Distances") {
        stringr::str_glue("Distance Trend (ΔÅ/{regressor})")
      } else if(trend == "Angles") {
        stringr::str_glue("Angle Trend (Δ°/{regressor})")
      } else {stop("Invalid trend input for trend visualization")}
    ) +
    ggplot2::coord_flip() # Flip coordinates back
}
light.trend <- function(trend, datasetType) {
  regressor <- ifelse(datasetType == "Pseudohelices", "MGy", "Wedge Number")
  
  ggplot2::ggplot(
    dplyr::filter(longData[[datasetType]], Estimate != "Contrast", Measurement == trend), # Plot only trends A and B (exclude contrast coefficients)
    ggplot2::aes(x = Residue, y = Coefficient) # Start off with inverted axes to allow for asterisk offset
  ) + 
    ggplot2::geom_hline(yintercept = 0, color = 'gray') + # Vertical line to show zero mark
    ggplot2::geom_segment( # Line portion of barbell
      data = dplyr::filter(wideData[[datasetType]], Measurement == trend),
      ggplot2::aes(
        x = Residue, 
        xend = Residue, 
        y = Coefficient_TrendA, 
        yend = Coefficient_TrendB, 
        color = PValue_Contrast <= 0.05
      ), 
      show.legend = FALSE, 
      inherit.aes = FALSE, 
      linewidth = 3,
      alpha = 0.5
    ) +
    ggplot2::scale_color_manual( # Significant contrasts get colored light blue
      breaks = c(FALSE, TRUE),
      values = c("gray", "skyblue2")
    ) +
    ggnewscale::new_scale_color() + # Need new color scale for dots
    ggplot2::scale_color_manual(
      "Chain", 
      labels = c("A", "B"), 
      breaks = c("TrendA", "TrendB"),
      values = c("#8100b6", "#5cb344"), 
    ) +
    ggtheme_light() +
    ggplot2::geom_point( # Blocks out ggplot2::geom_segment to allow for transparent dots
      color = "white", 
      size = 3,
    ) +
    ggplot2::geom_point( # Individual trend plotting
      ggplot2::aes(color = Estimate), 
      size = 3,
      alpha = 0.6
    ) +
    ggplot2::geom_text( # Asterisks for significance
      data = dplyr::filter(longData[[datasetType]], Estimate == "TrendA", Measurement == trend),
      ggplot2::aes(
        label = Significance,
        x = Residue, # Offset vertically for only trend A (proportional to number of data points)
        y = Coefficient,
      ),
      size = 3, 
      position = ggplot2::position_nudge(x = 0.25),
      color = '#8100b6',
      inherit.aes = FALSE
    ) +
    ggplot2::geom_text( # Asterisks for significance
      data = dplyr::filter(longData[[datasetType]], Estimate == "TrendB", Measurement == trend),
      ggplot2::aes(
        label = Significance,
        x = Residue, # Offset vertically for only trend A (proportional to number of data points)
        y = Coefficient,
      ),
      size = 3, 
      position = ggplot2::position_nudge(x = 0.3),
      color = '#4f9437',
      inherit.aes = FALSE
    ) +
    ggplot2::theme(
      legend.position.inside = c(0.87, 0.2), 
      panel.grid.major.y = ggplot2::element_blank()
    ) +
    ggplot2::scale_x_discrete(labels = c(
      "T1" = bquote(θ[1]), 
      "T2" = bquote(θ[2]), 
      "T3" = bquote(θ[3]), 
      "TH1" = bquote(θ[H1]), 
      "THN" = bquote(θ[HN]), 
      "THH" = bquote(θ[H-H]), 
      "TT" = bquote(θ[T]), 
      "CuNterm" = bquote(Cu-N[term]),
      "CuAx" = bquote(Cu-H[2]*O[ax]), 
      "CuEq" = bquote(Cu-H[2]*O[eq]), 
      "CuHis1ND1" = bquote(Cu-His1Nδ), 
      "CuHis84NE2" = bquote(Cu-His84Nε),
      "CuTyr" = bquote(Cu-Tyr168Oη),
      "Ax" = bquote(H[2]*O[ax]), 
      "Eq" = bquote(H[2]*O[eq]), 
      "CO2" = bquote(CO[2]), 
      "Glu" = bquote(Intact*" "*Glu30)
    )) +
    ggplot2::labs(
      x = if(trend %in% c("Occupancies", "BFactors")) {
        "Residue"
      } else if(trend == "Distances") {
        "Atom Pair"
      } else if(trend == "Angles") {
        "Angle ID"
      } else {stop("Invalid trend input for trend visualization")},
      y = if(trend == "Occupancies") {
        stringr::str_glue("Occupancy Trend (Δ/{regressor})")
      } else if(trend == "Distances") {
        stringr::str_glue("Distance Trend (ΔÅ/{regressor})")
      } else if(trend == "Angles") {
        stringr::str_glue("Angle Trend (Δ°/{regressor})")
      } else {stop("Invalid trend input for trend visualization")}
    ) +
    ggplot2::coord_flip()
}
darklighttrend <- function(trend) {
  
  ggplots$Light[[trend]]$Trends$Pseudohelices <<- light.trend(
    trend = trend,
    datasetType = "Pseudohelices"
  )
  
  ggplots$Dark[[trend]]$Trends$Pseudohelices <<- dark.trend(
    trend = trend,
    datasetType = "Pseudohelices"
  )
  
  ggplots$Light[[trend]]$Trends$Wedges <<- light.trend(
    trend = trend,
    datasetType = "Wedges"
  )
  
  ggplots$Dark[[trend]]$Trends$Wedges <<- dark.trend(
    trend = trend,
    datasetType = "Wedges"
  )
}

darklighttrend("Occupancies")
darklighttrend("Distances")
darklighttrend("Angles")
