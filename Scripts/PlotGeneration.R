# Function setup ----------------------------------------------------------

ggtheme_dark <- function() {
  list(
    scale_color_manual(
      "Molecule", 
      labels = c("A", "B"), 
      values = c("#5cb344", "#8100b6")
    ), 
    theme_dark(), 
    theme(
      #Overall elements
      rect = element_blank(), 
      text = element_text(color = "white"), 
      line = element_line(color = "black"), 
      
      #Legend positioning
      legend.position = "inside", 
      legend.position.inside = c(0.85, 0.25), 
      
      #Manual override of desired theme elements
      panel.background = element_rect(fill = "black", color = "white"), 
      legend.key = element_blank(), #panel.background automatically maps to legend.key; I want to override this
      legend.background = element_rect(fill = "black", color = "white"), 
      strip.background = element_rect(fill = "black", color = "white"),
      axis.text = element_text(color = 'gray')
    )
  )
} # Define a custom ggplot dark theme
ggtheme_light <- function() {
  list(
    scale_color_manual(
      "Molecule", 
      labels = c("A", "B"), 
      values = c("#5cb344", "#8100b6"), 
    ), 
    theme_bw(), 
    theme(
      #Overall elements
      text = element_text(color = "black"), 
      
      #Legend positioning
      legend.position = "inside", 
      legend.position.inside = c(0.85, 0.25), 
      
      #Manual override of desired theme elements
      legend.background = element_rect(fill = "white", color = "black"), 
      strip.background = element_rect(fill = "white"), 
      plot.background = element_blank()
    )
  )
} # Define a custom ggplot light theme
faceting <- function(facetVar, datasetType) {
  # Create a list containing all the necessary facet variables
  mapping_list <- list(
    Residue = c(
      "CO2" = "CO[2]", 
      "Ax"  = "H[2]*O[Ax]", 
      "Eq"  = "H[2]*O[Eq]", 
      "Glu" = "Glu[30]", 
      "Oxy" = "'Dioxygen'"
    ),
    AtomPair = c(
      "Cu-Tyr" = "Cu-Tyr[168]",
      "Cu-NTerm" = "Cu-N[term]",
      "Cu-His1ND1" = "Cu-His[1]*Nδ[1]",
      "Cu-His84NE2" = "Cu-His[84]*Nε[2]",
      "Cu-Eq" = "Cu-H[2]*O[Eq]",
      "Cu-Ax" = "Cu-H[2]*O[Ax]"
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
    facet_wrap(
      as.formula(paste("~", facetVar)), # Create a formula based on the input facetVar string
      scales = 'free',
      labeller = labeller(
        # Use a custom labeller from the mapping list; use tidy eval (!! to inject the user-supplied string and := to allow for LHS evaluation; R coerces the input string into a symbol)
        !!facetVar := as_labeller(mapping_list[[facetVar]], label_parsed), # 
        .default = label_parsed
      )),
    # Dynamically assign axis labels depending on the data type and facet used
    labs(
      x = if(datasetType == 'Pseudohelix') {
        "Dose (MGy)"
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
  
  facetted_pos_scales(
    y = if(facetVar == "AngleID") {
      list(
        scale_y_continuous(limits = c(94 - angle_window / 2, 94 + angle_window / 2)),
        scale_y_continuous(limits = c(94 - angle_window / 2, 94 + angle_window / 2)),
        scale_y_continuous(limits = c(171 - angle_window / 2, 171 + angle_window / 2)),
        scale_y_continuous(limits = c(180 - angle_window / 2, 180 + angle_window / 2)),
        scale_y_continuous(limits = c(68 - angle_window / 2, 68 + angle_window / 2)),
        scale_y_continuous(limits = c(176 - angle_window / 2, 176 + angle_window / 2)),
        scale_y_continuous(limits = c(5 - angle_window / 2, 5 + angle_window / 2))
    )} else if(facetVar == "Residue") {
      list(
        scale_y_continuous(limits = c(0.6 - occupancy_window / 2, 0.6 + occupancy_window / 2)),
        scale_y_continuous(limits = c(0.55 - occupancy_window / 2, 0.55 + occupancy_window / 2)),
        scale_y_continuous(limits = c(0.35 - occupancy_window / 2, 0.35 + occupancy_window / 2)),
        scale_y_continuous(limits = c(0.55 - occupancy_window / 2, 0.55 + occupancy_window / 2)),
        scale_y_continuous(limits = c(0.65 - occupancy_window / 2, 0.65 + occupancy_window / 2))
    )} else if(facetVar == "AtomPair") {
      list(
        scale_y_continuous(limits = c(2.375 - distance_window / 2, 2.375 + distance_window / 2)),
        scale_y_continuous(limits = c(1.99 - distance_window / 2, 1.99 + distance_window / 2)),
        scale_y_continuous(limits = c(1.95 - distance_window / 2, 1.95 + distance_window / 2)),
        scale_y_continuous(limits = c(1.98 - distance_window / 2, 1.98 + distance_window / 2)),
        scale_y_continuous(limits = c(2.18 - distance_window / 2, 2.18 + distance_window / 2)),
        scale_y_continuous(limits = c(2.66 - distance_window / 2, 2.66 + distance_window / 2))
    )} else {NULL}) # If input facetVar isn't known, just automatically scale the axes
}

scatter_dark <- function(data, mapping, facetVar, datasetType) {
  ggplot(data, mapping) +
    stat_smooth(
      method = "lm", 
      linewidth = 0, 
      fill = "gray", 
      show.legend = FALSE,
      na.rm = TRUE,
      fullrange = TRUE
    ) + #Standard error plotting
    stat_smooth(
      method = "lm", 
      linetype = 2, 
      se = FALSE,
      na.rm = TRUE,
      fullrange = TRUE
    ) + #Linear regression line
    geom_point(
      size = 0.3, 
      show.legend = FALSE
    ) + #Point for each occupancy value
    faceting(facetVar = facetVar, datasetType = datasetType) +
    ggtheme_dark() +
    coord_cartesian(expand = FALSE) +
    scales(facetVar = facetVar) +
    if(facetVar == "AngleID") {
      theme(legend.position.inside = c(0.85, 0.15))
    } else if(facetVar == "AtomPair") {
      theme(legend.position = "right")
    }
} #Make scatter plot with fitted linear regression
scatter_light <- function(data, mapping, facetVar, datasetType) {
  ggplot(data, mapping) +
    stat_smooth(
      method = "lm", 
      linewidth = 0, 
      fill = "gray", 
      show.legend = FALSE,
      na.rm = TRUE,
      fullrange = TRUE
    ) + #Standard error plotting
    stat_smooth(
      method = "lm", 
      linetype = 2, 
      se = FALSE,
      na.rm = TRUE,
      fullrange = TRUE
    ) + #Linear regression line
    geom_point(
      size = 0.3, 
      show.legend = FALSE
    ) + #Point for each occupancy value
    faceting(facetVar = facetVar, datasetType = datasetType) +
    ggtheme_light() +
    coord_cartesian(expand = FALSE) +
    scales(facetVar = facetVar) +
    if(facetVar == "AngleID") {
      theme(legend.position.inside = c(0.85, 0.15))
    } else if(facetVar == "AtomPair") {
      theme(legend.position = "right")
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
      data = data$Pseudohelices,
      mapping = aes(
        x = Dose,
        y = .data[[yvar]],
        color = Molecule
      ),
      facetVar = facetStr,
      datasetType = "Pseudohelix"
    ),
    Wedges = scatter_light(
      data = data$Wedges,
      mapping = aes(
        x = WedgeNumber,
        y = .data[[yvar]],
        color = Molecule
      ),
      facetVar = facetStr,
      datasetType = "Pseudohelix"
    )
  )
  
  darkplots <- list(
    Pseudohelices = scatter_dark(
      data = data$Pseudohelices,
      mapping = aes(
        x = Dose,
        y = .data[[yvar]],
        color = Molecule
      ),
      facetVar = facetStr,
      datasetType = "Pseudohelix"
    ),
    Wedges = scatter_dark(
      data = data$Wedges,
      mapping = aes(
        x = WedgeNumber,
        y = .data[[yvar]],
        color = Molecule
      ),
      facetVar = facetStr,
      datasetType = "Pseudohelix"
    )
  )
  
  ggplots$Light[[key]] <<- lightplots
  ggplots$Dark[[key]] <<- darkplots
} #Combination of scatter_dark and scatter_light to ease the plotting of multiple datasets with multiple themes

# Scatter plots -----------------------------------------------------------

ggdarklight(stackedOccupancies, "Occupancies")
ggdarklight(stackedDistances, "Distances")
ggdarklight(stackedAngles, "Angles")

ggplots$Dark$Dose$DWDs <- ggplot(dwds, aes(x = datasetNumber, y = dwd_MGy, color = datasetType)) +
  geom_point() +
  ggtheme_dark() +
  scale_color_manual(
    "Dataset Type", 
    labels = c("Wedges", "Pseudohelices"), 
    values = c("#0096c5", "#b8008c")
  ) +
  coord_cartesian(expand = FALSE) +
  labs(
    x = "Dataset Number",
    y = "Density-Weighted Dose (MGy)"
  )

ggplots$Light$Dose$DWDs <- ggplot(dwds, aes(x = datasetNumber, y = dwd_MGy, color = datasetType)) +
  geom_point() +
  ggtheme_light() +
  scale_color_manual(
    "Dataset Type", 
    labels = c("Wedges", "Pseudohelices"), 
    values = c("#0096c5", "#b8008c")
  ) +
  coord_cartesian(expand = FALSE) +
  labs(
    x = "Dataset Number",
    y = "Density-Weighted Dose (MGy)"
  )

ggplots$Dark$Stats$CrystalStats <- ggplot(crystal_stats$combined, aes(x = dataset_number, y = value, color = dataset_type)) +
  geom_point() +
  ggtheme_dark() +
  scale_color_manual(
    "Dataset Type",
    labels = c("Wedges", "Pseudohelices"), 
    values = c("#0096c5", "#b8008c")
  ) +
  facet_wrap(
    . ~ statistic,
    scales = "free",
    ncol = 3,
    labeller = as_labeller(
      c(
        cc1_2 = "CC[1/2]",
        multiplicity = "Multiplicity",
        r_free = "R[free]",
        r_work = "R[work]",
        completeness_percent = "'Completeness (%)'",
        mean_i_sigma_i = "'Mean '*I/σ[I]",
        wilson_b_factor = "'Wilson B-factor'",
        rms_angles = "RMS[angles]",
        rms_bonds = "RMS[bonds]"
      ),
      label_parsed
    )
  ) +
  theme(axis.title.y = element_blank(), legend.position.inside = c(0.85, 0.1)) +
  labs(x = "Dataset Number")
ggplots$Light$Stats$CrystalStats <- ggplot(crystal_stats$combined, aes(x = dataset_number, y = value, color = dataset_type)) +
  geom_point() +
  ggtheme_light() +
  scale_color_manual(
    "Dataset Type",
    labels = c("Wedges", "Pseudohelices"), 
    values = c("#0096c5", "#b8008c")
  ) +
  facet_wrap(
    . ~ statistic,
    scales = "free",
    ncol = 3,
    labeller = as_labeller(
      c(
        cc1_2 = "CC[1/2]",
        multiplicity = "Multiplicity",
        r_free = "R[free]",
        r_work = "R[work]",
        completeness_percent = "'Completeness (%)'",
        mean_i_sigma_i = "'Mean '*I/σ[I]",
        wilson_b_factor = "'Wilson B-factor'",
        rms_angles = "RMS[angles]",
        rms_bonds = "RMS[bonds]"
      ),
      label_parsed
    )
  ) +
  theme(axis.title.y = element_blank(), legend.position.inside = c(0.85, 0.1)) +
  labs(x = "Dataset Number")

# Trend plotting ----------------------------------------------------------

dark.trend <- function(trend, datasetType) {
  ggplot(
    filter(longData[[datasetType]], Estimate != "Contrast", Measurement == trend), # Plot only trends A and B (exclude contrast coefficients)
    aes(x = Residue, y = Coefficient) # Start off with inverted axes to allow for asterisk offset
  ) + 
    geom_hline(yintercept = 0, color = 'gray') + # Vertical line to show zero mark
    geom_segment( # Line portion of barbell
      data = filter(wideData[[datasetType]], Measurement == trend),
      aes(
        x = Residue, 
        xend = Residue, 
        y = Coefficient_TrendA, 
        yend = Coefficient_TrendB, 
        color = PValue_Contrast <= 0.05
      ), 
      show.legend = FALSE, 
      inherit.aes = FALSE, 
      linewidth = 5,
      alpha = 0.5
    ) +
    scale_color_manual( # Significant contrasts get colored light blue
      values = c("gray", "skyblue2"), 
    ) +
    new_scale_color() + # Need new color scale for dots
    ggtheme_dark() +
    geom_point( # Blocks out geom_segment to allow for transparent dots
      color = "black", 
      size = 5,
    ) +
    geom_point( # Individual trend plotting
      aes(color = Estimate), 
      size = 5,
      alpha = 0.6
    ) +
    geom_text( # Asterisks for significance
      data = filter(longData[[datasetType]], Estimate == "TrendA", Measurement == trend),
      aes(
        label = Significance,
        x = Residue, # Offset vertically for only trend A (proportional to number of data points)
        y = Coefficient,
      ),
      size = 6, 
      position = position_nudge(x = 0.1),
      color = '#B2DEAB',
      inherit.aes = FALSE
    ) +
    geom_text( # Asterisks for significance
      data = filter(longData[[datasetType]], Estimate == "TrendB", Measurement == trend),
      aes(
        label = Significance,
        x = Residue, # Offset vertically for only trend A (proportional to number of data points)
        y = Coefficient,
      ),
      size = 6, 
      position = position_nudge(x = 0.2),
      color = '#BD9ADB',
      inherit.aes = FALSE
    ) +
    theme(
      legend.position.inside = c(0.9, 0.15), 
      panel.grid.major.y = element_blank()
    ) +
    scale_x_discrete(labels = c(
      "T1" = bquote(θ[1]), 
      "T2" = bquote(θ[2]), 
      "T3" = bquote(θ[3]), 
      "TH1" = bquote(θ[H1]), 
      "THN" = bquote(θ[HN]), 
      "THH" = bquote(θ[H-H]), 
      "TT" = bquote(θ[T]), 
      "CuAx" = bquote(Cu-H[2]*O[Ax]), 
      "CuNterm" = bquote(Cu-N[term]),
      "CuEq" = bquote(Cu-H[2]*O[Eq]), 
      "CuHis1ND1" = bquote(Cu-His[1]*Nδ[1]), 
      "CuHis84NE2" = bquote(Cu-His[84]*Nε[2]),
      "CuTyr" = bquote(Cu-Tyr[168]), 
      "Ax" = bquote(H[2]*O[Ax]), 
      "Eq" = bquote(H[2]*O[Eq]), 
      "CO2" = bquote(CO[2]), 
      "Glu" = bquote(Glu[30])
    )) +
    labs(
      y = if(trend %in% c("Occupancies", "BFactors")) {
        "Residue"
      } else if(trend == "Distances") {
        "Atom Pair"
      } else if(trend == "Angles") {
        "Angle ID"
      } else {stop("Invalid trend input for trend visualization")},
      x = if(trend == "Occupancies") {
        "Occupancy Trend (Δ/MGy)"
      } else if(trend == "BFactors") {
        bquote("B-Factor Trend (" * Å^2 * "/MGy)")
      } else if(trend == "Distances") {
        "Distance Trend (ΔÅ/MGy)"
      } else if(trend == "Angles") {
        "Angle Trend (Δ°/MGy)"
      } else {stop("Invalid trend input for trend visualization")}
    ) +
    coord_flip() # Flip coordinates back
}
light.trend <- function(trend, datasetType) {
  ggplot(
    filter(longData[[datasetType]], Estimate != "Contrast", Measurement == trend), # Plot only trends A and B (exclude contrast coefficients)
    aes(x = Residue, y = Coefficient) # Start off with inverted axes to allow for asterisk offset
  ) + 
    geom_hline(yintercept = 0, color = 'gray') + # Vertical line to show zero mark
    geom_segment( # Line portion of barbell
      data = filter(wideData[[datasetType]], Measurement == trend),
      aes(
        x = Residue, 
        xend = Residue, 
        y = Coefficient_TrendA, 
        yend = Coefficient_TrendB, 
        color = PValue_Contrast <= 0.05
      ), 
      show.legend = FALSE, 
      inherit.aes = FALSE, 
      linewidth = 5,
      alpha = 0.5
    ) +
    scale_color_manual( # Significant contrasts get colored light blue
      values = c("gray", "skyblue2"), 
    ) +
    new_scale_color() + # Need new color scale for dots
    ggtheme_light() +
    geom_point( # Blocks out geom_segment to allow for transparent dots
      color = "white", 
      size = 5,
    ) +
    geom_point( # Individual trend plotting
      aes(color = Estimate), 
      size = 5,
      alpha = 0.6
    ) +
    geom_text( # Asterisks for significance
      data = filter(longData[[datasetType]], Estimate == "TrendA", Measurement == trend),
      aes(
        label = Significance,
        x = Residue, # Offset vertically for only trend A (proportional to number of data points)
        y = Coefficient,
      ),
      size = 6, 
      position = position_nudge(x = 0.15),
      color = '#4F9437',
      inherit.aes = FALSE
    ) +
    geom_text( # Asterisks for significance
      data = filter(longData[[datasetType]], Estimate == "TrendB", Measurement == trend),
      aes(
        label = Significance,
        x = Residue, # Offset vertically for only trend A (proportional to number of data points)
        y = Coefficient,
      ),
      size = 6, 
      position = position_nudge(x = 0.25),
      color = '#8100b6',
      inherit.aes = FALSE
    ) +
    theme(
      legend.position.inside = c(0.9, 0.15), 
      panel.grid.major.y = element_blank()
    ) +
    scale_x_discrete(labels = c(
      "T1" = bquote(θ[1]), 
      "T2" = bquote(θ[2]), 
      "T3" = bquote(θ[3]), 
      "TH1" = bquote(θ[H1]), 
      "THN" = bquote(θ[HN]), 
      "THH" = bquote(θ[H-H]), 
      "TT" = bquote(θ[T]), 
      "CuNterm" = bquote(Cu-N[term]),
      "CuAx" = bquote(Cu-H[2]*O[Ax]), 
      "CuEq" = bquote(Cu-H[2]*O[Eq]), 
      "CuHis1ND1" = bquote(Cu-His[1]*Nδ[1]), 
      "CuHis84NE2" = bquote(Cu-His[84]*Nε[2]),
      "CuTyr" = bquote(Cu-Tyr[168]),
      "Ax" = bquote(H[2]*O[Ax]), 
      "Eq" = bquote(H[2]*O[Eq]), 
      "CO2" = bquote(CO[2]), 
      "Glu" = bquote(Glu[30])
    )) +
    labs(
      y = if(trend %in% c("Occupancies", "BFactors")) {
        "Residue"
      } else if(trend == "Distances") {
        "Atom Pair"
      } else if(trend == "Angles") {
        "Angle ID"
      } else {stop("Invalid trend input for trend visualization")},
      x = if(trend == "Occupancies") {
        "Occupancy Trend (Δ/MGy)"
      } else if(trend == "BFactors") {
        bquote("B-Factor Trend (" * Å^2 * "/MGy)")
      } else if(trend == "Distances") {
        "Distance Trend (ΔÅ/MGy)"
      } else if(trend == "Angles") {
        "Angle Trend (Δ°/MGy)"
      } else {stop("Invalid trend input for trend visualization")}
    ) +
    coord_flip()
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
