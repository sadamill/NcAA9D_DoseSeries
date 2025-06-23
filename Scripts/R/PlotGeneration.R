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
} #Global ggplot dark theme
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
} #Global ggplot light theme
faceting <- function(facetVar, datasetType) {
  mapping_list <- list(
    Residue = c(
      "CO2" = "CO[2]", 
      "Ax"  = "H[2]*O[Ax]", 
      "Eq"  = "H[2]*O[Eq]", 
      "Glu" = "Glu[30]", 
      "Oxy" = "'Dioxygen'"
    ),
    bResidue = c(
      "CO2BF" = "CO[2]", 
      "AxBF"  = "H[2]*O[Ax]", 
      "EqBF"  = "H[2]*O[Eq]", 
      "GluBF" = "Glu[30]", 
      "OxyBF" = "'Dioxygen'"
    ),
    CuAtomPair = c(
      "Cu-Tyr" = "Cu-Tyr[168]",
      "Cu-NTerm" = "Cu-N[term]",
      "Cu-His1ND1" = "Cu-His[1]*Nδ[1]",
      "Cu-His84NE2" = "Cu-His[84]*Nε[2]",
      "Cu-Eq" = "Cu-H[2]*O[Eq]",
      "Cu-Ax" = "Cu-H[2]*O[Ax]"
    ),
    OAtomPair = c(
      "O1-Eq" = "O[prox]*-H[2]*O[Eq]", 
      "O2-Eq" = "O[dist]*-H[2]*O[Eq]", 
      "O1-His157NE2" = "O[prox]*-His[157]*N[ε2]", 
      "O2-His157NE2" = "O[dist]*-His[157]*N[ε2]", 
      "O1-GluOE1" = "O[prox]*-Glu[30]*O[ε1]", 
      "O2-GluOE1" = "O[dist]*-Glu[30]*O[ε1]", 
      "O1-GluOE2" = "O[prox]*-Glu[30]*O[ε2]", 
      "O2-GluOE2" = "O[dist]*-Glu[30]*O[ε2]"
    ),
    AngleID = c(
      "T1" = "θ[1]", 
      "T2" = "θ[2]", 
      "T3" = "θ[3]", 
      "TT" = "θ[T]", 
      "THH" = "θ[HH]", 
      "TH1HN1" = "θH[1]*HN[1]", 
      "TH1HN84" = "θH[1]*HN[84]", 
      "TOxy" = "θ[oxygen]"
    )
  )
  
  list(
    facet_wrap(
      as.formula(paste("~", facetVar)), 
      scales = 'free',
      labeller = labeller(
        !!facetVar := as_labeller(mapping_list[[facetVar]], label_parsed),
        .default = label_parsed
      )),
    labs(
      x = if(datasetType == 'Pseudohelix') {
        "Dose (MGy)"
      } else if(datasetType == "Wedge") {
        "Wedge Number"
      },
      y = if(facetVar == "Residue") {
        "Occupancy"
      } else if(facetVar == "bResidue") {
        bquote("B-Factor (" * Å^2 * ")")
      } else if(facetVar %in% c('CuAtomPair', 'OAtomPair')) {
        "Distance (Å)"
      } else if(facetVar == "AngleID") {
        "Angle (°)"
      }
    )
  )
}
scales <- function(facetVar) {
  
  angle_window <- 7
  occupancy_window <- 0.4
  bfactor_window <- 8
  cudistance_window <- 0.25
  
  facetted_pos_scales(y = if(facetVar == "AngleID") {
    list(
      scale_y_continuous(limits = c(94 - angle_window / 2, 94 + angle_window / 2)),
      scale_y_continuous(limits = c(94 - angle_window / 2, 94 + angle_window / 2)),
      scale_y_continuous(limits = c(171 - angle_window / 2, 171 + angle_window / 2)),
      scale_y_continuous(limits = c(180 - angle_window / 2, 180 + angle_window / 2)),
      scale_y_continuous(limits = c(176 - angle_window / 2, 176 + angle_window / 2)),
      scale_y_continuous(limits = c(68 - angle_window / 2, 68 + angle_window / 2)),
      scale_y_continuous(limits = c(120, 150)),
      scale_y_continuous(limits = c(5 - angle_window / 2, 5 + angle_window / 2))
    )} else if(facetVar == "Residue") {
      list(
        scale_y_continuous(limits = c(0.6 - occupancy_window / 2, 0.6 + occupancy_window / 2)),
        scale_y_continuous(limits = c(0.55 - occupancy_window / 2, 0.55 + occupancy_window / 2)),
        scale_y_continuous(limits = c(0.35 - occupancy_window / 2, 0.35 + occupancy_window / 2)),
        scale_y_continuous(limits = c(0.55 - occupancy_window / 2, 0.55 + occupancy_window / 2)),
        scale_y_continuous(limits = c(0.65 - occupancy_window / 2, 0.65 + occupancy_window / 2))
      )} else if(facetVar == "bResidue") {
        list(
          scale_y_continuous(limits = c(21 - bfactor_window / 2, 21 + bfactor_window / 2)),
          scale_y_continuous(limits = c(10.5 - bfactor_window / 2, 10.5 + bfactor_window / 2)),
          scale_y_continuous(limits = c(22.6 - bfactor_window / 2, 22.6 + bfactor_window / 2)),
          scale_y_continuous(limits = c(14 - bfactor_window / 2, 14 + bfactor_window / 2)),
          scale_y_continuous(limits = c(12.5 - bfactor_window / 2, 12.5 + bfactor_window / 2))
        )
      } else if(facetVar == "CuAtomPair") {
        list(
          scale_y_continuous(limits = c(2.375 - cudistance_window / 2, 2.375 + cudistance_window / 2)),
          scale_y_continuous(limits = c(1.99 - cudistance_window / 2, 1.99 + cudistance_window / 2)),
          scale_y_continuous(limits = c(1.95 - cudistance_window / 2, 1.95 + cudistance_window / 2)),
          scale_y_continuous(limits = c(1.98 - cudistance_window / 2, 1.98 + cudistance_window / 2)),
          scale_y_continuous(limits = c(2.18 - cudistance_window / 2, 2.18 + cudistance_window / 2)),
          scale_y_continuous(limits = c(2.66 - cudistance_window / 2, 2.66 + cudistance_window / 2))
        )
      } else {NULL})
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
    if(facetVar == "AngleID" | facetVar == "OAtomPair") {
      theme(legend.position.inside = c(0.85, 0.15))
    } else if(facetVar == "CuAtomPair") {
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
    if(facetVar == "AngleID" | facetVar == "OAtomPair") {
      theme(legend.position.inside = c(0.85, 0.15))
    } else if(facetVar == "CuAtomPair") {
      theme(legend.position = "right")
    }
} #Make scatter plot with fitted linear regression

ggdarklight <- function(data, key, key2 = NA) {
  yvar <- if(deparse(substitute(data)) == "stackedOccupancies") {
    'Occupancy'
  } else if(deparse(substitute(data)) == "stackedBFactors") {
    'bFactor'
  } else if(deparse(substitute(data)) %in% c("stackedDistances$Cu", "stackedDistances$Oxy")) {
    'Distance'
  } else if(deparse(substitute(data)) == "stackedAngles") {
    'Angle'
  } else {stop("Invalid data input")}
  
  facetStr <- if(deparse(substitute(data)) == "stackedOccupancies") {
    "Residue"
  } else if(deparse(substitute(data)) == "stackedBFactors") {
    "bResidue"
  } else if(deparse(substitute(data)) == "stackedDistances$Cu") {
    "CuAtomPair"
  } else if(deparse(substitute(data)) == "stackedDistances$Oxy") {
    "OAtomPair"
  } else if(deparse(substitute(data)) == "stackedAngles") {
    "AngleID"
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
  
  if(is.na(key2)) {
    ggplots$Light[[key]] <<- lightplots
    ggplots$Dark[[key]] <<- darkplots
  } else {
    ggplots$Light[[key]][[key2]] <<- lightplots
    ggplots$Dark[[key]][[key2]] <<- darkplots
  }

} #Combination of scatter_dark and scatter_light to ease the plotting of multiple datasets with multiple themes

# Scatter plots -----------------------------------------------------------

ggdarklight(stackedOccupancies, "Occupancies")
ggdarklight(stackedBFactors, "BFactors")
ggdarklight(stackedDistances$Cu, "Distances", "Cu")
ggdarklight(stackedDistances$Oxy, "Distances", "Oxy")
ggdarklight(stackedAngles, "Angles")

# Trend plotting ----------------------------------------------------------

dark.trend <- function(trend, datasetType, nudge_y_val = 0.2) {
  
  regressionSummaries[[trend]][[datasetType]] <- regressionSummaries[[trend]][[datasetType]] %>%
    mutate(Residue = factor(Residue, levels = levels(wideData[[trend]][[datasetType]]$Residue)))
  
  ggplot(
    subset(regressionSummaries[[trend]][[datasetType]], Estimate %in% c("TrendA", "TrendB")), 
    aes(x = Residue, y = Coefficient)
  ) + 
    geom_hline(yintercept = 0, color = 'gray') +
    geom_segment(
      data = wideData[[trend]][[datasetType]], 
      aes(
        x = Residue, 
        xend = Residue, 
        y = TrendA, 
        yend = TrendB, 
        color = Significance,
        alpha = Significance
      ), 
      show.legend = FALSE, 
      inherit.aes = FALSE, 
      linewidth = 5
    ) +
    scale_color_manual(
      values = c("gray", "skyblue2"), 
    ) +
    scale_alpha_manual(
      values = c(0.5, 0.5)
    ) +
    new_scale_color() +
    ggtheme_dark() +
    geom_point(
      color = "black", 
      size = 5,
    ) +
    geom_point(
      aes(color = Estimate), 
      size = 5,
      alpha = 0.6
    ) +
    geom_text(
      aes(
        label = Significance,
        x = as.numeric(Residue) + ifelse(Estimate == "TrendA", nrow(regressionSummaries[[trend]][[datasetType]])/3/72, nrow(regressionSummaries[[trend]])/3/38),
        y = Coefficient
      ),
      size = 6, 
      color = 'white',
      inherit.aes = FALSE
    ) +
    theme(
      legend.position.inside = c(0.94, 0.08), 
      panel.grid.major.y = element_blank()
    ) +
    scale_x_discrete(labels = c(
      "T1" = bquote(θ[1]), 
      "T2" = bquote(θ[2]), 
      "T3" = bquote(θ[3]), 
      "TH1HN1" = bquote(θH[1]*HN[1]), 
      "TH1HN84" = bquote(θH[1]*HN[84]), 
      "THH" = bquote(θ[HH]), 
      "TT" = bquote(θ[T]), 
      "CuAx" = bquote(Cu-H[2]*O[Ax]), 
      "CuEq" = bquote(Cu-H[2]*O[Eq]), 
      "CuHis1ND1" = bquote(Cu-His[1]*Nδ[1]), 
      "CuHis84NE2" = bquote(Cu-His[84]*Nε[2]), 
      "CuO1" = bquote(Cu-O[prox]), 
      "CuO2" = bquote(Cu-O[dist]), 
      "CuTyr" = bquote(Cu-Tyr[168]), 
      "O1Eq" = bquote(O[prox]*H[2]*O[Eq]), 
      "O2Eq" = bquote(O[dist]*H[2]*O[Eq]), 
      "O1GluOE1" = bquote(O[prox]*Glu[30]*Oε[1]), 
      "O2GluOE1" = bquote(O[dist]*Glu[30]*Oε[1]), 
      "O1GluOE2" = bquote(O[prox]*Glu[30]*Oε[2]), 
      "O2GluOE2" = bquote(O[dist]*Glu[30]*Oε[2]), 
      "O1His157NE2" = bquote(O[prox]*His[157]*Nε[2]), 
      "O2His157NE2" = bquote(O[dist]*His[157]*Nε[2]), 
      "Ax" = bquote(H[2]*O[Ax]), 
      "Eq" = bquote(H[2]*O[Eq]), 
      "CO2" = bquote(CO[2]), 
      "Glu" = bquote(Glu[30]), 
      "OxyBF" = bquote(Dioxygen), 
      "AxBF" = bquote(H[2]*O[Ax]), 
      "EqBF" = bquote(H[2]*O[Eq]), 
      "CO2BF" = bquote(CO[2]), 
      "GluBF" = bquote(Glu[30]), 
      "OxyBF" = bquote(Dioxygen)
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
light.trend <- function(trend, datasetType, nudge_y_val = 0.2) {
  
  regressionSummaries[[trend]][[datasetType]] <- regressionSummaries[[trend]][[datasetType]] %>%
    mutate(Residue = factor(Residue, levels = levels(wideData[[trend]][[datasetType]]$Residue)))
  
  ggplot(
    subset(regressionSummaries[[trend]][[datasetType]], Estimate %in% c("TrendA", "TrendB")), 
    aes(x = Residue, y = Coefficient)
  ) + 
    geom_hline(yintercept = 0, color = 'gray') +
    geom_segment(
      data = wideData[[trend]][[datasetType]], 
      aes(
        x = Residue, 
        xend = Residue, 
        y = TrendA, 
        yend = TrendB, 
        color = Significance,
        alpha = Significance
      ), 
      show.legend = FALSE, 
      inherit.aes = FALSE, 
      linewidth = 5
    ) +
    scale_color_manual(
      values = c("gray", "skyblue2"), 
    ) +
    scale_alpha_manual(
      values = c(0.5, 0.5)
    ) +
    new_scale_color() +
    ggtheme_light() +
    geom_point(
      color = "white", 
      size = 5,
    ) +
    geom_point(
      aes(color = Estimate), 
      size = 5,
      alpha = 0.6
    ) +
    geom_text(
      aes(
        label = Significance,
        x = as.numeric(Residue) + ifelse(Estimate == "TrendA", nrow(regressionSummaries[[trend]][[datasetType]])/3/72, nrow(regressionSummaries[[trend]])/3/38),
        y = Coefficient
      ),
      size = 6, 
      color = 'black',
      inherit.aes = FALSE
    ) +
    theme(
      legend.position.inside = c(0.94, 0.08), 
      panel.grid.major.y = element_blank()
    ) +
    scale_x_discrete(labels = c(
      "T1" = bquote(θ[1]), 
      "T2" = bquote(θ[2]), 
      "T3" = bquote(θ[3]), 
      "TH1HN1" = bquote(θH[1]*HN[1]), 
      "TH1HN84" = bquote(θH[1]*HN[84]), 
      "THH" = bquote(θ[HH]), 
      "TT" = bquote(θ[T]), 
      "CuAx" = bquote(Cu-H[2]*O[Ax]), 
      "CuEq" = bquote(Cu-H[2]*O[Eq]), 
      "CuHis1ND1" = bquote(Cu-His[1]*Nδ[1]), 
      "CuHis84NE2" = bquote(Cu-His[84]*Nε[2]), 
      "CuO1" = bquote(Cu-O[prox]), 
      "CuO2" = bquote(Cu-O[dist]), 
      "CuTyr" = bquote(Cu-Tyr[168]), 
      "O1Eq" = bquote(O[prox]*H[2]*O[Eq]), 
      "O2Eq" = bquote(O[dist]*H[2]*O[Eq]), 
      "O1GluOE1" = bquote(O[prox]*Glu[30]*Oε[1]), 
      "O2GluOE1" = bquote(O[dist]*Glu[30]*Oε[1]), 
      "O1GluOE2" = bquote(O[prox]*Glu[30]*Oε[2]), 
      "O2GluOE2" = bquote(O[dist]*Glu[30]*Oε[2]), 
      "O1His157NE2" = bquote(O[prox]*His[157]*Nε[2]), 
      "O2His157NE2" = bquote(O[dist]*His[157]*Nε[2]), 
      "Ax" = bquote(H[2]*O[Ax]), 
      "Eq" = bquote(H[2]*O[Eq]), 
      "CO2" = bquote(CO[2]), 
      "Glu" = bquote(Glu[30]), 
      "OxyBF" = bquote(Dioxygen), 
      "AxBF" = bquote(H[2]*O[Ax]), 
      "EqBF" = bquote(H[2]*O[Eq]), 
      "CO2BF" = bquote(CO[2]), 
      "GluBF" = bquote(Glu[30]), 
      "OxyBF" = bquote(Dioxygen)
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
darklighttrend <- function(trend, nudge_y_val = 0.2) {
  
  ggplots$Light[[trend]]$Trends$Pseudohelices <<- light.trend(
    trend = trend,
    datasetType = "Pseudohelices",
    nudge_y_val = nudge_y_val
  )
  
  ggplots$Dark[[trend]]$Trends$Pseudohelices <<- dark.trend(
    trend = trend,
    datasetType = "Pseudohelices",
    nudge_y_val = nudge_y_val
  )
  
  ggplots$Light[[trend]]$Trends$Wedges <<- light.trend(
    trend = trend,
    datasetType = "Wedges",
    nudge_y_val = nudge_y_val
  )
  
  ggplots$Dark[[trend]]$Trends$Wedges <<- dark.trend(
    trend = trend,
    datasetType = "Wedges",
    nudge_y_val = nudge_y_val
  )
}

darklighttrend("Occupancies")
darklighttrend("BFactors")
darklighttrend("Distances")
darklighttrend("Angles")
