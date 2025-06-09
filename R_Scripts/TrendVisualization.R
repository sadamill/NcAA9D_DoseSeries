library(ggnewscale)

###############################################################
#               LINEAR REGRESSION VISUALIZATION               #
###############################################################

make.wide <- function(trend,datasetType) {
  regressionSummaries[[trend]][[datasetType]] %>%
    filter(Estimate %in% c('TrendA', 'TrendB')) %>%
    pivot_wider(
      id_cols = Residue,
      names_from = Estimate,
      values_from = Coefficient
    ) %>%
    left_join(
      regressionSummaries[[trend]][[datasetType]] %>%
        filter(Estimate == 'Contrast') %>%
        select(Residue, PValue),
      by = 'Residue'
    ) %>%
    mutate(
      Significance = ifelse(PValue <= 0.05, '*', ' '),
      PValue = round(PValue, 3),
      MinTrend = pmin(TrendA, TrendB)
    ) %>%
    arrange(MinTrend) %>%
    mutate(Residue = factor(Residue, levels = Residue))
}

wideData <- list(
  Occupancies = list(
    Pseudohelices = make.wide('Occupancies', 'Pseudohelices'),
    Wedges = make.wide('Occupancies', 'Wedges')
  ),
  BFactors = list(
    Pseudohelices = make.wide('BFactors', 'Pseudohelices'),
    Wedges = make.wide('BFactors', 'Wedges')
  ),
  Distances = list(
    Pseudohelices = make.wide('Distances', 'Pseudohelices'),
    Wedges = make.wide('Distances', 'Wedges')
  ),
  Angles = list(
    Pseudohelices = make.wide('Angles', 'Pseudohelices'),
    Wedges = make.wide('Angles', 'Wedges')
  )
)

#Define function to plot different trends
plot.trend <- function(trend, datasetType, nudge_y_val = 0.2) {
  
  regressionSummaries[[trend]][[datasetType]] <- regressionSummaries[[trend]][[datasetType]] %>%
    mutate(Residue = factor(Residue, levels = levels(wideData[[trend]][[datasetType]]$Residue)))
  
  ggplot(
    subset(regressionSummaries[[trend]][[datasetType]], Estimate %in% c('TrendA', 'TrendB')),
    aes(y = Residue, x = Coefficient)
  ) +   
    geom_segment(
      data = wideData[[trend]][[datasetType]],
      aes(
        y = Residue,
        yend = Residue,
        x = TrendA,
        xend = TrendB,
        color = Significance
      ),
      show.legend = FALSE,
      inherit.aes = FALSE,
      linewidth = 5,
      alpha = 0.5
    ) +
    scale_color_manual(
      values = c('gray','cadetblue3'),
    ) +
    geom_vline(xintercept = 0) +
    new_scale_color() +
    ggtheme() +
    geom_point(
      aes(color = Estimate),
      size = 5
    ) +
    geom_text(
      aes(label = Significance),
      size = 6,
      position = position_nudge(y = nudge_y_val)
    ) +
    coord_cartesian(expand = TRUE)+ 
    theme(
      legend.position.inside = c(0.92,0.11),
      panel.grid.major.y = element_blank()
    ) +
    scale_y_discrete(labels = c(
      'T1' = bquote(θ[1]),
      'T2' = bquote(θ[2]),
      'T3' = bquote(θ[3]),
      'TH1HN1' = bquote(θH[1]*HN[1]),
      'TH1HN84' = bquote(θH[1]*HN[84]),
      'THH' = bquote(θ[HH]),
      'TT' = bquote(θ[T]),
      'CuAx' = bquote(Cu-H[2]*O[Ax]),
      'CuEq' = bquote(Cu-H[2]*O[Eq]),
      'CuHis1ND1' = bquote(Cu-His[1]*Nδ[1]),
      'CuHis84NE2' = bquote(Cu-His[84]*Nε[2]),
      'CuO1' = bquote(Cu-O[prox]),
      'CuO2' = bquote(Cu-O[dist]),
      'CuTyr' = bquote(Cu-Tyr[168]),
      'O1Eq' = bquote(O[prox]*H[2]*O[Eq]),
      'O2Eq' = bquote(O[dist]*H[2]*O[Eq]),
      'O1GluOE1' = bquote(O[prox]*Glu[30]*Oε[1]),
      'O2GluOE1' = bquote(O[dist]*Glu[30]*Oε[1]),
      'O1GluOE2' = bquote(O[prox]*Glu[30]*Oε[2]),
      'O2GluOE2' = bquote(O[dist]*Glu[30]*Oε[2]),
      'O1His157NE2' = bquote(O[prox]*His[157]*Nε[2]),
      'O2His157NE2' = bquote(O[dist]*His[157]*Nε[2]),
      'Ax' = bquote(H[2]*O[Ax]),
      'Eq' = bquote(H[2]*O[Eq]),
      'CO2' = bquote(CO[2]),
      'Glu' = bquote(Glu[30]),
      'OxyBF' = bquote(Dioxygen),
      'AxBF' = bquote(H[2]*O[Ax]),
      'EqBF' = bquote(H[2]*O[Eq]),
      'CO2BF' = bquote(CO[2]),
      'GluBF' = bquote(Glu[30]),
      'OxyBF' = bquote(Dioxygen)
    ))
}

ggplots$Occupancies$PseudohelixTrends <- plot.trend(
  trend = 'Occupancies',
  datasetType = 'Pseudohelices',
  nudge_y_val = 0.2
) +
  labs(x = 'Occupancy Trend (Δ/MGy)')
ggplots$Occupancies$WedgeTrends <- plot.trend(
  trend = 'Occupancies',
  datasetType = 'Wedges',
  nudge_y_val = 0.2
) +
  labs(x = 'Occupancy Trend (Δ/Wedge Number)')

ggplots$BFactors$PseudohelixTrends <- plot.trend(
  trend = 'BFactors',
  datasetType = 'Pseudohelices',
  nudge_y_val = 0.2
) +
  labs(x = bquote('B-Factor Trend ('*ΔÅ^2*'/MGy)'))
ggplots$BFactors$WedgeTrends <- plot.trend(
  trend = 'BFactors',
  datasetType = 'Wedges',
  nudge_y_val = 0.2
) +
  labs(x = bquote('B-Factor Trend ('*ΔÅ^2*'/Wedge Number)'))

ggplots$Distances$PseudohelixTrends <- plot.trend(
  trend = 'Distances',
  datasetType = 'Pseudohelices',
  nudge_y_val = 0.3
) +
  labs(x = 'Distance Trend (ΔÅ/MGy)')
ggplots$Distances$WedgeTrends <- plot.trend(
  trend = 'Distances',
  datasetType = 'Wedges',
  nudge_y_val = 0.3
) +
  labs(x = 'Distance Trend (ΔÅ/Wedge Number)')

ggplots$Angles$PseudohelixTrends <- plot.trend(
  trend = 'Angles',
  datasetType = 'Pseudohelices',
  nudge_y_val = 0.2
) +
  labs(x = 'Angle Trend (Δ°/MGy)')
ggplots$Angles$WedgeTrends <- plot.trend(
  trend = 'Angles',
  datasetType = 'Wedges',
  nudge_y_val = 0.2
) +
  labs(x = 'Angle Trend (Δ°/Wedge Number)')
