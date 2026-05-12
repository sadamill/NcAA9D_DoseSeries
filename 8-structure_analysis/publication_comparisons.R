library(tidyverse)
library(readxl)

# Import comparisons and get rid of any proteins treated with ligand or Asc
# since these we're focusing on resting-state, substrate free enzyme
lit_comparison <- read_csv(
  "8-structure_analysis/input/lpmo_comparisons.csv",
  col_types = list(protein = 
                     readr::col_factor(levels = c(
                       "AoAA13", "EfAA10A", 
                       "TaAA9A", "LsAA9A(Ec)", 
                       "LsAA9A(f)", "NcAA9D"
                     )))
) |> 
  filter_out(soak %in% c("cell3", "cell4", "ascorbate")) |> 
  mutate(
    min_fwd = min(fwd, na.rm = TRUE),
    max_fwd = max(fwd, na.rm = TRUE),
    .by = protein
  ) |> 
  pivot_longer(
    ends_with("_distance"),
    values_to = "distance",
    names_to = c("location", "in_out", NA),
    names_pattern = "cu-(.*)_(.*)_(.*)",
    names_ptypes = list("location" = factor(), "in_out" = factor())
  )

thresholds <- tibble(
  redox = factor(c("Cu(II)", "Mixed", "Cu(I)")),
  min_eq_distance = c(0, 2.2, 2.9),
  max_eq_distance = c(2.2, 2.9, 5),
  min_ax_distance = c(0, 2.7, 3.2),
  max_ax_distance = c(2.7, 3.2, 5)
)

scientific_10 <- function(x) {
  parse(text = ifelse(
    x == 0, 0,
    gsub("e\\+*", " %*% 10^", scales::label_scientific()(x))
  ))
}

ggplots$Light$publication_comparisons$dose <- lit_comparison |> 
  filter_out(is.na(fwd)) |> 
  ggplot(aes(
    x = fwd,
    xmin = min_fwd, 
    xmax = max_fwd,
    y = protein,
    fill = protein
  )) +
  geom_rect(color = "black", height = 0.3) +
  geom_point(alpha = 0.5) +
  scale_fill_discrete(palette = "Set2") +
  scale_y_discrete(
    labels = c(
      "NcAA9D" = expression(italic("Nc") * "AA9D"),
      "LsAA9A(f)" = expression(italic("Ls") * "AA9A(" * italic("f") * ")"),
      "LsAA9A(Ec)" = expression(italic("Ls") * "AA9A(" * italic("Ec") * ")"),
      "TaAA9A" = expression(italic("Ta") * "AA9A"),
      "EfAA10A" = expression(italic("Ef") * "AA10A"),
      "AoAA13" = expression(italic("Ao") * "AA13")
    )
  ) +
  scale_x_continuous(
    labels = scientific_10
  ) +
  labs(
    x = "Average FWD (Gy)",
    y = "Protein"
  ) +
  ggtheme_light(legend.position = "none")
  
ggplots$Light$publication_comparisons$distance <- lit_comparison |> 
  filter(authors == "Miller et al.",
         !is.na(distance)) |> 
  ggplot(aes(x = distance, y = location, shape = in_out)) +
  geom_rect(
    data = thresholds, 
    aes(xmin = min_eq_distance, 
        xmax = max_eq_distance, 
        y = "eq", 
        fill = redox), 
    inherit.aes = FALSE, 
    height = 1,
    alpha = 0.25,
    color = "black"
  ) +
  geom_rect(
    data = thresholds, 
    aes(xmin = min_ax_distance, 
        xmax = max_ax_distance, 
        y = "ax", 
        fill = redox), 
    inherit.aes = FALSE, 
    height = 1,
    alpha = 0.25,
    color = "black"
  ) +
  geom_point(alpha = 0.5) +
  scale_fill_manual(
    "Cu Redox State",
    breaks = c("Cu(II)", "Mixed", "Cu(I)"),
    values = c("green", "yellow", "red")
  ) +
  scale_y_discrete(labels = c(
    "ax" = bquote(H[2] * O[ax]),
    "eq" = bquote(H[2] * O[eq])
  )) +
  scale_shape_discrete(
    'Water Position',
    breaks = c("in", "out"),
    labels = c("In", "Out")
  ) +
  coord_cartesian(
    xlim = c(1.9, 3.6),
    expand = FALSE
  ) +
  ggtheme_light(legend.position = "right") +
  labs(x = "Distance (Å)",
       y = "Water Location")
