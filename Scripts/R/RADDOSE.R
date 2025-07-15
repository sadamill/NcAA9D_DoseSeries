allData <- list()

allData$wedgeDoseState <- lapply(1:38, function(i) {
  read.csv(paste0('/Users/sm9/Documents/17_Meilleur_Lab/NcLPMO9D_Dose_Series_Study/Analysis/Input/RADDOSE/Output/Wedges/Wedge', i, '/output-DoseState.csv'), header = FALSE) %>% 
    rename(x = V1, y = V2, z = V3, dose = V4, fluence = V5, elasticScattering = V6) %>% 
    mutate(number = paste(!!i), datatype = 'Wedge')
})

allData$pseudohelixDoseState <- lapply(1:36, function(i) {
  read.csv(paste0('/Users/sm9/Documents/17_Meilleur_Lab/NcLPMO9D_Dose_Series_Study/Analysis/Input/RADDOSE/Output/Pseudohelices/Pseudohelix', i, '/output-DoseState.csv'), header = FALSE) %>% 
    rename(x = V1, y = V2, z = V3, dose = V4, fluence = V5, elasticScattering = V6) %>% 
    mutate(number = paste(!!i), datatype = 'Pseudohelix')
})

allData$wedgeDoseDWDs <- lapply(1:38, function(i) {
  read.csv(paste0('/Users/sm9/Documents/17_Meilleur_Lab/NcLPMO9D_Dose_Series_Study/Analysis/Input/RADDOSE/Output/Wedges/Wedge', i, '/output-DWDs.csv'), header = TRUE) %>% 
    mutate(number = paste(!!i), datatype = 'Wedge')
})

allData$pseudohelixDWDs <- lapply(1:36, function(i) {
  read.csv(paste0('/Users/sm9/Documents/17_Meilleur_Lab/NcLPMO9D_Dose_Series_Study/Analysis/Input/RADDOSE/Output/Pseudohelices/Pseudohelix', i, '/output-DWDs.csv'), header = TRUE) %>% 
    mutate(number = paste(!!i), datatype = 'Pseudohelix')
})

# Define Crystal Box ------------------------------------------------------

# Define corners of prism
xrange <- range(allData$wedgeDoseState[[1]]$x)
yrange <- range(allData$wedgeDoseState[[1]]$y)
zrange <- range(allData$wedgeDoseState[[1]]$z)

corners <- expand.grid(
  x = xrange,
  y = yrange,
  z = zrange
)

# Define list of indices for valid prism edges
edgeIndices <- rbind(
  c(1,2), c(1,3), c(1,5),
  c(8,7), c(8,6), c(8,4),
  c(2,4), c(2,6),
  c(5,6), c(5,7),
  c(3,4), c(3,7)
)

# Convert indices into stop/end coordinates for edges Use NA spacers so each segment is separated
edges <- lapply(1:nrow(edgeIndices), function(i) {
  p1 <- corners[edgeIndices[i,1],]
  p2 <- corners[edgeIndices[i,2],]
  list(
    x = c(p1$x, p2$x, NA),
    y = c(p1$y, p2$y, NA),
    z = c(p1$z, p2$z, NA)
  )
})

edges <- tibble(
  x = sapply(edges,'[[',1) %>% 
    sapply(function(column) {
      x <- c(column)
    }),
  y = sapply(edges,'[[',2) %>% 
    sapply(function(column) {
      y <- c(column)
    }),
  z = sapply(edges,'[[',3) %>% 
    sapply(function(column) {
      z <- c(column)
    }),
)

# Plot Generation ---------------------------------------------------------

steps <- list()
fig <- plot_ly()
for (i in 1:length(wedgeDose)) {
  fig <- fig %>% add_trace(data = wedgeDose[[i]],
                           type = 'isosurface',
                           x = ~x,  y = ~y, z = ~z,
                           value = ~dwd,
                           visible = (i == 1), 
                           name = "\u200B", 
                           isomin = 30, isomax = Inf, 
                           opacity = 1, 
                           coloraxis = 'coloraxis') %>% 
    add_trace(
      data = wedgeDose[[i]],
      type = 'isosurface',
      x = ~x, y = ~y, z = ~z,
      value = ~dwd,
      visible = (i == 1), 
      name = "\u200B",
      isomin = 5, isomax = Inf,
      opacity = 0.25,
      coloraxis = 'coloraxis'
    ) %>%
    add_trace(
      data = wedgeDose[[i]],
      type = 'isosurface',
      x = ~x, y = ~y, z = ~z,
      value = ~dwd,
      visible = (i == 1), 
      name = "\u200B", 
      isomin = 1, isomax = Inf,
      opacity = 0.05,
      coloraxis = 'coloraxis'
    )
  
  numberTraces <- 3
  totalTraces <- length(wedgeDose) * numberTraces + 1  # +1 for box trace
  index <- ((i - 1) * numberTraces + 1):((i - 1) * numberTraces + numberTraces)
  step <- list(args = list('visible', rep(FALSE, length(wedgeDose) * numberTraces)),
               method = 'restyle',
               label = paste('Wedge', i))
  step$args[[2]][index] <- TRUE
  step$args[[2]][totalTraces] <- TRUE
  steps[[i]] <- step
}

fig <- fig %>%
  layout(
    scene = list(
      annotations = list(
        list( ###############################DFIX THIS RIGHT NOW
          showarrow = TRUE,
          arrowside = 'none', # No arrowhead
          x = 0, y = max(wedgeDose[[1]]$y), z = 0,
          ax = 100, ay = 0,
          text = "Rotation Axis",
          font = list(color = "black", size = 14),
          arrowcolor = "black"
        )
      ),
      aspectmode = "data",
      xaxis = list(
        title = "X (mm)",
        showgrid = FALSE,
        zeroline = FALSE,
        nticks = 3,
        ticks = 'outside',
        showspikes = FALSE
      ),
      yaxis = list(
        title = "Y (mm)", 
        showgrid = FALSE, 
        zeroline = FALSE, 
        ticks = 'outside', 
        showspikes = FALSE
      ),
      zaxis = list(
        title = "Z (mm)",
        showgrid = FALSE, 
        zeroline = FALSE, 
        nticks = 3, 
        ticks = 'outside', 
        showspikes = FALSE
      ),
      camera = list(
        center = list(x = 0, y = 0, z = 0),
        eye = list(x = 2.25, y = 0.01, z = 1.1)
      )
    ),
    coloraxis = list(
      cmax = 30,
      cmin = 0,
      colorbar = list(
        tickmode = 'array',
        tickvals = list(0, 1, 5, 30),
        title = list(
          text = 'Density-Weighted Dose (MGy)',
          side = 'right'
        )
      ),
      colorscale = list(
        list(0/30, 'white'),
        list(1/30, 'gray'),
        list(5/30, 'dodgerblue'),
        list(30/30, 'red')
      )
    ),
    hovermode = FALSE,
    sliders = list(
      list(active = 0,
           currentvalue = list(prefix = "Dataset: "),
           steps = steps)
    )
  )

fig <- fig %>% add_trace(
  type = 'scatter3d',
  mode = 'lines',
  x = edges$x,
  y = edges$y,
  z = edges$z,
  showlegend = FALSE,
  color = I('black')
)

fig
