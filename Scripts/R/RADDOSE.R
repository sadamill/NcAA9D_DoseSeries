
# Data import -------------------------------------------------------------

allData <- list() # Initalize a list to contain all the raw data

allData$wedgeDoseState <- lapply(1:38, function(i) {
  read.csv(paste0('Input/RADDOSE/Output/Wedges/Wedge', i, '/output-DoseState.csv'), header = FALSE) %>% 
    rename(x = V1, y = V2, z = V3, dose = V4, fluence = V5, elasticScattering = V6) %>% 
    mutate(number = paste(!!i), datatype = 'Wedge')
}) # Import all voxel dose state files for the wedges (to be used for dose state isosurfaces)

allData$pseudohelixDoseState <- lapply(1:36, function(i) {
  read.csv(paste0('Input/RADDOSE/Output/Pseudohelices/Pseudohelix', i, '/output-DoseState.csv'), header = FALSE) %>% 
    rename(x = V1, y = V2, z = V3, dose = V4, fluence = V5, elasticScattering = V6) %>% 
    mutate(number = paste(!!i), datatype = 'Pseudohelix')
}) # Import all voxel dose state files for the pseudohelices (to be used for dose state isosurfaces)

allData$wedgeDWDs <- lapply(1:38, function(i) {
  read.csv(paste0('Input/RADDOSE/Output/Wedges/Wedge', i, '/output-DWDs.csv'), header = TRUE) %>% 
    mutate(number = paste(!!i), datatype = 'Wedge')
}) # Import all wedge DWD traces for the wedges (for wedge and pseudohelix DWD calculations)

# DWD Analysis ------------------------------------------------------------

dwds <- tibble(
  datasetNumber = rep(1:36)
) # Initialize a data frame to put all the doses in

dwds$wedges <- sapply(2:37, function(wedgeNumber) {
  wedge <- allData$wedgeDWDs[[wedgeNumber]]
  mean(wedge$DWD)
}) # Wedges 2-37 were used, so calculate the average DWD for all these


# Pseudohelix average DWD calculations calculate the average DWD for 5 frames
# worth of each wedge (subwedge). Calculate the respective subwedge average DWDs
# and average these across wedges 2-37 to calculate the average DWD for a pseudhelix
dwds$pseudohelices <- sapply(2:37, function(pseudohelixNumber) {
  subwedgeAverages <- sapply(2:37, function(subwedgeNumber) {
    wedge <- allData$wedgeDWDs[[subwedgeNumber]]
    startAngle <- ((pseudohelixNumber - 1) * 5) + ((subwedgeNumber - 2) * 5) # Start angle depends on both the pseudohelix number and subwedge you are calculating
    subwedgeAverage <- filter(wedge, DWD.Angle >= startAngle & DWD.Angle <= startAngle + 4)$DWD %>% # Extract the DWD column for angles within a target subwedge (5 frames = 4° angular range) 
      mean() # Average the extracted DWDs
    return(subwedgeAverage)
  }) # Creates a vector of subwedge average DWDs
  pseudohelixAverage <- mean(subwedgeAverages) # Average the subwedge average DWDs to generate a single pseudohelix average DWD
  return(pseudohelixAverage)
})

# Reshape the dwd tibble so the wedge and pseudohelix dwds are stacked and labeled. This allows for proper ggplotting
dwds <- tibble(
  datasetNumber = rep(1:36, 2),
  dwd_MGy = c(dwds$wedges, dwds$pseudohelices),
  datasetType = c(rep("Wedges", 36), rep("Pseudohelices", 36))
)

# Dose state analysis -----------------------------------------------------
## Define Crystal Box ------------------------------------------------------

# Define corners of prism
xrange <- c(-125, 125)
yrange <- c(-593.5, 593.5)
zrange <- c(-125, 125)

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

## Set up plotly functions -------------------------------------------------

# Initialize empty list for use in plotly sliders
steps <- list()

# Create an empty plotly object
base_fig <- plot_ly()

add_dummy_points <- function(plot, linecolor) {
  plot <- plot %>% add_trace( # Add dummy scatter3d traces to act as placeholders for the legend
    p = plot,
    type = 'scatter3d',
    mode = 'markers',
    x = c(max(xrange) * 10), y = c(max(yrange) * 10), z = c(max(zrange) * 10), # Have to set points out of bounds so they are visible. Positions are based on the data so the aspect ratio stays correct.
    marker = list(
      color = 'gray',
      size = 20,
      symbol = "square",
      line = list(
        color = linecolor,
        width = 2
      )
    ),
    name = "1 MGy"
  ) %>%
    add_trace(
      type = 'scatter3d',
      mode = 'markers',
      x = c(max(xrange) * 10), y = c(max(yrange) * 10), z = c(max(zrange) * 10),
      marker = list(
        color = 'dodgerblue',
        size = 20,
        symbol = "square",
        line = list(
          color = linecolor,
          width = 2
        )
      )
      ,
      name = "5 MGy"
    ) %>%
    add_trace(
      type = 'scatter3d',
      mode = 'markers',
      x = c(max(xrange) * 10), y = c(max(yrange) * 10), z = c(max(zrange) * 10),
      marker = list(
        color = 'red',
        size = 20,
        symbol = "square",
        line = list(
          color = linecolor,
          width = 2
        )
      ),
      name = "20 MGy"
    ) %>% 
    add_trace(
      type = 'scatter3d',
      mode = 'lines',
      x = edges$x,
      y = edges$y,
      z = edges$z,
      showlegend = FALSE,
      line = list(
        color = linecolor
      )
    )
  
  return(plot)
}
add_isosurface_traces <- function(plot, dataset, dataset_type) {
  
  # Add necessary traces for all the wedges (isosurfaces for 1, 5, and 20 MGy)
  for (i in 1:length(dataset)) {
    plot <- plot %>% 
      add_trace(
        data = dataset[[i]],
        type = 'isosurface',
        x = ~x, y = ~y, z = ~z,
        value = ~dose,
        visible = (i == length(dataset)), # Only make this plot visible initially
        name = "\u200B", # Zero-width space to prevent the name showing on legend
        isomin = 20, isomax = Inf,
        opacity = 1,
        showscale = FALSE,
        colorscale = list(c(0, 'red'), c(1, 'red'))
      ) %>% 
      add_trace(
        data = dataset[[i]],
        type = 'isosurface',
        x = ~x, y = ~y, z = ~z,
        value = ~dose,
        visible = (i == length(dataset)), 
        name = "\u200B",
        isomin = 5, isomax = Inf,
        opacity = 0.25,
        showscale = FALSE,
        colorscale = list(c(0, 'dodgerblue'), c(1, 'dodgerblue'))
      ) %>%
      add_trace(
        data = dataset[[i]],
        type = 'isosurface',
        x = ~x, y = ~y, z = ~z,
        value = ~dose,
        visible = (i == length(dataset)), 
        name = "\u200B", 
        isomin = 1, isomax = Inf,
        opacity = 0.1,
        showscale = FALSE,
        colorscale = list(c(0, 'gray'), c(1, 'gray'))
      )
    
    numberBoxTraces <- 1 # Traces used to draw the box outline
    numberDummyTraces <- 3 # Dummy scatter3d traces used to make a discrete color legend
    isosPerWedge <- 3 # Number of isosurface traces per wedge
    totalTraces <- (length(dataset) * isosPerWedge) + numberBoxTraces + numberDummyTraces
    index <- ((i - 1) * isosPerWedge + 1):((i - 1) * isosPerWedge + isosPerWedge) # Generate index of traces of the current wedge being plotted
    step <- list(args = list('visible', rep(FALSE, totalTraces)), # Initialize a list of steps the length of the total number of traces
                 method = 'restyle',
                 label = paste(i))
    step$args[[2]][index + (numberDummyTraces + numberBoxTraces)] <- TRUE # Set traces for the current wedge iteration to be visible
    step$args[[2]][(1:(numberDummyTraces + numberBoxTraces))] <- TRUE # The static traces (for legend and outline box) are always visible
    steps[[i]] <<- step # Inject the current step parameters into the steps list
  }
  
  return(plot)
}
plotly_layout <- function(plot, dataset, dataset_type, linecolor, bgcolor) {
  layout(
    plot, 
    font = list(
      color = linecolor
    ),
    paper_bgcolor = bgcolor,
    scene = list(
      aspectmode = "data",
      bgcolor = bgcolor,
      xaxis = list(
        title = "X (mm)",
        showgrid = FALSE,
        zeroline = FALSE,
        nticks = 3,
        ticks = 'outside', 
        tickcolor = linecolor,
        showspikes = FALSE,
        autorange = FALSE, 
        range = xrange
      ),
      yaxis = list(
        title = "Y (mm)", 
        showgrid = FALSE, 
        zeroline = FALSE, 
        ticks = 'outside',  
        tickcolor = linecolor,
        showspikes = FALSE,
        autorange = FALSE, 
        range = yrange
      ),
      zaxis = list(
        title = "Z (mm)",
        showgrid = FALSE, 
        zeroline = FALSE, 
        nticks = 3, 
        ticks = 'outside', 
        tickcolor = linecolor,
        showspikes = FALSE,
        autorange = FALSE, 
        range = zrange
      ),
      camera = list(
        center = list(x = 0, y = 0, z = 0),
        eye = list(x = 2.25, y = 0, z = 0.75),
        projection = list(type = "perspective")
      )
    ),
    margin = list(
      b = 0,
      t = 0,
      l = 0,
      r = 0
    ),
    annotations = list(
      list(
        text = "Isosurface Contour Level",
        x = 0.5,
        xanchor = "center",
        y = 0.95,
        yanchor = "bottom",
        xref = "paper",
        yref = "paper",
        showarrow = FALSE,
        font = list(size = 14)
      )
    ),
    legend = list(
      orientation = "h",
      x = 0.5,
      xanchor = "center",
      y = 0.95,
      yanchor = "top",
      bgcolor = "rgba(0,0,0,0)",
      itemclick = FALSE,
      itemdoubleclick = FALSE
    ),
    hovermode = FALSE,
    sliders = list(
      list(
        active = length(dataset) - 1,
        currentvalue = list(prefix = paste0("Dataset: ", dataset_type, " ")),
        steps = steps,
        tickcolor = linecolor,
        pad = list(
          b = 20,
          t = 0,
          l = 20,
          r = 20
        )
      )
    )
  )
}

plotly_dose <- function(dataset, theme) {
  bgcolor <- if(theme == "light") {
    "white"
  } else if(theme == "dark") {
    "black"
  } else {stop("Improper theme input")}
  
  linecolor <- if(theme == "light") {
    "black"
  } else if(theme == "dark") {
    "white"
  } else {stop("Improper theme input")}
  
  dataset_type = if(grepl("wedge", deparse(substitute(dataset)), fixed = TRUE)) {
      "Wedge"
    } else if(grepl("pseudohelix", deparse(substitute(dataset)), fixed = TRUE)) {
      "Pseudohelix"
    } else {stop(print(deparse(substitute(dataset))))}
  
  fig <- base_fig %>% 
    add_dummy_points(linecolor = linecolor) %>% 
    add_isosurface_traces(dataset = dataset, dataset_type = dataset_type) %>% 
    plotly_layout(dataset = dataset, dataset_type = dataset_type, linecolor = linecolor, bgcolor = bgcolor)
  
  steps <<- list()  #  Clear the global steps list to prepare for next plot
  
  return(fig)
}

## Plotly generation -------------------------------------------------------

plotlys <- list()

plotlys$light$wedgeDoseState <- plotly_dose(dataset = allData$wedgeDoseState, theme = "light")
plotlys$light$pseudohelixDoseState <- plotly_dose(dataset = allData$pseudohelixDoseState, theme = "light")
plotlys$dark$wedgeDoseState <- plotly_dose(dataset = allData$wedgeDoseState, theme = "dark")
plotlys$dark$pseudohelixDoseState <- plotly_dose(dataset = allData$pseudohelixDoseState, theme = "dark")
