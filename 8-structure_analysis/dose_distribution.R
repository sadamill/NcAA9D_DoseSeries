# Data Import -------------------------------------------------------------

# Read in the start angles for each pseudohelix number, as calculated in
# section four.
pseudohelix_start_angles <- readr::read_csv("4-sampling/output/samples.csv") |> 
  dplyr::pull(start_angle)

pseudohelix_ddwd <- readr::read_csv("4-sampling/output/samples.csv") |> 
  dplyr::select(dataset_number, ddwd) |> 
  dplyr::mutate(pseudohelix_number = 1:36)
  
# Function to determine which RADDOSE frames need to be pulled for a given
# pseudohelix number, as dictated by the pseudohelix_start_angles object
# and a subwedge range of 5°.
determine_pseudohelix_frames <- function(pseudohelix_number) {
  
  wedge_info <- tibble::tibble(
    wedge = 2:37,
    start_angle = 1:36 * 5,
    end_angle = 1:36 * 5 + 180
  ) |> 
    dplyr::rowwise() |> 
    dplyr::mutate(frame_angles = list(
      seq(start_angle, 
          end_angle, 
          180 / 200)
    )) |> 
    tidyr::unnest_longer(
      frame_angles, 
      values_to = "frame_angle", 
      indices_to = "frame_number"
    )
  
  pseudohelix_start_angle <- pseudohelix_start_angles[pseudohelix_number]

  pseudohelix_subwedge_info <- tibble::tibble(
    subwedge = 1:36
  ) |> 
    dplyr::mutate(
      subwedge_start = pseudohelix_start_angle + (subwedge - 1) * 5,
      subwedge_end = subwedge_start + 5
    )
  
  purrr::map(1:36, \(subwedge) {
    subwedge_start <- pseudohelix_subwedge_info[subwedge, ]$subwedge_start
    subwedge_end <- pseudohelix_subwedge_info[subwedge, ]$subwedge_end
    
    wedge_info |> 
      dplyr::filter(
        wedge == subwedge + 1,
        frame_angle >= subwedge_start,
        frame_angle <= subwedge_end
      ) |> 
      dplyr::select(!start_angle:end_angle)
  }) |> 
    purrr::list_rbind()
  
}

# Function to read in voxel- and frame-wise data as output by RADDOSE-3D
# and transform it into proper Tidy format.
read_voxel_data <- function(filepath, measurement) {
  vroom::vroom(
    filepath,
    col_names = FALSE
  ) |> 
    dplyr::select(!where(anyNA)) |> 
    t() |> 
    tibble::as_tibble() |> 
    rlang::set_names(c("coordinates", stringr::str_glue("frame_{1:201}"))) |> 
    dplyr::mutate(coordinates = purrr::map(
      stringr::str_split(coordinates, "_"),
      \(coord) {
        c(x_coord = as.numeric(coord[1]),
          y_coord = as.numeric(coord[2]),
          z_coord = as.numeric(coord[3]))
      }
    )) |> 
    tidyr::unnest_wider(coordinates) |>
    tidyr::pivot_longer(
      starts_with("frame_"),
      names_to = "frame_number",
      names_prefix = "frame_",
      names_transform = as.numeric,
      values_to = as.character(stringr::str_glue("voxel_{measurement}"))
    ) |> 
    dplyr::mutate(across(everything(), as.numeric))
}

# Function to generate voxel- and frame-wise dose and fluence data from multiple
# subwedges. Uses the functions described above to do this.
generate_pseudohelix_voxel_data <- function(pseudohelix_number) {
  
  pseudohelix_frame_info <- determine_pseudohelix_frames(pseudohelix_number)
  
  purrr::map(1:36, \(subwedge_number) {
      wedge_number <- subwedge_number + 1
      this_subwedge <- pseudohelix_frame_info |> 
        dplyr::filter(wedge == wedge_number)
      frames <- this_subwedge$frame_number
      
      fst::read_fst(
        stringr::str_glue("8-structure_analysis/input/raddose_wedge_voxel_data/wedge_{wedge_number}_voxel_data.fst")
      ) |> 
        tibble::as_tibble() |> 
        dplyr::filter(frame_number %in% frames) |> 
        dplyr::left_join(
          dplyr::select(this_subwedge, frame_number, frame_angle),
          by = dplyr::join_by(frame_number)
        ) |> 
        dplyr::relocate(voxel_dose, voxel_fluence, .after = everything())
    }
  ) |> 
    purrr::list_rbind()
}

# Wedge voxel- and frame-wise data was saved in fst format. These
# files read drastically faster than CSV files, and are thus better for 
# the computations executed here, as the data is too large to store in RAM
# as Tibbles. Data.table also fails to read the wedge data, presumably because
# of its wide format.

if (all(!file.exists(stringr::str_glue("8-structure_analysis/input/raddose_wedge_voxel_data/wedge_{1:38}_voxel_data.fst")))) {
  for (i in 1:38) {
    dose_data <- read_voxel_data(
      filepath = stringr::str_glue("3-ddwd_calculation/output/raddose_output/wedge{i}/output-VoxelDose.csv"),
      measurement = "dose"
    )
    fluence_data <- read_voxel_data(
      filepath = stringr::str_glue("3-ddwd_calculation/output/raddose_output/wedge{i}/output-VoxelFluences.csv"),
      measurement = "fluence"
    )
    all_data <- dplyr::inner_join(
      dose_data, 
      fluence_data,
      by = dplyr::join_by(x_coord, y_coord, z_coord, frame_number)
    )
    fst::write_fst(all_data, stringr::str_glue("8-structure_analysis/input/raddose_wedge_voxel_data/wedge_{i}_voxel_data.fst"))
  }
  rm(dose_data, fluence_data, all_data)
}

# RDE Calculations ---------------------------------------------------------

# Define parameters for the Leal IDM, as described by Dickerson, et al., 2024
# in equation 4.
gamma <- 0.0257
b_naught <- 7.2160
beta <- 0.0307

# Import the BEST data, as this is required to calculate resolution-dependent
# intensity decay. Data was obtained from the source code of
# RADDOSE-3D (v5.01066).
best_data <- read.csv("8-structure_analysis/input/best_data.csv")

# Calculate midpoint values between BEST data points.
# The result is a tibble with three columns:
# h_squared contains the interpolated h^2 values from the BEST data
# i contains the interpolated intensity values from the BEST data
# dh contains the width of each h interval (NOTE: not h^2)
interpolated_values <- purrr::map(
  seq(nrow(best_data) - 1),
  \(idx) {
    tibble::tibble(
      h_squared = (best_data[[idx, 1]] + best_data[[idx + 1, 1]]) / 2,
      i = (best_data[[idx, 2]] + best_data[[idx + 1, 2]]) / 2,
      dh = sqrt(best_data[[idx + 1, 1]]) - sqrt(best_data[[idx, 1]])
    )
  }
) |> 
  purrr::list_rbind()

# For each interpolated value, calculate the log product and exponential product
# of the RADDOSE-implemented IDM (Dickerson et al., 2024, equation 4.
# Taking the natural log of the numerator of this equation, brings the
# exponentiated term down, allowing for faster calculation of η once the dose
# is known.
integral_products <- purrr::pmap(
  interpolated_values, \(h_squared, i, dh) {
    tibble::tibble(
      log_product = log(h_squared * i * dh * exp(-0.5 * b_naught * h_squared)),
      exponential_product = -0.5 * beta * h_squared
    )
  }) |>
  purrr::list_c()

# Function to calculate the integrated intensity given an input dose.
# If performs the integration across resolution shells.
get_integrated_intensity <- function(dose) {
  integral_sum <- purrr::reduce2(
    integral_products$log_product, 
    integral_products$exponential_product, 
    \(accumulated, log, expt) {
      accumulated + exp(expt * dose + log)
    }, 
    .init = 0
  )
  
  exp(-gamma^2 * dose^2) * integral_sum
}

# Need to calculate the intensity at zero dose to compare other values against.
zero_dose_integrated_intensity <- get_integrated_intensity(dose = 0)

# Contribution-Weighted Dose Distributions --------------------------------

# Calculate the relative diffraction efficiency (RDE) for each voxel using the
# voxel- and image-wise dose map as a template.
# Similar to how RADDOSE calculates DDWD, calculate each voxel's contribution
# by multiplying its fluence by its RDE. Normalize each image so its
# contributions all sum to one.
pseudohelix_voxel_data <- purrr::map(
  1:36, 
    \(pseudohelix_number) {
      generate_pseudohelix_voxel_data(pseudohelix_number) |> 
        # get rid of images where voxel wasn't exposed; these mess up the 
        # averaging and slow down compute time
        dplyr::filter_out(voxel_fluence == 0) |>
        dplyr::mutate(
          voxel_rde = get_integrated_intensity(dose = voxel_dose) / 
            zero_dose_integrated_intensity,
          voxel_contribution = voxel_fluence * voxel_rde
        ) |> 
        dplyr::mutate(pseudohelix_number = pseudohelix_number)
    },
  .progress = TRUE
) |> 
  purrr::list_rbind() |> 
  dplyr::arrange(pseudohelix_number, voxel_dose) |> 
  dplyr::mutate(
    total_contribution = sum(voxel_contribution),
    relative_contribution = voxel_contribution / total_contribution,
    cumulative_contribution = purrr::accumulate(relative_contribution, `+`),
    .by = pseudohelix_number
  ) |> 
  dplyr::relocate(pseudohelix_number, .before = x_coord)

# Binned Statistics -------------------------------------------------------

# Define bins for distribution visualization.
n_pseudohelices <- 36
max_dose <- max(pseudohelix_voxel_data$voxel_dose) |> 
  ceiling()
bin_width <- 0.5
bins <- seq(0, max_dose, bin_width)

# Make a tibble containing intensity contributions summed across dose bins.
# Ensure all frames and bins are filled in with zero values so that plots are
# completely filled in.
dose_bin_summary <- pseudohelix_voxel_data |> 
  dplyr::mutate(bin = cut(voxel_dose,
                   breaks = bins,
                   right = FALSE,
                   include.lowest = TRUE)) |> 
  dplyr::group_by(pseudohelix_number, bin) |> 
  dplyr::summarize(
    contribution = sum(relative_contribution),
    .groups = "drop"
  ) |> 
  tidyr::complete(pseudohelix_number, bin, fill = list(contribution = 0)) |> 
  dplyr::bind_cols(
    dose_min = rep(bins[-length(bins)], n_pseudohelices),
    dose_max = rep(bins[-1], n_pseudohelices),
    dose_middle = rep(bins[-1] - bin_width/2, n_pseudohelices)
  ) |> 
  dplyr::mutate(
    offset = (pseudohelix_number - 1) * 0.25, 
    total_offset = offset + contribution
  ) |> 
  dplyr::left_join(pseudohelix_ddwd)

# Plotting ----------------------------------------------------------------

# Make a stacked histogram of weighted dose distributions.
ggplots$Light$dose_distribution$dose_distribution <- dose_bin_summary |>
  ggplot2::ggplot(ggplot2::aes(
    x = dose_middle, 
    ymin = offset, 
    ymax = total_offset, 
    fill = pseudohelix_number, 
    color = pseudohelix_number, 
    group = pseudohelix_number
  )) +
  ggplot2::geom_ribbon(alpha = 0.5, linewidth = 0.2) +
  ggplot2::scale_y_continuous(breaks = NULL) +
  ggplot2::coord_cartesian(expand = FALSE, xlim = c(0, 30)) +
  ggplot2::labs(x = "Dose (MGy)", y = "Weighted Contribution") +
  ggtheme_light(
    legend.justification = c(1, 0),
    legend.position.inside = c(0.95, 0.05),
    legend.key.height = grid::unit(4, "mm")
  ) +
  ggplot2::scale_fill_gradient(
    "Pseudohelix\nNumber",
    breaks = c(1, 12, 24, 36),
    low = "orchid1",
    high = "darkorchid4",
    aesthetics = c("fill", "color")
  )

# Make an animated histogram to add in supplementary information.
anim <- dose_bin_summary |> 
  ggplot2::ggplot(ggplot2::aes(
    x = dose_middle, 
    ymin = 0, 
    ymax = contribution, 
    color = pseudohelix_number, 
    fill = pseudohelix_number
  )) +
  ggplot2::geom_ribbon(alpha = 0.5) +
  ggplot2::geom_vline(
    ggplot2::aes(xintercept = ddwd, color = pseudohelix_number),
    stat = "unique", 
    inherit.aes = FALSE
  ) +
  ggplot2::geom_text(
    ggplot2::aes(label = stringr::str_glue("Pseudohelix {pseudohelix_number}
                                Average DDWD = {ddwd}"), 
        x = ddwd + 1,
        color = pseudohelix_number
    ),
    inherit.aes = FALSE,
    stat = "unique",
    y = 0.175, 
    hjust = 0
  ) +
  ggplot2::coord_cartesian(ylim = c(0, 0.2), xlim = c(0, 40), expand = FALSE) +
  ggplot2::labs(x = "Dose (MGy)", y = "Weighted Contribution") +
  ggplot2::scale_fill_gradient(
    breaks = c(1, 12, 24, 36),
    low = "orchid1",
    high = "darkorchid4",
    aesthetics = c("fill", "color")
  ) +
  ggtheme_light(legend.position = "none") +
  gganimate::transition_states(
    pseudohelix_number,
    transition_length = 5,
    state_length = 0
  ) +
  gganimate::ease_aes("sine-in-out")
gganimate::anim_save(
  "test.gif",
  anim,
  nframes = 100,
  height = 5,
  width = 5,
  units = "in",
  res = 300
)
