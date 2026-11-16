#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(scales)
  library(patchwork)
  library(ragg)
  library(svglite)
})


# ============================================================
# Accepted versus rejected reads and read-length distributions
#
# Panel A:
#   Accepted and rejected percentages among mapped
#   mitochondrial reads for each species.
#
# Panel B:
#   Accepted and rejected read-length distributions.
#
# Both panels use identical numeric species-row positions.
# ============================================================


# ============================================================
# Local canonical-data directory
# ============================================================

project_dir <- local({
  workspace_args <- commandArgs(trailingOnly = TRUE)
  if (length(workspace_args) != 1L) {
    stop("Usage: Rscript SCRIPT.R WORKSPACE_DIRECTORY")
  }
  normalizePath(
    workspace_args[[1]],
    winslash = "/",
    mustWork = TRUE
  )
})

canonical_data_dir <- file.path(
  project_dir,
  "canonical_data"
)


analysis_dir <- file.path(
  canonical_data_dir,
  "analyses",
  "accept_reject_length"
)

read_file <- file.path(
  analysis_dir,
  "outputs",
  "full_nas_accept_reject_read_lengths.tsv.gz"
)

decision_summary_file <- file.path(
  analysis_dir,
  "outputs",
  "full_nas_accept_reject_summary.tsv"
)

length_summary_file <- file.path(
  analysis_dir,
  "outputs",
  "full_nas_accept_reject_length_summary.tsv"
)

output_dir <- file.path(
  analysis_dir,
  "outputs"
)

figure_dir <- file.path(
  analysis_dir,
  "figures"
)

log_dir <- file.path(
  analysis_dir,
  "logs"
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  figure_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  log_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ============================================================
# Output files
# ============================================================

figure_png <- file.path(
  figure_dir,
  "accept_reject_length_publication.png"
)

figure_pdf <- file.path(
  figure_dir,
  "accept_reject_length_publication.pdf"
)

figure_svg <- file.path(
  figure_dir,
  "accept_reject_length_publication.svg"
)

plot_summary_file <- file.path(
  output_dir,
  "full_nas_accept_reject_plot_summary_publication.tsv"
)

proportion_plot_data_file <- file.path(
  output_dir,
  "full_nas_accept_reject_proportion_plot_data.tsv"
)

length_plot_data_file <- file.path(
  output_dir,
  "full_nas_accept_reject_length_plot_data.tsv.gz"
)

display_summary_file <- file.path(
  output_dir,
  "full_nas_accept_reject_display_summary.tsv"
)

log_file <- file.path(
  log_dir,
  "plot_accept_reject_length_publication.log"
)

session_file <- file.path(
  log_dir,
  "R_sessionInfo_accept_reject_length_publication.txt"
)


# ============================================================
# Validate input files
# ============================================================

required_files <- c(
  read_file,
  decision_summary_file,
  length_summary_file
)

missing_files <- required_files[
  !file.exists(required_files)
]

if (length(missing_files) > 0) {
  stop(
    paste(
      "Missing required files:",
      paste(
        missing_files,
        collapse = ", "
      )
    )
  )
}


# ============================================================
# Read canonical inputs
# ============================================================

read_data <- fread(
  read_file,
  sep = "\t",
  header = TRUE,
  integer64 = "double"
) %>%
  as_tibble()

decision_summary <- fread(
  decision_summary_file,
  sep = "\t",
  header = TRUE,
  integer64 = "double"
) %>%
  as_tibble()

length_summary <- fread(
  length_summary_file,
  sep = "\t",
  header = TRUE,
  integer64 = "double"
) %>%
  as_tibble()


# ============================================================
# Required columns
# ============================================================

required_read_columns <- c(
  "dataset",
  "decision_class",
  "read_id",
  "reference",
  "species",
  "target_status",
  "read_length_bases",
  "cigar_aligned_query_bases",
  "aligned_fraction_of_read"
)

required_decision_columns <- c(
  "reference",
  "species",
  "target_status",
  "accepted_reads",
  "rejected_reads",
  "total_reads",
  "accepted_percent",
  "rejected_percent"
)

required_length_columns <- c(
  "reference",
  "species",
  "target_status",
  "decision_class",
  "n",
  "first_quartile_read_length_bases",
  "median_read_length_bases",
  "third_quartile_read_length_bases"
)

missing_read_columns <- setdiff(
  required_read_columns,
  names(read_data)
)

missing_decision_columns <- setdiff(
  required_decision_columns,
  names(decision_summary)
)

missing_length_columns <- setdiff(
  required_length_columns,
  names(length_summary)
)

if (length(missing_read_columns) > 0) {
  stop(
    paste(
      "Missing read-level columns:",
      paste(
        missing_read_columns,
        collapse = ", "
      )
    )
  )
}

if (length(missing_decision_columns) > 0) {
  stop(
    paste(
      "Missing decision-summary columns:",
      paste(
        missing_decision_columns,
        collapse = ", "
      )
    )
  )
}

if (length(missing_length_columns) > 0) {
  stop(
    paste(
      "Missing length-summary columns:",
      paste(
        missing_length_columns,
        collapse = ", "
      )
    )
  )
}


# ============================================================
# Validate canonical dimensions
# ============================================================

# Adapter-trimmed reanalysis counts
expected_counts <- c(
  "Accepted" = 76833,
  "Rejected" = 87054
)

observed_counts <- table(
  read_data$decision_class
)

for (decision in names(expected_counts)) {
  if (
    !decision %in% names(observed_counts) ||
    observed_counts[[decision]] != expected_counts[[decision]]
  ) {
    stop(
      paste(
        "Unexpected",
        decision,
        "record count"
      )
    )
  }
}

if (nrow(read_data) != 163887) {
  stop(
    paste(
      "Expected 163887 read-level rows, observed",
      nrow(read_data)
    )
  )
}

if (nrow(decision_summary) != 16) {
  stop(
    paste(
      "Expected 16 decision-summary rows, observed",
      nrow(decision_summary)
    )
  )
}

if (nrow(length_summary) != 32) {
  stop(
    paste(
      "Expected 32 length-summary rows, observed",
      nrow(length_summary)
    )
  )
}


# ============================================================
# Validate numerical values
# ============================================================

if (
  anyNA(read_data$read_length_bases) ||
  any(read_data$read_length_bases <= 0)
) {
  stop(
    "Missing or non-positive complete read lengths detected"
  )
}

if (
  anyNA(read_data$aligned_fraction_of_read) ||
  any(
    read_data$aligned_fraction_of_read <= 0 |
    read_data$aligned_fraction_of_read > 1
  )
) {
  stop(
    "Invalid aligned fractions detected"
  )
}

duplicate_keys <- duplicated(
  read_data %>%
    select(
      dataset,
      read_id
    )
)

if (any(duplicate_keys)) {
  stop(
    "Duplicate dataset and read_id records detected"
  )
}

percent_sums <- (
  decision_summary$accepted_percent +
    decision_summary$rejected_percent
)

if (
  !all(
    abs(
      percent_sums - 100
    ) < 1e-8
  )
) {
  stop(
    "Accepted and rejected percentages do not sum to 100"
  )
}


# ============================================================
# Species ordering
#
# Lowest accepted percentage appears at the bottom.
# Highest accepted percentage appears at the top.
# ============================================================

species_order <- decision_summary %>%
  arrange(
    accepted_percent,
    total_reads
  ) %>%
  pull(
    species
  ) %>%
  as.character()

if (length(unique(species_order)) != 16) {
  stop(
    "Expected 16 unique species"
  )
}

species_positions <- tibble(
  species = species_order,
  species_id = seq_along(
    species_order
  )
)


# ============================================================
# Abbreviated species display names
#
# Genus is abbreviated to its initial.
# Only genus and species epithet are displayed.
#
# Apis mellifera ligustica is displayed as A. mellifera.
# ============================================================

abbreviate_species <- function(species_name) {
  name_parts <- strsplit(
    species_name,
    split = " ",
    fixed = TRUE
  )[[1]]

  if (length(name_parts) < 2) {
    return(
      species_name
    )
  }

  paste0(
    substr(
      name_parts[1],
      1,
      1
    ),
    ". ",
    name_parts[2]
  )
}

species_abbreviations <- vapply(
  species_order,
  abbreviate_species,
  character(1)
)


# ============================================================
# Italic expression labels
# ============================================================

species_labels <- lapply(
  species_abbreviations,

  function(display_name) {
    bquote(
      italic(
        .(display_name)
      )
    )
  }
)

names(species_labels) <- as.character(
  species_positions$species_id
)

species_label_function <- function(x) {
  as.expression(
    unname(
      species_labels[
        as.character(
          x
        )
      ]
    )
  )
}


# ============================================================
# Shared decision-class levels
# ============================================================

decision_levels <- c(
  "Accepted",
  "Rejected"
)


# ============================================================
# Add shared numeric species positions
# ============================================================

decision_summary <- decision_summary %>%
  mutate(
    species = as.character(
      species
    )
  ) %>%
  left_join(
    species_positions,
    by = "species"
  )

read_data <- read_data %>%
  mutate(
    species = as.character(
      species
    ),

    decision_class = factor(
      decision_class,
      levels = decision_levels
    )
  ) %>%
  left_join(
    species_positions,
    by = "species"
  )

length_summary <- length_summary %>%
  mutate(
    species = as.character(
      species
    ),

    decision_class = factor(
      decision_class,
      levels = decision_levels
    )
  ) %>%
  left_join(
    species_positions,
    by = "species"
  )

if (
  anyNA(decision_summary$species_id) ||
  anyNA(read_data$species_id) ||
  anyNA(length_summary$species_id)
) {
  stop(
    "Missing species positions detected"
  )
}


# ============================================================
# Overall publication summary
# ============================================================

plot_summary <- read_data %>%
  group_by(
    decision_class
  ) %>%
  summarise(
    n = n(),

    minimum_read_length_bases =
      min(
        read_length_bases
      ),

    first_quartile_read_length_bases =
      as.numeric(
        quantile(
          read_length_bases,
          0.25
        )
      ),

    median_read_length_bases =
      median(
        read_length_bases
      ),

    mean_read_length_bases =
      mean(
        read_length_bases
      ),

    third_quartile_read_length_bases =
      as.numeric(
        quantile(
          read_length_bases,
          0.75
        )
      ),

    maximum_read_length_bases =
      max(
        read_length_bases
      ),

    standard_deviation_read_length_bases =
      sd(
        read_length_bases
      ),

    median_aligned_fraction_of_read =
      median(
        aligned_fraction_of_read
      ),

    .groups = "drop"
  )

fwrite(
  plot_summary,
  plot_summary_file,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)


# ============================================================
# Explicit stacked-bar rectangles
# ============================================================

accepted_rectangles <- decision_summary %>%
  transmute(
    species,
    species_id,

    decision_class = factor(
      "Accepted",
      levels = decision_levels
    ),

    xmin = 0,
    xmax = accepted_percent,
    ymin = species_id - 0.34,
    ymax = species_id + 0.34
  )

rejected_rectangles <- decision_summary %>%
  transmute(
    species,
    species_id,

    decision_class = factor(
      "Rejected",
      levels = decision_levels
    ),

    xmin = accepted_percent,
    xmax = 100,
    ymin = species_id - 0.34,
    ymax = species_id + 0.34
  )

proportion_data <- bind_rows(
  accepted_rectangles,
  rejected_rectangles
)

fwrite(
  proportion_data,
  proportion_plot_data_file,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)


# ============================================================
# Centred stacked annotations
#
# Each annotation is one centred text block:
#
#   97.6%
#   (n=288)
# ============================================================

count_labels <- decision_summary %>%
  mutate(
    annotation = paste0(
      sprintf(
        "%.1f%%",
        accepted_percent
      ),
      "\n(n=",
      format(
        total_reads,
        big.mark = ",",
        scientific = FALSE,
        trim = TRUE
      ),
      ")"
    )
  )


# ============================================================
# Fixed manuscript colours and common line thickness
# ============================================================

decision_colours <- c(
  "Accepted" = "#1B9E77",
  "Rejected" = "#7570B3"
)

common_line_width <- 0.30


# ============================================================
# Shared fill scale
# ============================================================

decision_fill_scale <- scale_fill_manual(
  values = decision_colours,
  breaks = decision_levels,
  labels = decision_levels,
  drop = FALSE
)


# ============================================================
# Shared legend guide
# ============================================================

decision_fill_guide <- guides(
  fill = guide_legend(
    nrow = 1,
    byrow = TRUE,
    title = NULL,

    override.aes = list(
      alpha = 1,
      colour = NA,
      linewidth = 0
    )
  )
)


# ============================================================
# Shared numeric species y-axis
# ============================================================

shared_y_limits <- c(
  0.50,
  length(species_order) + 0.50
)

shared_y_breaks <- species_positions$species_id


# ============================================================
# Shared publication theme
# ============================================================

theme_scientific_reports <- function() {
  theme_classic(
    base_size = 8,
    base_family = "sans"
  ) +
    theme(
      axis.line = element_line(
        colour = "black",
        linewidth = common_line_width
      ),

      axis.ticks = element_line(
        colour = "black",
        linewidth = common_line_width
      ),

      axis.ticks.length = unit(
        1.2,
        "mm"
      ),

      axis.text = element_text(
        colour = "black",
        size = 7.1
      ),

      axis.title = element_text(
        colour = "black",
        size = 7.75,
        face = "plain"
      ),

      axis.title.x = element_text(
        margin = margin(
          t = 3
        )
      ),

      panel.grid = element_blank(),

      panel.border = element_blank(),

      legend.position = "bottom",

      legend.direction = "horizontal",

      legend.title = element_blank(),

      legend.text = element_text(
        colour = "black",
        size = 7.25
      ),

      legend.key.width = unit(
        6,
        "mm"
      ),

      legend.key.height = unit(
        2.8,
        "mm"
      ),

      legend.spacing.x = unit(
        1.5,
        "mm"
      ),

      legend.margin = margin(
        t = 1,
        r = 0,
        b = 0,
        l = 0
      ),

      legend.box.margin = margin(
        0,
        0,
        0,
        0
      )
    )
}


# ============================================================
# Panel A: mapped mitochondrial-read percentages
#
# The plotting range extends to 118 for annotations.
#
# The automatic x-axis line is hidden and replaced with an
# explicit segment from 0 to 100 so the baseline does not
# continue beneath the annotations.
# ============================================================

panel_a <- ggplot() +

  geom_rect(
    data = proportion_data,

    aes(
      xmin = xmin,
      xmax = xmax,
      ymin = ymin,
      ymax = ymax,
      fill = decision_class
    ),

    colour = NA,

    show.legend = TRUE
  ) +

  geom_text(
    data = count_labels,

    aes(
      x = 109,
      y = species_id,
      label = annotation
    ),

    inherit.aes = FALSE,

    hjust = 0.5,
    vjust = 0.5,

    family = "sans",

    size = 1.72,

    lineheight = 0.88,

    colour = "black"
  ) +

  annotate(
    geom = "segment",

    x = 0,
    xend = 100,

    y = shared_y_limits[1],
    yend = shared_y_limits[1],

    colour = "black",

    linewidth = common_line_width,

    lineend = "butt"
  ) +

  decision_fill_scale +

  decision_fill_guide +

  scale_y_continuous(
    limits = shared_y_limits,

    breaks = shared_y_breaks,

    labels = species_label_function,

    expand = expansion(
      mult = c(
        0,
        0
      )
    )
  ) +

  scale_x_continuous(
    limits = c(
      0,
      118
    ),

    breaks = c(
      0,
      25,
      50,
      75,
      100
    ),

    labels = label_percent(
      scale = 1,
      accuracy = 1
    ),

    expand = expansion(
      mult = c(
        0,
        0
      )
    )
  ) +

  labs(
    tag = "A)",
    x = "Mapped mitochondrial reads (%)",
    y = NULL,
    fill = NULL
  ) +

  theme_scientific_reports() +

  theme(
    axis.text.y = element_text(
      colour = "black",
      size = 6.8,
      hjust = 1,

      margin = margin(
        r = 1.5
      )
    ),

    axis.line.x = element_blank(),

    axis.line.y = element_line(
      colour = "black",
      linewidth = common_line_width
    ),

    axis.ticks.x = element_line(
      colour = "black",
      linewidth = common_line_width
    ),

    axis.ticks.y = element_line(
      colour = "black",
      linewidth = common_line_width
    ),

    plot.tag = element_text(
      family = "sans",
      face = "bold",
      size = 9,
      colour = "black",
      hjust = 0,
      vjust = 1,

      margin = margin(
        b = 5
      )
    ),

    plot.tag.position = "topleft",

    plot.margin = margin(
      t = 3,
      r = 0,
      b = 2,
      l = 7
    )
  )


# ============================================================
# Panel B displayed read-length data
#
# The violin display is restricted to 200–20,000 bp.
#
# Complete canonical reads remain in read_data and in the
# exported length_plot_data_file.
# ============================================================

read_length_display_min <- 200
read_length_display_max <- 20000

read_length_plot_data <- read_data %>%
  filter(
    read_length_bases >=
      read_length_display_min,

    read_length_bases <=
      read_length_display_max
  ) %>%
  mutate(
    violin_y = case_when(
      decision_class == "Accepted" ~
        species_id - 0.20,

      decision_class == "Rejected" ~
        species_id + 0.20,

      TRUE ~
        as.numeric(
          species_id
        )
    )
  )


# ============================================================
# Export complete canonical read-level plotting data
# ============================================================

fwrite(
  read_data,
  length_plot_data_file,
  sep = "\t",
  quote = FALSE,
  na = "NA",
  compress = "gzip"
)


# ============================================================
# Display-range summary
# ============================================================

display_summary <- read_length_plot_data %>%
  count(
    decision_class,
    name = "displayed_reads"
  ) %>%
  mutate(
    canonical_reads = as.numeric(
      observed_counts[
        as.character(
          decision_class
        )
      ]
    ),

    excluded_from_display =
      canonical_reads -
      displayed_reads,

    displayed_percent =
      100 *
      displayed_reads /
      canonical_reads
  )

fwrite(
  display_summary,
  display_summary_file,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)


# ============================================================
# Panel B: read-length distributions
#
# Alpha = 1 keeps Accepted and Rejected colours identical to
# Panel A.
#
# Violin outline width matches the common axis-line width.
# ============================================================

read_length_breaks <- c(
  200,
  500,
  1000,
  3000,
  10000,
  20000
)

panel_b <- ggplot(
  read_length_plot_data,

  aes(
    x = read_length_bases,
    y = violin_y,
    fill = decision_class,

    group = interaction(
      species,
      decision_class,
      drop = TRUE
    )
  )
) +

  geom_violin(
    orientation = "y",

    width = 0.30,

    scale = "width",

    trim = FALSE,

    linewidth = common_line_width,

    colour = "grey25",

    alpha = 1,

    show.legend = FALSE
  ) +

  decision_fill_scale +

  scale_y_continuous(
    limits = shared_y_limits,

    breaks = shared_y_breaks,

    labels = NULL,

    expand = expansion(
      mult = c(
        0,
        0
      )
    )
  ) +

  scale_x_log10(
    limits = c(
      read_length_display_min,
      read_length_display_max
    ),

    breaks = read_length_breaks,

    labels = label_number(
      big.mark = ",",
      accuracy = 1
    ),

    minor_breaks = NULL,

    expand = expansion(
      mult = c(
        0.025,
        0.075
      )
    )
  ) +

  labs(
    tag = "B)",
    x = "Read length (bp)",
    y = NULL,
    fill = NULL
  ) +

  theme_scientific_reports() +

  theme(
    axis.text.y = element_blank(),

    axis.line.x = element_line(
      colour = "black",
      linewidth = common_line_width
    ),

    axis.line.y = element_line(
      colour = "black",
      linewidth = common_line_width
    ),

    axis.ticks.x = element_line(
      colour = "black",
      linewidth = common_line_width
    ),

    axis.ticks.y = element_line(
      colour = "black",
      linewidth = common_line_width
    ),

    plot.tag = element_text(
      family = "sans",
      face = "bold",
      size = 9,
      colour = "black",
      hjust = 0,
      vjust = 1,

      margin = margin(
        b = 5
      )
    ),

    plot.tag.position = "topleft",

    plot.margin = margin(
      t = 3,
      r = 8,
      b = 2,
      l = 7
    )
  )


# ============================================================
# Assemble panels using patchwork
#
# Panel B suppresses its copy of the legend, guaranteeing one
# Accepted/Rejected legend sourced from Panel A.
# ============================================================

publication_plot <- (
  panel_a +
    panel_b
) +

  patchwork::plot_layout(
    ncol = 2,

    widths = c(
      0.94,
      1.06
    ),

    guides = "collect"
  ) &

  theme(
    legend.position = "bottom",

    legend.title = element_blank(),

    legend.text = element_text(
      family = "sans",
      colour = "black",
      size = 7.25
    ),

    legend.key.width = unit(
      6,
      "mm"
    ),

    legend.key.height = unit(
      2.8,
      "mm"
    ),

    legend.spacing.x = unit(
      1.5,
      "mm"
    ),

    legend.margin = margin(
      t = 1,
      r = 0,
      b = 0,
      l = 0
    ),

    legend.box.margin = margin(
      0,
      0,
      0,
      0
    )
  )


# ============================================================
# Parse and build validation
#
# ggplot_build validates each panel.
# patchworkGrob validates the complete assembled figure.
# ============================================================

invisible(
  ggplot_build(
    panel_a
  )
)

invisible(
  ggplot_build(
    panel_b
  )
)

invisible(
  patchwork::patchworkGrob(
    publication_plot
  )
)

message(
  "PASS: panel plots and patchwork layout built successfully"
)


# ============================================================
# Display in RStudio
# ============================================================

print(
  publication_plot
)


# ============================================================
# Publication-quality PNG export
#
# Keep OneDrive paused while overwriting the figure files.
# ============================================================

ggsave(
  filename = figure_png,
  plot = publication_plot,
  width = 180,
  height = 132,
  units = "mm",
  dpi = 600,
  bg = "white",
  device = ragg::agg_png
)

if (
  !file.exists(figure_png) ||
  is.na(file.info(figure_png)$size) ||
  file.info(figure_png)$size <= 0
) {
  stop(
    "The PNG was not written successfully."
  )
}


# ============================================================
# PDF export
# ============================================================

pdf_saved <- tryCatch(
  {
    ggsave(
      filename = figure_pdf,
      plot = publication_plot,
      width = 180,
      height = 132,
      units = "mm",
      bg = "white",
      device = grDevices::cairo_pdf,
      family = "sans"
    )

    TRUE
  },

  error = function(error) {
    message(
      "cairo_pdf export failed; using standard PDF device."
    )

    ggsave(
      filename = figure_pdf,
      plot = publication_plot,
      width = 180,
      height = 132,
      units = "mm",
      bg = "white",
      device = grDevices::pdf,
      family = "sans",
      useDingbats = FALSE
    )

    TRUE
  }
)

if (
  !file.exists(figure_pdf) ||
  is.na(file.info(figure_pdf)$size) ||
  file.info(figure_pdf)$size <= 0
) {
  stop(
    "The PDF was not written successfully."
  )
}


# ============================================================
# SVG export
# ============================================================

ggsave(
  filename = figure_svg,
  plot = publication_plot,
  width = 180,
  height = 132,
  units = "mm",
  bg = "white",
  device = svglite::svglite
)

if (
  !file.exists(figure_svg) ||
  is.na(file.info(figure_svg)$size) ||
  file.info(figure_svg)$size <= 0
) {
  stop(
    "The SVG was not written successfully."
  )
}


# ============================================================
# Run log
# ============================================================

label_map <- tibble(
  canonical_species = species_order,
  displayed_species = species_abbreviations
)

log_lines <- c(
  paste(
    "Started:",
    format(
      Sys.time(),
      "%Y-%m-%dT%H:%M:%S%z"
    )
  ),

  paste(
    "Canonical-data directory:",
    canonical_data_dir
  ),

  paste(
    "Read-level input:",
    read_file
  ),

  paste(
    "Decision summary:",
    decision_summary_file
  ),

  paste(
    "Length summary:",
    length_summary_file
  ),

  "Sampling: none in canonical read table",

  paste(
    "Accepted canonical reads:",
    observed_counts[["Accepted"]]
  ),

  paste(
    "Rejected canonical reads:",
    observed_counts[["Rejected"]]
  ),

  "Decision colour Accepted: #1B9E77",

  "Decision colour Rejected: #7570B3",

  "Accepted and Rejected colours use alpha 1 in both panels",

  "Common line width: 0.30",

  "Axis lines, ticks and violin outlines use common line width",

  "Figure architecture: two aligned ggplots",

  "Layout package: patchwork",

  "Figure width: 180 mm",

  "Figure height: 132 mm",

  "PNG resolution: 600 dpi",

  "Panel tags: A), B)",

  "Panel-tag size: 9 pt",

  "Panel-tag bottom margin: 5 pt",

  "Panel-tag position: topleft",

  "Panel widths: 0.94 and 1.06",

  "Panels use identical numeric species y-axis limits",

  "Species display names use genus initial and species epithet",

  "Apis mellifera ligustica displayed as A. mellifera",

  "Accepted violin offset: -0.20",

  "Rejected violin offset: +0.20",

  "Violin width: 0.30",

  "Boxplots removed",

  "Accepted bar extends from 0 to accepted percentage",

  "Rejected bar extends from accepted percentage to 100",

  "Panel A x-axis baseline manually limited to 0-100",

  "Annotation format: centred two-line text block",

  "Panel A x-axis title: Mapped mitochondrial reads (%)",

  "Panel B x-axis title: Read length (bp)",

  "Panel B x-axis: log10",

  "Displayed read-length range: 200 to 20000 bp",

  paste0(
    "Displayed read-length labels: ",
    "200, 500, 1,000, 3,000, 10,000, 20,000"
  ),

  "Single figure legend sourced from Panel A",

  "Complete canonical read table retained in exported plotting data",

  paste(
    "Plot summary:",
    plot_summary_file
  ),

  paste(
    "Proportion plot data:",
    proportion_plot_data_file
  ),

  paste(
    "Read-level plot data:",
    length_plot_data_file
  ),

  paste(
    "Display summary:",
    display_summary_file
  ),

  paste(
    "Figure PNG:",
    figure_png
  ),

  paste(
    "Figure PDF:",
    figure_pdf
  ),

  paste(
    "Figure SVG:",
    figure_svg
  )
)

writeLines(
  log_lines,
  log_file
)

log_connection <- file(
  log_file,
  open = "at"
)

writeLines(
  "",
  log_connection
)

writeLines(
  "Overall canonical read-length summaries:",
  log_connection
)

writeLines(
  capture.output(
    print(
      plot_summary
    )
  ),
  log_connection
)

writeLines(
  "",
  log_connection
)

writeLines(
  "Read counts represented in displayed 200-20000 bp range:",
  log_connection
)

writeLines(
  capture.output(
    print(
      display_summary
    )
  ),
  log_connection
)

writeLines(
  "",
  log_connection
)

writeLines(
  "Canonical-to-display species labels:",
  log_connection
)

writeLines(
  capture.output(
    print(
      label_map
    )
  ),
  log_connection
)

writeLines(
  "",
  log_connection
)

writeLines(
  paste(
    "Finished:",
    format(
      Sys.time(),
      "%Y-%m-%dT%H:%M:%S%z"
    )
  ),
  log_connection
)

close(
  log_connection
)


# ============================================================
# Session information
# ============================================================

writeLines(
  capture.output(
    sessionInfo()
  ),
  session_file
)


# ============================================================
# Completion messages
# ============================================================

print(
  plot_summary
)

print(
  display_summary
)

print(
  label_map
)

message(
  "Wrote: ",
  plot_summary_file
)

message(
  "Wrote: ",
  proportion_plot_data_file
)

message(
  "Wrote: ",
  length_plot_data_file
)

message(
  "Wrote: ",
  display_summary_file
)

message(
  "Wrote: ",
  figure_png
)

message(
  "Wrote: ",
  figure_pdf
)

message(
  "Wrote: ",
  figure_svg
)

message(
  "Wrote: ",
  log_file
)

message(
  "Wrote: ",
  session_file
)
