#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(scales)
  library(cowplot)
  library(ragg)
  library(svglite)
})

# ============================================================
# Figure 2: Run performance
#
# A): Whole-flowcell active channels
# B): Half-flowcell active channels
# C): Whole-flowcell cumulative pass bases
# D): Half-flowcell cumulative pass bases
#
# Each panel retains its own ticks and numerical labels.
# Y-axis titles are shown only in the left column.
# One x-axis title is shared across the complete figure.
# ============================================================


# ============================================================
# Current project run-performance directory
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

analysis_dir <- file.path(
  project_dir,
  "analysis",
  "run_performance"
)

input_file <- file.path(
  analysis_dir,
  "outputs",
  "run_performance_binned.tsv"
)

summary_input_file <- file.path(
  analysis_dir,
  "outputs",
  "run_performance_summary.tsv"
)

validation_input_file <- file.path(
  analysis_dir,
  "outputs",
  "run_performance_validation.tsv"
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

plot_data_file <- file.path(
  output_dir,
  "run_performance_plot_data_publication.tsv"
)

plot_summary_file <- file.path(
  output_dir,
  "run_performance_plot_summary_publication.tsv"
)

figure_png <- file.path(
  figure_dir,
  "figure_2_run_performance_publication.png"
)

figure_pdf <- file.path(
  figure_dir,
  "figure_2_run_performance_publication.pdf"
)

figure_svg <- file.path(
  figure_dir,
  "figure_2_run_performance_publication.svg"
)

log_file <- file.path(
  log_dir,
  "plot_run_performance_publication.log"
)

session_file <- file.path(
  log_dir,
  "R_sessionInfo_run_performance_publication.txt"
)


# ============================================================
# Confirm required files
# ============================================================

required_files <- c(
  input_file,
  summary_input_file,
  validation_input_file
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
# Read validated run-performance outputs
# ============================================================

binned <- fread(
  input_file,
  sep = "\t",
  header = TRUE,
  integer64 = "double"
) %>%
  as_tibble()

run_summary <- fread(
  summary_input_file,
  sep = "\t",
  header = TRUE,
  integer64 = "double"
) %>%
  as_tibble()

validation <- fread(
  validation_input_file,
  sep = "\t",
  header = TRUE,
  integer64 = "double"
) %>%
  as_tibble()


# ============================================================
# Structural validation
# ============================================================

required_binned_columns <- c(
  "dataset",
  "method",
  "run_type",
  "time_bin_sec",
  "time_bin_hours",
  "productive_channels",
  "starting_productive_channels",
  "productive_channels_pct_of_first_nonzero_bin",
  "cumulative_pass_reads",
  "cumulative_pass_bases",
  "cumulative_pass_bases_Gb"
)

missing_binned_columns <- setdiff(
  required_binned_columns,
  names(binned)
)

if (length(missing_binned_columns) > 0) {
  stop(
    paste(
      "Missing binned columns:",
      paste(
        missing_binned_columns,
        collapse = ", "
      )
    )
  )
}

if (nrow(binned) != 576) {
  stop(
    paste(
      "Expected 576 binned rows, observed",
      nrow(binned)
    )
  )
}

if (nrow(run_summary) != 4) {
  stop(
    paste(
      "Expected 4 summary rows, observed",
      nrow(run_summary)
    )
  )
}

if (nrow(validation) != 4) {
  stop(
    paste(
      "Expected 4 validation rows, observed",
      nrow(validation)
    )
  )
}

validation_fields <- c(
  "row_partition_passed",
  "read_sum_passed",
  "base_sum_passed",
  "bin_count_passed"
)

if (
  !all(
    unlist(
      validation[
        validation_fields
      ]
    )
  )
) {
  stop(
    "One or more run-performance validations failed"
  )
}

expected_datasets <- c(
  "full_meta",
  "full_nas",
  "half_meta",
  "half_nas"
)

if (!setequal(
  binned$dataset,
  expected_datasets
)) {
  stop(
    "Unexpected analytical datasets detected"
  )
}

bin_counts <- binned %>%
  count(
    dataset
  )

if (
  nrow(bin_counts) != 4 ||
  any(bin_counts$n != 144)
) {
  stop(
    "Expected 144 ten-minute bins for each dataset"
  )
}


# ============================================================
# Validate cumulative series
# ============================================================

monotonicity_check <- binned %>%
  arrange(
    dataset,
    time_bin_sec
  ) %>%
  group_by(
    dataset
  ) %>%
  summarise(
    cumulative_reads_monotonic = all(
      diff(
        cumulative_pass_reads
      ) >= 0
    ),

    cumulative_bases_monotonic = all(
      diff(
        cumulative_pass_bases
      ) >= 0
    ),

    cumulative_Gb_monotonic = all(
      diff(
        cumulative_pass_bases_Gb
      ) >= -1e-12
    ),

    .groups = "drop"
  )

if (
  any(
    monotonicity_check$
    cumulative_reads_monotonic != TRUE
  ) ||
  any(
    monotonicity_check$
    cumulative_bases_monotonic != TRUE
  ) ||
  any(
    monotonicity_check$
    cumulative_Gb_monotonic != TRUE
  )
) {
  stop(
    "One or more cumulative series are not monotonic"
  )
}


# ============================================================
# Factor ordering and plotting data
# ============================================================

method_levels <- c(
  "NAS",
  "Metagenomic"
)

run_type_levels <- c(
  "Whole-flowcell",
  "Half-flowcell"
)

plot_data <- binned %>%
  mutate(
    method = factor(
      method,
      levels = method_levels
    ),

    run_type = factor(
      run_type,
      levels = run_type_levels
    )
  ) %>%
  arrange(
    run_type,
    method,
    time_bin_sec
  ) %>%
  select(
    dataset,
    method,
    run_type,
    time_bin_sec,
    time_bin_hours,
    productive_channels,
    starting_productive_channels,
    productive_channels_pct_of_first_nonzero_bin,
    cumulative_pass_bases_Gb
  )

fwrite(
  plot_data,
  plot_data_file,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

plot_summary <- run_summary %>%
  mutate(
    method = factor(
      method,
      levels = method_levels
    ),

    run_type = factor(
      run_type,
      levels = run_type_levels
    ),

    final_active_channels_pct_of_starting =
      100 *
      final_productive_channels /
      starting_productive_channels
  ) %>%
  arrange(
    run_type,
    method
  )

fwrite(
  plot_summary,
  plot_summary_file,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)


# ============================================================
# Consistent manuscript colour system
# ============================================================

method_colours <- c(
  "NAS" = "#D55E00",
  "Metagenomic" = "#0072B2"
)

method_linetypes <- c(
  "NAS" = "solid",
  "Metagenomic" = "longdash"
)


# ============================================================
# Publication theme
# ============================================================

theme_scientific_reports <- function() {
  theme_classic(
    base_size = 8,
    base_family = "sans"
  ) +
    theme(
      axis.line = element_line(
        colour = "black",
        linewidth = 0.3
      ),

      axis.ticks = element_line(
        colour = "black",
        linewidth = 0.3
      ),

      axis.ticks.length = unit(
        1.25,
        "mm"
      ),

      axis.text = element_text(
        colour = "black",
        size = 7.25
      ),

      axis.title = element_text(
        colour = "black",
        size = 7.75,
        face = "plain"
      ),

      plot.title = element_text(
        colour = "black",
        size = 9,
        face = "bold",
        hjust = 0.5,

        margin = margin(
          b = 2.5
        )
      ),

      legend.position = "bottom",

      legend.direction = "horizontal",

      legend.title = element_blank(),

      legend.text = element_text(
        colour = "black",
        size = 7.5
      ),

      legend.key.width = unit(
        9,
        "mm"
      ),

      legend.key.height = unit(
        3.2,
        "mm"
      ),

      legend.spacing.x = unit(
        1.2,
        "mm"
      ),

      panel.grid = element_blank(),

      panel.border = element_blank(),

      plot.margin = margin(
        1.5,
        2,
        1.5,
        2
      )
    )
}


# ============================================================
# Shared scales
# ============================================================

method_colour_scale <- scale_colour_manual(
  values = method_colours,
  breaks = method_levels,
  drop = FALSE
)

method_linetype_scale <- scale_linetype_manual(
  values = method_linetypes,
  breaks = method_levels,
  drop = FALSE
)

# The right-side expansion leaves room after the 24-hour tick
# so the "24" label is not cropped.
x_scale_all <- scale_x_continuous(
  limits = c(
    0,
    24
  ),

  breaks = seq(
    0,
    24,
    by = 2
  ),

  labels = label_number(
    accuracy = 1
  ),

  expand = expansion(
    add = c(
      0,
      0.45
    )
  )
)


# ============================================================
# Panel A: Whole-flowcell active channels
# ============================================================

panel_a <- ggplot(
  plot_data %>%
    filter(
      run_type == "Whole-flowcell"
    ),

  aes(
    x = time_bin_hours,
    y =
      productive_channels_pct_of_first_nonzero_bin,
    colour = method,
    linetype = method,
    group = method
  )
) +
  geom_line(
    linewidth = 0.58,
    lineend = "round"
  ) +
  method_colour_scale +
  method_linetype_scale +
  x_scale_all +
  scale_y_continuous(
    limits = c(
      30,
      101
    ),

    breaks = seq(
      30,
      100,
      by = 10
    ),

    labels = label_number(
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
    title = "Whole-flowcell",
    x = NULL,
    y = "% of starting active channels",
    colour = NULL,
    linetype = NULL
  ) +
  theme_scientific_reports() +
  theme(
    legend.position = "none",

    axis.title.y = element_text(
      margin = margin(
        r = 2.5
      )
    ),

    plot.margin = margin(
      3.5,
      1.5,
      0.5,
      2
    )
  )


# ============================================================
# Panel B: Half-flowcell active channels
# ============================================================

panel_b <- ggplot(
  plot_data %>%
    filter(
      run_type == "Half-flowcell"
    ),

  aes(
    x = time_bin_hours,
    y =
      productive_channels_pct_of_first_nonzero_bin,
    colour = method,
    linetype = method,
    group = method
  )
) +
  geom_line(
    linewidth = 0.58,
    lineend = "round"
  ) +
  method_colour_scale +
  method_linetype_scale +
  x_scale_all +
  scale_y_continuous(
    limits = c(
      30,
      101
    ),

    breaks = seq(
      30,
      100,
      by = 10
    ),

    labels = label_number(
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
    title = "Half-flowcell",
    x = NULL,
    y = NULL,
    colour = NULL,
    linetype = NULL
  ) +
  theme_scientific_reports() +
  theme(
    legend.position = "none",

    plot.margin = margin(
      3.5,
      2,
      0.5,
      1.5
    )
  )


# ============================================================
# Panel C: Whole-flowcell cumulative pass bases
# ============================================================

panel_c <- ggplot(
  plot_data %>%
    filter(
      run_type == "Whole-flowcell"
    ),

  aes(
    x = time_bin_hours,
    y = cumulative_pass_bases_Gb,
    colour = method,
    linetype = method,
    group = method
  )
) +
  geom_line(
    linewidth = 0.58,
    lineend = "round"
  ) +
  method_colour_scale +
  method_linetype_scale +
  x_scale_all +
  scale_y_continuous(
    limits = c(
      0,
      12
    ),

    breaks = seq(
      0,
      12,
      by = 2
    ),

    labels = label_number(
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
    x = NULL,
    y = "Cumulative pass bases (Gb)",
    colour = NULL,
    linetype = NULL
  ) +
  theme_scientific_reports() +
  theme(
    axis.title.y = element_text(
      margin = margin(
        r = 2.5
      )
    ),

    plot.margin = margin(
      0.5,
      1.5,
      0.5,
      2
    )
  )


# ============================================================
# Panel D: Half-flowcell cumulative pass bases
# ============================================================

panel_d <- ggplot(
  plot_data %>%
    filter(
      run_type == "Half-flowcell"
    ),

  aes(
    x = time_bin_hours,
    y = cumulative_pass_bases_Gb,
    colour = method,
    linetype = method,
    group = method
  )
) +
  geom_line(
    linewidth = 0.58,
    lineend = "round"
  ) +
  method_colour_scale +
  method_linetype_scale +
  x_scale_all +
  scale_y_continuous(
    limits = c(
      0,
      5
    ),

    breaks = seq(
      0,
      5,
      by = 1
    ),

    labels = label_number(
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
    x = NULL,
    y = NULL,
    colour = NULL,
    linetype = NULL
  ) +
  theme_scientific_reports() +
  theme(
    plot.margin = margin(
      0.5,
      2,
      0.5,
      1.5
    )
  )


# ============================================================
# Extract one shared legend
# ============================================================

legend_source <- panel_c +
  theme(
    legend.position = "bottom",

    legend.margin = margin(
      0,
      0,
      0,
      0
    ),

    legend.box.margin = margin(
      0,
      0,
      0,
      0
    )
  )

shared_legend <- cowplot::get_legend(
  legend_source
)

panel_c <- panel_c +
  theme(
    legend.position = "none"
  )

panel_d <- panel_d +
  theme(
    legend.position = "none"
  )


# ============================================================
# Assemble aligned 2 x 2 panel grid
#
# A closing parenthesis is included in each panel label.
# ============================================================

panel_grid <- cowplot::plot_grid(
  panel_a,
  panel_b,
  panel_c,
  panel_d,

  ncol = 2,

  rel_widths = c(
    1,
    1
  ),

  rel_heights = c(
    1,
    1
  ),

  align = "hv",

  axis = "tblr",

  labels = c(
    "A)",
    "B)",
    "C)",
    "D)"
  ),

  label_size = 9,

  label_fontface = "bold",

  label_x = c(
    0.008,
    0.008,
    0.008,
    0.008
  ),

  label_y = c(
    0.994,
    0.994,
    0.994,
    0.994
  ),

  hjust = 0,

  vjust = 1
)


# ============================================================
# Shared x-axis title
# ============================================================

shared_x_title <- cowplot::ggdraw() +
  cowplot::draw_label(
    "Time since run start (hours)",

    fontfamily = "sans",

    fontface = "plain",

    size = 7.75,

    x = 0.5,

    y = 0.5,

    hjust = 0.5,

    vjust = 0.5
  )


# ============================================================
# Final Figure 2 assembly
# ============================================================

performance_plot <- cowplot::plot_grid(
  panel_grid,
  shared_x_title,
  shared_legend,

  ncol = 1,

  rel_heights = c(
    1,
    0.036,
    0.052
  )
)

print(
  performance_plot
)


# ============================================================
# Publication-quality PNG export
#
# Render to a temporary file first, then copy into OneDrive.
# ============================================================

temporary_png <- tempfile(
  pattern = "figure_2_run_performance_",
  fileext = ".png"
)

ragg::agg_png(
  filename = temporary_png,
  width = 180,
  height = 116,
  units = "mm",
  res = 600,
  background = "white"
)

print(
  performance_plot
)

grDevices::dev.off()

if (
  !file.exists(temporary_png) ||
  is.na(file.info(temporary_png)$size) ||
  file.info(temporary_png)$size <= 0
) {
  stop(
    "The temporary PNG was not rendered successfully."
  )
}

if (file.exists(figure_png)) {
  unlink(
    figure_png,
    force = TRUE
  )
}

if (file.exists(figure_png)) {
  stop(
    paste0(
      "The previous PNG is locked and could not be removed:\n",
      figure_png
    )
  )
}

png_copy_succeeded <- file.copy(
  from = temporary_png,
  to = figure_png,
  overwrite = TRUE,
  copy.date = FALSE
)

unlink(
  temporary_png,
  force = TRUE
)

if (
  !png_copy_succeeded ||
  !file.exists(figure_png) ||
  is.na(file.info(figure_png)$size) ||
  file.info(figure_png)$size <= 0
) {
  stop(
    paste0(
      "The PNG rendered successfully but could not be copied ",
      "to the figure directory:\n",
      figure_png
    )
  )
}

message(
  "Wrote and validated PNG: ",
  figure_png,
  " (",
  format(
    file.info(figure_png)$size,
    big.mark = ",",
    scientific = FALSE
  ),
  " bytes)"
)


# ============================================================
# PDF export
# ============================================================

pdf_saved <- tryCatch(
  {
    ggsave(
      filename = figure_pdf,
      plot = performance_plot,
      width = 180,
      height = 116,
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
      plot = performance_plot,
      width = 180,
      height = 116,
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

if (requireNamespace(
  "svglite",
  quietly = TRUE
)) {
  ggsave(
    filename = figure_svg,
    plot = performance_plot,
    width = 180,
    height = 116,
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
}


# ============================================================
# Run log
# ============================================================

log_lines <- c(
  paste(
    "Started:",
    format(
      Sys.time(),
      "%Y-%m-%dT%H:%M:%S%z"
    )
  ),

  paste(
    "Project directory:",
    project_dir
  ),

  paste(
    "Binned input:",
    input_file
  ),

  paste(
    "Summary input:",
    summary_input_file
  ),

  paste(
    "Validation input:",
    validation_input_file
  ),

  "Figure width: 180 mm",

  "Figure height: 116 mm",

  "PNG device: ragg::agg_png",

  "Raster resolution: 600 dpi",

  "Font family: sans",

  "Panel labels: A), B), C), D)",

  "Base font size: 8 pt",

  "Axis tick-label size: 7.25 pt",

  "Axis-title size: 7.75 pt",

  "Panel-title size: 9 pt",

  "Panel-label size: 9 pt",

  "Legend-text size: 7.5 pt",

  "Line width: 0.58",

  "NAS line type: solid",

  "Metagenomic line type: longdash",

  "Time resolution: 10-minute bins",

  "Analysis window: 24 hours",

  "Time-axis interval: 2 hours",

  "Right-side x-axis expansion: 0.45 hours",

  paste0(
    "Active-channel definition: unique channel values ",
    "producing at least one pass read within the bin"
  ),

  paste0(
    "% of starting active channels: active channels in ",
    "the bin divided by active channels in the first ",
    "non-zero bin for the same analytical dataset"
  ),

  paste0(
    "Interpretation warning: active channels are defined ",
    "from pass-read production and are not direct ",
    "measurements of pores available for immediate sequencing"
  ),

  paste0(
    "Pass-base definition: cumulative sum of ",
    "sequence_length_template from pre-demultiplexing ",
    "offline Dorado SUP pass-read summaries"
  ),

  paste0(
    "Whole-flowcell scope: all offline SUP pass reads ",
    "from the complete physical run"
  ),

  paste0(
    "Half-flowcell allocation: NAS channels 1-256; ",
    "Metagenomic channels 257-512"
  ),

  "Every panel retains axis ticks and numerical values",

  "Y-axis titles appear only in the left column",

  "One shared x-axis title is used beneath the panel grid",

  "Active-channel y-axis: 30 to 101 percent",

  "Whole-flowcell cumulative pass-base y-axis: 0 to 12 Gb",

  "Half-flowcell cumulative pass-base y-axis: 0 to 5 Gb",

  paste(
    "Plot data:",
    plot_data_file
  ),

  paste(
    "Plot summary:",
    plot_summary_file
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
  "Run-performance summary:",
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
  "Cumulative monotonicity validation:",
  log_connection
)

writeLines(
  capture.output(
    print(
      monotonicity_check
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

writeLines(
  capture.output(
    sessionInfo()
  ),
  session_file
)

message(
  "Wrote: ",
  plot_data_file
)

message(
  "Wrote: ",
  plot_summary_file
)

message(
  "Wrote: ",
  log_file
)

message(
  "Wrote: ",
  session_file
)
