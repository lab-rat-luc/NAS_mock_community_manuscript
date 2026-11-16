#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(scales)
  library(patchwork)
})

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
  "canonical_data",
  "analyses",
  "coverage_over_time"
)

time_series_file <- file.path(
  analysis_dir,
  "outputs",
  "target_coverage_over_time.tsv"
)

threshold_file <- file.path(
  analysis_dir,
  "outputs",
  "target_coverage_threshold_times.tsv"
)

final_summary_file <- file.path(
  analysis_dir,
  "outputs",
  "target_coverage_final_summary.tsv"
)

join_validation_file <- file.path(
  analysis_dir,
  "outputs",
  "target_timing_join_validation.tsv"
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
  figure_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  log_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

breadth_png <- file.path(
  figure_dir,
  "target_coverage_breadth_over_time.png"
)

breadth_pdf <- file.path(
  figure_dir,
  "target_coverage_breadth_over_time.pdf"
)

depth_png <- file.path(
  figure_dir,
  "target_mean_depth_over_time.png"
)

depth_pdf <- file.path(
  figure_dir,
  "target_mean_depth_over_time.pdf"
)

combined_png <- file.path(
  figure_dir,
  "target_coverage_breadth_and_depth.png"
)

combined_pdf <- file.path(
  figure_dir,
  "target_coverage_breadth_and_depth.pdf"
)

threshold_png <- file.path(
  figure_dir,
  "target_coverage_threshold_times.png"
)

threshold_pdf <- file.path(
  figure_dir,
  "target_coverage_threshold_times.pdf"
)

plot_summary_file <- file.path(
  analysis_dir,
  "outputs",
  "coverage_plot_summary.tsv"
)

log_file <- file.path(
  log_dir,
  "plot_coverage_over_time_run.log"
)

session_file <- file.path(
  log_dir,
  "R_sessionInfo.txt"
)

required_files <- c(
  time_series_file,
  threshold_file,
  final_summary_file,
  join_validation_file
)

missing_files <- required_files[
  !file.exists(required_files)
]

if (length(missing_files) > 0) {
  stop(
    paste(
      "Missing required coverage files:",
      paste(
        missing_files,
        collapse = ", "
      )
    )
  )
}

coverage_data <- fread(
  time_series_file,
  sep = "\t",
  header = TRUE
) %>%
  as_tibble()

threshold_data <- fread(
  threshold_file,
  sep = "\t",
  header = TRUE
) %>%
  as_tibble()

final_summary <- fread(
  final_summary_file,
  sep = "\t",
  header = TRUE
) %>%
  as_tibble()

join_validation <- fread(
  join_validation_file,
  sep = "\t",
  header = TRUE
) %>%
  as_tibble()

required_coverage_columns <- c(
  "dataset",
  "method",
  "reference",
  "species",
  "reference_length",
  "time_hours",
  "cumulative_reads",
  "cumulative_cigar_aligned_query_bases",
  "covered_reference_bases_1x",
  "coverage_breadth_percent_1x",
  "cumulative_mean_depth"
)

required_threshold_columns <- c(
  "dataset",
  "method",
  "reference",
  "species",
  "breadth_threshold_percent",
  "first_hour_reached",
  "reached_within_24_hours"
)

required_final_columns <- c(
  "dataset",
  "method",
  "reference",
  "species",
  "final_hour",
  "final_cumulative_reads",
  "final_coverage_breadth_percent_1x",
  "final_cumulative_mean_depth"
)

missing_coverage_columns <- setdiff(
  required_coverage_columns,
  names(coverage_data)
)

missing_threshold_columns <- setdiff(
  required_threshold_columns,
  names(threshold_data)
)

missing_final_columns <- setdiff(
  required_final_columns,
  names(final_summary)
)

if (length(missing_coverage_columns) > 0) {
  stop(
    paste(
      "Missing time-series columns:",
      paste(
        missing_coverage_columns,
        collapse = ", "
      )
    )
  )
}

if (length(missing_threshold_columns) > 0) {
  stop(
    paste(
      "Missing threshold columns:",
      paste(
        missing_threshold_columns,
        collapse = ", "
      )
    )
  )
}

if (length(missing_final_columns) > 0) {
  stop(
    paste(
      "Missing final-summary columns:",
      paste(
        missing_final_columns,
        collapse = ", "
      )
    )
  )
}

if (nrow(coverage_data) != 144) {
  stop(
    paste(
      "Expected 144 coverage rows, observed",
      nrow(coverage_data)
    )
  )
}

if (nrow(threshold_data) != 30) {
  stop(
    paste(
      "Expected 30 threshold rows, observed",
      nrow(threshold_data)
    )
  )
}

if (nrow(final_summary) != 6) {
  stop(
    paste(
      "Expected 6 final-summary rows, observed",
      nrow(final_summary)
    )
  )
}

if (
  any(
    join_validation$
      timing_join_passed != TRUE
  )
) {
  stop(
    "One or more timing joins did not pass validation"
  )
}

if (
  any(
    coverage_data$time_hours < 1 |
      coverage_data$time_hours > 24
  )
) {
  stop(
    "Coverage time values fall outside hours 1 to 24"
  )
}

if (
  any(
    coverage_data$
      coverage_breadth_percent_1x < 0 |
      coverage_data$
        coverage_breadth_percent_1x > 100
  )
) {
  stop(
    "Coverage breadth falls outside 0 to 100 percent"
  )
}

if (
  any(
    coverage_data$cumulative_mean_depth < 0
  )
) {
  stop(
    "Negative cumulative mean depth detected"
  )
}

group_counts <- coverage_data %>%
  count(
    dataset,
    reference
  )

if (
  nrow(group_counts) != 6 ||
  any(group_counts$n != 24)
) {
  stop(
    paste(
      "Expected 24 hourly records for each of",
      "six dataset-reference combinations"
    )
  )
}

monotonicity_check <- coverage_data %>%
  arrange(
    dataset,
    reference,
    time_hours
  ) %>%
  group_by(
    dataset,
    reference
  ) %>%
  summarise(
    reads_monotonic =
      all(
        diff(cumulative_reads) >= 0
      ),
    bases_monotonic =
      all(
        diff(
          cumulative_cigar_aligned_query_bases
        ) >= 0
      ),
    breadth_monotonic =
      all(
        diff(
          coverage_breadth_percent_1x
        ) >= -1e-12
      ),
    depth_monotonic =
      all(
        diff(
          cumulative_mean_depth
        ) >= -1e-12
      ),
    .groups = "drop"
  )

if (
  any(
    monotonicity_check$
      reads_monotonic != TRUE
  ) ||
  any(
    monotonicity_check$
      bases_monotonic != TRUE
  ) ||
  any(
    monotonicity_check$
      breadth_monotonic != TRUE
  ) ||
  any(
    monotonicity_check$
      depth_monotonic != TRUE
  )
) {
  stop(
    "One or more cumulative time series are not monotonic"
  )
}

species_levels <- c(
  "Diuraphis noxia",
  "Trichogramma pretiosum",
  "Myzus persicae"
)

method_levels <- c(
  "Metagenomic",
  "NAS"
)

coverage_data <- coverage_data %>%
  mutate(
    species = factor(
      species,
      levels = species_levels
    ),
    method = factor(
      method,
      levels = method_levels
    )
  ) %>%
  arrange(
    species,
    method,
    time_hours
  )

threshold_data <- threshold_data %>%
  mutate(
    species = factor(
      species,
      levels = species_levels
    ),
    method = factor(
      method,
      levels = method_levels
    ),
    threshold_label = paste0(
      format(
        breadth_threshold_percent,
        trim = TRUE,
        scientific = FALSE
      ),
      "%"
    ),
    threshold_label = factor(
      threshold_label,
      levels = c(
        "50%",
        "75%",
        "90%",
        "95%",
        "99%"
      )
    ),
    first_hour_plot = if_else(
      reached_within_24_hours,
      as.numeric(first_hour_reached),
      NA_real_
    )
  )

final_summary <- final_summary %>%
  mutate(
    species = factor(
      species,
      levels = species_levels
    ),
    method = factor(
      method,
      levels = method_levels
    )
  ) %>%
  arrange(
    species,
    method
  )

method_colours <- c(
  "Metagenomic" = "#0072B2",
  "NAS" = "#D55E00"
)

breadth_plot <- ggplot(
  coverage_data,
  aes(
    x = time_hours,
    y = coverage_breadth_percent_1x,
    colour = method,
    group = method
  )
) +
  geom_line(
    linewidth = 0.9,
    lineend = "round"
  ) +
  geom_point(
    size = 1.4,
    alpha = 0.8
  ) +
  facet_wrap(
    ~ species,
    nrow = 1
  ) +
  scale_colour_manual(
    values = method_colours
  ) +
  scale_x_continuous(
    limits = c(
      1,
      24
    ),
    breaks = c(
      1,
      4,
      8,
      12,
      16,
      20,
      24
    )
  ) +
  scale_y_continuous(
    limits = c(
      0,
      100
    ),
    breaks = seq(
      0,
      100,
      by = 20
    ),
    labels = label_percent(
      scale = 1
    )
  ) +
  labs(
    x = "Time since run start (hours)",
    y = "Mitogenome breadth at >=1x (%)",
    colour = NULL
  ) +
  theme_bw(
    base_size = 11
  ) +
  theme(
    strip.text = element_text(
      face = "italic"
    ),
    panel.grid.minor =
      element_blank(),
    legend.position =
      "bottom"
  )

depth_plot <- ggplot(
  coverage_data,
  aes(
    x = time_hours,
    y = cumulative_mean_depth,
    colour = method,
    group = method
  )
) +
  geom_line(
    linewidth = 0.9,
    lineend = "round"
  ) +
  geom_point(
    size = 1.4,
    alpha = 0.8
  ) +
  facet_wrap(
    ~ species,
    nrow = 1,
    scales = "free_y"
  ) +
  scale_colour_manual(
    values = method_colours
  ) +
  scale_x_continuous(
    limits = c(
      1,
      24
    ),
    breaks = c(
      1,
      4,
      8,
      12,
      16,
      20,
      24
    )
  ) +
  scale_y_continuous(
    labels = label_number(
      accuracy = 0.1
    ),
    expand = expansion(
      mult = c(
        0,
        0.05
      )
    )
  ) +
  labs(
    x = "Time since run start (hours)",
    y = "Mean mitogenome coverage depth (x)",
    colour = NULL
  ) +
  theme_bw(
    base_size = 11
  ) +
  theme(
    strip.text = element_text(
      face = "italic"
    ),
    panel.grid.minor =
      element_blank(),
    legend.position =
      "bottom"
  )

combined_plot <- (
  breadth_plot /
    depth_plot
) +
  plot_layout(
    guides = "collect"
  ) +
  plot_annotation(
    tag_levels = "A"
  ) &
  theme(
    legend.position = "bottom"
  )

threshold_plot <- ggplot(
  threshold_data %>%
    filter(
      reached_within_24_hours
    ),
  aes(
    x = threshold_label,
    y = first_hour_plot,
    colour = method,
    group = method
  )
) +
  geom_line(
    linewidth = 0.75
  ) +
  geom_point(
    size = 2.3
  ) +
  facet_wrap(
    ~ species,
    nrow = 1
  ) +
  scale_colour_manual(
    values = method_colours
  ) +
  scale_y_continuous(
    limits = c(
      0,
      24
    ),
    breaks = seq(
      0,
      24,
      by = 4
    )
  ) +
  labs(
    x = "Coverage breadth threshold",
    y = "First hourly bin reaching threshold",
    colour = NULL
  ) +
  theme_bw(
    base_size = 11
  ) +
  theme(
    strip.text = element_text(
      face = "italic"
    ),
    panel.grid.minor =
      element_blank(),
    legend.position =
      "bottom"
  )

plot_summary <- final_summary %>%
  select(
    dataset,
    method,
    species,
    final_hour,
    final_cumulative_reads,
    final_coverage_breadth_percent_1x,
    final_cumulative_mean_depth,
    total_target_alignments_all_times,
    target_alignments_within_24_hours,
    target_alignments_after_24_hours
  )

fwrite(
  plot_summary,
  plot_summary_file,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

ggsave(
  breadth_png,
  breadth_plot,
  width = 11,
  height = 4,
  dpi = 300
)

ggsave(
  breadth_pdf,
  breadth_plot,
  width = 11,
  height = 4
)

ggsave(
  depth_png,
  depth_plot,
  width = 11,
  height = 4,
  dpi = 300
)

ggsave(
  depth_pdf,
  depth_plot,
  width = 11,
  height = 4
)

ggsave(
  combined_png,
  combined_plot,
  width = 11,
  height = 8,
  dpi = 300
)

ggsave(
  combined_pdf,
  combined_plot,
  width = 11,
  height = 8
)

ggsave(
  threshold_png,
  threshold_plot,
  width = 11,
  height = 4,
  dpi = 300
)

ggsave(
  threshold_pdf,
  threshold_plot,
  width = 11,
  height = 4
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
    "Coverage time series:",
    time_series_file
  ),
  paste(
    "Threshold table:",
    threshold_file
  ),
  paste(
    "Final summary:",
    final_summary_file
  ),
  paste(
    "Timing join validation:",
    join_validation_file
  ),
  paste0(
    "Breadth definition: reference positions ",
    "with depth >= 1 divided by reference length"
  ),
  paste0(
    "Mean-depth definition: cumulative M + = + X ",
    "reference coverage divided by reference length"
  ),
  "Timing resolution: cumulative hourly bins",
  "Analysis window: hours 1 to 24",
  "Reference topology: linear",
  "Datasets: full_meta and full_nas_accepted",
  "Targets: Diuraphis noxia, Myzus persicae, Trichogramma pretiosum",
  paste(
    "Breadth figure:",
    breadth_png
  ),
  paste(
    "Depth figure:",
    depth_png
  ),
  paste(
    "Combined figure:",
    combined_png
  ),
  paste(
    "Threshold figure:",
    threshold_png
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
  "Final 24-hour summary:",
  log_connection
)

writeLines(
  capture.output(
    print(plot_summary)
  ),
  log_connection
)

writeLines(
  "",
  log_connection
)

writeLines(
  "Cumulative-series monotonicity validation:",
  log_connection
)

writeLines(
  capture.output(
    print(monotonicity_check)
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

close(log_connection)

writeLines(
  capture.output(
    sessionInfo()
  ),
  session_file
)

print(plot_summary)

message(
  "Wrote: ",
  plot_summary_file
)

message(
  "Wrote: ",
  breadth_png
)

message(
  "Wrote: ",
  breadth_pdf
)

message(
  "Wrote: ",
  depth_png
)

message(
  "Wrote: ",
  depth_pdf
)

message(
  "Wrote: ",
  combined_png
)

message(
  "Wrote: ",
  combined_pdf
)

message(
  "Wrote: ",
  threshold_png
)

message(
  "Wrote: ",
  threshold_pdf
)
