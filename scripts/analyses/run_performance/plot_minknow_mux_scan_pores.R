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
  "analysis",
  "run_performance"
)

observations_file <- file.path(
  analysis_dir,
  "outputs",
  "minknow_mux_scan_pore_observations.tsv"
)

summary_file <- file.path(
  analysis_dir,
  "outputs",
  "minknow_mux_scan_pore_summary.tsv"
)

figure_dir <- file.path(
  analysis_dir,
  "figures"
)

log_dir <- file.path(
  analysis_dir,
  "logs"
)

dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

figure_png <- file.path(
  figure_dir,
  "supplementary_minknow_mux_scan_pores.png"
)

figure_pdf <- file.path(
  figure_dir,
  "supplementary_minknow_mux_scan_pores.pdf"
)

plot_data_file <- file.path(
  analysis_dir,
  "outputs",
  "minknow_mux_scan_pore_plot_data.tsv"
)

log_file <- file.path(
  log_dir,
  "plot_minknow_mux_scan_pores_run.log"
)

session_file <- file.path(
  log_dir,
  "R_sessionInfo_minknow_pores.txt"
)

if (!file.exists(observations_file)) {
  stop(paste("Missing observations table:", observations_file))
}

if (!file.exists(summary_file)) {
  stop(paste("Missing summary table:", summary_file))
}

observations <- fread(
  observations_file,
  sep = "\t",
  header = TRUE
) %>%
  as_tibble()

summary_data <- fread(
  summary_file,
  sep = "\t",
  header = TRUE
) %>%
  as_tibble()

required_columns <- c(
  "run",
  "run_label",
  "observation_index",
  "elapsed_hours_from_sequencing_start",
  "total_pores_found",
  "total_pores_percent_of_initial",
  "pores_available_for_immediate_sequencing",
  "available_pores_percent_of_initial_available",
  "available_fraction_of_total_pores"
)

missing_columns <- setdiff(required_columns, names(observations))
if (length(missing_columns) > 0) {
  stop(
    paste(
      "Missing observation columns:",
      paste(missing_columns, collapse = ", ")
    )
  )
}

if (nrow(observations) != 48) {
  stop(paste("Expected 48 observations, observed", nrow(observations)))
}

if (nrow(summary_data) != 3) {
  stop(paste("Expected 3 summary rows, observed", nrow(summary_data)))
}

run_levels <- c(
  "Whole-flowcell Metagenomic",
  "Whole-flowcell NAS",
  "Half-flowcell physical run"
)

run_colours <- c(
  "Whole-flowcell Metagenomic" = "#0072B2",
  "Whole-flowcell NAS" = "#D55E00",
  "Half-flowcell physical run" = "#009E73"
)

run_linetypes <- c(
  "Whole-flowcell Metagenomic" = "22",
  "Whole-flowcell NAS" = "solid",
  "Half-flowcell physical run" = "dotdash"
)

observations <- observations %>%
  mutate(
    run_label = factor(run_label, levels = run_levels)
  ) %>%
  arrange(run_label, elapsed_hours_from_sequencing_start)

series_check <- observations %>%
  group_by(run, run_label) %>%
  summarise(
    observations = n(),
    times_strictly_increasing = all(
      diff(elapsed_hours_from_sequencing_start) > 0
    ),
    available_not_greater_than_total = all(
      pores_available_for_immediate_sequencing <= total_pores_found
    ),
    .groups = "drop"
  )

if (
  any(series_check$observations != 16) ||
  any(series_check$times_strictly_increasing != TRUE) ||
  any(series_check$available_not_greater_than_total != TRUE)
) {
  stop("One or more mux-scan series validations failed")
}

plot_data <- observations %>%
  select(
    run,
    run_label,
    observation_index,
    elapsed_hours_from_sequencing_start,
    total_pores_found,
    pores_available_for_immediate_sequencing,
    total_pores_percent_of_initial,
    available_pores_percent_of_initial_available,
    available_fraction_of_total_pores
  ) %>%
  pivot_longer(
    cols = c(
      total_pores_found,
      pores_available_for_immediate_sequencing
    ),
    names_to = "metric_key",
    values_to = "pore_count"
  ) %>%
  mutate(
    metric = recode(
      metric_key,
      "total_pores_found" = "Total pores found during mux scan",
      "pores_available_for_immediate_sequencing" =
        "Pores available for immediate sequencing"
    ),
    metric = factor(
      metric,
      levels = c(
        "Pores available for immediate sequencing",
        "Total pores found during mux scan"
      )
    )
  )

fwrite(
  plot_data,
  plot_data_file,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

p_available <- ggplot(
  observations,
  aes(
    x = elapsed_hours_from_sequencing_start,
    y = pores_available_for_immediate_sequencing,
    colour = run_label,
    linetype = run_label,
    group = run_label
  )
) +
  geom_line(linewidth = 0.9, lineend = "round") +
  geom_point(size = 2.0) +
  scale_colour_manual(values = run_colours, breaks = run_levels) +
  scale_linetype_manual(values = run_linetypes, breaks = run_levels) +
  scale_x_continuous(
    limits = c(0, 24),
    breaks = seq(0, 24, by = 4),
    expand = expansion(mult = c(0, 0.01))
  ) +
  scale_y_continuous(
    labels = label_number(big.mark = ","),
    expand = expansion(mult = c(0.03, 0.06))
  ) +
  labs(
    x = NULL,
    y = "Pores available for\nimmediate sequencing",
    colour = NULL,
    linetype = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "bottom"
  )

p_total <- ggplot(
  observations,
  aes(
    x = elapsed_hours_from_sequencing_start,
    y = total_pores_found,
    colour = run_label,
    linetype = run_label,
    group = run_label
  )
) +
  geom_line(linewidth = 0.9, lineend = "round") +
  geom_point(size = 2.0) +
  scale_colour_manual(values = run_colours, breaks = run_levels) +
  scale_linetype_manual(values = run_linetypes, breaks = run_levels) +
  scale_x_continuous(
    limits = c(0, 24),
    breaks = seq(0, 24, by = 4),
    expand = expansion(mult = c(0, 0.01))
  ) +
  scale_y_continuous(
    labels = label_number(big.mark = ","),
    expand = expansion(mult = c(0.03, 0.06))
  ) +
  labs(
    x = "Time since sequencing start (hours)",
    y = "Total pores found\nduring mux scan",
    colour = NULL,
    linetype = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "bottom"
  )

combined_plot <- (
  p_available /
  p_total
) +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "A") &
  theme(legend.position = "bottom")

ggsave(
  figure_png,
  combined_plot,
  width = 10,
  height = 8,
  dpi = 300
)

ggsave(
  figure_pdf,
  combined_plot,
  width = 10,
  height = 8
)

log_lines <- c(
  paste("Started:", format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")),
  paste("Observations input:", observations_file),
  paste("Summary input:", summary_file),
  "Figure role: supplementary run-performance figure",
  "Panel A: pores available for immediate sequencing at each mux scan",
  "Panel B: total pores found during each mux scan",
  "Time origin: sequencing.start message for each physical run",
  "Lines: whole-flowcell Metagenomic, whole-flowcell NAS, half-flowcell physical run",
  "Half-flowcell scope: whole physical run; no channel-region allocation",
  paste("Plot data:", plot_data_file),
  paste("Figure PNG:", figure_png),
  paste("Figure PDF:", figure_pdf)
)

writeLines(log_lines, log_file)
log_connection <- file(log_file, open = "at")
writeLines("", log_connection)
writeLines("Series validation:", log_connection)
writeLines(capture.output(print(series_check)), log_connection)
writeLines("", log_connection)
writeLines("Mux-scan summary:", log_connection)
writeLines(capture.output(print(summary_data)), log_connection)
writeLines("", log_connection)
writeLines(
  paste("Finished:", format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")),
  log_connection
)
close(log_connection)

writeLines(capture.output(sessionInfo()), session_file)

print(summary_data)
message("Wrote: ", plot_data_file)
message("Wrote: ", figure_png)
message("Wrote: ", figure_pdf)
