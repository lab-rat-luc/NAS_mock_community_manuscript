#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(scales)
})

started_time <- Sys.time()

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

canonical_dir <- file.path(
  project_dir,
  "canonical_data"
)


analysis_dir <- file.path(project_dir, "analysis", "ra_full_vs_half")
input_file <- file.path(canonical_dir, "outputs", "canonical_mito_mapping_summary.tsv")
reference_file <- file.path(canonical_dir, "inputs", "reference_manifest.tsv")
output_dir <- file.path(analysis_dir, "outputs")
figure_dir <- file.path(analysis_dir, "figures")
log_dir <- file.path(analysis_dir, "logs")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

regression_input_file <- file.path(output_dir, "ra_full_vs_half_regression_input.tsv")
statistics_file <- file.path(output_dir, "ra_full_vs_half_regression_statistics.tsv")
coefficients_file <- file.path(output_dir, "ra_full_vs_half_model_coefficients.tsv")
validation_file <- file.path(output_dir, "ra_full_vs_half_validation.tsv")
colour_file <- file.path(output_dir, "ra_full_vs_half_species_colour_map.tsv")
export_file <- file.path(output_dir, "ra_full_vs_half_figure_export_status.tsv")
figure_png <- file.path(figure_dir, "ra_full_vs_half_two_panel.png")
figure_pdf <- file.path(figure_dir, "ra_full_vs_half_two_panel.pdf")
figure_svg <- file.path(figure_dir, "ra_full_vs_half_two_panel.svg")
log_file <- file.path(log_dir, "ra_full_vs_half_run.log")
session_file <- file.path(log_dir, "sessionInfo.txt")
provenance_file <- file.path(log_dir, "ra_full_vs_half_provenance.tsv")

write_tsv <- function(x, path) {
  data.table::fwrite(x, path, sep = "\t", quote = FALSE, na = "NA")
}

abbreviate_species <- function(x) {
  paste0(substr(x, 1, 1), ". ", sub("^[^ ]+ ", "", x))
}

fit_panel <- function(d) {
  panel_name <- unique(as.character(d$panel))
  if (length(panel_name) != 1L) stop("fit_panel received more than one panel")

  model <- lm(log10_half_RA ~ log10_full_RA, data = d)
  model_summary <- summary(model)
  coefficient_matrix <- model_summary$coefficients
  confidence_intervals <- confint(model, level = 0.95)
  f_statistic <- model_summary$fstatistic
  model_p <- pf(
    unname(f_statistic[1]),
    unname(f_statistic[2]),
    unname(f_statistic[3]),
    lower.tail = FALSE
  )
  pearson_test <- cor.test(d$log10_full_RA, d$log10_half_RA, method = "pearson")

  statistics <- data.frame(
    panel = panel_name,
    n = nrow(d),
    predictor = "log10 full-flowcell mitochondrial RA percent",
    response = "log10 half-flowcell mitochondrial RA percent",
    transformation = "log10",
    zero_handling = "none required",
    intercept = unname(coefficient_matrix["(Intercept)", "Estimate"]),
    slope = unname(coefficient_matrix["log10_full_RA", "Estimate"]),
    intercept_standard_error = unname(coefficient_matrix["(Intercept)", "Std. Error"]),
    slope_standard_error = unname(coefficient_matrix["log10_full_RA", "Std. Error"]),
    intercept_CI_lower = unname(confidence_intervals["(Intercept)", 1]),
    intercept_CI_upper = unname(confidence_intervals["(Intercept)", 2]),
    slope_CI_lower = unname(confidence_intervals["log10_full_RA", 1]),
    slope_CI_upper = unname(confidence_intervals["log10_full_RA", 2]),
    R_squared = model_summary$r.squared,
    adjusted_R_squared = model_summary$adj.r.squared,
    Pearson_r = unname(pearson_test$estimate),
    Pearson_p_value = pearson_test$p.value,
    model_p_value = unname(model_p),
    residual_standard_error = model_summary$sigma,
    stringsAsFactors = FALSE
  )

  coefficients <- data.frame(
    panel = panel_name,
    term = rownames(coefficient_matrix),
    estimate = coefficient_matrix[, "Estimate"],
    standard_error = coefficient_matrix[, "Std. Error"],
    t_value = coefficient_matrix[, "t value"],
    p_value = coefficient_matrix[, "Pr(>|t|)"],
    CI_lower = confidence_intervals[, 1],
    CI_upper = confidence_intervals[, 2],
    row.names = NULL,
    stringsAsFactors = FALSE
  )

  list(statistics = statistics, coefficients = coefficients)
}

export_plot <- function(plot_object, destination, device, dpi = NULL) {
  temporary_file <- tempfile(fileext = paste0(".", tools::file_ext(destination)))
  error_message <- ""

  tryCatch(
    {
      if (is.null(dpi)) {
        ggsave(
          temporary_file,
          plot_object,
          device = device,
          width = 14,
          height = 7.5,
          units = "in",
          bg = "white"
        )
      } else {
        ggsave(
          temporary_file,
          plot_object,
          device = device,
          width = 14,
          height = 7.5,
          units = "in",
          dpi = dpi,
          bg = "white"
        )
      }
    },
    error = function(e) {
      error_message <<- conditionMessage(e)
    }
  )

  temporary_ok <- file.exists(temporary_file) && file.info(temporary_file)$size > 0
  copied <- FALSE

  if (temporary_ok && error_message == "") {
    if (file.exists(destination)) file.remove(destination)
    copied <- file.copy(temporary_file, destination, overwrite = TRUE)
  }

  if (file.exists(temporary_file)) unlink(temporary_file)

  destination_ok <- file.exists(destination) && file.info(destination)$size > 0
  success <- temporary_ok && copied && destination_ok && error_message == ""

  data.frame(
    format = toupper(device),
    destination = destination,
    success = success,
    file_size_bytes = if (destination_ok) file.info(destination)$size else 0,
    error = error_message,
    stringsAsFactors = FALSE
  )
}

if (!file.exists(input_file)) stop(paste("Missing canonical input:", input_file))
if (!file.exists(reference_file)) stop(paste("Missing reference manifest:", reference_file))

mito <- as.data.frame(data.table::fread(input_file), stringsAsFactors = FALSE)
references <- as.data.frame(data.table::fread(reference_file), stringsAsFactors = FALSE)

required_mito <- c("dataset", "reference", "mitochondrial_RA_percent")
required_reference <- c("reference", "species", "target_status")
missing_mito <- setdiff(required_mito, names(mito))
missing_reference <- setdiff(required_reference, names(references))

if (length(missing_mito) > 0) stop(paste("Missing canonical columns:", paste(missing_mito, collapse = ", ")))
if (length(missing_reference) > 0) stop(paste("Missing manifest columns:", paste(missing_reference, collapse = ", ")))
if (nrow(references) != 16) stop(paste("Expected 16 manifest rows, observed", nrow(references)))
if (anyDuplicated(references$reference) > 0) stop("Duplicate references in manifest")

dataset_map <- data.frame(
  dataset = c("full_meta", "half_meta", "full_nas_accepted", "half_nas_accepted"),
  panel = c("Metagenomic", "Metagenomic", "NAS", "NAS"),
  flowcell_scope = c("full", "half", "full", "half"),
  stringsAsFactors = FALSE
)

selected <- mito %>%
  filter(dataset %in% dataset_map$dataset) %>%
  select(dataset, reference, mitochondrial_RA_percent) %>%
  left_join(dataset_map, by = "dataset") %>%
  left_join(
    references %>% select(reference, species, target_status),
    by = "reference"
  ) %>%
  as.data.frame(stringsAsFactors = FALSE)

if (nrow(selected) != 64) stop(paste("Expected 64 selected rows, observed", nrow(selected)))
if (anyDuplicated(selected[c("dataset", "reference")]) > 0) stop("Duplicate dataset-reference records")
if (anyNA(selected$species)) stop("Missing species metadata after manifest join")
if (anyNA(selected$mitochondrial_RA_percent)) stop("Missing mitochondrial RA values")
if (any(selected$mitochondrial_RA_percent <= 0)) stop("Non-positive mitochondrial RA values detected")

dataset_counts <- selected %>% count(dataset, name = "observed_rows")
if (nrow(dataset_counts) != 4 || any(dataset_counts$observed_rows != 16)) {
  stop("Expected 16 rows in each of the four selected datasets")
}

ra_sums <- selected %>%
  group_by(dataset) %>%
  summarise(
    mitochondrial_RA_percent_sum = sum(mitochondrial_RA_percent),
    .groups = "drop"
  )

ra_tolerance <- 0.000001
if (any(abs(ra_sums$mitochondrial_RA_percent_sum - 100) > ra_tolerance)) {
  stop("At least one dataset mitochondrial RA sum differs from 100 percent")
}

species_labels <- references %>%
  transmute(
    reference,
    species_display = ifelse(grepl("^Apis mellifera", species), "Apis mellifera", species),
    species_abbreviation = abbreviate_species(species_display)
  )

wide <- selected %>%
  select(panel, reference, target_status, flowcell_scope, mitochondrial_RA_percent) %>%
  pivot_wider(names_from = flowcell_scope, values_from = mitochondrial_RA_percent) %>%
  left_join(species_labels, by = "reference") %>%
  arrange(panel, reference) %>%
  mutate(
    log10_full_RA = log10(full),
    log10_half_RA = log10(half)
  ) %>%
  as.data.frame(stringsAsFactors = FALSE)

if (anyNA(wide$full) || anyNA(wide$half)) stop("Missing paired full-flowcell or half-flowcell values")
if (anyDuplicated(wide[c("panel", "reference")]) > 0) stop("Duplicate panel-reference records")
if (any(!is.finite(wide$log10_full_RA)) || any(!is.finite(wide$log10_half_RA))) {
  stop("Non-finite log10 values detected")
}

panel_counts <- wide %>% count(panel, name = "observed_species")
if (nrow(panel_counts) != 2 || any(panel_counts$observed_species != 16)) {
  stop("Expected 16 species in each panel")
}

panel_levels <- c("Metagenomic", "NAS")
wide$panel <- factor(wide$panel, levels = panel_levels)
write_tsv(wide, regression_input_file)

meta_result <- fit_panel(wide %>% filter(panel == "Metagenomic"))
nas_result <- fit_panel(wide %>% filter(panel == "NAS"))

statistics <- bind_rows(meta_result$statistics, nas_result$statistics) %>%
  mutate(panel = factor(panel, levels = panel_levels)) %>%
  arrange(panel) %>%
  mutate(panel = as.character(panel))

coefficients <- bind_rows(meta_result$coefficients, nas_result$coefficients) %>%
  mutate(panel = factor(panel, levels = panel_levels)) %>%
  arrange(panel, term) %>%
  mutate(panel = as.character(panel))

write_tsv(statistics, statistics_file)
write_tsv(coefficients, coefficients_file)

species_palette <- c(
  "#0072B2", "#D55E00", "#009E73", "#CC79A7",
  "#E69F00", "#56B4E9", "#7B2CBF", "#8C564B",
  "#00A6D6", "#E7298A", "#66A61E", "#7570B3",
  "#A6761D", "#1B9E77", "#4D4D4D", "#B2182B"
)

species_levels <- wide %>%
  filter(panel == "Metagenomic") %>%
  mutate(mean_meta_RA = (full + half) / 2) %>%
  arrange(mean_meta_RA) %>%
  pull(species_display) %>%
  as.character()

if (length(species_levels) != 16 || anyDuplicated(species_levels) > 0) {
  stop("Species colour mapping requires 16 unique species")
}

colour_map <- setNames(species_palette, species_levels)
abbreviation_map <- setNames(abbreviate_species(species_levels), species_levels)
wide$species_display <- factor(wide$species_display, levels = species_levels)

species_colour_map <- data.frame(
  species_display = species_levels,
  species_abbreviation = unname(abbreviation_map[species_levels]),
  colour = unname(colour_map[species_levels]),
  stringsAsFactors = FALSE
)
write_tsv(species_colour_map, colour_file)

minimum_ra <- min(c(wide$full, wide$half), na.rm = TRUE)
maximum_ra <- max(c(wide$full, wide$half), na.rm = TRUE)

axis_lower <- 0.03
axis_upper <- 50
axis_breaks <- c(0.03, 0.1, 0.3, 1, 3, 10, 30)
axis_labels <- c("0.03", "0.1", "0.3", "1", "3", "10", "30")

if (minimum_ra < axis_lower || maximum_ra > axis_upper) {
  stop(
    paste0(
      "Observed RA values fall outside the configured axis range. ",
      "Observed range: ",
      sprintf("%.6f", minimum_ra),
      " to ",
      sprintf("%.6f", maximum_ra),
      "; configured range: ",
      axis_lower,
      " to ",
      axis_upper
    )
  )
}

annotation_data <- statistics %>%
  mutate(
    panel = factor(panel, levels = panel_levels),
    x = 0.22,
    y = 6.5,
    label = sprintf("R\u00B2 = %.3f", R_squared)
  )

finite_logs <- all(is.finite(wide$log10_full_RA)) && all(is.finite(wide$log10_half_RA))
duplicate_pairs <- sum(duplicated(wide[c("panel", "reference")]))
ra_sums_pass <- all(abs(ra_sums$mitochondrial_RA_percent_sum - 100) <= ra_tolerance)

validation <- data.frame(
  validation = c(
    "selected_rows",
    "rows_per_dataset",
    "species_per_panel",
    "positive_RA_values",
    "finite_log10_values",
    "duplicate_panel_reference_records",
    "dataset_RA_sums",
    "R_squared_range"
  ),
  passed = c(
    nrow(selected) == 64,
    all(dataset_counts$observed_rows == 16),
    all(panel_counts$observed_species == 16),
    minimum_ra > 0,
    finite_logs,
    duplicate_pairs == 0,
    ra_sums_pass,
    all(statistics$R_squared >= 0 & statistics$R_squared <= 1)
  ),
  stringsAsFactors = FALSE
)
write_tsv(validation, validation_file)
if (any(!validation$passed)) stop("One or more validation checks failed")

plot_object <- ggplot(
  wide,
  aes(x = full, y = half, fill = species_display)
) +
  geom_abline(
    intercept = 0,
    slope = 1,
    linetype = "22",
    colour = "#D73027",
    linewidth = 0.65,
    lineend = "butt"
  ) +
  geom_point(
    shape = 21,
    size = 3.1,
    stroke = 0.35,
    colour = "black",
    alpha = 0.95
  ) +
  facet_wrap(
    ~ panel,
    nrow = 1,
    labeller = as_labeller(
      c(
        "Metagenomic" = "A) Metagenomic",
        "NAS" = "B) NAS"
      )
    )
  ) +
  scale_fill_manual(
    values = colour_map,
    breaks = species_levels,
    labels = unname(abbreviation_map[species_levels]),
    drop = FALSE
  ) +
  scale_x_log10(
    limits = c(axis_lower, axis_upper),
    breaks = axis_breaks,
    labels = axis_labels,
    expand = expansion(mult = c(0, 0))
  ) +
  scale_y_log10(
    limits = c(axis_lower, axis_upper),
    breaks = axis_breaks,
    labels = axis_labels,
    expand = expansion(mult = c(0, 0))
  ) +
  coord_equal(clip = "off") +
  geom_text(
    data = annotation_data,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    hjust = 0,
    colour = "black",
    size = 4.3
  ) +
  labs(
    x = "Full-flowcell RA (%)",
    y = "Half-flowcell RA (%)",
    fill = NULL
  ) +
  guides(
    fill = guide_legend(
      nrow = 4,
      byrow = FALSE,
      override.aes = list(size = 3, alpha = 1)
    )
  ) +
  theme_bw(base_size = 13, base_family = "sans") +
  theme(
    panel.border = element_rect(colour = "#666666", fill = NA, linewidth = 0.6),
    panel.grid.major = element_line(colour = "#EFEFEF", linewidth = 0.4),
    panel.grid.minor = element_blank(),
    panel.spacing = grid::unit(0.85, "cm"),
    strip.background = element_blank(),
    strip.text = element_text(
      face = "bold",
      colour = "black",
      size = 13,
      margin = margin(t = 2, r = 0, b = 10, l = 0)
    ),
    axis.title = element_text(colour = "black", size = 12),
    axis.title.x = element_text(
      margin = margin(t = 8, r = 0, b = 0, l = 0)
    ),
    axis.title.y = element_text(
      margin = margin(t = 0, r = 8, b = 0, l = 0)
    ),
    axis.text = element_text(colour = "black", size = 9.5),
    legend.position = "bottom",
    legend.justification = "center",
    legend.key = element_blank(),
    legend.text = element_text(face = "italic", colour = "black", size = 9),
    legend.box.spacing = grid::unit(0.08, "cm"),
    legend.margin = margin(t = -2, r = 0, b = 0, l = 0),
    plot.margin = margin(t = 16, r = 10, b = 2, l = 10)
  )

print(plot_object)

figure_export_status <- bind_rows(
  export_plot(plot_object, figure_png, "png", dpi = 300),
  export_plot(plot_object, figure_pdf, "pdf"),
  export_plot(plot_object, figure_svg, "svg")
)
write_tsv(figure_export_status, export_file)

if (any(!figure_export_status$success)) {
  print(figure_export_status)
  stop("One or more figure exports failed")
}

finished_time <- Sys.time()

provenance <- data.frame(
  field = c(
    "analysis_name",
    "started",
    "finished",
    "canonical_input",
    "reference_manifest",
    "RA_definition",
    "metagenomic_comparison",
    "accepted_NAS_comparison",
    "regression_model",
    "zero_handling",
    "species_per_panel"
  ),
  value = c(
    "ra_full_vs_half",
    format(started_time, "%Y-%m-%dT%H:%M:%S%z"),
    format(finished_time, "%Y-%m-%dT%H:%M:%S%z"),
    input_file,
    reference_file,
    "Species CIGAR M, = and X aligned query bases divided by the corresponding total across all 16 mitochondrial references",
    "half_meta versus full_meta",
    "half_nas_accepted versus full_nas_accepted",
    "lm(log10(half mitochondrial RA percent) ~ log10(full mitochondrial RA percent))",
    "none required",
    "16"
  ),
  stringsAsFactors = FALSE
)
write_tsv(provenance, provenance_file)

log_lines <- c(
  paste("Started:", format(started_time, "%Y-%m-%dT%H:%M:%S%z")),
  paste("Finished:", format(finished_time, "%Y-%m-%dT%H:%M:%S%z")),
  paste("Canonical input:", input_file),
  paste("Reference manifest:", reference_file),
  "RA definition: species CIGAR M + = + X aligned query bases divided by the corresponding total across all 16 mitochondrial references",
  "Model: lm(log10(half mitochondrial RA percent) ~ log10(full mitochondrial RA percent))",
  "Zero handling: none required",
  "Species included per panel: 16",
  "",
  "Regression statistics:",
  capture.output(print(statistics)),
  "",
  "Model coefficients:",
  capture.output(print(coefficients)),
  "",
  "Validation results:",
  capture.output(print(validation)),
  "",
  "Figure export status:",
  capture.output(print(figure_export_status))
)
writeLines(log_lines, log_file, useBytes = TRUE)
writeLines(capture.output(sessionInfo()), session_file, useBytes = TRUE)

print(statistics)
print(figure_export_status)
message("Wrote figures to: ", figure_dir)
message("Wrote outputs to: ", output_dir)
message("Wrote logs to: ", log_dir)
message("PASS: analysis and exports completed successfully")
