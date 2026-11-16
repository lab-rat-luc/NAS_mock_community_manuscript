#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
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

analysis_dir <- file.path(
  project_dir,
  "analysis",
  "ra_deviation"
)

output_dir <- file.path(analysis_dir, "outputs")
figure_dir <- file.path(analysis_dir, "figures")
log_dir <- file.path(analysis_dir, "logs")

mapping_file <- file.path(
  canonical_dir,
  "outputs",
  "canonical_mito_mapping_summary.tsv"
)

reference_file <- file.path(
  canonical_dir,
  "inputs",
  "reference_manifest.tsv"
)

colour_map_file <- file.path(
  project_dir,
  "analysis",
  "ra_full_vs_half",
  "outputs",
  "ra_full_vs_half_species_colour_map.tsv"
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

analysis_output <- file.path(
  output_dir,
  "expected_vs_observed_linear_regression_input.tsv"
)

statistics_output <- file.path(
  output_dir,
  "expected_vs_observed_linear_regression_statistics.tsv"
)

coefficients_output <- file.path(
  output_dir,
  "expected_vs_observed_linear_regression_model_coefficients.tsv"
)

diagnostics_output <- file.path(
  output_dir,
  "expected_vs_observed_linear_regression_diagnostics.tsv"
)

validation_output <- file.path(
  output_dir,
  "expected_vs_observed_linear_regression_validation.tsv"
)

figure_export_output <- file.path(
  output_dir,
  "expected_vs_observed_linear_regression_figure_export_status.tsv"
)

meta_png <- file.path(
  figure_dir,
  "expected_vs_observed_linear_regression_meta.png"
)

meta_pdf <- file.path(
  figure_dir,
  "expected_vs_observed_linear_regression_meta.pdf"
)

meta_svg <- file.path(
  figure_dir,
  "expected_vs_observed_linear_regression_meta.svg"
)

nas_png <- file.path(
  figure_dir,
  "expected_vs_observed_linear_regression_nas.png"
)

nas_pdf <- file.path(
  figure_dir,
  "expected_vs_observed_linear_regression_nas.pdf"
)

nas_svg <- file.path(
  figure_dir,
  "expected_vs_observed_linear_regression_nas.svg"
)

log_file <- file.path(
  log_dir,
  "ra_deviation_linear_regression_run.log"
)

session_file <- file.path(
  log_dir,
  "ra_deviation_linear_regression_sessionInfo.txt"
)

provenance_file <- file.path(
  log_dir,
  "ra_deviation_linear_regression_provenance.tsv"
)

console_file <- file.path(
  log_dir,
  "ra_deviation_linear_regression_console.log"
)

write_tsv <- function(x, path) {
  data.table::fwrite(
    x,
    file = path,
    sep = "\t",
    quote = FALSE,
    na = "NA"
  )
}

format_plotmath_p <- function(p) {
  if (is.na(p)) {
    return("italic(P) == NA")
  }

  if (p < 0.001) {
    return("italic(P) < 0.001")
  }

  paste0(
    "italic(P) == ",
    formatC(
      p,
      format = "f",
      digits = 3
    )
  )
}

export_plot <- function(
    plot_object,
    destination,
    device,
    dpi = NULL
) {
  temporary_file <- tempfile(
    fileext = paste0(
      ".",
      tools::file_ext(destination)
    )
  )

  error_message <- ""

  tryCatch(
    {
      if (is.null(dpi)) {
        ggsave(
          filename = temporary_file,
          plot = plot_object,
          device = device,
          width = 14,
          height = 7.5,
          units = "in",
          bg = "white"
        )
      } else {
        ggsave(
          filename = temporary_file,
          plot = plot_object,
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

  temporary_ok <- (
    file.exists(temporary_file) &&
      file.info(temporary_file)$size > 0
  )

  copied <- FALSE

  if (temporary_ok && error_message == "") {
    if (file.exists(destination)) {
      file.remove(destination)
    }

    copied <- file.copy(
      temporary_file,
      destination,
      overwrite = TRUE
    )
  }

  if (file.exists(temporary_file)) {
    unlink(temporary_file)
  }

  destination_ok <- (
    file.exists(destination) &&
      file.info(destination)$size > 0
  )

  success <- (
    temporary_ok &&
      copied &&
      destination_ok &&
      error_message == ""
  )

  data.frame(
    format = toupper(device),
    destination = destination,
    success = success,
    file_size_bytes = if (
      destination_ok
    ) {
      file.info(destination)$size
    } else {
      0
    },
    error = error_message,
    stringsAsFactors = FALSE
  )
}

calculate_model_outputs <- function(group_data) {
  dataset_values <- unique(group_data$dataset)
  method_values <- unique(group_data$method)
  scope_values <- unique(group_data$flowcell_scope)
  panel_values <- unique(group_data$panel)

  if (
    length(dataset_values) != 1 ||
    length(method_values) != 1 ||
    length(scope_values) != 1 ||
    length(panel_values) != 1
  ) {
    stop("Model group contains inconsistent dataset metadata")
  }

  model <- lm(
    log10_observed_RA ~ log10_expected_RA,
    data = group_data
  )

  model_summary <- summary(model)
  model_coefficients <- coef(model_summary)
  coefficient_intervals <- confint(
    model,
    level = 0.95
  )

  model_f_statistic <- model_summary$fstatistic

  model_p_value <- pf(
    model_f_statistic[["value"]],
    model_f_statistic[["numdf"]],
    model_f_statistic[["dendf"]],
    lower.tail = FALSE
  )

  pearson_test <- cor.test(
    group_data$log10_expected_RA,
    group_data$log10_observed_RA,
    method = "pearson"
  )

  influence <- influence.measures(model)
  leverage_values <- hatvalues(model)
  cooks_values <- cooks.distance(model)
  standardised_residuals <- rstandard(model)
  studentised_residuals <- rstudent(model)
  fitted_values <- fitted(model)
  raw_residuals <- residuals(model)

  cooks_threshold <- 4 / nrow(group_data)

  statistics <- data.frame(
    dataset = dataset_values,
    method = method_values,
    flowcell_scope = scope_values,
    panel = panel_values,
    n = nrow(group_data),
    predictor = "log10 expected DNA RA percent",
    response = "log10 observed mitochondrial RA percent",
    transformation = "log10",
    zero_handling = "none required",
    intercept = unname(coef(model)[["(Intercept)"]]),
    slope = unname(coef(model)[["log10_expected_RA"]]),
    intercept_standard_error =
      model_coefficients["(Intercept)", "Std. Error"],
    slope_standard_error =
      model_coefficients["log10_expected_RA", "Std. Error"],
    intercept_CI_lower =
      coefficient_intervals["(Intercept)", 1],
    intercept_CI_upper =
      coefficient_intervals["(Intercept)", 2],
    slope_CI_lower =
      coefficient_intervals["log10_expected_RA", 1],
    slope_CI_upper =
      coefficient_intervals["log10_expected_RA", 2],
    R_squared = unname(model_summary$r.squared),
    adjusted_R_squared =
      unname(model_summary$adj.r.squared),
    Pearson_r = unname(pearson_test$estimate),
    Pearson_p_value = pearson_test$p.value,
    model_p_value = unname(model_p_value),
    residual_standard_error =
      unname(model_summary$sigma),
    mean_signed_log10_deviation = mean(
      group_data$log10_signed_deviation
    ),
    median_signed_log10_deviation = median(
      group_data$log10_signed_deviation
    ),
    mean_absolute_log10_deviation = mean(
      group_data$log10_absolute_deviation
    ),
    median_absolute_log10_deviation = median(
      group_data$log10_absolute_deviation
    ),
    cooks_distance_threshold = cooks_threshold,
    influential_observations = sum(
      cooks_values > cooks_threshold
    ),
    stringsAsFactors = FALSE
  )

  coefficients <- data.frame(
    dataset = dataset_values,
    method = method_values,
    flowcell_scope = scope_values,
    panel = panel_values,
    term = rownames(model_coefficients),
    estimate = model_coefficients[, "Estimate"],
    standard_error =
      model_coefficients[, "Std. Error"],
    t_value = model_coefficients[, "t value"],
    p_value = model_coefficients[, "Pr(>|t|)"],
    CI_lower = coefficient_intervals[, 1],
    CI_upper = coefficient_intervals[, 2],
    stringsAsFactors = FALSE,
    row.names = NULL
  )

  diagnostics <- group_data %>%
    mutate(
      fitted_log10_observed_RA =
        as.numeric(fitted_values),
      raw_residual = as.numeric(raw_residuals),
      standardised_residual =
        as.numeric(standardised_residuals),
      studentised_residual =
        as.numeric(studentised_residuals),
      leverage = as.numeric(leverage_values),
      cooks_distance = as.numeric(cooks_values),
      cooks_distance_threshold = cooks_threshold,
      influential_flag =
        cooks_distance > cooks_distance_threshold,
      dfbetas_intercept =
        as.numeric(influence$infmat[, "dfb.1_"]),
      dfbetas_slope =
        as.numeric(influence$infmat[, "dfb.l10_"])
    ) %>%
    select(
      dataset,
      method,
      flowcell_scope,
      panel,
      reference,
      species,
      species_display,
      target_status,
      expected_DNA_RA_percent,
      observed_mitochondrial_RA_percent,
      log10_expected_RA,
      log10_observed_RA,
      fitted_log10_observed_RA,
      raw_residual,
      standardised_residual,
      studentised_residual,
      leverage,
      cooks_distance,
      cooks_distance_threshold,
      influential_flag,
      dfbetas_intercept,
      dfbetas_slope
    ) %>%
    as.data.frame(
      stringsAsFactors = FALSE
    )

  prediction_grid <- data.frame(
    log10_expected_RA = seq(
      min(group_data$log10_expected_RA),
      max(group_data$log10_expected_RA),
      length.out = 200
    )
  )

  predictions <- predict(
    model,
    newdata = prediction_grid,
    interval = "confidence",
    level = 0.95
  )

  prediction_data <- data.frame(
    dataset = dataset_values,
    method = method_values,
    flowcell_scope = scope_values,
    panel = panel_values,
    log10_expected_RA =
      prediction_grid$log10_expected_RA,
    fitted_log10_observed_RA =
      predictions[, "fit"],
    lower_confidence_limit =
      predictions[, "lwr"],
    upper_confidence_limit =
      predictions[, "upr"],
    stringsAsFactors = FALSE
  )

  list(
    model = model,
    statistics = statistics,
    coefficients = coefficients,
    diagnostics = diagnostics,
    predictions = prediction_data
  )
}

required_files <- c(
  mapping_file,
  reference_file,
  colour_map_file
)

missing_files <- required_files[
  !file.exists(required_files)
]

if (length(missing_files) > 0) {
  stop(
    paste(
      "Missing required files:",
      paste(missing_files, collapse = "; ")
    )
  )
}

mapping <- as.data.frame(
  data.table::fread(mapping_file),
  stringsAsFactors = FALSE
)

references <- as.data.frame(
  data.table::fread(reference_file),
  stringsAsFactors = FALSE
)

species_colours <- as.data.frame(
  data.table::fread(colour_map_file),
  stringsAsFactors = FALSE
)

required_mapping_columns <- c(
  "dataset",
  "reference",
  "mitochondrial_RA_percent"
)

required_reference_columns <- c(
  "reference",
  "species",
  "target_status",
  "expected_DNA_RA_percent"
)

required_colour_columns <- c(
  "species_display",
  "species_abbreviation",
  "colour"
)

missing_mapping_columns <- setdiff(
  required_mapping_columns,
  names(mapping)
)

missing_reference_columns <- setdiff(
  required_reference_columns,
  names(references)
)

missing_colour_columns <- setdiff(
  required_colour_columns,
  names(species_colours)
)

if (length(missing_mapping_columns) > 0) {
  stop(
    paste(
      "Missing mapping columns:",
      paste(missing_mapping_columns, collapse = ", ")
    )
  )
}

if (length(missing_reference_columns) > 0) {
  stop(
    paste(
      "Missing reference columns:",
      paste(missing_reference_columns, collapse = ", ")
    )
  )
}

if (length(missing_colour_columns) > 0) {
  stop(
    paste(
      "Missing colour-map columns:",
      paste(missing_colour_columns, collapse = ", ")
    )
  )
}

if (nrow(references) != 16) {
  stop(
    paste(
      "Expected 16 reference rows, observed",
      nrow(references)
    )
  )
}

if (nrow(species_colours) != 16) {
  stop(
    paste(
      "Expected 16 species colours, observed",
      nrow(species_colours)
    )
  )
}

if (anyDuplicated(references$reference) > 0) {
  stop("Duplicate references detected in reference manifest")
}

if (anyDuplicated(species_colours$species_display) > 0) {
  stop("Duplicate species detected in species colour map")
}

dataset_map <- data.frame(
  dataset = c(
    "full_meta",
    "half_meta",
    "full_nas_accepted",
    "half_nas_accepted"
  ),
  method = c(
    "Metagenomic",
    "Metagenomic",
    "NAS",
    "NAS"
  ),
  flowcell_scope = c(
    "Full-flowcell",
    "Half-flowcell",
    "Full-flowcell",
    "Half-flowcell"
  ),
  panel = c(
    "A) Full-flowcell",
    "B) Half-flowcell",
    "A) Full-flowcell",
    "B) Half-flowcell"
  ),
  stringsAsFactors = FALSE
)

species_labels <- references %>%
  transmute(
    reference,
    species_display = ifelse(
      grepl("^Apis mellifera", species),
      "Apis mellifera",
      species
    )
  )

analysis_data <- mapping %>%
  filter(dataset %in% dataset_map$dataset) %>%
  select(
    dataset,
    reference,
    mitochondrial_RA_percent
  ) %>%
  left_join(
    dataset_map,
    by = "dataset"
  ) %>%
  left_join(
    references %>%
      select(
        reference,
        species,
        target_status,
        expected_DNA_RA_percent
      ),
    by = "reference"
  ) %>%
  left_join(
    species_labels,
    by = "reference"
  ) %>%
  rename(
    observed_mitochondrial_RA_percent =
      mitochondrial_RA_percent
  ) %>%
  mutate(
    expected_DNA_RA_percent = as.numeric(
      expected_DNA_RA_percent
    ),
    observed_mitochondrial_RA_percent =
      as.numeric(
        observed_mitochondrial_RA_percent
      ),
    log10_expected_RA = log10(
      expected_DNA_RA_percent
    ),
    log10_observed_RA = log10(
      observed_mitochondrial_RA_percent
    ),
    log10_signed_deviation =
      log10_observed_RA - log10_expected_RA,
    log10_absolute_deviation = abs(
      log10_signed_deviation
    ),
    observed_to_expected_ratio =
      observed_mitochondrial_RA_percent /
      expected_DNA_RA_percent,
    direction = case_when(
      log10_signed_deviation > 0 ~
        "observed_above_expected",
      log10_signed_deviation < 0 ~
        "observed_below_expected",
      TRUE ~
        "observed_equals_expected"
    )
  ) %>%
  arrange(
    method,
    flowcell_scope,
    reference
  ) %>%
  as.data.frame(
    stringsAsFactors = FALSE
  )

if (nrow(analysis_data) != 64) {
  stop(
    paste(
      "Expected 64 analysis rows, observed",
      nrow(analysis_data)
    )
  )
}

if (anyNA(analysis_data$species_display)) {
  stop("Missing species labels after metadata join")
}

if (anyNA(analysis_data$expected_DNA_RA_percent)) {
  stop("Missing expected DNA RA values detected")
}

if (anyNA(analysis_data$observed_mitochondrial_RA_percent)) {
  stop("Missing observed mitochondrial RA values detected")
}

if (any(analysis_data$expected_DNA_RA_percent <= 0)) {
  stop("Zero or negative expected DNA RA values detected")
}

if (any(analysis_data$observed_mitochondrial_RA_percent <= 0)) {
  stop("Zero or negative observed mitochondrial RA values detected")
}

if (
  any(!is.finite(analysis_data$log10_expected_RA)) ||
  any(!is.finite(analysis_data$log10_observed_RA))
) {
  stop("Non-finite log10 RA values detected")
}

duplicate_count <- sum(
  duplicated(
    analysis_data[
      c("dataset", "reference")
    ]
  )
)

if (duplicate_count > 0) {
  stop("Duplicate dataset-reference records detected")
}

group_counts <- analysis_data %>%
  count(
    dataset,
    method,
    flowcell_scope,
    name = "observed_species"
  )

if (
  nrow(group_counts) != 4 ||
  any(group_counts$observed_species != 16)
) {
  stop("Expected 16 species in each of four datasets")
}

write_tsv(
  analysis_data,
  analysis_output
)

model_results <- split(
  analysis_data,
  analysis_data$dataset
) %>%
  lapply(calculate_model_outputs)

statistics <- bind_rows(
  lapply(
    model_results,
    function(x) x$statistics
  )
) %>%
  mutate(
    method = factor(
      method,
      levels = c("Metagenomic", "NAS")
    ),
    flowcell_scope = factor(
      flowcell_scope,
      levels = c("Full-flowcell", "Half-flowcell")
    )
  ) %>%
  arrange(
    method,
    flowcell_scope
  ) %>%
  mutate(
    method = as.character(method),
    flowcell_scope = as.character(flowcell_scope)
  ) %>%
  as.data.frame(
    stringsAsFactors = FALSE
  )

coefficients <- bind_rows(
  lapply(
    model_results,
    function(x) x$coefficients
  )
) %>%
  mutate(
    method = factor(
      method,
      levels = c("Metagenomic", "NAS")
    ),
    flowcell_scope = factor(
      flowcell_scope,
      levels = c("Full-flowcell", "Half-flowcell")
    )
  ) %>%
  arrange(
    method,
    flowcell_scope,
    term
  ) %>%
  mutate(
    method = as.character(method),
    flowcell_scope = as.character(flowcell_scope)
  ) %>%
  as.data.frame(
    stringsAsFactors = FALSE
  )

diagnostics <- bind_rows(
  lapply(
    model_results,
    function(x) x$diagnostics
  )
) %>%
  mutate(
    method = factor(
      method,
      levels = c("Metagenomic", "NAS")
    ),
    flowcell_scope = factor(
      flowcell_scope,
      levels = c("Full-flowcell", "Half-flowcell")
    )
  ) %>%
  arrange(
    method,
    flowcell_scope,
    reference
  ) %>%
  mutate(
    method = as.character(method),
    flowcell_scope = as.character(flowcell_scope)
  ) %>%
  as.data.frame(
    stringsAsFactors = FALSE
  )

prediction_data <- bind_rows(
  lapply(
    model_results,
    function(x) x$predictions
  )
) %>%
  as.data.frame(
    stringsAsFactors = FALSE
  )

if (nrow(statistics) != 4) {
  stop(
    paste(
      "Expected four regression-statistics rows, observed",
      nrow(statistics)
    )
  )
}

if (nrow(coefficients) != 8) {
  stop(
    paste(
      "Expected eight model-coefficient rows, observed",
      nrow(coefficients)
    )
  )
}

if (nrow(diagnostics) != 64) {
  stop(
    paste(
      "Expected 64 model-diagnostic rows, observed",
      nrow(diagnostics)
    )
  )
}

if (any(!is.finite(statistics$intercept))) {
  stop("Non-finite regression intercept detected")
}

if (any(!is.finite(statistics$slope))) {
  stop("Non-finite regression slope detected")
}

if (
  any(statistics$R_squared < 0) ||
  any(statistics$R_squared > 1)
) {
  stop("One or more R-squared values fall outside 0 to 1")
}

if (
  any(!is.finite(
    diagnostics$fitted_log10_observed_RA
  )) ||
  any(!is.finite(diagnostics$raw_residual)) ||
  any(!is.finite(diagnostics$leverage)) ||
  any(!is.finite(diagnostics$cooks_distance))
) {
  stop("Non-finite regression diagnostic detected")
}

write_tsv(
  statistics,
  statistics_output
)

write_tsv(
  coefficients,
  coefficients_output
)

write_tsv(
  diagnostics,
  diagnostics_output
)

species_levels <- species_colours$species_display

colour_map <- setNames(
  species_colours$colour,
  species_colours$species_display
)

abbreviation_map <- setNames(
  species_colours$species_abbreviation,
  species_colours$species_display
)

if (
  !setequal(
    unique(analysis_data$species_display),
    species_levels
  )
) {
  stop(
    paste0(
      "Species in RA-deviation analysis do not match ",
      "the Supplementary Figure 1 colour map"
    )
  )
}

analysis_data$species_display <- factor(
  analysis_data$species_display,
  levels = species_levels
)

panel_levels <- c(
  "A) Full-flowcell",
  "B) Half-flowcell"
)

analysis_data$panel <- factor(
  analysis_data$panel,
  levels = panel_levels
)

prediction_data$panel <- factor(
  prediction_data$panel,
  levels = panel_levels
)

observed_axis_minimum <- min(
  analysis_data$expected_DNA_RA_percent,
  analysis_data$observed_mitochondrial_RA_percent,
  na.rm = TRUE
)

observed_axis_maximum <- max(
  analysis_data$expected_DNA_RA_percent,
  analysis_data$observed_mitochondrial_RA_percent,
  na.rm = TRUE
)

axis_lower <- 0.03
axis_upper <- 50

axis_breaks_raw <- c(
  0.03,
  0.1,
  0.3,
  1,
  3,
  10,
  30
)

axis_breaks_log10 <- log10(axis_breaks_raw)

axis_labels <- c(
  "0.03",
  "0.1",
  "0.3",
  "1",
  "3",
  "10",
  "30"
)

if (
  observed_axis_minimum < axis_lower ||
  observed_axis_maximum > axis_upper
) {
  stop(
    paste0(
      "Observed values fall outside the configured axis range. ",
      "Observed range: ",
      sprintf("%.6f", observed_axis_minimum),
      " to ",
      sprintf("%.6f", observed_axis_maximum),
      "; configured range: ",
      axis_lower,
      " to ",
      axis_upper
    )
  )
}

make_plot <- function(
    plot_data,
    plot_statistics,
    plot_predictions,
    method_name
) {
  annotation_data <- plot_statistics %>%
    mutate(
      panel = factor(
        panel,
        levels = panel_levels
      ),
      x = log10(0.06),
      y = log10(18),
      label = paste0(
        "italic(R)^2 == ",
        sprintf("%.3f", R_squared)
      )
    )

  ggplot(
    plot_data,
    aes(
      x = log10_expected_RA,
      y = log10_observed_RA,
      fill = species_display
    )
  ) +
    geom_line(
      data = plot_predictions,
      aes(
        x = log10_expected_RA,
        y = fitted_log10_observed_RA,
        group = panel
      ),
      inherit.aes = FALSE,
      colour = "black",
      linewidth = 0.75
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
      drop = FALSE
    ) +
    geom_text(
      data = annotation_data,
      aes(
        x = x,
        y = y,
        label = label
      ),
      inherit.aes = FALSE,
      hjust = 0,
      vjust = 1,
      colour = "black",
      size = 3.8,
      lineheight = 1.05,
      parse = TRUE
    ) +
    scale_fill_manual(
      values = colour_map,
      breaks = species_levels,
      labels = unname(
        abbreviation_map[species_levels]
      ),
      drop = FALSE
    ) +
    scale_x_continuous(
      limits = log10(
        c(axis_lower, axis_upper)
      ),
      breaks = axis_breaks_log10,
      labels = axis_labels,
      expand = expansion(
        mult = c(0, 0)
      )
    ) +
    scale_y_continuous(
      limits = log10(
        c(axis_lower, axis_upper)
      ),
      breaks = axis_breaks_log10,
      labels = axis_labels,
      expand = expansion(
        mult = c(0, 0)
      )
    ) +
    coord_equal(
      clip = "off"
    ) +
    labs(
      title = method_name,
      x = "Expected DNA RA (%)",
      y = "Observed mitochondrial RA (%)",
      fill = NULL
    ) +
    guides(
      fill = guide_legend(
        nrow = 4,
        byrow = FALSE,
        override.aes = list(
          size = 3,
          alpha = 1
        )
      )
    ) +
    theme_bw(
      base_size = 13,
      base_family = "sans"
    ) +
    theme(
      plot.title = element_text(
        face = "bold",
        colour = "black",
        size = 13,
        hjust = 0.5,
        margin = margin(
          t = 2,
          r = 0,
          b = 10,
          l = 0
        )
      ),
      strip.background = element_blank(),
      strip.text = element_text(
        face = "bold",
        colour = "black",
        size = 13,
        margin = margin(
          t = 2,
          r = 0,
          b = 10,
          l = 0
        )
      ),
      panel.border = element_rect(
        colour = "#666666",
        fill = NA,
        linewidth = 0.6
      ),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.spacing = grid::unit(
        0.85,
        "cm"
      ),
      axis.title = element_text(
        colour = "black",
        size = 12
      ),
      axis.title.x = element_text(
        margin = margin(
          t = 8,
          r = 0,
          b = 0,
          l = 0
        )
      ),
      axis.title.y = element_text(
        margin = margin(
          t = 0,
          r = 8,
          b = 0,
          l = 0
        )
      ),
      axis.text = element_text(
        colour = "black",
        size = 9.5
      ),
      axis.ticks = element_line(
        colour = "black",
        linewidth = 0.45
      ),
      legend.position = "bottom",
      legend.justification = "center",
      legend.key = element_blank(),
      legend.text = element_text(
        face = "italic",
        colour = "black",
        size = 9
      ),
      legend.box.spacing = grid::unit(
        0.08,
        "cm"
      ),
      legend.margin = margin(
        t = -2,
        r = 0,
        b = 0,
        l = 0
      ),
      plot.margin = margin(
        t = 16,
        r = 10,
        b = 2,
        l = 10
      )
    )
}

meta_data <- analysis_data %>%
  filter(method == "Metagenomic")

meta_statistics <- statistics %>%
  filter(method == "Metagenomic")

meta_predictions <- prediction_data %>%
  filter(method == "Metagenomic")

nas_data <- analysis_data %>%
  filter(method == "NAS")

nas_statistics <- statistics %>%
  filter(method == "NAS")

nas_predictions <- prediction_data %>%
  filter(method == "NAS")

if (
  nrow(meta_data) != 32 ||
  nrow(meta_statistics) != 2
) {
  stop(
    "Unexpected Metagenomic data or statistics row count"
  )
}

if (
  nrow(nas_data) != 32 ||
  nrow(nas_statistics) != 2
) {
  stop(
    "Unexpected NAS data or statistics row count"
  )
}

meta_plot <- make_plot(
  meta_data,
  meta_statistics,
  meta_predictions,
  "Metagenomic"
)

nas_plot <- make_plot(
  nas_data,
  nas_statistics,
  nas_predictions,
  "NAS"
)

validation <- data.frame(
  validation = c(
    "selected_rows",
    "rows_per_dataset",
    "positive_expected_RA",
    "positive_observed_RA",
    "finite_log10_values",
    "duplicate_dataset_reference_records",
    "regression_statistics_rows",
    "model_coefficient_rows",
    "diagnostic_rows",
    "finite_regression_coefficients",
    "finite_regression_diagnostics",
    "R_squared_range",
    "confidence_interval_order",
    "species_colour_map_match",
    "log_axis_range_contains_all_values"
  ),
  passed = c(
    nrow(analysis_data) == 64,
    all(group_counts$observed_species == 16),
    all(analysis_data$expected_DNA_RA_percent > 0),
    all(
      analysis_data$observed_mitochondrial_RA_percent > 0
    ),
    all(is.finite(analysis_data$log10_expected_RA)) &&
      all(is.finite(analysis_data$log10_observed_RA)),
    duplicate_count == 0,
    nrow(statistics) == 4,
    nrow(coefficients) == 8,
    nrow(diagnostics) == 64,
    all(is.finite(statistics$intercept)) &&
      all(is.finite(statistics$slope)),
    all(
      is.finite(
        diagnostics$fitted_log10_observed_RA
      )
    ) &&
      all(is.finite(diagnostics$raw_residual)) &&
      all(is.finite(diagnostics$leverage)) &&
      all(is.finite(diagnostics$cooks_distance)),
    all(
      statistics$R_squared >= 0 &
        statistics$R_squared <= 1
    ),
    all(
      statistics$intercept_CI_lower <=
        statistics$intercept
    ) &&
      all(
        statistics$intercept <=
          statistics$intercept_CI_upper
      ) &&
      all(
        statistics$slope_CI_lower <=
          statistics$slope
      ) &&
      all(
        statistics$slope <=
          statistics$slope_CI_upper
      ),
    setequal(
      as.character(
        unique(
          analysis_data$species_display
        )
      ),
      species_levels
    ),
    observed_axis_minimum >= axis_lower &&
      observed_axis_maximum <= axis_upper
  ),
  stringsAsFactors = FALSE
)

write_tsv(
  validation,
  validation_output
)

if (any(!validation$passed)) {
  stop(
    paste(
      "Validation failure:",
      paste(
        validation$validation[
          !validation$passed
        ],
        collapse = ", "
      )
    )
  )
}

figure_export_status <- bind_rows(
  export_plot(
    meta_plot,
    meta_png,
    "png",
    dpi = 300
  ),
  export_plot(
    meta_plot,
    meta_pdf,
    "pdf"
  ),
  export_plot(
    meta_plot,
    meta_svg,
    "svg"
  ),
  export_plot(
    nas_plot,
    nas_png,
    "png",
    dpi = 300
  ),
  export_plot(
    nas_plot,
    nas_pdf,
    "pdf"
  ),
  export_plot(
    nas_plot,
    nas_svg,
    "svg"
  )
)

write_tsv(
  figure_export_status,
  figure_export_output
)

if (any(!figure_export_status$success)) {
  print(
    figure_export_status %>%
      filter(!success)
  )

  stop("One or more figure exports failed")
}

finished_time <- Sys.time()

provenance <- data.frame(
  field = c(
    "analysis_name",
    "analysis_status",
    "started",
    "finished",
    "canonical_mapping_summary",
    "reference_manifest",
    "species_colour_map",
    "observed_RA_definition",
    "deviation_definition",
    "regression_model",
    "regression_values",
    "figure_annotation",
    "confidence_band",
    "identity_line",
    "zero_handling",
    "species_per_dataset",
    "cooks_distance_threshold",
    "influential_observation_handling",
    "axis_transformation",
    "shared_axis_limits"
  ),
  value = c(
    "ra_deviation_linear_regression",
    "manuscript analysis",
    format(
      started_time,
      "%Y-%m-%dT%H:%M:%S%z"
    ),
    format(
      finished_time,
      "%Y-%m-%dT%H:%M:%S%z"
    ),
    mapping_file,
    reference_file,
    colour_map_file,
    paste0(
      "Species CIGAR M, = and X aligned query bases ",
      "divided by total corresponding aligned query bases ",
      "across all 16 mitochondrial references"
    ),
    paste0(
      "log10 observed mitochondrial RA minus ",
      "log10 expected DNA RA"
    ),
    paste0(
      "lm(log10 observed mitochondrial RA percent ~ ",
      "log10 expected DNA RA percent), fitted separately ",
      "for each dataset"
    ),
    "log10-transformed RA percent",
    "R-squared only",
    "not displayed",
    "slope 1 and intercept 0 in log10 space",
    "none required",
    "16",
    "4/n within each dataset",
    paste0(
      "Influential observations flagged in diagnostics; ",
      "no observations excluded"
    ),
    "log10 values plotted with original-percent labels",
    paste(
      axis_lower,
      axis_upper,
      sep = " to "
    )
  ),
  stringsAsFactors = FALSE
)

write_tsv(
  provenance,
  provenance_file
)

log_lines <- c(
  paste(
    "Started:",
    format(
      started_time,
      "%Y-%m-%dT%H:%M:%S%z"
    )
  ),
  paste(
    "Finished:",
    format(
      finished_time,
      "%Y-%m-%dT%H:%M:%S%z"
    )
  ),
  paste(
    "Canonical mapping summary:",
    mapping_file
  ),
  paste(
    "Reference manifest:",
    reference_file
  ),
  paste(
    "Species colour map:",
    colour_map_file
  ),
  paste0(
    "Observed RA definition: species CIGAR M + = + X ",
    "aligned query bases divided by total corresponding ",
    "aligned query bases across all 16 mitochondrial references"
  ),
  paste0(
    "Deviation definition: log10 observed mitochondrial RA ",
    "minus log10 expected DNA RA"
  ),
  paste0(
    "Regression model: lm(log10 observed mitochondrial RA ",
    "percent ~ log10 expected DNA RA percent), fitted ",
    "separately for each dataset"
  ),
  "Figure annotation: R-squared only",
  "Confidence band displayed: FALSE",
  paste0(
    "Full model coefficients, confidence intervals and ",
    "diagnostics retained in tabular outputs"
  ),
  "Zero handling: none required",
  "Species included per dataset: 16",
  "Cook's distance threshold: 4/n within each dataset",
  paste0(
    "Influential observations: flagged in diagnostics; ",
    "no observations excluded"
  ),
  paste0(
    "Identity line: slope 1 and intercept 0 in log10 space"
  ),
  paste0(
    "Axis display: log10 values with original-percent labels"
  ),
  paste(
    "Shared axis limits:",
    axis_lower,
    "to",
    axis_upper
  ),
  "",
  "Regression statistics:",
  capture.output(
    print(statistics)
  ),
  "",
  "Model coefficients:",
  capture.output(
    print(coefficients)
  ),
  "",
  "Influential observations:",
  capture.output(
    print(
      diagnostics %>%
        filter(influential_flag) %>%
        select(
          dataset,
          reference,
          species,
          cooks_distance,
          cooks_distance_threshold,
          standardised_residual,
          leverage
        )
    )
  ),
  "",
  "Validation results:",
  capture.output(
    print(validation)
  ),
  "",
  "Figure export status:",
  capture.output(
    print(figure_export_status)
  )
)

writeLines(
  log_lines,
  log_file,
  useBytes = TRUE
)

writeLines(
  capture.output(
    sessionInfo()
  ),
  session_file,
  useBytes = TRUE
)

console_lines <- c(
  capture.output(
    print(statistics)
  ),
  "",
  capture.output(
    print(coefficients)
  ),
  "",
  capture.output(
    print(
      diagnostics %>%
        filter(influential_flag) %>%
        select(
          dataset,
          reference,
          species,
          cooks_distance,
          cooks_distance_threshold,
          standardised_residual,
          leverage
        )
    )
  ),
  "",
  capture.output(
    print(validation)
  ),
  "",
  capture.output(
    print(figure_export_status)
  ),
  "",
  paste(
    "Wrote analysis outputs to:",
    output_dir
  ),
  paste(
    "Wrote figures to:",
    figure_dir
  ),
  paste(
    "Wrote logs to:",
    log_dir
  ),
  paste0(
    "PASS: RA-deviation linear-regression analysis ",
    "and exports completed successfully"
  )
)

writeLines(
  console_lines,
  console_file,
  useBytes = TRUE
)

print(statistics)
print(coefficients)
print(
  diagnostics %>%
    filter(influential_flag) %>%
    select(
      dataset,
      reference,
      species,
      cooks_distance,
      cooks_distance_threshold,
      standardised_residual,
      leverage
    )
)
print(validation)
print(figure_export_status)

message(
  "Wrote analysis outputs to: ",
  output_dir
)

message(
  "Wrote figures to: ",
  figure_dir
)

message(
  "Wrote logs to: ",
  log_dir
)

message(
  paste0(
    "PASS: RA-deviation linear-regression analysis ",
    "and exports completed successfully"
  )
)
