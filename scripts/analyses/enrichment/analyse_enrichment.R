#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(ragg)
  library(svglite)
  library(grid)
})


# ============================================================
# Enrichment by yield and composition
#
# Figure outputs:
#
#   1. Whole-flowcell manuscript figure.
#   2. Half-flowcell supplementary figure.
#
# Each figure contains:
#
#   A) Enrichment by yield.
#   B) Enrichment by composition.
#
# Statistical definitions
# -----------------------
#
# Enrichment by yield:
#
#   species-specific accepted NAS yield per hour
#   ------------------------------------------------
#   species-specific Metagenomic yield per hour
#
# Both sequencing runs covered the same 24-hour duration.
# Therefore, the common time denominator cancels:
#
#   (NAS aligned bases / 24 h)
#   --------------------------
#   (Meta aligned bases / 24 h)
#
# equals:
#
#   NAS aligned bases
#   -----------------
#   Meta aligned bases
#
# Enrichment by composition:
#
#   species mitochondrial relative-abundance fraction in NAS
#   ---------------------------------------------------------
#   species mitochondrial relative-abundance fraction in Meta
#
# Both ratios are expressed as:
#
#   log2(
#     accepted NAS /
#     Metagenomic
#   )
#
# Interpretation:
#
#   negative log2 fold change:
#     greater in Metagenomic sequencing;
#
#   zero:
#     equal between methods;
#
#   positive log2 fold change:
#     greater in accepted NAS.
#
# Examples:
#
#   -2 = Metagenomic is fourfold greater;
#   -1 = Metagenomic is twofold greater;
#    0 = equal;
#   +1 = accepted NAS is twofold greater;
#   +2 = accepted NAS is fourfold greater.
#
# The shared x-axis range is calculated from all yield and
# composition values in both whole-flowcell and half-flowcell
# comparisons.
#
# The smallest symmetric integer limit containing every value
# is used in both figures.
#
# Provenance outputs retained:
#
#   complete enrichment results;
#   descriptive enrichment summary;
#   plotting data;
#   species-label map;
#   calculation validation;
#   capture-class validation;
#   axis summary;
#   figure export status;
#   run log;
#   R session information.
#
# Outputs not generated:
#
#   combined whole- and half-flowcell figure;
#   directional reciprocal-fold figure;
#   signed raw-fold figure.
# ============================================================


# ============================================================
# Fixed analysis settings
# ============================================================

sequencing_duration_hours <- 24

calculation_tolerance <- 1e-10


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
  project_dir,
  "analysis",
  "enrichment"
)

mapping_file <- file.path(
  canonical_data_dir,
  "outputs",
  "canonical_mito_mapping_summary.tsv"
)

reference_file <- file.path(
  canonical_data_dir,
  "inputs",
  "reference_manifest.tsv"
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
# Provenance output files
# ============================================================

results_file <- file.path(
  output_dir,
  "enrichment_results.tsv"
)

summary_file <- file.path(
  output_dir,
  "enrichment_summary.tsv"
)

plot_data_file <- file.path(
  output_dir,
  "enrichment_log2_plot_data_publication.tsv"
)

species_label_file <- file.path(
  output_dir,
  "enrichment_species_label_map.tsv"
)

calculation_validation_file <- file.path(
  output_dir,
  "enrichment_log2_calculation_validation.tsv"
)

class_validation_file <- file.path(
  output_dir,
  "enrichment_capture_class_validation.tsv"
)

axis_summary_file <- file.path(
  output_dir,
  "enrichment_log2_axis_summary_publication.tsv"
)

export_status_file <- file.path(
  output_dir,
  "enrichment_figure_export_status.tsv"
)


# ============================================================
# Whole-flowcell manuscript figure outputs
# ============================================================

whole_figure_png <- file.path(
  figure_dir,
  paste0(
    "enrichment_yield_and_composition_",
    "whole_flowcell_publication.png"
  )
)

whole_figure_pdf <- file.path(
  figure_dir,
  paste0(
    "enrichment_yield_and_composition_",
    "whole_flowcell_publication.pdf"
  )
)

whole_figure_svg <- file.path(
  figure_dir,
  paste0(
    "enrichment_yield_and_composition_",
    "whole_flowcell_publication.svg"
  )
)


# ============================================================
# Half-flowcell supplementary figure outputs
# ============================================================

half_figure_png <- file.path(
  figure_dir,
  paste0(
    "enrichment_yield_and_composition_",
    "half_flowcell_supplementary.png"
  )
)

half_figure_pdf <- file.path(
  figure_dir,
  paste0(
    "enrichment_yield_and_composition_",
    "half_flowcell_supplementary.pdf"
  )
)

half_figure_svg <- file.path(
  figure_dir,
  paste0(
    "enrichment_yield_and_composition_",
    "half_flowcell_supplementary.svg"
  )
)


# ============================================================
# Log and session-information files
# ============================================================

log_file <- file.path(
  log_dir,
  "enrichment_publication_run.log"
)

session_file <- file.path(
  log_dir,
  "R_sessionInfo_enrichment_publication.txt"
)


# ============================================================
# Validate required input files
# ============================================================

required_files <- c(
  mapping_file,
  reference_file
)

missing_files <- required_files[
  !file.exists(
    required_files
  )
]

if (length(missing_files) > 0) {
  stop(
    paste(
      "Missing required enrichment input files:",
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

mapping <- fread(
  mapping_file,
  sep = "\t",
  header = TRUE,
  integer64 = "double"
) %>%
  as_tibble()

references <- fread(
  reference_file,
  sep = "\t",
  header = TRUE,
  integer64 = "double"
) %>%
  as_tibble()


# ============================================================
# Required input columns
# ============================================================

required_mapping_columns <- c(
  "dataset",
  "reference",
  "cigar_aligned_query_bases",
  "mitochondrial_RA_fraction"
)

required_reference_columns <- c(
  "reference",
  "species",
  "target_status",
  "expected_DNA_RA_percent"
)


# ============================================================
# Column-validation helper
# ============================================================

check_columns <- function(
    data,
    required_columns,
    data_name
) {
  missing_columns <- setdiff(
    required_columns,
    names(
      data
    )
  )

  if (length(missing_columns) > 0) {
    stop(
      paste0(
        "Missing required columns in ",
        data_name,
        ": ",
        paste(
          missing_columns,
          collapse = ", "
        )
      )
    )
  }
}


check_columns(
  mapping,
  required_mapping_columns,
  "canonical mitochondrial mapping summary"
)

check_columns(
  references,
  required_reference_columns,
  "reference manifest"
)


# ============================================================
# Validate reference metadata
# ============================================================

if (
  any(
    duplicated(
      references$reference
    )
  )
) {
  stop(
    "Duplicate references detected in reference manifest"
  )
}

if (
  anyNA(
    references$reference
  ) ||
  anyNA(
    references$species
  ) ||
  anyNA(
    references$target_status
  )
) {
  stop(
    "Missing required reference metadata detected"
  )
}


# ============================================================
# Validate required analytical datasets
# ============================================================

required_datasets <- c(
  "full_meta",
  "full_nas_accepted",
  "half_meta",
  "half_nas_accepted"
)

missing_datasets <- setdiff(
  required_datasets,
  unique(
    mapping$dataset
  )
)

if (length(missing_datasets) > 0) {
  stop(
    paste(
      "Missing required mapping datasets:",
      paste(
        missing_datasets,
        collapse = ", "
      )
    )
  )
}


# ============================================================
# Construct one enrichment comparison
# ============================================================

build_comparison <- function(
    mapping_data,
    reference_data,
    meta_dataset,
    nas_dataset,
    comparison_label,
    duration_hours
) {
  meta_data <- mapping_data %>%
    filter(
      dataset ==
        meta_dataset
    ) %>%
    transmute(
      reference,

      meta_cigar_aligned_query_bases =
        as.numeric(
          cigar_aligned_query_bases
        ),

      meta_mitochondrial_RA_fraction =
        as.numeric(
          mitochondrial_RA_fraction
        )
    )

  nas_data <- mapping_data %>%
    filter(
      dataset ==
        nas_dataset
    ) %>%
    transmute(
      reference,

      nas_cigar_aligned_query_bases =
        as.numeric(
          cigar_aligned_query_bases
        ),

      nas_mitochondrial_RA_fraction =
        as.numeric(
          mitochondrial_RA_fraction
        )
    )

  if (nrow(meta_data) != 16) {
    stop(
      paste(
        "Expected 16 Metagenomic rows for",
        meta_dataset,
        "but observed",
        nrow(
          meta_data
        )
      )
    )
  }

  if (nrow(nas_data) != 16) {
    stop(
      paste(
        "Expected 16 NAS rows for",
        nas_dataset,
        "but observed",
        nrow(
          nas_data
        )
      )
    )
  }

  if (
    any(
      duplicated(
        meta_data$reference
      )
    )
  ) {
    stop(
      paste(
        "Duplicate Metagenomic references in",
        comparison_label
      )
    )
  }

  if (
    any(
      duplicated(
        nas_data$reference
      )
    )
  ) {
    stop(
      paste(
        "Duplicate NAS references in",
        comparison_label
      )
    )
  }

  comparison_data <- meta_data %>%
    inner_join(
      nas_data,
      by = "reference"
    ) %>%
    left_join(
      reference_data %>%
        select(
          reference,
          species,
          target_status,
          expected_DNA_RA_percent
        ),
      by = "reference"
    ) %>%
    mutate(
      comparison =
        comparison_label,

      meta_dataset =
        meta_dataset,

      nas_dataset =
        nas_dataset,

      nas_duration_hours =
        duration_hours,

      meta_duration_hours =
        duration_hours
    )

  if (nrow(comparison_data) != 16) {
    stop(
      paste(
        "Expected 16 joined rows for",
        comparison_label,
        "but observed",
        nrow(
          comparison_data
        )
      )
    )
  }

  if (
    anyNA(
      comparison_data$species
    ) ||
    anyNA(
      comparison_data$target_status
    )
  ) {
    stop(
      paste(
        "Missing reference metadata in",
        comparison_label
      )
    )
  }

  input_values <- comparison_data %>%
    select(
      meta_cigar_aligned_query_bases,
      nas_cigar_aligned_query_bases,
      meta_mitochondrial_RA_fraction,
      nas_mitochondrial_RA_fraction
    )

  if (
    anyNA(
      input_values
    )
  ) {
    stop(
      paste(
        "Missing enrichment input values in",
        comparison_label
      )
    )
  }

  if (
    any(
      input_values <= 0
    )
  ) {
    stop(
      paste(
        "Zero or negative enrichment input values in",
        comparison_label
      )
    )
  }

  comparison_data %>%
    mutate(
      nas_target_yield_per_hour =
        nas_cigar_aligned_query_bases /
        nas_duration_hours,

      meta_target_yield_per_hour =
        meta_cigar_aligned_query_bases /
        meta_duration_hours,

      enrichment_by_yield =
        nas_target_yield_per_hour /
        meta_target_yield_per_hour,

      enrichment_by_yield_simplified =
        nas_cigar_aligned_query_bases /
        meta_cigar_aligned_query_bases,

      enrichment_by_composition =
        nas_mitochondrial_RA_fraction /
        meta_mitochondrial_RA_fraction,

      log2_enrichment_by_yield =
        log2(
          enrichment_by_yield
        ),

      log2_enrichment_by_composition =
        log2(
          enrichment_by_composition
        ),

      enrichment_direction_by_yield =
        case_when(
          log2_enrichment_by_yield > 0 ~
            "Greater in accepted NAS",

          log2_enrichment_by_yield < 0 ~
            "Greater in Metagenomic",

          TRUE ~
            "Equal"
        ),

      enrichment_direction_by_composition =
        case_when(
          log2_enrichment_by_composition > 0 ~
            "Greater in accepted NAS",

          log2_enrichment_by_composition < 0 ~
            "Greater in Metagenomic",

          TRUE ~
            "Equal"
        )
    )
}


# ============================================================
# Build requested comparisons
# ============================================================

whole_flowcell <- build_comparison(
  mapping_data = mapping,
  reference_data = references,
  meta_dataset = "full_meta",
  nas_dataset = "full_nas_accepted",
  comparison_label = "Whole-flowcell",
  duration_hours = sequencing_duration_hours
)

half_flowcell <- build_comparison(
  mapping_data = mapping,
  reference_data = references,
  meta_dataset = "half_meta",
  nas_dataset = "half_nas_accepted",
  comparison_label = "Half-flowcell",
  duration_hours = sequencing_duration_hours
)


# ============================================================
# Combine comparison results
# ============================================================

comparison_levels <- c(
  "Whole-flowcell",
  "Half-flowcell"
)

capture_class_levels <- c(
  "NAS target",
  "Related non-target",
  "Other non-target"
)

results <- bind_rows(
  whole_flowcell,
  half_flowcell
) %>%
  mutate(
    comparison = factor(
      comparison,
      levels =
        comparison_levels
    ),

    capture_class = recode(
      target_status,

      "NAS_target" =
        "NAS target",

      "related_non_target" =
        "Related non-target",

      "other_non_target" =
        "Other non-target",

      .default =
        NA_character_
    ),

    capture_class = factor(
      capture_class,
      levels =
        capture_class_levels
    )
  ) %>%
  arrange(
    comparison,
    desc(
      log2_enrichment_by_yield
    )
  )


# ============================================================
# Validate complete results
# ============================================================

if (nrow(results) != 32) {
  stop(
    paste(
      "Expected 32 enrichment rows, observed",
      nrow(
        results
      )
    )
  )
}

if (
  anyNA(
    results$capture_class
  )
) {
  stop(
    "Unexpected or missing capture classes detected"
  )
}

duplicate_keys <- duplicated(
  results %>%
    select(
      comparison,
      reference
    )
)

if (
  any(
    duplicate_keys
  )
) {
  stop(
    "Duplicate comparison-reference records detected"
  )
}

if (
  anyNA(
    results$log2_enrichment_by_yield
  ) ||
  anyNA(
    results$log2_enrichment_by_composition
  ) ||
  any(
    !is.finite(
      results$log2_enrichment_by_yield
    )
  ) ||
  any(
    !is.finite(
      results$log2_enrichment_by_composition
    )
  )
) {
  stop(
    "Invalid log2 enrichment values detected"
  )
}


# ============================================================
# Calculation validation
#
# Validates:
#
#   1. NAS and Metagenomic durations are identical.
#
#   2. The explicit per-hour yield ratio equals the simplified
#      aligned-base ratio because both durations are 24 hours.
#
#   3. Exponentiating the log2 value reconstructs the original
#      enrichment ratio.
#
#   4. Fold-change signs match the source ratios.
# ============================================================

calculation_validation <- results %>%
  transmute(
    comparison,
    reference,
    species,

    nas_duration_hours,
    meta_duration_hours,

    duration_match_passed =
      abs(
        nas_duration_hours -
          meta_duration_hours
      ) <= calculation_tolerance,

    yield_ratio_per_hour =
      enrichment_by_yield,

    yield_ratio_simplified =
      enrichment_by_yield_simplified,

    yield_time_cancellation_error =
      abs(
        enrichment_by_yield -
          enrichment_by_yield_simplified
      ),

    yield_time_cancellation_passed =
      yield_time_cancellation_error <=
      calculation_tolerance,

    yield_ratio_reconstructed =
      2 ^
      log2_enrichment_by_yield,

    yield_log2_absolute_error =
      abs(
        enrichment_by_yield -
          yield_ratio_reconstructed
      ),

    yield_log2_calculation_passed =
      yield_log2_absolute_error <=
      calculation_tolerance,

    yield_direction_passed =
      case_when(
        enrichment_by_yield > 1 ~
          log2_enrichment_by_yield > 0,

        enrichment_by_yield < 1 ~
          log2_enrichment_by_yield < 0,

        TRUE ~
          abs(
            log2_enrichment_by_yield
          ) <= calculation_tolerance
      ),

    composition_ratio =
      enrichment_by_composition,

    composition_ratio_reconstructed =
      2 ^
      log2_enrichment_by_composition,

    composition_log2_absolute_error =
      abs(
        enrichment_by_composition -
          composition_ratio_reconstructed
      ),

    composition_log2_calculation_passed =
      composition_log2_absolute_error <=
      calculation_tolerance,

    composition_direction_passed =
      case_when(
        enrichment_by_composition > 1 ~
          log2_enrichment_by_composition > 0,

        enrichment_by_composition < 1 ~
          log2_enrichment_by_composition < 0,

        TRUE ~
          abs(
            log2_enrichment_by_composition
          ) <= calculation_tolerance
      )
  )


if (
  any(
    calculation_validation$
    duration_match_passed != TRUE
  )
) {
  stop(
    "NAS and Metagenomic sequencing durations do not match"
  )
}

if (
  any(
    calculation_validation$
    yield_time_cancellation_passed != TRUE
  )
) {
  stop(
    paste0(
      "The explicit per-hour yield ratio does not equal ",
      "the simplified aligned-base ratio"
    )
  )
}

if (
  any(
    calculation_validation$
    yield_log2_calculation_passed != TRUE
  ) ||
  any(
    calculation_validation$
    composition_log2_calculation_passed != TRUE
  )
) {
  stop(
    "One or more log2 fold-change calculations failed validation"
  )
}

if (
  any(
    calculation_validation$
    yield_direction_passed != TRUE
  ) ||
  any(
    calculation_validation$
    composition_direction_passed != TRUE
  )
) {
  stop(
    "One or more log2 fold-change directions failed validation"
  )
}


fwrite(
  calculation_validation,
  calculation_validation_file,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)


# ============================================================
# Capture-class validation
# ============================================================

expected_class_counts <- tidyr::expand_grid(
  comparison =
    comparison_levels,

  capture_class =
    capture_class_levels
) %>%
  mutate(
    expected_n = case_when(
      capture_class ==
        "NAS target" ~ 3L,

      capture_class ==
        "Related non-target" ~ 1L,

      capture_class ==
        "Other non-target" ~ 12L,

      TRUE ~
        NA_integer_
    )
  )

observed_class_counts <- results %>%
  transmute(
    comparison =
      as.character(
        comparison
      ),

    capture_class =
      as.character(
        capture_class
      )
  ) %>%
  count(
    comparison,
    capture_class,
    name =
      "observed_n"
  )

class_count_validation <- expected_class_counts %>%
  left_join(
    observed_class_counts,
    by = c(
      "comparison",
      "capture_class"
    )
  ) %>%
  mutate(
    observed_n = replace_na(
      observed_n,
      0L
    ),

    passed =
      expected_n ==
      observed_n
  )

if (
  any(
    class_count_validation$passed != TRUE
  )
) {
  stop(
    paste0(
      "Capture-class counts differ from the expected ",
      "3 NAS targets, 1 related non-target and 12 other ",
      "non-targets per comparison"
    )
  )
}


fwrite(
  class_count_validation,
  class_validation_file,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)


# ============================================================
# Write complete enrichment results
# ============================================================

fwrite(
  results,
  results_file,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)


# ============================================================
# Geometric-mean helper
# ============================================================

geometric_mean <- function(x) {
  if (
    anyNA(
      x
    ) ||
    any(
      x <= 0
    )
  ) {
    stop(
      "Geometric mean requires positive, non-missing values"
    )
  }

  exp(
    mean(
      log(
        x
      )
    )
  )
}


# ============================================================
# Descriptive summary by capture class
# ============================================================

summary_by_class <- results %>%
  group_by(
    comparison,
    capture_class
  ) %>%
  summarise(
    n = n(),

    geometric_mean_yield_enrichment =
      geometric_mean(
        enrichment_by_yield
      ),

    median_yield_enrichment =
      median(
        enrichment_by_yield
      ),

    minimum_yield_enrichment =
      min(
        enrichment_by_yield
      ),

    maximum_yield_enrichment =
      max(
        enrichment_by_yield
      ),

    geometric_mean_composition_enrichment =
      geometric_mean(
        enrichment_by_composition
      ),

    median_composition_enrichment =
      median(
        enrichment_by_composition
      ),

    minimum_composition_enrichment =
      min(
        enrichment_by_composition
      ),

    maximum_composition_enrichment =
      max(
        enrichment_by_composition
      ),

    mean_log2_yield_enrichment =
      mean(
        log2_enrichment_by_yield
      ),

    median_log2_yield_enrichment =
      median(
        log2_enrichment_by_yield
      ),

    minimum_log2_yield_enrichment =
      min(
        log2_enrichment_by_yield
      ),

    maximum_log2_yield_enrichment =
      max(
        log2_enrichment_by_yield
      ),

    mean_log2_composition_enrichment =
      mean(
        log2_enrichment_by_composition
      ),

    median_log2_composition_enrichment =
      median(
        log2_enrichment_by_composition
      ),

    minimum_log2_composition_enrichment =
      min(
        log2_enrichment_by_composition
      ),

    maximum_log2_composition_enrichment =
      max(
        log2_enrichment_by_composition
      ),

    .groups =
      "drop"
  ) %>%
  mutate(
    summary_scope =
      "capture_class"
  )


# ============================================================
# Overall descriptive summary
# ============================================================

summary_overall <- results %>%
  group_by(
    comparison
  ) %>%
  summarise(
    n = n(),

    geometric_mean_yield_enrichment =
      geometric_mean(
        enrichment_by_yield
      ),

    median_yield_enrichment =
      median(
        enrichment_by_yield
      ),

    minimum_yield_enrichment =
      min(
        enrichment_by_yield
      ),

    maximum_yield_enrichment =
      max(
        enrichment_by_yield
      ),

    geometric_mean_composition_enrichment =
      geometric_mean(
        enrichment_by_composition
      ),

    median_composition_enrichment =
      median(
        enrichment_by_composition
      ),

    minimum_composition_enrichment =
      min(
        enrichment_by_composition
      ),

    maximum_composition_enrichment =
      max(
        enrichment_by_composition
      ),

    mean_log2_yield_enrichment =
      mean(
        log2_enrichment_by_yield
      ),

    median_log2_yield_enrichment =
      median(
        log2_enrichment_by_yield
      ),

    minimum_log2_yield_enrichment =
      min(
        log2_enrichment_by_yield
      ),

    maximum_log2_yield_enrichment =
      max(
        log2_enrichment_by_yield
      ),

    mean_log2_composition_enrichment =
      mean(
        log2_enrichment_by_composition
      ),

    median_log2_composition_enrichment =
      median(
        log2_enrichment_by_composition
      ),

    minimum_log2_composition_enrichment =
      min(
        log2_enrichment_by_composition
      ),

    maximum_log2_composition_enrichment =
      max(
        log2_enrichment_by_composition
      ),

    .groups =
      "drop"
  ) %>%
  mutate(
    capture_class =
      "All species",

    summary_scope =
      "overall"
  )


# ============================================================
# Combined descriptive summary
# ============================================================

summary_table <- bind_rows(
  summary_by_class,
  summary_overall
) %>%
  select(
    comparison,
    summary_scope,
    capture_class,
    everything()
  ) %>%
  arrange(
    comparison,
    summary_scope,
    capture_class
  )


fwrite(
  summary_table,
  summary_file,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)


# ============================================================
# Species-label helper
#
# Genus is abbreviated to its initial.
#
# Only genus and species epithet are displayed.
#
# Additional subspecies or strain terms are omitted from the
# plotting label.
# ============================================================

abbreviate_species <- function(
    species_name
) {
  name_parts <- strsplit(
    as.character(
      species_name
    ),
    split = " ",
    fixed = TRUE
  )[[1]]

  if (length(name_parts) < 2) {
    return(
      as.character(
        species_name
      )
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


species_label_map <- results %>%
  distinct(
    reference,
    species
  ) %>%
  mutate(
    species_display = vapply(
      species,
      abbreviate_species,
      character(1)
    )
  )


if (
  nrow(
    species_label_map
  ) != 16
) {
  stop(
    paste(
      "Expected 16 species labels, observed",
      nrow(
        species_label_map
      )
    )
  )
}

if (
  any(
    duplicated(
      species_label_map$species_display
    )
  )
) {
  stop(
    paste0(
      "Duplicate abbreviated species labels detected. ",
      "Species labels must be unique."
    )
  )
}


fwrite(
  species_label_map,
  species_label_file,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)


# ============================================================
# Shared species order
#
# Species are ordered using whole-flowcell enrichment by yield.
#
# Greatest positive accepted-NAS enrichment appears at the top.
#
# Greatest negative Metagenomic enrichment appears at the
# bottom.
# ============================================================

species_order <- results %>%
  filter(
    comparison ==
      "Whole-flowcell"
  ) %>%
  arrange(
    log2_enrichment_by_yield
  ) %>%
  pull(
    species
  ) %>%
  as.character()


if (
  length(
    unique(
      species_order
    )
  ) != 16
) {
  stop(
    "Expected 16 unique species in plotting order"
  )
}


species_display_order <- vapply(
  species_order,
  abbreviate_species,
  character(1)
)


# ============================================================
# Italic species labels
# ============================================================

species_expression_labels <- lapply(
  species_display_order,

  function(
    display_name
  ) {
    bquote(
      italic(
        .(display_name)
      )
    )
  }
)

names(
  species_expression_labels
) <- species_display_order


species_label_function <- function(x) {
  as.expression(
    unname(
      species_expression_labels[
        as.character(
          x
        )
      ]
    )
  )
}


# ============================================================
# Long-format plotting data
# ============================================================

plot_data <- results %>%
  left_join(
    species_label_map,
    by = c(
      "reference",
      "species"
    )
  ) %>%
  select(
    reference,
    species,
    species_display,
    capture_class,
    comparison,
    enrichment_by_yield,
    enrichment_by_composition,
    log2_enrichment_by_yield,
    log2_enrichment_by_composition
  ) %>%
  pivot_longer(
    cols = c(
      log2_enrichment_by_yield,
      log2_enrichment_by_composition
    ),

    names_to =
      "metric",

    values_to =
      "log2_fold_change"
  ) %>%
  mutate(
    metric = recode(
      metric,

      "log2_enrichment_by_yield" =
        "Enrichment by yield",

      "log2_enrichment_by_composition" =
        "Enrichment by composition"
    ),

    metric = factor(
      metric,
      levels = c(
        "Enrichment by yield",
        "Enrichment by composition"
      )
    ),

    species_display = factor(
      species_display,
      levels =
        species_display_order
    ),

    direction = case_when(
      log2_fold_change < 0 ~
        "Greater in Metagenomic",

      log2_fold_change > 0 ~
        "Greater in accepted NAS",

      TRUE ~
        "Equal"
    )
  ) %>%
  arrange(
    comparison,
    metric,
    species_display
  )


if (
  nrow(
    plot_data
  ) != 64
) {
  stop(
    paste(
      "Expected 64 plotting rows, observed",
      nrow(
        plot_data
      )
    )
  )
}

if (
  anyNA(
    plot_data$log2_fold_change
  ) ||
  any(
    !is.finite(
      plot_data$log2_fold_change
    )
  )
) {
  stop(
    "Invalid log2 fold-change values detected"
  )
}


fwrite(
  plot_data,
  plot_data_file,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)


# ============================================================
# Shared data-derived symmetric log2 scale
#
# The maximum absolute value is calculated from:
#
#   whole-flowcell enrichment by yield;
#   whole-flowcell enrichment by composition;
#   half-flowcell enrichment by yield;
#   half-flowcell enrichment by composition.
#
# The next integer is used as the common symmetric bound.
#
# For example:
#
#   maximum absolute observed value = 2.31
#
# gives:
#
#   scale = -3 to +3
# ============================================================

observed_log2_minimum <- min(
  plot_data$log2_fold_change,
  na.rm = TRUE
)

observed_log2_maximum <- max(
  plot_data$log2_fold_change,
  na.rm = TRUE
)

maximum_absolute_log2 <- max(
  abs(
    plot_data$log2_fold_change
  ),
  na.rm = TRUE
)

log2_axis_limit <- ceiling(
  maximum_absolute_log2
)

if (
  !is.finite(
    log2_axis_limit
  ) ||
  log2_axis_limit <= 0
) {
  stop(
    "Unable to calculate a valid shared log2 axis limit"
  )
}

log2_axis_lower <-
  -log2_axis_limit

log2_axis_upper <-
  log2_axis_limit

log2_axis_breaks <- seq(
  log2_axis_lower,
  log2_axis_upper,
  by = 1
)


if (
  observed_log2_minimum <
  log2_axis_lower ||
  observed_log2_maximum >
  log2_axis_upper
) {
  stop(
    paste0(
      "Observed log2 fold changes exceed the calculated ",
      "shared plotting range. Observed range: ",
      format(
        observed_log2_minimum,
        digits = 8,
        trim = TRUE
      ),
      " to ",
      format(
        observed_log2_maximum,
        digits = 8,
        trim = TRUE
      ),
      "; plotting range: ",
      log2_axis_lower,
      " to ",
      log2_axis_upper
    )
  )
}


# ============================================================
# Signed log2 tick labels
#
# Negative values use a Unicode minus sign.
#
# Positive values use an explicit plus sign.
#
# Zero is shown without a sign.
# ============================================================

signed_log2_labeller <- function(x) {
  vapply(
    x,

    function(value) {
      if (
        abs(
          value
        ) < calculation_tolerance
      ) {
        return(
          "0"
        )
      }

      if (value < 0) {
        paste0(
          "\u2212",
          abs(
            value
          )
        )
      } else {
        paste0(
          "+",
          value
        )
      }
    },

    character(1)
  )
}


# ============================================================
# Shared x scale
#
# Scale limits are applied with coord_cartesian() rather than
# the scale limits argument.
#
# This avoids dropping values during statistical processing.
# ============================================================

shared_log2_x_scale <- scale_x_continuous(
  breaks =
    log2_axis_breaks,

  labels =
    signed_log2_labeller,

  expand = expansion(
    mult = c(
      0.015,
      0.015
    )
  )
)


# ============================================================
# Axis summary
# ============================================================

axis_summary <- tibble(
  transformation =
    "log2",

  plotted_ratio =
    "accepted NAS / Metagenomic",

  sequencing_duration_hours =
    sequencing_duration_hours,

  yield_time_denominator =
    paste0(
      sequencing_duration_hours,
      " hours for both methods; common denominator cancels"
    ),

  scaling_scope =
    paste0(
      "Whole-flowcell and half-flowcell comparisons; ",
      "yield and composition"
    ),

  scaling_rule =
    paste0(
      "Ceiling of maximum absolute observed log2 fold change"
    ),

  minimum_observed_log2 =
    observed_log2_minimum,

  maximum_observed_log2 =
    observed_log2_maximum,

  maximum_absolute_observed_log2 =
    maximum_absolute_log2,

  plotted_lower_log2_limit =
    log2_axis_lower,

  plotted_upper_log2_limit =
    log2_axis_upper,

  tick_interval =
    1,

  negative_values_mean =
    "greater in Metagenomic",

  positive_values_mean =
    "greater in accepted NAS"
)


fwrite(
  axis_summary,
  axis_summary_file,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)


# ============================================================
# Fixed manuscript colours
# ============================================================

capture_colours <- c(
  "NAS target" =
    "#E7298A",

  "Related non-target" =
    "#00A6D6",

  "Other non-target" =
    "#707070"
)

capture_fill_scale <- scale_fill_manual(
  values =
    capture_colours,

  breaks =
    capture_class_levels,

  limits =
    capture_class_levels,

  drop =
    FALSE
)


# ============================================================
# Publication geometry
# ============================================================

axis_line_width <- 0.30
bar_outline_width <- 0.18
zero_line_width <- 0.30
bar_width <- 0.68


# ============================================================
# Shared publication theme
# ============================================================

theme_enrichment <- function() {
  theme_classic(
    base_size = 8,
    base_family = "sans"
  ) +
    theme(
      axis.line = element_line(
        colour = "black",
        linewidth =
          axis_line_width
      ),

      axis.ticks = element_line(
        colour = "black",
        linewidth =
          axis_line_width
      ),

      axis.ticks.length = unit(
        1.25,
        "mm"
      ),

      axis.text = element_text(
        colour = "black",
        size = 7.1
      ),

      axis.text.y = element_text(
        colour = "black",
        size = 6.7,
        hjust = 1,

        margin = margin(
          r = 1.5
        )
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

      panel.grid.major.x = element_line(
        colour = "grey85",
        linewidth = 0.22
      ),

      panel.grid.major.y =
        element_blank(),

      panel.grid.minor =
        element_blank(),

      panel.border =
        element_blank(),

      legend.position =
        "bottom",

      legend.direction =
        "horizontal",

      legend.title =
        element_blank(),

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

      plot.tag.position =
        "topleft"
    )
}


# ============================================================
# Build one analytical panel
#
# Negative and positive directions are communicated by the
# signed x-axis values.
#
# The ratio direction is defined in the figure legend.
# ============================================================

build_enrichment_panel <- function(
    panel_data,
    panel_tag,
    panel_title,
    show_species_labels,
    show_legend
) {
  panel_plot <- ggplot(
    panel_data,

    aes(
      x = log2_fold_change,
      y = species_display,
      fill = capture_class
    )
  ) +

    geom_col(
      width =
        bar_width,

      colour =
        "grey25",

      linewidth =
        bar_outline_width,

      show.legend =
        show_legend
    ) +

    geom_vline(
      xintercept = 0,
      colour = "black",
      linewidth = zero_line_width,
      linetype = "31",
      lineend = "butt"
    ) +

    capture_fill_scale +

    shared_log2_x_scale +

    scale_y_discrete(
      limits =
        species_display_order,

      labels =
        species_label_function,

      drop =
        FALSE,

      expand = expansion(
        add = c(
          0.35,
          0.35
        )
      )
    ) +

    coord_cartesian(
      xlim = c(
        log2_axis_lower,
        log2_axis_upper
      ),

      clip =
        "on"
    ) +

    labs(
      tag =
        panel_tag,

      title =
        panel_title,

      x = expression(
        log[2] ~
          "fold change"
      ),

      y =
        NULL,

      fill =
        NULL
    ) +

    theme_enrichment() +

    theme(
      plot.title = element_text(
        family = "sans",
        colour = "black",
        size = 9,
        face = "bold",
        hjust = 0.5,
        vjust = 1,

        margin = margin(
          b = 6
        )
      ),

      plot.margin = margin(
        t = 3,
        r = 2,
        b = 2,
        l = 2
      )
    )

  if (
    !show_species_labels
  ) {
    panel_plot <- panel_plot +

      theme(
        axis.text.y =
          element_blank(),

        legend.position =
          "none",

        plot.margin = margin(
          t = 3,
          r = 2,
          b = 2,
          l = 3
        )
      )
  }

  panel_plot
}


# ============================================================
# Build one two-panel comparison figure
# ============================================================

build_enrichment_figure <- function(
    comparison_name
) {
  comparison_data <- plot_data %>%
    filter(
      comparison ==
        comparison_name
    )

  if (
    nrow(
      comparison_data
    ) != 32
  ) {
    stop(
      paste(
        "Expected 32 plotting rows for",
        comparison_name,
        "but observed",
        nrow(
          comparison_data
        )
      )
    )
  }

  yield_data <- comparison_data %>%
    filter(
      metric ==
        "Enrichment by yield"
    )

  composition_data <- comparison_data %>%
    filter(
      metric ==
        "Enrichment by composition"
    )

  if (
    nrow(
      yield_data
    ) != 16
  ) {
    stop(
      paste(
        "Expected 16 yield rows for",
        comparison_name
      )
    )
  }

  if (
    nrow(
      composition_data
    ) != 16
  ) {
    stop(
      paste(
        "Expected 16 composition rows for",
        comparison_name
      )
    )
  }

  yield_panel <- build_enrichment_panel(
    panel_data =
      yield_data,

    panel_tag =
      "A)",

    panel_title =
      "Enrichment by yield",

    show_species_labels =
      TRUE,

    show_legend =
      TRUE
  )

  composition_panel <- build_enrichment_panel(
    panel_data =
      composition_data,

    panel_tag =
      "B)",

    panel_title =
      "Enrichment by composition",

    show_species_labels =
      FALSE,

    show_legend =
      FALSE
  )

  comparison_figure <- (
    yield_panel +
      composition_panel
  ) +

    plot_layout(
      ncol = 2,

      widths = c(
        1,
        1
      ),

      guides =
        "collect"
    ) &

    theme(
      legend.position =
        "bottom"
    )

  comparison_figure
}


# ============================================================
# Create requested figures
# ============================================================

whole_flowcell_plot <- build_enrichment_figure(
  comparison_name =
    "Whole-flowcell"
)

half_flowcell_plot <- build_enrichment_figure(
  comparison_name =
    "Half-flowcell"
)


# ============================================================
# Plot-build validation
# ============================================================

invisible(
  patchwork::patchworkGrob(
    whole_flowcell_plot
  )
)

invisible(
  patchwork::patchworkGrob(
    half_flowcell_plot
  )
)

message(
  paste0(
    "PASS: whole-flowcell manuscript and ",
    "half-flowcell supplementary plots built successfully"
  )
)


# ============================================================
# Display manuscript figure in RStudio
# ============================================================

print(
  whole_flowcell_plot
)


# ============================================================
# Replace a destination using a completed temporary file
# ============================================================

replace_destination_file <- function(
    temporary_file,
    destination_file
) {
  if (
    !file.exists(
      temporary_file
    )
  ) {
    stop(
      "Temporary figure file was not created"
    )
  }

  temporary_info <- file.info(
    temporary_file
  )

  if (
    is.na(
      temporary_info$size
    ) ||
    temporary_info$size <= 0
  ) {
    stop(
      "Temporary figure file is empty"
    )
  }

  if (
    file.exists(
      destination_file
    )
  ) {
    removed <- file.remove(
      destination_file
    )

    if (
      !removed &&
      file.exists(
        destination_file
      )
    ) {
      stop(
        paste0(
          "Unable to remove the existing destination file. ",
          "Close any application using the file and pause ",
          "OneDrive synchronisation: ",
          destination_file
        )
      )
    }
  }

  copied <- file.copy(
    from =
      temporary_file,

    to =
      destination_file,

    overwrite =
      FALSE,

    copy.mode =
      TRUE,

    copy.date =
      TRUE
  )

  if (!copied) {
    stop(
      paste0(
        "Unable to copy the completed figure to: ",
        destination_file
      )
    )
  }

  destination_info <- file.info(
    destination_file
  )

  if (
    !file.exists(
      destination_file
    ) ||
    is.na(
      destination_info$size
    ) ||
    destination_info$size <= 0
  ) {
    stop(
      "Copied destination file is missing or empty"
    )
  }

  if (
    destination_info$size !=
    temporary_info$size
  ) {
    stop(
      paste0(
        "Destination file size differs from the temporary ",
        "source. Temporary size: ",
        temporary_info$size,
        " bytes; destination size: ",
        destination_info$size,
        " bytes."
      )
    )
  }

  list(
    size_bytes =
      destination_info$size,

    modified_time =
      format(
        destination_info$mtime,
        "%Y-%m-%dT%H:%M:%S%z"
      )
  )
}


# ============================================================
# PNG export
# ============================================================

export_png <- function(
    plot_object,
    destination_file
) {
  temporary_file <- tempfile(
    pattern =
      "enrichment_",

    fileext =
      ".png"
  )

  on.exit(
    {
      if (
        file.exists(
          temporary_file
        )
      ) {
        unlink(
          temporary_file,
          force = TRUE
        )
      }
    },
    add = TRUE
  )

  tryCatch(
    {
      warning_messages <- character(
        0
      )

      withCallingHandlers(
        {
          ggsave(
            filename =
              temporary_file,

            plot =
              plot_object,

            width =
              180,

            height =
              126,

            units =
              "mm",

            dpi =
              600,

            bg =
              "white",

            device =
              ragg::agg_png
          )
        },

        warning = function(
    warning_condition
        ) {
          warning_messages <<- c(
            warning_messages,
            conditionMessage(
              warning_condition
            )
          )

          invokeRestart(
            "muffleWarning"
          )
        }
      )

      if (
        length(
          warning_messages
        ) > 0
      ) {
        stop(
          paste(
            warning_messages,
            collapse = "; "
          )
        )
      }

      copied_info <- replace_destination_file(
        temporary_file =
          temporary_file,

        destination_file =
          destination_file
      )

      list(
        success =
          TRUE,

        message =
          "Rendered, copied and verified",

        size_bytes =
          copied_info$size_bytes,

        modified_time =
          copied_info$modified_time
      )
    },

    error = function(
    error_condition
    ) {
      list(
        success =
          FALSE,

        message =
          conditionMessage(
            error_condition
          ),

        size_bytes =
          NA_real_,

        modified_time =
          NA_character_
      )
    }
  )
}


# ============================================================
# PDF export
# ============================================================

export_pdf <- function(
    plot_object,
    destination_file
) {
  temporary_file <- tempfile(
    pattern =
      "enrichment_",

    fileext =
      ".pdf"
  )

  on.exit(
    {
      if (
        file.exists(
          temporary_file
        )
      ) {
        unlink(
          temporary_file,
          force = TRUE
        )
      }
    },
    add = TRUE
  )

  tryCatch(
    {
      warning_messages <- character(
        0
      )

      cairo_success <- tryCatch(
        {
          withCallingHandlers(
            {
              ggsave(
                filename =
                  temporary_file,

                plot =
                  plot_object,

                width =
                  180,

                height =
                  126,

                units =
                  "mm",

                bg =
                  "white",

                device =
                  grDevices::cairo_pdf,

                family =
                  "sans"
              )
            },

            warning = function(
    warning_condition
            ) {
              warning_messages <<- c(
                warning_messages,
                conditionMessage(
                  warning_condition
                )
              )

              invokeRestart(
                "muffleWarning"
              )
            }
          )

          file.exists(
            temporary_file
          ) &&
            !is.na(
              file.info(
                temporary_file
              )$size
            ) &&
            file.info(
              temporary_file
            )$size > 0
        },

    error = function(
    error_condition
    ) {
      FALSE
    }
      )

      if (
        !cairo_success ||
        length(
          warning_messages
        ) > 0
      ) {
        if (
          file.exists(
            temporary_file
          )
        ) {
          unlink(
            temporary_file,
            force = TRUE
          )
        }

        warning_messages <- character(
          0
        )

        withCallingHandlers(
          {
            ggsave(
              filename =
                temporary_file,

              plot =
                plot_object,

              width =
                180,

              height =
                126,

              units =
                "mm",

              bg =
                "white",

              device =
                grDevices::pdf,

              family =
                "sans",

              useDingbats =
                FALSE
            )
          },

          warning = function(
    warning_condition
          ) {
            warning_messages <<- c(
              warning_messages,
              conditionMessage(
                warning_condition
              )
            )

            invokeRestart(
              "muffleWarning"
            )
          }
        )
      }

      if (
        length(
          warning_messages
        ) > 0
      ) {
        stop(
          paste(
            warning_messages,
            collapse = "; "
          )
        )
      }

      copied_info <- replace_destination_file(
        temporary_file =
          temporary_file,

        destination_file =
          destination_file
      )

      list(
        success =
          TRUE,

        message = if (
          cairo_success
        ) {
          "Rendered with cairo_pdf, copied and verified"
        } else {
          "Rendered with standard PDF device, copied and verified"
        },

        size_bytes =
          copied_info$size_bytes,

        modified_time =
          copied_info$modified_time
      )
    },

    error = function(
    error_condition
    ) {
      list(
        success =
          FALSE,

        message =
          conditionMessage(
            error_condition
          ),

        size_bytes =
          NA_real_,

        modified_time =
          NA_character_
      )
    }
  )
}


# ============================================================
# SVG export
# ============================================================

export_svg <- function(
    plot_object,
    destination_file
) {
  temporary_file <- tempfile(
    pattern =
      "enrichment_",

    fileext =
      ".svg"
  )

  on.exit(
    {
      if (
        file.exists(
          temporary_file
        )
      ) {
        unlink(
          temporary_file,
          force = TRUE
        )
      }
    },
    add = TRUE
  )

  tryCatch(
    {
      warning_messages <- character(
        0
      )

      withCallingHandlers(
        {
          ggsave(
            filename =
              temporary_file,

            plot =
              plot_object,

            width =
              180,

            height =
              126,

            units =
              "mm",

            bg =
              "white",

            device =
              svglite::svglite
          )
        },

        warning = function(
    warning_condition
        ) {
          warning_messages <<- c(
            warning_messages,
            conditionMessage(
              warning_condition
            )
          )

          invokeRestart(
            "muffleWarning"
          )
        }
      )

      if (
        length(
          warning_messages
        ) > 0
      ) {
        stop(
          paste(
            warning_messages,
            collapse = "; "
          )
        )
      }

      copied_info <- replace_destination_file(
        temporary_file =
          temporary_file,

        destination_file =
          destination_file
      )

      list(
        success =
          TRUE,

        message =
          "Rendered, copied and verified",

        size_bytes =
          copied_info$size_bytes,

        modified_time =
          copied_info$modified_time
      )
    },

    error = function(
    error_condition
    ) {
      list(
        success =
          FALSE,

        message =
          conditionMessage(
            error_condition
          ),

        size_bytes =
          NA_real_,

        modified_time =
          NA_character_
      )
    }
  )
}


# ============================================================
# Export all six requested figure files
# ============================================================

whole_png_result <- export_png(
  whole_flowcell_plot,
  whole_figure_png
)

whole_pdf_result <- export_pdf(
  whole_flowcell_plot,
  whole_figure_pdf
)

whole_svg_result <- export_svg(
  whole_flowcell_plot,
  whole_figure_svg
)

half_png_result <- export_png(
  half_flowcell_plot,
  half_figure_png
)

half_pdf_result <- export_pdf(
  half_flowcell_plot,
  half_figure_pdf
)

half_svg_result <- export_svg(
  half_flowcell_plot,
  half_figure_svg
)


# ============================================================
# Export-status table
# ============================================================

export_status <- tibble(
  comparison = c(
    "Whole-flowcell",
    "Whole-flowcell",
    "Whole-flowcell",
    "Half-flowcell",
    "Half-flowcell",
    "Half-flowcell"
  ),

  format = c(
    "PNG",
    "PDF",
    "SVG",
    "PNG",
    "PDF",
    "SVG"
  ),

  file = c(
    whole_figure_png,
    whole_figure_pdf,
    whole_figure_svg,
    half_figure_png,
    half_figure_pdf,
    half_figure_svg
  ),

  success = c(
    whole_png_result$success,
    whole_pdf_result$success,
    whole_svg_result$success,
    half_png_result$success,
    half_pdf_result$success,
    half_svg_result$success
  ),

  message = c(
    whole_png_result$message,
    whole_pdf_result$message,
    whole_svg_result$message,
    half_png_result$message,
    half_pdf_result$message,
    half_svg_result$message
  ),

  size_bytes = c(
    whole_png_result$size_bytes,
    whole_pdf_result$size_bytes,
    whole_svg_result$size_bytes,
    half_png_result$size_bytes,
    half_pdf_result$size_bytes,
    half_svg_result$size_bytes
  ),

  modified_time = c(
    whole_png_result$modified_time,
    whole_pdf_result$modified_time,
    whole_svg_result$modified_time,
    half_png_result$modified_time,
    half_pdf_result$modified_time,
    half_svg_result$modified_time
  )
)


fwrite(
  export_status,
  export_status_file,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)


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
    "Canonical-data directory:",
    canonical_data_dir
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
    "Sequencing duration per method:",
    sequencing_duration_hours,
    "hours"
  ),

  paste0(
    "Enrichment by yield: species-specific aligned ",
    "mitochondrial query-base yield per hour in accepted NAS ",
    "divided by the corresponding yield per hour in ",
    "Metagenomic sequencing"
  ),

  paste0(
    "Yield time denominator: both methods were analysed over ",
    sequencing_duration_hours,
    " hours, so the common time denominator cancels"
  ),

  paste0(
    "Enrichment by composition: species-specific ",
    "mitochondrial relative-abundance fraction in accepted NAS ",
    "divided by the corresponding fraction in Metagenomic ",
    "sequencing"
  ),

  "Plotting transformation: log2 ratio",

  "Plotted ratio: accepted NAS / Metagenomic",

  "Negative values: greater in Metagenomic sequencing",

  "Positive values: greater in accepted NAS",

  "Zero: equal between methods",

  "Whole-flowcell: full_nas_accepted divided by full_meta",

  "Half-flowcell: half_nas_accepted divided by half_meta",

  "Main manuscript figure: whole-flowcell comparison",

  "Supplementary figure: half-flowcell comparison",

  "Combined four-panel figure: not generated",

  "Signed raw-fold figure: not generated",

  "Panel A title: Enrichment by yield",

  "Panel B title: Enrichment by composition",

  "Species labels shown in Panel A only",

  "Zero handling: none required",

  "Pseudocount: none",

  "Species included per comparison: 16",

  "Species order: whole-flowcell yield enrichment",

  "Species order shared across both figures",

  "NAS target colour: #E7298A",

  "Related non-target colour: #00A6D6",

  "Other non-target colour: #707070",

  "Axis-line and tick width: 0.30",

  "Bar-outline width: 0.18",

  "Zero-reference line width: 0.30",

  "Zero-reference line: dashed",

  "Panel labels: A), B)",

  "Panel-label size: 9 pt",

  paste(
    "Minimum observed log2 fold change:",
    format(
      observed_log2_minimum,
      digits = 8,
      trim = TRUE
    )
  ),

  paste(
    "Maximum observed log2 fold change:",
    format(
      observed_log2_maximum,
      digits = 8,
      trim = TRUE
    )
  ),

  paste(
    "Maximum absolute observed log2 fold change:",
    format(
      maximum_absolute_log2,
      digits = 8,
      trim = TRUE
    )
  ),

  paste(
    "Shared lower log2 limit:",
    log2_axis_lower
  ),

  paste(
    "Shared upper log2 limit:",
    log2_axis_upper
  ),

  paste(
    "Shared tick interval:",
    1
  ),

  paste0(
    "Scale-selection rule: ceiling of maximum absolute ",
    "observed log2 fold change across both comparisons and ",
    "both enrichment metrics"
  ),

  "Figure width: 180 mm",

  "Figure height: 126 mm",

  "PNG resolution: 600 dpi",

  paste(
    "Results table:",
    results_file
  ),

  paste(
    "Summary table:",
    summary_file
  ),

  paste(
    "Plot-data table:",
    plot_data_file
  ),

  paste(
    "Species-label map:",
    species_label_file
  ),

  paste(
    "Calculation validation:",
    calculation_validation_file
  ),

  paste(
    "Capture-class validation:",
    class_validation_file
  ),

  paste(
    "Axis summary:",
    axis_summary_file
  ),

  paste(
    "Export-status table:",
    export_status_file
  ),

  paste(
    "Whole-flowcell manuscript PNG:",
    whole_figure_png
  ),

  paste(
    "Whole-flowcell manuscript PDF:",
    whole_figure_pdf
  ),

  paste(
    "Whole-flowcell manuscript SVG:",
    whole_figure_svg
  ),

  paste(
    "Half-flowcell supplementary PNG:",
    half_figure_png
  ),

  paste(
    "Half-flowcell supplementary PDF:",
    half_figure_pdf
  ),

  paste(
    "Half-flowcell supplementary SVG:",
    half_figure_svg
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
  "Log2 calculation validation:",
  log_connection
)

writeLines(
  capture.output(
    print(
      calculation_validation
    )
  ),
  log_connection
)

writeLines(
  "",
  log_connection
)

writeLines(
  "Capture-class count validation:",
  log_connection
)

writeLines(
  capture.output(
    print(
      class_count_validation
    )
  ),
  log_connection
)

writeLines(
  "",
  log_connection
)

writeLines(
  "Descriptive enrichment summary:",
  log_connection
)

writeLines(
  capture.output(
    print(
      summary_table
    )
  ),
  log_connection
)

writeLines(
  "",
  log_connection
)

writeLines(
  "Axis summary:",
  log_connection
)

writeLines(
  capture.output(
    print(
      axis_summary
    )
  ),
  log_connection
)

writeLines(
  "",
  log_connection
)

writeLines(
  "Figure export status:",
  log_connection
)

writeLines(
  capture.output(
    print(
      export_status
    )
  ),
  log_connection
)

writeLines(
  "",
  log_connection
)

writeLines(
  "Species display order from bottom to top:",
  log_connection
)

writeLines(
  species_display_order,
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
# Final completion output
# ============================================================

print(
  summary_table
)

print(
  axis_summary
)

print(
  export_status
)


if (
  any(
    export_status$success != TRUE
  )
) {
  failed_exports <- export_status %>%
    filter(
      success != TRUE
    )

  stop(
    paste0(
      "One or more requested figure exports failed:\n",
      paste(
        paste0(
          failed_exports$comparison,
          " ",
          failed_exports$format,
          ": ",
          failed_exports$message
        ),
        collapse = "\n"
      )
    )
  )
}


message(
  "PASS: all six figure files were newly rendered, copied and verified"
)

message(
  "Shared log2 plotting range: ",
  log2_axis_lower,
  " to +",
  log2_axis_upper
)

message(
  "Wrote: ",
  results_file
)

message(
  "Wrote: ",
  summary_file
)

message(
  "Wrote: ",
  plot_data_file
)

message(
  "Wrote: ",
  species_label_file
)

message(
  "Wrote: ",
  calculation_validation_file
)

message(
  "Wrote: ",
  class_validation_file
)

message(
  "Wrote: ",
  axis_summary_file
)

message(
  "Wrote: ",
  export_status_file
)

message(
  "Wrote: ",
  whole_figure_png
)

message(
  "Wrote: ",
  whole_figure_pdf
)

message(
  "Wrote: ",
  whole_figure_svg
)

message(
  "Wrote: ",
  half_figure_png
)

message(
  "Wrote: ",
  half_figure_pdf
)

message(
  "Wrote: ",
  half_figure_svg
)

message(
  "Wrote: ",
  log_file
)

message(
  "Wrote: ",
  session_file
)
