#!/usr/bin/env Rscript
# r_analysis/main_analysis_pipeline.R
# Main pipeline script that orchestrates all R analysis steps

# Load required functions
source("functions/utils.R")

# Main pipeline function
main_r_analysis_pipeline <- function(config_file) {

  cat("Starting Complete R Analysis Pipeline\n")
  cat("=====================================\n\n")

  # Load required packages
  load_required_packages()

  # Read configuration
  config <- read_config(config_file)

  # Validate configuration has all required sections
  required_sections <- c("nanotel_analysis", "mapping_analysis", "methylation_analysis")
  missing_sections <- required_sections[!(required_sections %in% names(config))]
  if (length(missing_sections) > 0) {
    stop("Missing required configuration sections: ", paste(missing_sections, collapse = ", "))
  }

  pipeline_start_time <- Sys.time()
  results <- list()

  log_message("Starting complete R analysis pipeline")
  log_message(paste("Base output directory:", config$base_output_dir))

  # Step 1: NanoTel Analysis
  if (config$run_nanotel_analysis %||% TRUE) {
    log_message("=" * 50)
    log_message("STEP 1: NANOTEL ANALYSIS")
    log_message("=" * 50)

    tryCatch({
      # Create temporary config file for NanoTel analysis
      nanotel_config_file <- create_temp_config(config$nanotel_analysis, "nanotel_temp_config.json")

      # Source and run NanoTel analysis
      source("batch_nanotel_analysis.R")
      nanotel_result <- main_nanotel_analysis(nanotel_config_file)

      results$nanotel <- nanotel_result
      log_message("✓ NanoTel analysis completed successfully")

      # Clean up temp file
      file.remove(nanotel_config_file)

    }, error = function(e) {
      log_message(paste("✗ NanoTel analysis failed:", e$message), "ERROR")
      if (config$stop_on_error %||% TRUE) {
        stop("Pipeline stopped due to NanoTel analysis failure")
      }
    })
  } else {
    log_message("Skipping NanoTel analysis (disabled in config)")
  }

  # Step 2: Mapping Analysis
  if (config$run_mapping_analysis %||% TRUE) {
    log_message("=" * 50)
    log_message("STEP 2: MAPPING ANALYSIS")
    log_message("=" * 50)

    tryCatch({
      # Update mapping config with NanoTel results if available
      mapping_config <- config$mapping_analysis
      if (!is.null(results$nanotel)) {
        mapping_config$filtered_nanotel_dir <- config$nanotel_analysis$output_dir
      }

      # Create temporary config file for mapping analysis
      mapping_config_file <- create_temp_config(mapping_config, "mapping_temp_config.json")

      # Source and run mapping analysis
      source("batch_mapping_analysis.R")
      mapping_result <- main_mapping_analysis(mapping_config_file)

      results$mapping <- mapping_result
      log_message("✓ Mapping analysis completed successfully")

      # Clean up temp file
      file.remove(mapping_config_file)

    }, error = function(e) {
      log_message(paste("✗ Mapping analysis failed:", e$message), "ERROR")
      if (config$stop_on_error %||% TRUE) {
        stop("Pipeline stopped due to mapping analysis failure")
      }
    })
  } else {
    log_message("Skipping mapping analysis (disabled in config)")
  }

  # Step 3: Methylation Analysis
  if (config$run_methylation_analysis %||% TRUE) {
    log_message("=" * 50)
    log_message("STEP 3: METHYLATION ANALYSIS")
    log_message("=" * 50)

    tryCatch({
      # Update methylation config with mapping results if available
      methylation_config <- config$methylation_analysis
      if (!is.null(results$mapping)) {
        methylation_config$pileup_bed_dir <- config$mapping_analysis$output_dir
      }

      # Create temporary config file for methylation analysis
      methylation_config_file <- create_temp_config(methylation_config, "methylation_temp_config.json")

      # Source and run methylation analysis
      source("batch_methylation_prep.R")
      methylation_result <- main_methylation_analysis(methylation_config_file)

      results$methylation <- methylation_result
      log_message("✓ Methylation analysis completed successfully")

      # Clean up temp file
      file.remove(methylation_config_file)

    }, error = function(e) {
      log_message(paste("✗ Methylation analysis failed:", e$message), "ERROR")
      if (config$stop_on_error %||% TRUE) {
        stop("Pipeline stopped due to methylation analysis failure")
      }
    })
  } else {
    log_message("Skipping methylation analysis (disabled in config)")
  }

  pipeline_end_time <- Sys.time()
  pipeline_duration <- as.numeric(difftime(pipeline_end_time, pipeline_start_time, units = "mins"))

  # Generate final pipeline report
  generate_pipeline_report(results, config, pipeline_duration)

  log_message("=" * 50)
  log_message("PIPELINE COMPLETED SUCCESSFULLY")
  log_message("=" * 50)
  log_message(paste("Total duration:", round(pipeline_duration, 1), "minutes"))

  return(results)
}

# Create temporary configuration file
create_temp_config <- function(config_section, filename) {
  temp_file <- file.path(tempdir(), filename)
  jsonlite::write_json(config_section, temp_file, auto_unbox = TRUE, pretty = TRUE)
  return(temp_file)
}

# Generate comprehensive pipeline report
generate_pipeline_report <- function(results, config, duration) {

  log_message("Generating comprehensive pipeline report")

  report_file <- file.path(config$base_output_dir, "complete_pipeline_report.txt")

  report_lines <- c(
    "=" * 80,
    "COMPLETE R ANALYSIS PIPELINE REPORT",
    "=" * 80,
    paste("Pipeline date:", Sys.Date()),
    paste("Pipeline start time:", format(Sys.time() - duration * 60, "%H:%M:%S")),
    paste("Pipeline end time:", format(Sys.time(), "%H:%M:%S")),
    paste("Total duration:", round(duration, 1), "minutes"),
    "",
    "PIPELINE CONFIGURATION:",
    paste("  Base output directory:", config$base_output_dir),
    paste("  NanoTel analysis:", ifelse(config$run_nanotel_analysis %||% TRUE, "ENABLED", "DISABLED")),
    paste("  Mapping analysis:", ifelse(config$run_mapping_analysis %||% TRUE, "ENABLED", "DISABLED")),
    paste("  Methylation analysis:", ifelse(config$run_methylation_analysis %||% TRUE, "ENABLED", "DISABLED")),
    paste("  Stop on error:", ifelse(config$stop_on_error %||% TRUE, "YES", "NO")),
    "",
    "ANALYSIS RESULTS:"
  )

  # NanoTel results
  if (!is.null(results$nanotel)) {
    nanotel <- results$nanotel
    report_lines <- c(report_lines,
                      "",
                      "1. NANOTEL ANALYSIS:",
                      paste("   Status: SUCCESS"),
                      paste("   Barcodes processed:", nanotel$barcodes_processed),
                      paste("   Files processed:", nanotel$files_processed),
                      paste("   Output directory:", config$nanotel_analysis$output_dir)
    )
  } else {
    report_lines <- c(report_lines,
                      "",
                      "1. NANOTEL ANALYSIS:",
                      "   Status: SKIPPED OR FAILED"
    )
  }

  # Mapping results
  if (!is.null(results$mapping)) {
    mapping <- results$mapping
    report_lines <- c(report_lines,
                      "",
                      "2. MAPPING ANALYSIS:",
                      paste("   Status: SUCCESS"),
                      paste("   Total processed:", mapping$total_processed),
                      paste("   Successful barcodes:", length(mapping$successful_results)),
                      paste("   Failed barcodes:", length(mapping$failed_barcodes)),
                      paste("   Output directory:", config$mapping_analysis$output_dir)
    )

    if (length(mapping$failed_barcodes) > 0) {
      report_lines <- c(report_lines,
                        paste("   Failed barcode list:", paste(mapping$failed_barcodes, collapse = ", "))
      )
    }
  } else {
    report_lines <- c(report_lines,
                      "",
                      "2. MAPPING ANALYSIS:",
                      "   Status: SKIPPED OR FAILED"
    )
  }

  # Methylation results
  if (!is.null(results$methylation)) {
    methylation <- results$methylation
    report_lines <- c(report_lines,
                      "",
                      "3. METHYLATION ANALYSIS:",
                      paste("   Status: SUCCESS"),
                      paste("   BED files processed:", methylation$processed_files),
                      paste("   Total methylation sites:", methylation$total_sites),
                      paste("   Output directory:", config$methylation_analysis$output_dir)
    )
  } else {
    report_lines <- c(report_lines,
                      "",
                      "3. METHYLATION ANALYSIS:",
                      "   Status: SKIPPED OR FAILED"
    )
  }

  report_lines <- c(report_lines,
                    "",
                    "OUTPUT STRUCTURE:",
                    paste("  ", config$base_output_dir, "/"),
                    "    ├── nanotel_output/",
                    "    │   ├── filtered_summary*.csv",
                    "    │   ├── nanotel_summary_statistics.csv",
                    "    │   └── nanotel_analysis_report.txt",
                    "    ├── mapping_output/",
                    "    │   ├── mapped*.csv",
                    "    │   ├── filtered_*.bam",
                    "    │   ├── pileup-*.bed",
                    "    │   └── mapping_analysis_report.txt",
                    "    ├── methylation_output/",
                    "    │   ├── plots/",
                    "    │   ├── processed_data/",
                    "    │   ├── shiny_app/ (if enabled)",
                    "    │   ├── methylation_summary_statistics.csv",
                    "    │   └── methylation_analysis_report.txt",
                    "    └── complete_pipeline_report.txt",
                    "",
                    "NEXT STEPS:",
                    "  1. Review individual analysis reports for detailed results",
                    "  2. Check any failed barcodes and investigate issues",
                    "  3. Use interactive Shiny app for methylation visualization",
                    "  4. Proceed with downstream analysis using processed data",
                    "",
                    "=" * 80
  )

  # Write report
  writeLines(report_lines, report_file)

  log_message(paste("Complete pipeline report saved to:", basename(report_file)))

  return(report_file)
}

# Command line interface
if (!interactive()) {
  # Parse command line arguments
  args <- commandArgs(trailingOnly = TRUE)

  if (length(args) != 1) {
    cat("Usage: Rscript main_analysis_pipeline.R <config_file>\n")
    cat("\n")
    cat("Arguments:\n")
    cat("  config_file: JSON configuration file with complete pipeline parameters\n")
    cat("\n")
    cat("Example:\n")
    cat("  Rscript main_analysis_pipeline.R config/pipeline_config.json\n")
    cat("\n")
    cat("The configuration file should contain sections for:\n")
    cat("  - nanotel_analysis\n")
    cat("  - mapping_analysis\n")
    cat("  - methylation_analysis\n")
    quit(status = 1)
  }

  config_file <- args[1]

  # Validate config file exists
  if (!file.exists(config_file)) {
    cat("Error: Configuration file not found:", config_file, "\n")
    quit(status = 1)
  }

  # Run complete pipeline
  tryCatch({
    results <- main_r_analysis_pipeline(config_file)
    cat("\n🎉 Complete R analysis pipeline finished successfully!\n")

    # Print summary
    if (!is.null(results$nanotel)) {
      cat("✓ NanoTel: Processed", results$nanotel$barcodes_processed, "barcodes\n")
    }
    if (!is.null(results$mapping)) {
      cat("✓ Mapping: Processed", length(results$mapping$successful_results), "barcodes\n")
    }
    if (!is.null(results$methylation)) {
      cat("✓ Methylation:", results$methylation$total_sites, "sites analyzed\n")
    }

  }, error = function(e) {
    cat("❌ Pipeline failed:", e$message, "\n")
    quit(status = 1)
  })
}
