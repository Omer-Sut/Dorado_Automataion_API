# r_analysis/functions/mapping_functions.R
# Functions for alignment mapping and BAM processing

source("functions/utils.R")

# Process alignment and NanoTel data for a single barcode
process_barcode_mapping <- function(alignment_summary_path,
                                    filtered_nanotel_path,
                                    bam_file_path,
                                    output_dir,
                                    barcode_name,
                                    min_mapq = 10,
                                    head_max_start = 5000,
                                    tail_min_end = 35000) {

  log_message(paste("Processing mapping for barcode:", barcode_name))

  # Validate input files
  validate_file_exists(alignment_summary_path, "Alignment summary")
  validate_file_exists(filtered_nanotel_path, "Filtered NanoTel data")
  validate_file_exists(bam_file_path, "BAM file")

  ensure_directory_exists(output_dir)

  # Read alignment summary
  minimap <- safe_read_csv(alignment_summary_path, col_types = cols())

  # Remove unwanted columns if they exist
  unwanted_cols <- c("filename", "alignment_bed_hits")
  minimap <- minimap[, !(names(minimap) %in% unwanted_cols)]

  log_message(paste("Loaded alignment data:", nrow(minimap), "alignments"))

  # Read filtered NanoTel data
  nanotel_data <- safe_read_csv(filtered_nanotel_path, col_types = cols())

  # Remove unwanted columns
  unwanted_nanotel_cols <- c("...1", "telo_density", "Telomere_start",
                             "Telomere_end", "Telomere_length",
                             "running_median", "seqLen_runningMED")
  nanotel_data <- nanotel_data[, !(names(nanotel_data) %in% unwanted_nanotel_cols)]

  log_message(paste("Loaded NanoTel data:", nrow(nanotel_data), "reads"))

  # Filter minimap to only reads that have matching values in NanoTel read_id
  minimap_filtered <- minimap %>%
    filter(sapply(.data[["read_id"]], function(x) {
      any(grepl(x, nanotel_data[["read_id"]], fixed = TRUE))
    }))

  log_message(paste("After NanoTel matching:", nrow(minimap_filtered), "alignments"))

  # Remove unmapped reads (those with "*" in alignment_direction)
  minimap_filtered <- minimap_filtered %>%
    filter(!grepl("\\*", alignment_direction))

  log_message(paste("After removing unmapped:", nrow(minimap_filtered), "alignments"))

  # Apply quality and position filters
  minimap_filtered <- minimap_filtered %>%
    filter(
      # Direction and genome position rules
      (grepl("Tail", alignment_genome) & alignment_direction == "-") |
        (grepl("Head", alignment_genome) & alignment_direction == "+"),
      # Quality filter
      alignment_mapq >= min_mapq,
      # Position-specific filters
      (grepl("Tail", alignment_genome) & alignment_genome_end >= tail_min_end) |
        (grepl("Head", alignment_genome) & alignment_genome_start <= head_max_start)
    )

  log_message(paste("After quality/position filtering:", nrow(minimap_filtered), "alignments"))

  if (nrow(minimap_filtered) == 0) {
    warning("No alignments passed filters for barcode: ", barcode_name)
    return(list(mapped_data = data.frame(), filtered_ids = character()))
  }

  # Merge with NanoTel data
  merged_data <- merge(minimap_filtered, nanotel_data, by = "read_id")

  log_message(paste("Final merged data:", nrow(merged_data), "reads"))

  # Save mapped data
  mapped_output_file <- file.path(output_dir, paste0("mapped", barcode_name, ".csv"))
  safe_write_csv(merged_data, mapped_output_file)

  # Extract read IDs for BAM filtering
  filtered_read_ids <- merged_data$read_id

  return(list(
    mapped_data = merged_data,
    filtered_ids = filtered_read_ids,
    output_file = mapped_output_file
  ))
}

# Extract unique read IDs from BAM file
get_bam_read_ids <- function(bam_file_path) {

  log_message(paste("Extracting read IDs from BAM:", basename(bam_file_path)))

  validate_file_exists(bam_file_path, "BAM file")

  # Build samtools command
  command <- paste(
    "samtools view", shQuote(bam_file_path),
    "| cut -f1 | sort | uniq"
  )

  tryCatch({
    read_ids <- system(command, intern = TRUE)
    log_message(paste("Extracted", length(read_ids), "unique read IDs from BAM"))
    return(read_ids)
  }, error = function(e) {
    stop("Failed to extract read IDs from BAM: ", e$message)
  })
}

# Filter BAM file based on read ID list
create_filtered_bam <- function(filtered_read_ids,
                                input_bam_path,
                                output_bam_path,
                                temp_ids_file = NULL) {

  log_message(paste("Creating filtered BAM with", length(filtered_read_ids), "reads"))

  validate_file_exists(input_bam_path, "Input BAM file")
  ensure_directory_exists(dirname(output_bam_path))

  # Create temporary file for read IDs if not provided
  if (is.null(temp_ids_file)) {
    temp_ids_file <- paste0(output_bam_path, "_temp_ids.txt")
  }

  # Write read IDs to temporary file
  writeLines(filtered_read_ids, temp_ids_file)

  tryCatch({
    # Filter BAM using samtools
    filter_cmd <- paste(
      "samtools view -N", shQuote(temp_ids_file),
      "-b", shQuote(input_bam_path),
      ">", shQuote(output_bam_path)
    )
    run_system_command(filter_cmd, "BAM filtering")

    # Index the filtered BAM file
    index_cmd <- paste("samtools index", shQuote(output_bam_path))
    run_system_command(index_cmd, "BAM indexing")

    # Clean up temporary file
    if (file.exists(temp_ids_file)) {
      file.remove(temp_ids_file)
    }

    log_message(paste("Filtered BAM created:", basename(output_bam_path)))
    return(output_bam_path)

  }, error = function(e) {
    # Clean up on error
    if (file.exists(temp_ids_file)) {
      file.remove(temp_ids_file)
    }
    stop("Failed to create filtered BAM: ", e$message)
  })
}

# Run modkit pileup on filtered BAM
run_modkit_pileup <- function(filtered_bam_path,
                              output_bed_path,
                              filter_threshold_c = 0.75,
                              filter_threshold_g = 0.75,
                              mod_threshold_m = 0.75) {

  log_message(paste("Running modkit pileup on:", basename(filtered_bam_path)))

  validate_file_exists(filtered_bam_path, "Filtered BAM file")
  ensure_directory_exists(dirname(output_bed_path))

  # Build modkit command
  modkit_cmd <- paste(
    "modkit pileup",
    paste0("--filter-threshold C:", filter_threshold_c),
    paste0("--filter-threshold G:", filter_threshold_g),
    paste0("--mod-threshold m:", mod_threshold_m),
    shQuote(filtered_bam_path),
    shQuote(output_bed_path)
  )

  tryCatch({
    run_system_command(modkit_cmd, "Modkit pileup")
    log_message(paste("Modkit pileup completed:", basename(output_bed_path)))
    return(output_bed_path)
  }, error = function(e) {
    stop("Modkit pileup failed: ", e$message)
  })
}

# Complete mapping workflow for a single barcode
process_complete_barcode_workflow <- function(barcode_config) {

  barcode_name <- barcode_config$barcode_name
  log_message(paste("Starting complete workflow for:", barcode_name))

  # Step 1: Process mapping
  mapping_result <- process_barcode_mapping(
    alignment_summary_path = barcode_config$alignment_summary_path,
    filtered_nanotel_path = barcode_config$filtered_nanotel_path,
    bam_file_path = barcode_config$bam_file_path,
    output_dir = barcode_config$output_dir,
    barcode_name = barcode_name,
    min_mapq = barcode_config$min_mapq,
    head_max_start = barcode_config$head_max_start,
    tail_min_end = barcode_config$tail_min_end
  )

  if (length(mapping_result$filtered_ids) == 0) {
    warning("No reads to process for ", barcode_name)
    return(NULL)
  }

  # Step 2: Create filtered BAM
  filtered_bam_path <- file.path(barcode_config$output_dir,
                                 paste0("filtered_", barcode_name, ".bam"))

  create_filtered_bam(
    filtered_read_ids = mapping_result$filtered_ids,
    input_bam_path = barcode_config$bam_file_path,
    output_bam_path = filtered_bam_path
  )

  # Step 3: Run modkit pileup
  pileup_bed_path <- file.path(barcode_config$output_dir,
                               paste0("pileup-", tolower(barcode_name), ".bed"))

  run_modkit_pileup(
    filtered_bam_path = filtered_bam_path,
    output_bed_path = pileup_bed_path,
    filter_threshold_c = barcode_config$filter_threshold_c,
    filter_threshold_g = barcode_config$filter_threshold_g,
    mod_threshold_m = barcode_config$mod_threshold_m
  )

  log_message(paste("Complete workflow finished for:", barcode_name))

  return(list(
    barcode = barcode_name,
    mapped_file = mapping_result$output_file,
    filtered_bam = filtered_bam_path,
    pileup_bed = pileup_bed_path,
    read_count = length(mapping_result$filtered_ids)
  ))
}

# Find BAM files for barcodes
find_barcode_bam_files <- function(bam_directory, barcode_pattern = "barcode\\d+") {

  log_message(paste("Looking for BAM files in:", bam_directory))

  bam_files <- find_files_by_pattern(bam_directory, "\\.bam$", recursive = TRUE)

  # Filter for barcode files
  barcode_bams <- bam_files[grepl(barcode_pattern, basename(bam_files), ignore.case = TRUE)]

  if (length(barcode_bams) == 0) {
    warning("No barcode BAM files found in: ", bam_directory)
    return(data.frame())
  }

  # Extract barcode names
  bam_info <- data.frame(
    bam_path = barcode_bams,
    barcode = sapply(barcode_bams, extract_barcode_from_path),
    stringsAsFactors = FALSE
  )

  log_message(paste("Found", nrow(bam_info), "barcode BAM files"))

  return(bam_info)
}
