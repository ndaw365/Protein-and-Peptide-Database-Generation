# Fixed MOLT4 Mutation Analysis Script

# Load required libraries
library(tidyverse)  # For data manipulation
library(httr)       # For API requests
library(stringr)    # For string manipulation
library(readr)      # For reading/writing CSVs
library(knitr)      # For tables

# Read the input CSV file - adjust path as needed
molt4_data <- read.csv("~/Downloads/MOLT4 mutations.csv", stringsAsFactors = FALSE)

# Display the first few rows to understand the data structure
head(molt4_data)

# Check column names to confirm we have "Uniprot ID" and "Variant Info" columns
colnames(molt4_data)

# Count total rows in the dataset
total_entries <- nrow(molt4_data)
cat("Total entries in the dataset:", total_entries, "\n")

# Check if columns exist (case insensitive check)
uniprot_col <- grep("(?i)uniprot.*id", colnames(molt4_data), perl = TRUE)
variant_col <- grep("(?i)variant.*info", colnames(molt4_data), perl = TRUE)

if (length(uniprot_col) == 0 || length(variant_col) == 0) {
  stop("Required columns not found. Please check column names in the CSV file.")
}

# Assign proper column names
uniprot_col_name <- colnames(molt4_data)[uniprot_col]
variant_col_name <- colnames(molt4_data)[variant_col]

# Count entries with non-empty UniProt IDs
entries_with_uniprot <- molt4_data %>%
  filter(!is.na(!!sym(uniprot_col_name)) & !!sym(uniprot_col_name) != "") %>%
  nrow()

# Count entries with both UniProt IDs and missense variants
entries_with_missense <- molt4_data %>%
  filter(!is.na(!!sym(uniprot_col_name)) & !!sym(uniprot_col_name) != "") %>%
  filter(grepl("missense_variant", !!sym(variant_col_name), ignore.case = TRUE)) %>%
  nrow()

# Create a filtered dataframe with only the entries we want
filtered_data <- molt4_data %>%
  filter(!is.na(!!sym(uniprot_col_name)) & !!sym(uniprot_col_name) != "") %>%
  filter(grepl("missense_variant", !!sym(variant_col_name), ignore.case = TRUE))

# Document the counts
cat("Proteins with UniProt IDs:", entries_with_uniprot, "\n")
cat("Proteins with missense variants:", entries_with_missense, "\n")

# Function to fetch canonical sequence from UniProt using the current REST API
get_canonical_sequence <- function(uniprot_id) {
  if (is.na(uniprot_id) || uniprot_id == "") {
    return(NA_character_)
  }
  
  # The current UniProt REST API URL for retrieving FASTA sequence
  url <- paste0("https://rest.uniprot.org/uniprotkb/", uniprot_id, ".fasta")
  
  tryCatch({
    response <- GET(url)
    
    if (status_code(response) == 200) {
      # Parse FASTA content
      fasta_content <- content(response, "text")
      
      # Skip the header line and combine remaining lines
      sequence_lines <- strsplit(fasta_content, "\n")[[1]]
      
      if (length(sequence_lines) < 2) {
        message(paste("No sequence data found for", uniprot_id))
        return(NA_character_)
      }
      
      # Remove the header (first line that starts with ">")
      sequence_lines <- sequence_lines[!grepl("^>", sequence_lines)]
      
      # Join sequence lines and remove any whitespace
      sequence <- gsub("\\s+", "", paste(sequence_lines, collapse = ""))
      
      if (nchar(sequence) == 0) {
        message(paste("Empty sequence for", uniprot_id))
        return(NA_character_)
      }
      
      return(sequence)
    } else {
      message(paste("Failed to retrieve data for", uniprot_id, 
                    "- Status code:", status_code(response)))
      return(NA_character_)
    }
  }, error = function(e) {
    message(paste("Error retrieving data for", uniprot_id, "-", e$message))
    return(NA_character_)
  })
}

# Add a progress indicator since API calls can take time
cat("Retrieving canonical sequences from UniProt...\n")

# Create a new dataframe with results
result_data <- filtered_data %>%
  mutate(canonical_sequence = NA_character_,
         sequence_retrieved = FALSE)

# Function to add delay between API calls
add_delay <- function(delay_seconds = 0.5) {
  Sys.sleep(delay_seconds)
}

# Process in smaller batches to avoid overwhelming the API and implement rate limiting
batch_size <- 5  # Smaller batch size to be more conservative
total_batches <- ceiling(nrow(result_data) / batch_size)

for (i in 1:total_batches) {
  start_idx <- ((i-1) * batch_size) + 1
  end_idx <- min(i * batch_size, nrow(result_data))
  
  if (start_idx <= nrow(result_data)) {
    cat(paste0("Processing batch ", i, " of ", total_batches, " (rows ", start_idx, "-", end_idx, ")\n"))
    
    for (j in start_idx:end_idx) {
      if (j <= nrow(result_data)) {  # Ensure we don't go out of bounds
        uniprot_id <- result_data[j, uniprot_col_name]
        
        # Skip empty or NA UniProt IDs
        if (!is.na(uniprot_id) && nchar(uniprot_id) > 0) {
          cat(paste0("  Retrieving sequence for ", uniprot_id, "... "))
          
          seq_result <- get_canonical_sequence(uniprot_id)
          
          if (!is.na(seq_result)) {
            result_data$canonical_sequence[j] <- seq_result
            result_data$sequence_retrieved[j] <- TRUE
            cat("Success! Length:", nchar(seq_result), "characters\n")
          } else {
            cat("Failed to retrieve sequence.\n")
          }
          
          # Add a delay between API calls to respect rate limits
          add_delay(1)  # 1 second delay between calls
        } else {
          cat(paste0("  Skipping empty UniProt ID at row ", j, "\n"))
        }
      }
    }
    
    # Add a longer delay between batches
    if(i < total_batches) {
      cat("Pausing between batches to respect API rate limits...\n")
      add_delay(3)  # 3 seconds delay between batches
    }
  }
}

# Count proteins with retrieved sequences
proteins_with_sequences <- sum(result_data$sequence_retrieved, na.rm = TRUE)

# Document the count
cat("\nSummary:\n")
cat("Proteins with retrieved canonical sequences:", proteins_with_sequences, "\n")

# Check if we have any successful retrievals
if(proteins_with_sequences == 0) {
  cat("\nWARNING: No sequences were successfully retrieved. Please check your UniProt IDs and network connection.\n")
}

# Save the filtered data with canonical sequences
write.csv(result_data, "MOLT4_mutations_with_sequences.csv", row.names = FALSE)

# Create a summary dataframe for documentation
summary_data <- data.frame(
  Metric = c("Total entries", 
             "Proteins with UniProt IDs", 
             "Proteins with missense variants",
             "Proteins with retrieved canonical sequences"),
  Count = c(total_entries, 
            entries_with_uniprot, 
            entries_with_missense,
            proteins_with_sequences)
)

# Save the summary
write.csv(summary_data, "MOLT4_analysis_summary.csv", row.names = FALSE)

# Display the summary table
knitr::kable(summary_data, caption = "MOLT4 Mutations Analysis Summary")

# =========================================================================
# PART 2: Process the sequences and apply mutations
# =========================================================================

# Load the data (in case you're running this part separately)
# Comment this out if running the script continuously
# data <- read.csv("MOLT4_mutations_with_sequences.csv", stringsAsFactors = FALSE)

# Function to apply a mutation to a sequence
apply_mutation <- function(sequence, mutation) {
  # Remove 'p.' prefix if present
  mutation <- gsub("^p\\.", "", mutation)
  
  # Parse mutation (e.g., P750Q)
  pattern <- "([A-Z])(\\d+)([A-Z])"
  match_result <- stringr::str_match(mutation, pattern)
  
  if (is.null(match_result) || nrow(match_result) == 0 || any(is.na(match_result))) {
    return(NA)
  }
  
  original_aa <- match_result[1, 2]
  position <- as.numeric(match_result[1, 3])
  new_aa <- match_result[1, 4]
  
  # Check for NA values in sequence or position
  if (is.na(sequence) || is.na(position)) {
    return(NA)
  }
  
  # Validate position and sequence length
  if (position > nchar(sequence) || position < 1) {
    return(NA)
  }
  
  # Validate that the character at the position matches the expected original amino acid
  current_aa <- substring(sequence, position, position)
  if (current_aa != original_aa) {
    warning(paste(
      "Mismatch at position", position,
      "Expected:", original_aa,
      "Found:", current_aa,
      "in mutation:", mutation
    ))
    return(NA)
  }
  
  # Apply mutation
  result <- sequence
  substr(result, position, position) <- new_aa
  return(result)
}

# Function to reverse a sequence
reverse_sequence <- function(sequence) {
  if(is.na(sequence)) return(NA)
  chars <- strsplit(sequence, "")[[1]]
  reversed <- paste0(rev(chars), collapse = "")
  return(reversed)
}

# First, create a clean version of Protein.Change
processed_data <- result_data %>%
  filter(
    !!sym(variant_col_name) == "missense_variant",
    !is.na(canonical_sequence)
  ) %>%
  mutate(
    Clean_Protein_Change = gsub("^p\\.", "", Protein.Change)
  )

# Next, apply mutations one by one to avoid errors
processed_data$mutated_sequence <- NA
for(i in 1:nrow(processed_data)) {
  processed_data$mutated_sequence[i] <- apply_mutation(
    processed_data$canonical_sequence[i], 
    processed_data$Protein.Change[i]
  )
}

# Then continue with the rest of the processing
processed_data <- processed_data %>%
  mutate(
    # Reverse original and mutated sequences
    wt_reversed = sapply(canonical_sequence, reverse_sequence),
    mutated_reversed = sapply(mutated_sequence, reverse_sequence),
    # Create tags and headers
    Variant_Tag = paste0(Gene, "_", Clean_Protein_Change),
    GN = paste0("GN=", Gene, Clean_Protein_Change),
    Organism = "OS=Homo sapiens",
    Header_Fwd = paste0(">Fwd_sp", Clean_Protein_Change, "|", !!sym(uniprot_col_name), "|", Variant_Tag, " ", Organism, " ", GN),
    Header_Rev = paste0(">Rev_sp", Clean_Protein_Change, "|", !!sym(uniprot_col_name), "|", Variant_Tag, " ", Organism, " ", GN)
  )

# Show summary of processed data
summary_data <- processed_data %>%
  select(Gene, Clean_Protein_Change, canonical_sequence, mutated_sequence)
head(summary_data)

# Count successful mutations
successful_mutations <- sum(!is.na(processed_data$mutated_sequence))
cat("Successfully processed", successful_mutations, "out of", nrow(processed_data), "missense variants\n")

# Save output to a CSV file
output_file <- "MOLT4_mutated_protein_output.csv"
write_csv(processed_data, output_file)
cat("Output saved to", output_file, "\n")

# Show final output columns
final_output <- processed_data %>%
  select(Gene, Clean_Protein_Change, canonical_sequence, mutated_sequence, 
         wt_reversed, mutated_reversed, !!sym(uniprot_col_name), Organism, GN)
head(final_output)

# Function to write in FASTA format
write_combined_fasta <- function(data, output_file) {
  # Open file connection
  file_conn <- file(output_file, "w")
  
  # Process each row
  for (i in 1:nrow(data)) {
    # Skip entries with missing sequences
    if(is.na(data$mutated_sequence[i])) next
    
    # Forward mutated sequence
    writeLines(data$Header_Fwd[i], file_conn)
    writeLines(data$mutated_sequence[i], file_conn)
    
    # Reverse mutated sequence
    writeLines(data$Header_Rev[i], file_conn)
    writeLines(data$mutated_reversed[i], file_conn)
  }
  
  # Close file connection
  close(file_conn)
}

# Output FASTA file
output_fasta <- "MOLT4_mutated_protein_database.fasta"

# Generate FASTA file
write_combined_fasta(processed_data, output_fasta)
cat("FASTA file saved to", output_fasta, "\n")



# =========================================================================
# PART 3: Create Peptide Database from Mutated Proteins - CORRECTED VERSION
# =========================================================================

# CORRECTED Function to find peptide sequence around mutation site
# This version INCLUDES the K/R cleavage sites in the peptide
extract_peptide_around_mutation <- function(sequence, mutation, protein_name = "", is_reverse = FALSE, original_seq_length = NULL) {
  # Check for valid inputs
  if (is.na(sequence) || is.na(mutation) || sequence == "" || mutation == "") {
    return(list(peptide = NA, start_pos = NA, end_pos = NA, 
                mutation_pos_in_peptide = NA, original_mutation_pos = NA))
  }
  
  # Ensure sequence is character
  sequence <- as.character(sequence)
  mutation <- as.character(mutation)
  
  # Remove 'p.' prefix if present
  mutation <- gsub("^p\\.", "", mutation)
  
  # Parse mutation to get position
  pattern <- "([A-Z])(\\d+)([A-Z])"
  match_result <- stringr::str_match(mutation, pattern)
  
  if (is.null(match_result) || nrow(match_result) == 0 || any(is.na(match_result))) {
    warning(paste("Could not parse mutation:", mutation, "for protein:", protein_name))
    return(list(peptide = NA, start_pos = NA, end_pos = NA, 
                mutation_pos_in_peptide = NA, original_mutation_pos = NA))
  }
  
  original_mutation_position <- as.numeric(match_result[1, 3])
  
  # For reverse sequences, calculate the mutation position in the reversed sequence
  if (is_reverse && !is.null(original_seq_length)) {
    mutation_position <- original_seq_length - original_mutation_position + 1
  } else {
    mutation_position <- original_mutation_position
  }
  
  seq_length <- nchar(sequence)
  
  if (is.na(mutation_position) || mutation_position < 1 || mutation_position > seq_length) {
    warning(paste("Invalid mutation position:", mutation_position, "for protein:", protein_name, 
                  "Sequence length:", seq_length, "Is reverse:", is_reverse))
    return(list(peptide = NA, start_pos = NA, end_pos = NA, 
                mutation_pos_in_peptide = NA, original_mutation_pos = original_mutation_position))
  }
  
  # Initialize boundaries to sequence limits
  left_start <- 1
  right_end <- seq_length
  k_r_count_left <- 0
  k_r_count_right <- 0
  
  # Find 2nd K or R to the left of mutation site
  # We want to INCLUDE this K/R in our peptide
  if (mutation_position > 1) {
    search_range <- seq(mutation_position - 1, 1, by = -1)
    for (i in search_range) {
      current_aa <- substr(sequence, i, i)
      if (current_aa %in% c("K", "R")) {
        k_r_count_left <- k_r_count_left + 1
        if (k_r_count_left == 2) {
          left_start <- i  # INCLUDE the 2nd K/R (changed from i + 1)
          break
        }
      }
    }
  }
  
  # Find 2nd K or R to the right of mutation site  
  # We want to INCLUDE this K/R in our peptide
  if (mutation_position < seq_length) {
    search_range <- seq(mutation_position + 1, seq_length, by = 1)
    for (i in search_range) {
      current_aa <- substr(sequence, i, i)
      if (current_aa %in% c("K", "R")) {
        k_r_count_right <- k_r_count_right + 1
        if (k_r_count_right == 2) {
          right_end <- i  # INCLUDE the 2nd K/R (changed from i - 1)
          break
        }
      }
    }
  }
  
  # Validate boundaries
  if (left_start > right_end || left_start < 1 || right_end > seq_length) {
    warning(paste("Invalid peptide boundaries for protein:", protein_name, 
                  "Left start:", left_start, "Right end:", right_end, "Seq length:", seq_length,
                  "Is reverse:", is_reverse))
    return(list(peptide = NA, start_pos = NA, end_pos = NA, 
                mutation_pos_in_peptide = NA, original_mutation_pos = original_mutation_position))
  }
  
  # Extract peptide sequence
  peptide <- substr(sequence, left_start, right_end)
  
  # Calculate mutation position within peptide
  adjusted_mutation_pos <- mutation_position - left_start + 1
  
  # Validate that the mutation position is within the peptide
  if (adjusted_mutation_pos < 1 || adjusted_mutation_pos > nchar(peptide)) {
    warning(paste("Mutation position not within extracted peptide for protein:", protein_name,
                  "Adjusted pos:", adjusted_mutation_pos, "Peptide length:", nchar(peptide),
                  "Is reverse:", is_reverse))
    return(list(peptide = NA, start_pos = NA, end_pos = NA, 
                mutation_pos_in_peptide = NA, original_mutation_pos = original_mutation_position))
  }
  
  # Return successful result
  return(list(
    peptide = peptide,
    start_pos = left_start,
    end_pos = right_end,
    mutation_pos_in_peptide = adjusted_mutation_pos,
    original_mutation_pos = original_mutation_position,
    mutation_pos_in_sequence = mutation_position
  ))
}

# Apply peptide extraction to all mutated sequences (both forward and reverse)
cat("Extracting CORRECTED peptides around mutation sites (including K/R cleavage sites)...\n")

# Create peptide data for both forward and reverse sequences
peptide_data <- processed_data %>%
  filter(!is.na(mutated_sequence), !is.na(mutated_reversed)) %>%
  mutate(
    # Get original sequence length for reverse calculations
    original_seq_length = nchar(mutated_sequence),
    # Extract peptides from forward sequences
    peptide_info_fwd = map2(mutated_sequence, Protein.Change, 
                            ~extract_peptide_around_mutation(.x, .y, paste0(Gene, "_forward"), 
                                                             is_reverse = FALSE)),
    # Extract peptides from reverse sequences with correct position calculation
    peptide_info_rev = pmap(list(mutated_reversed, Protein.Change, original_seq_length), 
                            ~extract_peptide_around_mutation(..1, ..2, paste0(Gene, "_reverse"), 
                                                             is_reverse = TRUE, 
                                                             original_seq_length = ..3))
  )

# Extract forward peptide information
peptide_data <- peptide_data %>%
  mutate(
    # Forward peptide info
    has_valid_peptide_fwd = map_lgl(peptide_info_fwd, ~!is.na(.x$peptide)),
    peptide_sequence_fwd = map_chr(peptide_info_fwd, ~ifelse(is.na(.x$peptide), NA_character_, .x$peptide)),
    peptide_start_fwd = map_dbl(peptide_info_fwd, ~ifelse(is.na(.x$start_pos), NA_real_, .x$start_pos)),
    peptide_end_fwd = map_dbl(peptide_info_fwd, ~ifelse(is.na(.x$end_pos), NA_real_, .x$end_pos)),
    mutation_pos_in_peptide_fwd = map_dbl(peptide_info_fwd, ~ifelse(is.na(.x$mutation_pos_in_peptide), NA_real_, .x$mutation_pos_in_peptide)),
    
    # Reverse peptide info
    has_valid_peptide_rev = map_lgl(peptide_info_rev, ~!is.na(.x$peptide)),
    peptide_sequence_rev = map_chr(peptide_info_rev, ~ifelse(is.na(.x$peptide), NA_character_, .x$peptide)),
    peptide_start_rev = map_dbl(peptide_info_rev, ~ifelse(is.na(.x$start_pos), NA_real_, .x$start_pos)),
    peptide_end_rev = map_dbl(peptide_info_rev, ~ifelse(is.na(.x$end_pos), NA_real_, .x$end_pos)),
    mutation_pos_in_peptide_rev = map_dbl(peptide_info_rev, ~ifelse(is.na(.x$mutation_pos_in_peptide), NA_real_, .x$mutation_pos_in_peptide))
  )

# Create headers with the specified format: Fwd/Rev_spVariantTag | Uniprot ID | Organism (OS) Gene Name (GN)
peptide_data <- peptide_data %>%
  mutate(
    # Forward peptide header
    Peptide_Header_Fwd = paste0("Fwd_sp", Clean_Protein_Change, "|", !!sym(uniprot_col_name), "|", Variant_Tag, " ", Organism, " GN=", Gene),
    # Reverse peptide header  
    Peptide_Header_Rev = paste0("Rev_sp", Clean_Protein_Change, "|", !!sym(uniprot_col_name), "|", Variant_Tag, " ", Organism, " GN=", Gene)
  )

# Filter for entries with at least one valid peptide (forward or reverse)
valid_peptides <- peptide_data %>%
  filter(has_valid_peptide_fwd == TRUE | has_valid_peptide_rev == TRUE)

# Count successful peptide extractions
successful_peptides_fwd <- sum(peptide_data$has_valid_peptide_fwd, na.rm = TRUE)
successful_peptides_rev <- sum(peptide_data$has_valid_peptide_rev, na.rm = TRUE)
total_peptide_entries <- successful_peptides_fwd + successful_peptides_rev

cat("Successfully extracted", successful_peptides_fwd, "forward peptides\n")
cat("Successfully extracted", successful_peptides_rev, "reverse peptides\n")
cat("Total peptide entries:", total_peptide_entries, "\n")

# Function to write peptide FASTA file with both forward and reverse sequences
write_peptide_fasta <- function(data, output_file) {
  # Open file connection
  file_conn <- file(output_file, "w")
  
  # Process each row
  for (i in 1:nrow(data)) {
    # Write forward peptide if valid
    if(!is.na(data$peptide_sequence_fwd[i]) && data$has_valid_peptide_fwd[i]) {
      writeLines(paste0(">", data$Peptide_Header_Fwd[i]), file_conn)
      writeLines(data$peptide_sequence_fwd[i], file_conn)
    }
    
    # Write reverse peptide if valid
    if(!is.na(data$peptide_sequence_rev[i]) && data$has_valid_peptide_rev[i]) {
      writeLines(paste0(">", data$Peptide_Header_Rev[i]), file_conn)
      writeLines(data$peptide_sequence_rev[i], file_conn)
    }
  }
  
  # Close file connection
  close(file_conn)
}

# Output peptide FASTA file
output_peptide_fasta <- "MOLT4_mutated_peptide_database.fasta"

# Generate peptide FASTA file
write_peptide_fasta(valid_peptides, output_peptide_fasta)
cat("CORRECTED Peptide FASTA file saved to", output_peptide_fasta, "\n")

# Save peptide data to CSV for reference
peptide_output_csv <- "MOLT4_peptide_data.csv"
peptide_summary <- valid_peptides %>%
  select(Gene, Clean_Protein_Change, !!sym(uniprot_col_name), 
         # Forward peptide info
         peptide_sequence_fwd, peptide_start_fwd, peptide_end_fwd, 
         mutation_pos_in_peptide_fwd, Peptide_Header_Fwd,
         # Reverse peptide info
         peptide_sequence_rev, peptide_start_rev, peptide_end_rev, 
         mutation_pos_in_peptide_rev, Peptide_Header_Rev,
         # Validity flags
         has_valid_peptide_fwd, has_valid_peptide_rev)

write_csv(peptide_summary, peptide_output_csv)
cat("Peptide data saved to", peptide_output_csv, "\n")

# Display summary statistics
cat("\n=== PEPTIDE DATABASE SUMMARY ===\n")
cat("Total mutated proteins:", nrow(processed_data), "\n")
cat("Forward peptide extractions:", successful_peptides_fwd, "\n")
cat("Reverse peptide extractions:", successful_peptides_rev, "\n")
cat("Total peptide entries in FASTA:", total_peptide_entries, "\n")
cat("Proteins with at least one valid peptide:", nrow(valid_peptides), "\n")

# Show some example peptides
if(nrow(valid_peptides) > 0) {
  cat("\nExample peptide headers and sequences (first 3):\n")
  
  for(i in 1:min(3, nrow(valid_peptides))) {
    if(valid_peptides$has_valid_peptide_fwd[i]) {
      cat(paste0("Forward: >", valid_peptides$Peptide_Header_Fwd[i], "\n"))
      cat(paste0("Sequence: ", valid_peptides$peptide_sequence_fwd[i], "\n"))
      cat(paste0("Mutation at position: ", valid_peptides$mutation_pos_in_peptide_fwd[i], "\n\n"))
    }
    
    if(valid_peptides$has_valid_peptide_rev[i]) {
      cat(paste0("Reverse: >", valid_peptides$Peptide_Header_Rev[i], "\n"))
      cat(paste0("Sequence: ", valid_peptides$peptide_sequence_rev[i], "\n"))
      cat(paste0("Mutation at position: ", valid_peptides$mutation_pos_in_peptide_rev[i], "\n\n"))
    }
  }
}

# Final summary
cat("\n=== FILES GENERATED ===\n")
cat("1. Protein sequences:", output_fasta, "\n")
cat("2. Peptide sequences:", output_peptide_fasta, "\n")
cat("3. Protein data CSV:", output_file, "\n")
cat("4. Peptide data CSV:", peptide_output_csv, "\n")
cat("5. Analysis summary:", "MOLT4_analysis_summary.csv", "\n")

# =========================================================================
# PART 4: VERIFICATION OF PEPTIDE TRUNCATIONS
# =========================================================================

cat("\n=== VERIFICATION OF PEPTIDE TRUNCATIONS ===\n")

# Fixed verification function that properly checks K/R boundaries
# CORRECTED verification function that matches the extraction logic
verify_peptide_extraction <- function(original_seq, peptide_seq, start_pos, end_pos, 
                                      mutation_pos, mutation_pos_in_peptide, 
                                      protein_change, sequence_type = "forward") {
  
  verification_results <- list(
    protein_change = protein_change,
    sequence_type = sequence_type,
    original_length = nchar(original_seq),
    peptide_length = nchar(peptide_seq),
    start_pos = start_pos,
    end_pos = end_pos,
    mutation_pos = mutation_pos,
    mutation_pos_in_peptide = mutation_pos_in_peptide,
    tests_passed = 0,
    total_tests = 0,
    issues = character(0)
  )
  
  # Test 1: Check if peptide matches extracted region from original sequence
  verification_results$total_tests <- verification_results$total_tests + 1
  extracted_from_original <- substr(original_seq, start_pos, end_pos)
  if (identical(extracted_from_original, peptide_seq)) {
    verification_results$tests_passed <- verification_results$tests_passed + 1
  } else {
    verification_results$issues <- c(verification_results$issues, 
                                     "Peptide doesn't match extracted region from original sequence")
  }
  
  # Test 2: Check if mutation position is within peptide bounds
  verification_results$total_tests <- verification_results$total_tests + 1
  if (mutation_pos_in_peptide >= 1 && mutation_pos_in_peptide <= nchar(peptide_seq)) {
    verification_results$tests_passed <- verification_results$tests_passed + 1
  } else {
    verification_results$issues <- c(verification_results$issues, 
                                     "Mutation position is outside peptide bounds")
  }
  
  # Test 3: Check if mutation position calculation is correct
  verification_results$total_tests <- verification_results$total_tests + 1
  expected_mut_pos_in_peptide <- mutation_pos - start_pos + 1
  if (mutation_pos_in_peptide == expected_mut_pos_in_peptide) {
    verification_results$tests_passed <- verification_results$tests_passed + 1
  } else {
    verification_results$issues <- c(verification_results$issues, 
                                     paste("Incorrect mutation position calculation. Expected:", 
                                           expected_mut_pos_in_peptide, "Got:", mutation_pos_in_peptide))
  }
  
  # Test 4: CORRECTED - Check for K/R boundaries according to ACTUAL extraction logic
  # The extraction function INCLUDES the 2nd K/R in the peptide, so we need to match that logic
  verification_results$total_tests <- verification_results$total_tests + 1
  k_r_boundary_check_passed <- TRUE
  
  # Check left boundary: verify that start_pos is correctly positioned relative to K/R
  if (start_pos > 1) {
    # Count K/R residues to the left of mutation position
    k_r_count_left <- 0
    expected_start <- 1  # Default if we don't find 2 K/R residues
    
    # Search from mutation position leftward
    for (i in seq(mutation_pos - 1, 1, by = -1)) {
      current_aa <- substr(original_seq, i, i)
      if (current_aa %in% c("K", "R")) {
        k_r_count_left <- k_r_count_left + 1
        if (k_r_count_left == 2) {
          expected_start <- i  # INCLUDE the 2nd K/R (matches extraction function)
          break
        }
      }
    }
    
    # Check if actual start position matches expected
    if (start_pos != expected_start) {
      verification_results$issues <- c(verification_results$issues, 
                                       paste("Left boundary mismatch. Expected start:", expected_start, 
                                             "Actual start:", start_pos))
      k_r_boundary_check_passed <- FALSE
    }
  }
  
  # Check right boundary: verify that end_pos is correctly positioned relative to K/R
  if (end_pos < nchar(original_seq)) {
    # Count K/R residues to the right of mutation position
    k_r_count_right <- 0
    expected_end <- nchar(original_seq)  # Default if we don't find 2 K/R residues
    
    # Search from mutation position rightward
    for (i in seq(mutation_pos + 1, nchar(original_seq), by = 1)) {
      current_aa <- substr(original_seq, i, i)
      if (current_aa %in% c("K", "R")) {
        k_r_count_right <- k_r_count_right + 1
        if (k_r_count_right == 2) {
          expected_end <- i  # INCLUDE the 2nd K/R (matches extraction function)
          break
        }
      }
    }
    
    # Check if actual end position matches expected
    if (end_pos != expected_end) {
      verification_results$issues <- c(verification_results$issues, 
                                       paste("Right boundary mismatch. Expected end:", expected_end, 
                                             "Actual end:", end_pos))
      k_r_boundary_check_passed <- FALSE
    }
  }
  
  if (k_r_boundary_check_passed) {
    verification_results$tests_passed <- verification_results$tests_passed + 1
  }
  
  # Test 5: Verify the mutation is actually present in the peptide
  verification_results$total_tests <- verification_results$total_tests + 1
  
  # Parse the mutation
  mutation_clean <- gsub("^p\\.", "", protein_change)
  pattern <- "([A-Z])(\\d+)([A-Z])"
  match_result <- stringr::str_match(mutation_clean, pattern)
  
  if (!is.na(match_result[1,4])) {
    new_aa <- match_result[1,4]
    actual_aa_in_peptide <- substr(peptide_seq, mutation_pos_in_peptide, mutation_pos_in_peptide)
    
    if (actual_aa_in_peptide == new_aa) {
      verification_results$tests_passed <- verification_results$tests_passed + 1
    } else {
      verification_results$issues <- c(verification_results$issues, 
                                       paste("Mutation not found in peptide. Expected:", new_aa, 
                                             "Found:", actual_aa_in_peptide, "at position", mutation_pos_in_peptide))
    }
  } else {
    verification_results$issues <- c(verification_results$issues, "Could not parse mutation")
  }
  
  return(verification_results)
}

# CORRECTED diagnostic function to show K/R positions around mutation
show_kr_diagnostic <- function(sequence, mutation_pos, protein_change, sequence_type = "forward") {
  cat(paste0("\n=== K/R DIAGNOSTIC for ", protein_change, " (", sequence_type, ") ===\n"))
  cat(paste0("Sequence length: ", nchar(sequence), "\n"))
  cat(paste0("Mutation position: ", mutation_pos, "\n"))
  
  # Show sequence around mutation (±20 positions)
  context_start <- max(1, mutation_pos - 20)
  context_end <- min(nchar(sequence), mutation_pos + 20)
  context_seq <- substr(sequence, context_start, context_end)
  
  cat(paste0("Context sequence (pos ", context_start, "-", context_end, "): ", context_seq, "\n"))
  
  # Find all K/R positions in the sequence
  kr_positions <- gregexpr("[KR]", sequence)[[1]]
  if (kr_positions[1] != -1) {
    cat("K/R positions in sequence: ")
    for (pos in kr_positions) {
      aa <- substr(sequence, pos, pos)
      cat(paste0(pos, "(", aa, ") "))
    }
    cat("\n")
    
    # Count K/R to the left and right of mutation
    left_kr <- kr_positions[kr_positions < mutation_pos]
    right_kr <- kr_positions[kr_positions > mutation_pos]
    
    cat(paste0("K/R positions left of mutation: ", length(left_kr), " - "))
    if (length(left_kr) > 0) {
      for (pos in left_kr) {
        aa <- substr(sequence, pos, pos)
        cat(paste0(pos, "(", aa, ") "))
      }
    }
    cat("\n")
    
    cat(paste0("K/R positions right of mutation: ", length(right_kr), " - "))
    if (length(right_kr) > 0) {
      for (pos in right_kr) {
        aa <- substr(sequence, pos, pos)
        cat(paste0(pos, "(", aa, ") "))
      }
    }
    cat("\n")
    
    # Determine expected boundaries based on CORRECTED extraction algorithm
    expected_start <- 1
    expected_end <- nchar(sequence)
    
    if (length(left_kr) >= 2) {
      expected_start <- left_kr[length(left_kr) - 1]  # INCLUDE 2nd K/R from right
    }
    
    if (length(right_kr) >= 2) {
      expected_end <- right_kr[2]  # INCLUDE 2nd K/R from left
    }
    
    cat(paste0("Expected peptide boundaries: ", expected_start, "-", expected_end, "\n"))
    expected_peptide <- substr(sequence, expected_start, expected_end)
    cat(paste0("Expected peptide: ", expected_peptide, "\n"))
    cat(paste0("Expected peptide length: ", nchar(expected_peptide), "\n"))
    
  } else {
    cat("No K/R residues found in sequence\n")
  }
}

# Run enhanced verification with diagnostics for problem cases
cat("\n=== ENHANCED VERIFICATION WITH DIAGNOSTICS ===\n")

# Re-run verification on the same sample
sample_size <- min(5, nrow(valid_peptides))  # Smaller sample for detailed diagnostics
enhanced_verification_results <- list()

set.seed(123)  # Same seed for reproducible results 
sample_indices <- sample(1:nrow(valid_peptides), sample_size)

for (i in sample_indices) {
  row_data <- valid_peptides[i, ]
  
  # Verify forward peptide if it exists
  if (row_data$has_valid_peptide_fwd && !is.na(row_data$peptide_sequence_fwd)) {
    original_mut_pos <- as.numeric(stringr::str_extract(gsub("^p\\.", "", row_data$Protein.Change), "\\d+"))
    
    fwd_result <- verify_peptide_extraction(
      original_seq = row_data$mutated_sequence,
      peptide_seq = row_data$peptide_sequence_fwd,
      start_pos = row_data$peptide_start_fwd,
      end_pos = row_data$peptide_end_fwd,
      mutation_pos = original_mut_pos,
      mutation_pos_in_peptide = row_data$mutation_pos_in_peptide_fwd,
      protein_change = row_data$Protein.Change,
      sequence_type = "forward"
    )
    enhanced_verification_results[[length(enhanced_verification_results) + 1]] <- fwd_result
    
    # Show diagnostic if there are issues
    if (length(fwd_result$issues) > 0) {
      show_kr_diagnostic(row_data$mutated_sequence, original_mut_pos, 
                         row_data$Protein.Change, "forward")
    }
  }
  
  # Verify reverse peptide if it exists
  if (row_data$has_valid_peptide_rev && !is.na(row_data$peptide_sequence_rev)) {
    original_seq_len <- nchar(row_data$mutated_sequence)
    original_mut_pos <- as.numeric(stringr::str_extract(gsub("^p\\.", "", row_data$Protein.Change), "\\d+"))
    reverse_mut_pos <- original_seq_len - original_mut_pos + 1
    
    rev_result <- verify_peptide_extraction(
      original_seq = row_data$mutated_reversed,
      peptide_seq = row_data$peptide_sequence_rev,
      start_pos = row_data$peptide_start_rev,
      end_pos = row_data$peptide_end_rev,
      mutation_pos = reverse_mut_pos,
      mutation_pos_in_peptide = row_data$mutation_pos_in_peptide_rev,
      protein_change = row_data$Protein.Change,
      sequence_type = "reverse"
    )
    enhanced_verification_results[[length(enhanced_verification_results) + 1]] <- rev_result
    
    # Show diagnostic if there are issues
    if (length(rev_result$issues) > 0) {
      show_kr_diagnostic(row_data$mutated_reversed, reverse_mut_pos, 
                         row_data$Protein.Change, "reverse")
    }
  }
}

# Summarize enhanced verification results
total_enhanced_verifications <- length(enhanced_verification_results)
total_enhanced_tests <- sum(sapply(enhanced_verification_results, function(x) x$total_tests))
total_enhanced_passed <- sum(sapply(enhanced_verification_results, function(x) x$tests_passed))

cat("\n=== ENHANCED VERIFICATION SUMMARY ===\n")
cat("Total peptides verified:", total_enhanced_verifications, "\n")
cat("Total tests performed:", total_enhanced_tests, "\n")
cat("Tests passed:", total_enhanced_passed, "\n")
cat("Overall success rate:", round(total_enhanced_passed/total_enhanced_tests * 100, 1), "%\n")

# Show remaining issues after fix
remaining_issues <- enhanced_verification_results[sapply(enhanced_verification_results, function(x) length(x$issues) > 0)]

if (length(remaining_issues) > 0) {
  cat("\n=== REMAINING ISSUES AFTER FIX ===\n")
  for (i in 1:length(remaining_issues)) {
    result <- remaining_issues[[i]]
    cat(paste0("Protein: ", result$protein_change, " (", result$sequence_type, ")\n"))
    cat(paste0("Tests passed: ", result$tests_passed, "/", result$total_tests, "\n"))
    for (issue in result$issues) {
      cat(paste0("  - ", issue, "\n"))
    }
    cat("\n")
  }
} else {
  cat("\n✓ All verified peptide extractions passed validation!\n")
}

# =========================================================================
# PART 5: VARIANT PEPTIDE LENGTH DISTRIBUTION ANALYSIS
# =========================================================================

# Load additional required libraries for visualization
library(ggplot2)
library(gridExtra)
library(RColorBrewer)
cat("\n=== VARIANT PEPTIDE LENGTH DISTRIBUTION ANALYSIS ===\n")
# Function to extract peptide from 2nd K/R to 3rd K/R around mutation site
extract_2nd_to_3rd_kr_peptide <- function(sequence, mutation, protein_name = "", is_reverse = FALSE, original_seq_length = NULL) {
  # Check for valid inputs
  if (is.na(sequence) || is.na(mutation) || sequence == "" || mutation == "") {
    return(list(peptide = NA, start_pos = NA, end_pos = NA, 
                mutation_pos_in_peptide = NA, original_mutation_pos = NA, length = NA))
  }
  
  # Ensure sequence is character
  sequence <- as.character(sequence)
  mutation <- as.character(mutation)
  
  # Remove 'p.' prefix if present
  mutation <- gsub("^p\\.", "", mutation)
  
  # Parse mutation to get position
  pattern <- "([A-Z])(\\d+)([A-Z])"
  match_result <- stringr::str_match(mutation, pattern)
  
  if (is.null(match_result) || nrow(match_result) == 0 || any(is.na(match_result))) {
    return(list(peptide = NA, start_pos = NA, end_pos = NA, 
                mutation_pos_in_peptide = NA, original_mutation_pos = NA, length = NA))
  }
  
  original_mutation_position <- as.numeric(match_result[1, 3])
  
  # For reverse sequences, calculate the mutation position in the reversed sequence
  if (is_reverse && !is.null(original_seq_length)) {
    mutation_position <- original_seq_length - original_mutation_position + 1
  } else {
    mutation_position <- original_mutation_position
  }
  
  seq_length <- nchar(sequence)
  
  if (is.na(mutation_position) || mutation_position < 1 || mutation_position > seq_length) {
    return(list(peptide = NA, start_pos = NA, end_pos = NA, 
                mutation_pos_in_peptide = NA, original_mutation_pos = original_mutation_position, length = NA))
  }
  
  # Find all K/R positions in the sequence
  kr_positions <- gregexpr("[KR]", sequence)[[1]]
  
  if (kr_positions[1] == -1 || length(kr_positions) < 4) {
    # Not enough K/R residues (need at least 4 total)
    return(list(peptide = NA, start_pos = NA, end_pos = NA, 
                mutation_pos_in_peptide = NA, original_mutation_pos = original_mutation_position, length = NA))
  }
  
  # Split K/R positions into those before and after mutation
  left_kr <- kr_positions[kr_positions < mutation_position]
  right_kr <- kr_positions[kr_positions > mutation_position]
  
  # We need at least 2 K/R on each side to get 2nd to 3rd
  if (length(left_kr) < 2 || length(right_kr) < 2) {
    return(list(peptide = NA, start_pos = NA, end_pos = NA, 
                mutation_pos_in_peptide = NA, original_mutation_pos = original_mutation_position, length = NA))
  }
  
  # Get 2nd K/R from the right (going leftward from mutation)
  second_kr_left <- left_kr[length(left_kr) - 1]
  
  # Get 2nd K/R from the left (going rightward from mutation)  
  second_kr_right <- right_kr[2]
  
  # Extract peptide from 2nd K/R on left to 2nd K/R on right (inclusive)
  peptide_start <- second_kr_left
  peptide_end <- second_kr_right
  
  # Extract peptide sequence
  peptide <- substr(sequence, peptide_start, peptide_end)
  
  # Calculate mutation position within peptide
  adjusted_mutation_pos <- mutation_position - peptide_start + 1
  
  # Validate that the mutation position is within the peptide
  if (adjusted_mutation_pos < 1 || adjusted_mutation_pos > nchar(peptide)) {
    return(list(peptide = NA, start_pos = NA, end_pos = NA, 
                mutation_pos_in_peptide = NA, original_mutation_pos = original_mutation_position, length = NA))
  }
  
  # Return successful result
  return(list(
    peptide = peptide,
    start_pos = peptide_start,
    end_pos = peptide_end,
    mutation_pos_in_peptide = adjusted_mutation_pos,
    original_mutation_pos = original_mutation_position,
    mutation_pos_in_sequence = mutation_position,
    length = nchar(peptide)
  ))
}
# Apply the 2nd-to-3rd K/R peptide extraction to all sequences
cat("Extracting peptides from 2nd K/R to 3rd K/R around mutation sites...\n")
# Create peptide data for both forward and reverse sequences
variant_peptide_data <- processed_data %>%
  filter(!is.na(mutated_sequence), !is.na(mutated_reversed)) %>%
  mutate(
    # Get original sequence length for reverse calculations
    original_seq_length = nchar(mutated_sequence),
    # Extract variant peptides from forward sequences (2nd to 3rd K/R)
    variant_peptide_info_fwd = map2(mutated_sequence, Protein.Change, 
                                    ~extract_2nd_to_3rd_kr_peptide(.x, .y, paste0(Gene, "_forward"), 
                                                                   is_reverse = FALSE)),
    # Extract variant peptides from reverse sequences
    variant_peptide_info_rev = pmap(list(mutated_reversed, Protein.Change, original_seq_length), 
                                    ~extract_2nd_to_3rd_kr_peptide(..1, ..2, paste0(Gene, "_reverse"), 
                                                                   is_reverse = TRUE, 
                                                                   original_seq_length = ..3))
  )
# Extract peptide length information
variant_peptide_data <- variant_peptide_data %>%
  mutate(
    # Forward peptide info
    variant_peptide_length_fwd = map_dbl(variant_peptide_info_fwd, ~ifelse(is.na(.x$length), NA_real_, .x$length)),
    variant_peptide_seq_fwd = map_chr(variant_peptide_info_fwd, ~ifelse(is.na(.x$peptide), NA_character_, .x$peptide)),
    has_valid_variant_fwd = !is.na(variant_peptide_length_fwd),
    
    # Reverse peptide info
    variant_peptide_length_rev = map_dbl(variant_peptide_info_rev, ~ifelse(is.na(.x$length), NA_real_, .x$length)),
    variant_peptide_seq_rev = map_chr(variant_peptide_info_rev, ~ifelse(is.na(.x$peptide), NA_character_, .x$peptide)),
    has_valid_variant_rev = !is.na(variant_peptide_length_rev)
  )
# Collect all valid peptide lengths for visualization
all_peptide_lengths <- c()
# Add forward peptide lengths
forward_lengths <- variant_peptide_data$variant_peptide_length_fwd[!is.na(variant_peptide_data$variant_peptide_length_fwd)]
all_peptide_lengths <- c(all_peptide_lengths, forward_lengths)
# Add reverse peptide lengths  
reverse_lengths <- variant_peptide_data$variant_peptide_length_rev[!is.na(variant_peptide_data$variant_peptide_length_rev)]
all_peptide_lengths <- c(all_peptide_lengths, reverse_lengths)
# Create data frame for plotting with three categories
peptide_length_df <- data.frame(
  Length = all_peptide_lengths,
  Category = ifelse(all_peptide_lengths < 7, "Below_7", 
                    ifelse(all_peptide_lengths <= 50, "7_to_50", "Above_50"))
)
# Count statistics
total_variant_peptides <- length(all_peptide_lengths)
below_threshold <- sum(all_peptide_lengths < 7)
optimal_range <- sum(all_peptide_lengths >= 7 & all_peptide_lengths <= 50)
above_optimal <- sum(all_peptide_lengths > 50)
successful_forward <- sum(variant_peptide_data$has_valid_variant_fwd, na.rm = TRUE)
successful_reverse <- sum(variant_peptide_data$has_valid_variant_rev, na.rm = TRUE)
cat("Successfully extracted", successful_forward, "forward variant peptides\n")
cat("Successfully extracted", successful_reverse, "reverse variant peptides\n")
cat("Total variant peptide entries:", total_variant_peptides, "\n")
cat("Peptides below 7 AA threshold:", below_threshold, "(", round(below_threshold/total_variant_peptides*100, 1), "%)\n")
cat("Peptides in optimal range (7-50 AA):", optimal_range, "(", round(optimal_range/total_variant_peptides*100, 1), "%)\n")
cat("Peptides above optimal range (>50 AA):", above_optimal, "(", round(above_optimal/total_variant_peptides*100, 1), "%)\n")
# Create histogram visualization
if (total_variant_peptides > 0) {
  
  # Create the main histogram
  p1 <- ggplot(peptide_length_df, aes(x = Length, fill = Category)) +
    geom_histogram(binwidth = 1, color = "black", alpha = 0.7) +
    geom_vline(xintercept = 7, color = "red", linetype = "dashed", size = 1.2) +
    geom_vline(xintercept = 50, color = "orange", linetype = "dashed", size = 1.2) +
    annotate("rect", xmin = 7, xmax = 50, ymin = -Inf, ymax = Inf, 
             alpha = 0.1, fill = "green") +
    scale_fill_manual(values = c("Below_7" = "#FF6B6B", 
                                 "7_to_50" = "#4ECDC4", 
                                 "Above_50" = "#FFB347"),
                      labels = c("Below_7" = "< 7 AA (Below Threshold)", 
                                 "7_to_50" = "7-50 AA (Optimal Range)",
                                 "Above_50" = "> 50 AA (Above Optimal)"),
                      name = "Peptide Category") +
    # Keep x-axis range from 0 to 100
    xlim(0, 100) +
    labs(title = "Distribution of Variant Peptide Lengths",
         subtitle = paste0("Length of Tryptic Peptide Flanking the Variant Site in Custom Database\n",
                           "Total peptides: ", total_variant_peptides, 
                           " | Below threshold: ", below_threshold, " (", round(below_threshold/total_variant_peptides*100, 1), "%)",
                           " | Optimal range: ", optimal_range, " (", round(optimal_range/total_variant_peptides*100, 1), "%)",
                           " | Above optimal: ", above_optimal, " (", round(above_optimal/total_variant_peptides*100, 1), "%)"),
         x = "Peptide Length (amino acids)",
         y = "Count") +
    theme_minimal() +
    theme(plot.title = element_text(size = 14, face = "bold"),
          plot.subtitle = element_text(size = 11),
          axis.title = element_text(size = 12),
          legend.title = element_text(size = 11),
          legend.position = "bottom")
  
  
  
  # Display the plots
  print(p1)
  cat("\n")

  
  # Save the plots
  ggsave("MOLT4_variant_peptide_length_histogram.png", plot = p1, width = 12, height = 8, dpi = 300)
}
