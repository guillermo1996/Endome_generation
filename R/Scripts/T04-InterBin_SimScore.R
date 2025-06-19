## _________________________________________________
##
## Inter-bin Similarity Score
##
## Aim:
##
## Author: Guillermo Rocamora Pérez
##
## Contributors:
##
## Date Created: 2025-05-20
##
## Copyright (c) Guillermo Rocamora Pérez, year
##
## Email: guillermorocamora@gmail.com
##
## Latest version: v1.0 (2025-05-20)
##
## _________________________________________________
##
## - Notes:
##
## - Changelog:
##
## - Please contact guillermorocamora@gmail.com for further assistance.
## _________________________________________________

# ORF_Similarity.R

# Required libraries
library(Biostrings)        # For sequence alignment
library(BiocParallel)      # For parallel processing
library(data.table)        # For efficient data handling
library(RSQLite)           # For database storage

# Configuration
N_CORES <- min(parallel::detectCores() - 1, 32)  # Limit to 32 cores max
CHUNK_SIZE <- 100            # Process 1M comparisons per chunk
TEMP_DIR <- "temp_results"    # Directory for temporary results
OUTPUT_DB <- "orf_similarity.sqlite"  # Final output database

# Create temp directory if it doesn't exist
if (!dir.exists(TEMP_DIR)) {
  dir.create(TEMP_DIR)
}


# Set up parallel backend
register(MulticoreParam(workers = N_CORES,
                        tasks = 0,  # 0 means one task per worker
                        progressbar = TRUE,
                        stop.on.error = FALSE))

# Function to generate pairs in chunks
generate_pairs_chunk <- function(n, chunk_size, chunk_num) {
  # Calculate total number of pairs
  total_pairs <- choose(n, 2)

  # Calculate start and end indices for this chunk (1-based)
  start_idx <- (chunk_num - 1) * chunk_size + 1
  end_idx <- min(chunk_num * chunk_size, total_pairs)
  chunk_size_actual <- end_idx - start_idx + 1

  # Pre-allocate matrix for results
  pairs <- matrix(0L, nrow = chunk_size_actual, ncol = 2)

  # Fill the matrix with (i,j) pairs
  pos <- 1
  current_k <- 1  # Current linear index

  # Iterate through all possible i,j pairs
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      # Check if this pair is in our chunk
      if (current_k >= start_idx && current_k <= end_idx) {
        pairs[pos, ] <- c(i, j)
        pos <- pos + 1
      }
      current_k <- current_k + 1

      # Early exit if we've passed our chunk
      if (current_k > end_idx) break
    }
    if (current_k > end_idx) break
  }

  return(pairs)
}

calculate_transcript_similarity <- function(pair, orf_list, transcripts) {
  i <- pair[1]
  j <- pair[2]
  tx1 <- transcripts[i]
  tx2 <- transcripts[j]

  orfs1 <- orf_list[[tx1]]
  orfs2 <- orf_list[[tx2]]

  # Skip if either transcript has no ORFs
  if (length(orfs1) == 0 || length(orfs2) == 0) {
    return(data.table(
      transcript1 = tx1,
      transcript2 = tx2,
      similarity = NA_real_
    ))
  }

  seqsA <- orf_list[[tx1]]
  seqsB <- orf_list[[tx2]]
  sim <- measureInterSimilarity(seqsA, seqsB)

  # Return result as data.table
  data.table(
    transcript1 = tx1,
    transcript2 = tx2,
    similarity = sim
  )
}


# Function to process a chunk of comparisons
process_chunk <- function(chunk_num, n, chunk_size, transcripts, orf_list) {
  # Generate pairs for this chunk
  pairs <- generate_pairs_chunk(n, chunk_size, chunk_num)

  # Process all pairs in this chunk
  results <- rbindlist(
    lapply(1:nrow(pairs), function(k) {
      calculate_transcript_similarity(pairs[k, ], orf_list, transcripts)
    })
  )

  # Save results
  chunk_file <- file.path(TEMP_DIR,
                          sprintf("chunk_%09d.csv", chunk_num))
  fwrite(results, file = chunk_file)
  return(chunk_file)
}

# Main function to process all comparisons
calculate_all_similarities <- function(orf_list) {
  # Get all transcript names
  transcripts <- names(orf_list)
  n <- length(transcripts)

  # Calculate total number of comparisons
  total_pairs <- choose(n, 2)
  message("Total pairs to process: ", format(total_pairs, big.mark = ","))

  # Calculate number of chunks needed
  n_chunks <- ceiling(total_pairs / CHUNK_SIZE)
  message("Processing in ", n_chunks, " chunks of size ", CHUNK_SIZE)

  # Process chunks in parallel
  chunk_files <- bplapply(1:n_chunks, function(chunk_num) {
    process_chunk(chunk_num, n, CHUNK_SIZE, transcripts, orf_list)
  })

  # Combine results
  combine_results(unlist(chunk_files), OUTPUT_DB)

  return(OUTPUT_DB)
}

# Function to combine chunk results into final database
combine_results <- function(chunk_files, output_db) {
  # Initialize database connection
  con <- dbConnect(RSQLite::SQLite(), output_db)

  # Create results table if it doesn't exist
  dbExecute(con, "
    CREATE TABLE IF NOT EXISTS orf_similarity (
      transcript1 TEXT,
      transcript2 TEXT,
      similarity REAL,
      PRIMARY KEY (transcript1, transcript2)
    ) WITHOUT ROWID;
  ")

  # Process each chunk file
  for (chunk_file in chunk_files) {
    # Import data using SQLite's import functionality
    dbExecute(con, sprintf("
      CREATE TEMP TABLE temp_import(
        transcript1 TEXT,
        transcript2 TEXT,
        similarity REAL
      );

      .mode csv
      .import '%s' temp_import

      INSERT OR IGNORE INTO orf_similarity
      SELECT * FROM temp_import;

      DROP TABLE temp_import;
    ", normalizePath(chunk_file)))

    # Clean up
    if (file.exists(chunk_file)) {
      file.remove(chunk_file)
    }
  }

  # Create index for faster queries
  dbExecute(con, "CREATE INDEX IF NOT EXISTS idx_transcript1 ON orf_similarity(transcript1)")
  dbExecute(con, "CREATE INDEX IF NOT EXISTS idx_transcript2 ON orf_similarity(transcript2)")

  # Close connection
  dbDisconnect(con)

  # Clean up temp directory if empty
  if (dir.exists(TEMP_DIR) && length(list.files(TEMP_DIR)) == 0) {
    unlink(TEMP_DIR, recursive = TRUE)
  }
}

# Example usage:
# 1. Load your ORF data into a named list where names are transcript IDs
#    and values are lists of ORF sequences (AAStringSet objects)
# orf_list <- list(
#   "tx1" = list(AAString("MSTRSG..."), AAString("MAG...")),
#   "tx2" = list(AAString("MSTRSG...")),
#   ...
# )
#
# 2. Run the analysis
# result_db <- calculate_all_similarities(orf_list)
#
# 3. Query results
# con <- dbConnect(RSQLite::SQLite(), result_db)
# results <- dbGetQuery(con, "SELECT * FROM orf_similarity WHERE similarity > 0.8")
# dbDisconnect(con)
