library(dplyr)
library(tidyr)
library(stringr)
library(purrr)

# Set the parent directory containing the sample folders
parent_dir <- "/home/grocamora/RytenLab-Research/Data/Lee2023"  # Change this to your actual path

# Toggle to include I2 files
include_I2 <- F

# Get all fastq.gz files recursively
fastq_files <- list.files(parent_dir, pattern = "fastq.gz$", full.names = TRUE, recursive = TRUE)

# Extract sample IDs from file paths
file_info <- tibble(
  file_path = fastq_files,
  folder = dirname(file_path),
  file_name = basename(file_path)
) %>%
  mutate(
    sample_id = str_extract(file_name, "SRR[0-9]+"),
    file_type = "fastq"
  )

# Filter out I2 files if needed
if (!include_I2) {
  file_info <- file_info %>% filter(!str_detect(file_name, "_I2.fastq.gz$"))
}

# Aggregate files per sample_id
sample_df <- file_info %>%
  group_by(sample_id, file_type) %>%
  summarise(files = paste(file_path, collapse = ";"), .groups = "drop")

# Print result
glimpse(sample_df)

# Save output
write.csv(sample_df, file.path(parent_dir, "sample_fastq_files.csv"), row.names = FALSE, quote = F)
