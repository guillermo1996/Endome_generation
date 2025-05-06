list.files("/home/grocamora/RytenLab-Research/37-UTRome_pipeline/References/seqkit-subsample/seqkit_subsample_p1.0", full.names = T) %>%
  tibble::as_tibble() %>%
  dplyr::mutate(sample_id = gsub(".fastq.gz", "", gsub("subsample_", "", basename(value)))) %>%
  dplyr::mutate(file_type = "fastq") %>%
  dplyr::select(sample_id, file_type, files = value) %>%
  vroom::vroom_write("/home/grocamora/RytenLab-Research/38-Endome_generation/workflow/data/sample_sheet.csv", delim = ",")
