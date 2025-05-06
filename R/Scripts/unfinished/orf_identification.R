sq3_corrected <- rtracklayer::import.gff("/home/grocamora/RytenLab-Research/38-Endome_generation/results/Sqanti3_QC/sq3.annotated_corrected.gtf")
sq3_corrected_cds <- rtracklayer::import.gff("/home/grocamora/RytenLab-Research/38-Endome_generation/results/Sqanti3_QC/sq3.annotated_corrected.gtf.cds.gff")
ensembl_gtf <- rtracklayer::import.gff("/home/grocamora/RytenLab-Research/Additional_files/ENSEMBL/Homo_sapiens.GRCh38.113.gtf")
agat_test <- rtracklayer::import.gff("/home/grocamora/RytenLab-Research/38-Endome_generation/results/Sqanti3_Rescue/test_agat.gtf")
agat_test2 <- rtracklayer::import.gff("/home/grocamora/RytenLab-Research/38-Endome_generation/results/Sqanti3_Rescue/test_agat2.gtf")

rescued_gtf <- rtracklayer::import.gff("/home/grocamora/RytenLab-Research/38-Endome_generation/results/Sqanti3_Rescue/sq3.annotated_rescued.gtf")
sq3_classification <- vroom::vroom("/home/grocamora/RytenLab-Research/38-Endome_generation/results/Sqanti3_Filter/sq3.annotated_RulesFilter_result_classification.txt")



sq3_corrected[mcols(sq3_corrected)[["transcript_id"]] == "ENST00000641515.2"]
sq3_corrected_cds[mcols(sq3_corrected_cds)[["transcript_id"]] == "ENST00000641515.2"]
rescued_gtf[mcols(rescued_gtf)[["transcript_id"]] == "ENST00000641515.2"]
extractTranscriptSeqs(rescued_gtf[1:10])


txs_of_interest <- c("ENST00000641515.2", "ENST00000611770", "ENST00000403646", "ENST00000503018")
tx_study <- txs_of_interest[4]

agat_test %>% plyranges::filter(grepl(tx_study, transcript_id))
rescued_gtf %>% plyranges::filter(grepl(tx_study, transcript_id))
ensembl_gtf %>% plyranges::filter(grepl(tx_study, transcript_id))


a <- ensembl_gtf %>% plyranges::filter(grepl(txs_of_interest[2], transcript_id))
b <- ensembl_gtf %>% plyranges::filter(grepl(txs_of_interest[3], transcript_id))

agat_test %>% plyranges::filter(gene_id == "ENSG00000124343.14")
