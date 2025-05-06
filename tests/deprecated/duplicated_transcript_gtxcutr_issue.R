gtf <- tibble::as_tibble(rtracklayer::import.gff("/home/grocamora/RytenLab-Research/38-Endome_generation/results/gtxcutr/sq3.annotated.gtxcutr.w500.gtf.gz"))

test_tx <- "ENST00000771700.1"
test0 <- tibble::as_tibble(rtracklayer::import.gff("/home/grocamora/RytenLab-Research/38-Endome_generation/results/test.gtf.gz"))
test1 <- tibble::as_tibble(rtracklayer::import.gff("/home/grocamora/RytenLab-Research/38-Endome_generation/results/Sqanti3_Rescue/sq3.annotated_rescued.gtf"))
test2 <- tibble::as_tibble(rtracklayer::import.gff("/home/grocamora/RytenLab-Research/37-UTRome_pipeline/homo_sapiens/sq3.annotated.txcutr.w500.gtf.gz"))
test3 <- tibble::as_tibble(rtracklayer::import.gff("/home/grocamora/RytenLab-Research/38-Endome_generation/results/AGAT/sq3.annotated.gtf"))
test4 <- tibble::as_tibble(rtracklayer::import.gff("/home/grocamora/RytenLab-Research/38-Endome_generation/results/ORF_Category/sq3.annotated_orf.gtf"))
test5 <- tibble::as_tibble(rtracklayer::import.gff("/home/grocamora/RytenLab-Research/38-Endome_generation/results/ORF_Filter/sq3.annotated_orf.filter.gtf"))

gtf %>% dplyr::filter(transcript_id == "ENST00000771700.1")
test0 %>% dplyr::filter(transcript_id == "ENST00000771700.1")
test1 %>% dplyr::filter(transcript_id == "ENST00000771700.1")
test2 %>% dplyr::filter(transcript_id == "ENST00000771700.1")
test3 %>% dplyr::filter(transcript_id == "ENST00000771700.1")
test4 %>% dplyr::filter(transcript_id == "ENST00000771700.1")
test5 %>% dplyr::filter(transcript_id == "ENST00000771700.1")

# My gtxcutr approach generates a phantom transcript assigned to another gene...
# And I don't know why. -> When transcripts are generated, the gene_id is also
# inferred. This is extracted from the mapTxToGene dictionary! However, this
# gene id might not be the same as the one from the exon lists, which kepts the
# gene_id from the original data. Thus, exons and transcripts does not have the
# same gene ID, so the program generates the redundancy. Suggested solution is
# maybe to remove gene id information from the exon, and let the program fill
# the voids. Default txcutr keeps the old gene ID, not sure how or why
#
# Please research why some transcripts are being pruned to 500 width (without
# reason!) -> Trasncripts are generated from the exons, so it is likely that
# monoexons only have 500 of width
#
# Why are there duplicated transcript_ids? ENST00000850701.1 -> Maybe because
# they do not represent transcripts, but rather exons.

test <- grTxs %>% plyranges::filter(grepl("ENST0000077170", transcript_id))

mcols(test)["type"] <- "transcript"
mcols(test)["gene_id"] <- mapTxToGene[as.character(test$transcript_id)]
mapTxToGene["ENST00000771700.1"]
