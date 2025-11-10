################################################################################
## 05 - ENDome truncation
################################################################################

## Variables
################################################################################
truncation_path = lambda x: Path(results_output_path) / "05-Truncation" / x

## Functions
################################################################################

## Rules
################################################################################
rule gtxcutr_truncation:
    message: """--- gTxcutr Truncation ----"""
    input: 
        gtf = rules.ORF_filtration.output.gtf
    output: 
        gtf = truncation_path("gtxcutr/{prefix}.{orf_filter}.gtxcutr.w{width}.{txEnd}.gtf"),
        fa = truncation_path("gtxcutr/{prefix}.{orf_filter}.gtxcutr.w{width}.{txEnd}.fa.gz"),
        transcript_overlap = truncation_path("gtxcutr/{prefix}.{orf_filter}.gtxcutr.w{width}.{txEnd}.overlaps.tsv"),
        merge_table = truncation_path("gtxcutr/{prefix}.{orf_filter}.gtxcutr.w{width}.{txEnd}.merge.tsv")
    log: log_path("gtxcutr/{prefix}.{orf_filter}.gtxcutr.w{width}.{txEnd}.log")
    benchmark: benchmark_path("gtxcutr/{prefix}.{orf_filter}.gtxcutr.w{width}.{txEnd}.tsv")
    wildcard_constraints: 
        txEnd = "\\dp"
    params:
        mergeDist = 200,
        genome = "hg38"
    conda: "../envs/gtxcutr.yaml"
    threads: 12
    script: "../scripts/05a-gtxcutr.R" # Modified `txcutr.R` script to include my version of the package


# To Do: Check Tama Collapse approach to solve the cascade binning issue in gtxcutr: https://github.com/GenomeRIK/tama/wiki/Tama-Collapse