################################################################################
## 05 - ENDome truncation
################################################################################

## Variables
################################################################################

### Configurations
gtxcutr_settings = config["gtxcutr_settings"][config["gtxcutr_profile"]]

step05_params = {
    **gtxcutr_settings,
    **({"toy_data": True} if config.get("use_toy_data", False) else {})
}

### Compute the step hash
step05_hash = compute_hash(get_cumulative_params(5))

### Paths
truncation_path = lambda x: Path(results_path) / f"05-Truncation_{step04_hash}" / x
truncation_logs_path = lambda x: log_path(x, step04_hash)
truncation_benchmark_path = lambda x: benchmark_path(x, step04_hash)

## Functions
################################################################################
create_save_params_rule(step_num=5, step_name="05-Truncation", step_dir=truncation_path(""), step_params=step05_params)

## Rules
################################################################################
rule gtxcutr_truncation:
    message: """--- gTxcutr Truncation ---"""
    input: 
        gtf = rules.ORF_filtration.output.gtf
    output: 
        gtf = truncation_path("gtxcutr/{prefix}.{orf_filter}.gtxcutr.w{width}.{txEnd}.gtf"),
        fa = truncation_path("gtxcutr/{prefix}.{orf_filter}.gtxcutr.w{width}.{txEnd}.fa.gz"),
        transcript_overlap = truncation_path("gtxcutr/{prefix}.{orf_filter}.gtxcutr.w{width}.{txEnd}.overlaps.tsv"),
        merge_table = truncation_path("gtxcutr/{prefix}.{orf_filter}.gtxcutr.w{width}.{txEnd}.merge.tsv")
    log: truncation_logs_path("gtxcutr/{prefix}.{orf_filter}.gtxcutr.w{width}.{txEnd}.log")
    benchmark: truncation_benchmark_path("gtxcutr/{prefix}.{orf_filter}.gtxcutr.w{width}.{txEnd}.tsv")
    params:
        mergeDist = step05_params["merge_distance"],
        genome = step05_params["genome"]
    conda: "../envs/gtxcutr.yaml"
    threads: 12
    script: "../scripts/05a-gtxcutr.R" # Modified `txcutr.R` script to include my version of the package


# To Do: Check Tama Collapse approach to solve the cascade binning issue in gtxcutr: https://github.com/GenomeRIK/tama/wiki/Tama-Collapse