################################################################################
## 05 - ENDome truncation
################################################################################

## Variables
################################################################################
step05_name = "05-Truncation"

### Configurations
gtxcutr_settings = config["gtxcutr_settings"][config["gtxcutr_profile"]]

step05_params = {
    "gtxcutr_settings": gtxcutr_settings,
    **({"toy_data": True} if config.get("use_toy_data", False) else {})
}

### Compute hash and set global parameters
global_params.update({step05_name: step05_params})
step05_hash = compute_hash(global_params)
global_params.update({f"{step05_name}_{step05_hash}": global_params.pop(step05_name)})

### Paths
truncation_path = lambda x: Path(results_path) / f"{step05_name}-{step05_hash}" / x
truncation_log_path = lambda x: Path(results_path) / f"{step05_name}-{step05_hash}" / log_path / x
truncation_benchmark_path = lambda x: Path(results_path) / f"{step05_name}-{step05_hash}" / benchmark_path / x

## Functions
################################################################################
create_save_params_rule(step05_name, truncation_path, step05_params, global_params)

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
    log: truncation_log_path("gtxcutr/{prefix}.{orf_filter}.gtxcutr.w{width}.{txEnd}.log")
    benchmark: truncation_benchmark_path("gtxcutr/{prefix}.{orf_filter}.gtxcutr.w{width}.{txEnd}.tsv")
    params:
        mergeDist = gtxcutr_settings["merge_distance"],
        genome = gtxcutr_settings["genome"]
    conda: "../envs/gtxcutr.yaml"
    threads: 12
    script: "../scripts/05a-gtxcutr.R" # Modified `txcutr.R` script to include my version of the package


# To Do: Check Tama Collapse approach to solve the cascade binning issue in gtxcutr: https://github.com/GenomeRIK/tama/wiki/Tama-Collapse