################################################################################
## 05 - ENDome truncation
################################################################################

## Debug: benchmark loading time
_start_time = time.perf_counter()

## Variables
################################################################################
### Register the step: generates the hash suffix, creates the parameter rule and returns the output_path helpers.
step05 = register_step(
    name="05-Truncation",
    params={
        "gtxcutr_presets": build_tool_settings(config, "gtxcutr_settings", "gtxcutr_preset"),
    },
)

### Resolved settings used directly in the rule bodies below
gtxcutr_settings = resolve_preset(config, "gtxcutr_settings", "gtxcutr_preset")

## Rules
################################################################################
rule gtxcutr_truncation:
    message: """--- gTxcutr Truncation ---"""
    input:
        gtf = rules.ORF_filtration.output.gtf
    output:
        gtf = step05.path("gtxcutr/{prefix}.{orf_filter}.gtxcutr.w{width}.{txEnd}.gtf"),
        fa = step05.path("gtxcutr/{prefix}.{orf_filter}.gtxcutr.w{width}.{txEnd}.fa.gz"),
        transcript_overlap = step05.path("gtxcutr/{prefix}.{orf_filter}.gtxcutr.w{width}.{txEnd}.overlaps.tsv"),
        merge_table = step05.path("gtxcutr/{prefix}.{orf_filter}.gtxcutr.w{width}.{txEnd}.merge.tsv")
    log: step05.log("gtxcutr/{prefix}.{orf_filter}.gtxcutr.w{width}.{txEnd}.log")
    benchmark: step05.benchmark("gtxcutr/{prefix}.{orf_filter}.gtxcutr.w{width}.{txEnd}.tsv")
    params:
        mergeDist = gtxcutr_settings["merge_distance"],
        genome = gtxcutr_settings["genome"]
    conda: "../envs/r.yaml"
    threads: 12
    script: "../scripts/05a-gtxcutr.R" # Modified `txcutr.R` script to include my version of the package

rule kallisto_index:
    message: "--- Building kallisto index for the ENDome ---"
    input:
        fa = step05.path("gtxcutr/{prefix}.{orf_filter}.gtxcutr.w{width}.{txEnd}.fa.gz")
    output:
        kdx = step05.path("kallisto_index/{prefix}.{orf_filter}.w{width}.{txEnd}.kdx")
    log: step05.log("kallisto_index/{prefix}.{orf_filter}.w{width}.{txEnd}.log")
    benchmark: step05.benchmark("kallisto_index/{prefix}.{orf_filter}.w{width}.{txEnd}.tsv")
    conda: "../envs/scutrquant-kallisto-bustools.yaml"
    shell:
        "kallisto index -i {output.kdx} {input.fa} 2>&1 | tee {log}"

# To Do: Check Tama Collapse approach to solve the cascade binning issue in gtxcutr: https://github.com/GenomeRIK/tama/wiki/Tama-Collapse

## Debug: benchmark loading time
_log(f"\t+ {step05.name} imported in {time.perf_counter() - _start_time:.3f}s")
