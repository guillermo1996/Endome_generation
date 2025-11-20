################################################################################
## ORF Identification
################################################################################

## Variables
################################################################################

### Configurations
orfannotate_settings = config["orfannotate_settings"][config["orfannotate_profile"]]
orf_filter_settings = config["orf_filter_settings"][config["orf_filter_profile"]]

step04_params = {
    **orfannotate_settings,
    **orf_filter_settings,
    **({"toy_data": True} if config.get("use_toy_data", False) else {})
}

### Compute the step hash
step04_hash = compute_hash(get_cumulative_params(4))

### Paths
orf_path = lambda x: Path(results_path) / f"04-ORF_Identification_{step04_hash}" / x
orf_logs_path = lambda x: log_path(x, step04_hash)
orf_benchmark_path = lambda x: benchmark_path(x, step04_hash)

### ORF Download path
orfannotate_version = step04_params["version"]
orfannotate_tar_url = f"https://github.com/egustavsson/ORFannotate/archive/refs/tags/{orfannotate_version}.tar.gz"
orfannotate_env_url = f"https://raw.githubusercontent.com/egustavsson/ORFannotate/refs/tags/{orfannotate_version}/ORFannotate.conda_env.yml"

## Functions
################################################################################
create_save_params_rule(step_num=4, step_name="04-ORF_Identification", step_dir=orf_path(""), step_params=step04_params)

def orfannotate_input(wildcards):
    if step04_params.get("gffread_clean", False):
        return rules.pigeon_filter.output.filtered_gtf
    else:
        return rules.gffread.output.gtf

## Rules
################################################################################
rule download_ORF_annotate:
    message: "--- Downloading and Extracting ORF Annotate ---"
    output:
        orfanntoate_dir = directory("tools/ORFannotate"),
        script = "tools/ORFannotate/ORFannotate.py"
    params:
        url = orfannotate_tar_url
    shell:
        """
        mkdir -p {output.orfanntoate_dir}
        curl -L {params.url} | tar -xz -C {output.orfanntoate_dir} --strip-components=1
        """

rule gffread: 
    message: """--- Cleaning gtf with gffread ---"""
    input: 
        gtf = artifact_path("Classify_Filter/{prefix}.sorted.filtered.gtf")
    output: 
        gtf = orf_path("gffread/{prefix}.gtf")
    log: orf_logs_path("ORF_Identification/gffread_{prefix}.log")
    benchmark: orf_benchmark_path("ORFannotate/gffread_{prefix}.tsv")
    conda: "../envs/gffread.yaml"
    shell:
        "gffread -E {input.gtf} -T -o {output.gtf} 2>&1 | tee {log}"

rule ORF_annotate:
    message: """--- ORF Identification with ORFannotate.py ---"""
    input:
        script = rules.download_ORF_annotate.output.script,
        gtf = orfannotate_input,
        ref_genome = ref_genome
    output: 
        summary_tsv = orf_path("ORFannotate/{prefix}/ORFannotate_summary.tsv"),
        gtf = orf_path("ORFannotate/{prefix}/ORFannotate_annotated.gtf"),
        protein_fa = orf_path("ORFannotate/{prefix}/protein.fa")
    log: orf_logs_path("ORFannotate/{prefix}.log")
    benchmark: orf_benchmark_path("ORFannotate/{prefix}.tsv")
    params:
        out_dir = lambda w, output: Path(output.summary_tsv).parent
    threads: 1
    conda: orfannotate_env_url
    shell:
        "python {input.script} --gtf {input.gtf} --fa {input.ref_genome} --outdir {params.out_dir} 2>&1 | tee {log}"

rule ORF_annotate_clean:
    message: """--- Clean ORFannotate output with gffread ---"""
    input:
        gtf = orf_path("ORFannotate/{prefix}/ORFannotate_annotated.gtf")
    output:
        gtf = orf_path("ORFannotate/{prefix}/ORFannotate_annotated_clean.gtf"),
    log: orf_logs_path("ORF_Identification/gffread_clean_{prefix}.log")
    benchmark: orf_benchmark_path("ORFannotate/gffread_clean_{prefix}.tsv")
    conda:  "../envs/gffread.yaml"
    shell:
        "gffread -E {input.gtf} -T -o {output.gtf} 2>&1 | tee {log}"

rule ORF_categorization:
    message: """--- Transcript categorization ---"""
    input:
        gtf = orf_path("ORFannotate/{prefix}/ORFannotate_annotated_clean.gtf"),
        ref_annotation = ref_annotation,
        protein_fa = orf_path("ORFannotate/{prefix}/protein.fa"),
        orf_summary = orf_path("ORFannotate/{prefix}/ORFannotate_summary.tsv"),
        pigeon_summary = artifact_path("Classify_Filter/{prefix}.pigeon_classification.filtered_lite_classification.txt")
    output:
        gtf = orf_path("ORF_Category/{prefix}_orf.gtf"),
        isoform_summary = orf_path("ORF_Category/{prefix}_isoform_summary.tsv"),
    log: orf_logs_path("ORF_Category/{prefix}_orf.log")
    benchmark: orf_benchmark_path("ORF_Category/{prefix}_orf.tsv")
    conda: "../envs/r.yaml"
    script: "../scripts/04a-ORF_Categorization.R"

rule ORF_filtration:
    message: """--- Transcript filtration ---"""
    input:
        gtf = rules.ORF_categorization.output.gtf
    output:
        gtf = orf_path("ORF_Filtration/{prefix}.{orf_filter}.orf_filter.gtf")
    log: orf_logs_path("ORF_Category/{prefix}.{orf_filter}_orf.log")
    benchmark: orf_benchmark_path("ORF_Category/{prefix}.{orf_filter}_orf.tsv")
    params:
        main_config = lambda wc: wc.orf_filter,
        valid_ref_gene_type = step04_params["valid_ref_gene_type"],
        valid_ref_tx_type = step04_params["valid_ref_tx_type"],
        valid_orfannotate_type = step04_params["valid_orfannotate_type"],
        in_ref_filter = step04_params["in_ref_filter"],
    conda: "../envs/r.yaml"
    script: "../scripts/04b-ORF_Filtration.R"