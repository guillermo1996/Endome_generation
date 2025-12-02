################################################################################
## ORF Identification
################################################################################

## Variables
################################################################################
step04_name = "04-ORF_Identification"

### Configurations
orfannotate_settings = config["orfannotate_settings"][config["orfannotate_profile"]]
orf_filter_settings = config["orf_filter_settings"][config["orf_filter_profile"]]

step04_params = {
    "orfannotate_settings": orfannotate_settings,
    "orf_filter_settings": orf_filter_settings,
    **({"toy_data": True} if config.get("use_toy_data", False) else {})
}

### Compute hash and set global parameters
global_params.update({step04_name: step04_params})
step04_hash = compute_hash(global_params)
global_params.update({f"{step04_name}_{step04_hash}": global_params.pop(step04_name)})

### Paths
orf_path = lambda x: Path(results_path) / f"{step04_name}_{step04_hash}" / x
orf_log_path = lambda x: Path(results_path) / f"{step04_name}_{step04_hash}" / log_path / x
orf_benchmark_path = lambda x: Path(results_path) / f"{step04_name}_{step04_hash}" / benchmark_path / x

### ORF Download path
orfannotate_version = orfannotate_settings["version"]
orfannotate_tar_url = f"https://github.com/egustavsson/ORFannotate/archive/refs/tags/{orfannotate_version}.tar.gz"
orfannotate_env_url = f"https://raw.githubusercontent.com/egustavsson/ORFannotate/refs/tags/{orfannotate_version}/ORFannotate.conda_env.yml"

## Functions
################################################################################
create_save_params_rule(step04_name, orf_path, step04_params, global_params)

def orfannotate_input(wildcards):
    if orfannotate_settings.get("gffread_clean", False):
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
    log: orf_log_path("ORF_Identification/gffread_{prefix}.log")
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
    log: orf_log_path("ORFannotate/{prefix}.log")
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
    log: orf_log_path("ORF_Identification/gffread_clean_{prefix}.log")
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
    log: orf_log_path("ORF_Category/{prefix}_orf.log")
    benchmark: orf_benchmark_path("ORF_Category/{prefix}_orf.tsv")
    conda: "../envs/r.yaml"
    script: "../scripts/04a-ORF_Categorization.R"


def select_source_dataset(wildcards):
    parts = wildcards.prefix.split(".")

    if parts[0] == "gencode":
        return ref_annotation
    else:
        return rules.ORF_categorization.output.gtf

rule ORF_filtration:
    message: """--- Transcript filtration ---"""
    input:
        gtf = select_source_dataset
    output:
        gtf = orf_path("ORF_Filtration/{prefix}.{orf_filter}.orf_filter.gtf")
    log: orf_log_path("ORF_Category/{prefix}.{orf_filter}_orf.log")
    benchmark: orf_benchmark_path("ORF_Category/{prefix}.{orf_filter}_orf.tsv")
    params:
        main_config = lambda wc: wc.orf_filter,
        valid_ref_gene_type = orf_filter_settings["valid_ref_gene_type"],
        valid_ref_tx_type = orf_filter_settings["valid_ref_tx_type"],
        valid_orfannotate_type = orf_filter_settings["valid_orfannotate_type"],
        in_ref_filter = orf_filter_settings["in_ref_filter"],
    conda: "../envs/r.yaml"
    script: "../scripts/04b-ORF_Filtration.R"