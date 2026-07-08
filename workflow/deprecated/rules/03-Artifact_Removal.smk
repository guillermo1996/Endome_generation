################################################################################
## Artifact Removal
################################################################################

## Variables
################################################################################
step03_name = "03-Artifact_Removal"

### Configurations
pigeon_settings = config["pigeon_settings"][config["pigeon_profile"]]

step03_params = {
    "pigeon_settings": pigeon_settings,
    **({"toy_data": True} if config.get("use_toy_data", False) else {})
}

### Compute hash and set global parameters
global_params.update({step03_name: step03_params})
step03_hash = compute_hash(global_params)
global_params.update({f"{step03_name}_{step03_hash}": global_params.pop(step03_name)})

### Paths
artifact_path = lambda x: Path(results_path) / f"{step03_name}-{step03_hash}" / x
artifact_log_path = lambda x: Path(results_path) / f"{step03_name}-{step03_hash}" / log_path / x
artifact_benchmark_path = lambda x: Path(results_path) / f"{step03_name}-{step03_hash}" / benchmark_path / x

## Functions
################################################################################
create_save_params_rule(step03_name, artifact_path, step03_params, global_params)

def add_suffix_to_filename(path, suffix):
    base, ext = os.path.splitext(path)
    return base + suffix + ext

def pigeon_prepare_input(wildcards):
    parts = wildcards.prefix.split(".")

    if wildcards.prefix == "ref":
        return ref_annotation
    else:
        return expand(transcriptome_path("gffcompare/{dataset}.{group}.gffcompare.annotated.clean.gtf"), dataset = parts[0], group = parts[1])[0]

## Rules
################################################################################
ref_annotation_sorted = add_suffix_to_filename(ref_annotation, ".sorted")
rule pigeon_prepare_reference:
    message: """--- Pigeon Prepare Reference ---"""
    input: ref_annotation
    output: ref_annotation_sorted
    conda: "../envs/pigeon.yaml"
    shell: "pigeon prepare {input}"

rule pigeon_prepare:
    message: """--- Pigeon Prepare ---"""
    input:
        annotation = pigeon_prepare_input,
        genome = ref_genome
    output: 
        sorted_gtf = artifact_path("Prepare/{prefix}.pigeon.sorted.gtf")
    log: artifact_log_path("Pigeon_prepare/{prefix}.log")
    benchmark: artifact_benchmark_path("Pigeon_prepare/{prefix}.tsv")
    params:
        pigeon_output = lambda w, input: add_suffix_to_filename(input.annotation, ".sorted")
    conda: "../envs/pigeon.yaml"
    shell: 
        "pigeon prepare {input.annotation} {input.genome} 2>&1 | tee {log}; "
        "mv {params.pigeon_output} {output.sorted_gtf}; "
        "mv {params.pigeon_output}.pgi {output.sorted_gtf}.pgi; "

rule pigeon_classify:
    message: """--- Pigeon Classify ---"""
    input:
        sorted_gtf = artifact_path("Prepare/{prefix}.pigeon.sorted.gtf"),
        sorted_annotation = ref_annotation_sorted,
        genome = ref_genome
    output: 
        classification_txt = artifact_path("Classify_Filter/{prefix}.pigeon_classification.txt")
    log: artifact_log_path("Pigeon/classify_{prefix}.log")
    benchmark: artifact_benchmark_path("Pigeon/classify_{prefix}.pigeon.tsv")
    params:
        pigeon_output = lambda w, output: Path(output.classification_txt).parent
    threads: 16
    conda: "../envs/pigeon.yaml"
    shell:
        "pigeon classify -d {params.pigeon_output} -o {wildcards.prefix}.pigeon -j {threads} --log-level INFO "
        "{input.sorted_gtf} {input.sorted_annotation} {input.genome} 2>&1 | tee {log}"

rule pigeon_filter:
    message: """--- Pigeon Filter ---"""
    input:
        classification_txt = artifact_path("Classify_Filter/{prefix}.pigeon_classification.txt"),
        sorted_annotation = artifact_path("Prepare/{prefix}.pigeon.sorted.gtf")
    output:
        filtered_gtf = artifact_path("Classify_Filter/{prefix}.pigeon.sorted.filtered.gtf"),
        pigeon_summary = artifact_path("Classify_Filter/{prefix}.pigeon_classification.filtered_lite_classification.txt")
    log: artifact_log_path("Pigeon/filter_{prefix}.log")
    benchmark: artifact_benchmark_path("Pigeon/filter_{prefix}.tsv")
    params:
        pigeon_output = lambda w, input: add_suffix_to_filename(input.sorted_annotation, ".filtered_lite")
    threads: 16
    conda: "../envs/pigeon.yaml"
    shell:
        "pigeon filter -j {threads} --log-level INFO "
        "{input.classification_txt} --isoforms {input.sorted_annotation} 2>&1 | tee {log}; "
        "mv {params.pigeon_output} {output.filtered_gtf}"