################################################################################
## Artifact Removal
################################################################################

## Debug: benchmark loading time
_start_time = time.perf_counter()

## Variables
################################################################################
### Register the step: generates the hash suffix, creates the parameter rule and returns the output_path helpers.
step03 = register_step(
    name="03-Artifact_Removal",
    params={
        "pigeon_presets": build_tool_settings(config, "pigeon_settings", "pigeon_preset"),
    },
)

## Functions
################################################################################
def add_suffix_to_filename(path, suffix):
    base, ext = os.path.splitext(path)
    return base + suffix + ext

def pigeon_prepare_input(wildcards):
    parts = wildcards.prefix.split(".")

    if wildcards.prefix == "ref":
        return ref_annotation
    else:
        # return rules.manual_filter.output.gtf
        return expand(step02.path("transcriptome_assembly/{dataset}.{group}.{merge_method}.annotated.clean.gtf"), dataset = parts[0], group = parts[1], merge_method = parts[2])[0]

def pigeon_prepare_input(wildcards):
    parts = wildcards.prefix.split(".")

    if wildcards.prefix == "ref":
        return ref_annotation
    elif wildcards.prefix == "gencode":
        return expand(rules.manual_filter.output.gtf, dataset = wildcards.prefix, group = "none", merge_method = "ref")[0]
    else:
        return expand(rules.manual_filter.output.gtf, dataset = parts[0], group = parts[1], merge_method = parts[2])[0]

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
        sorted_gtf = step03.path("Prepare/{prefix}.pigeon.sorted.gtf")
    log: step03.log("Pigeon_prepare/{prefix}.log")
    benchmark: step03.benchmark("Pigeon_prepare/{prefix}.tsv")
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
        sorted_gtf = step03.path("Prepare/{prefix}.pigeon.sorted.gtf"),
        sorted_annotation = ref_annotation_sorted,
        genome = ref_genome
    output:
        classification_txt = step03.path("Classify_Filter/{prefix}.pigeon_classification.txt")
    log: step03.log("Pigeon/classify_{prefix}.log")
    benchmark: step03.benchmark("Pigeon/classify_{prefix}.pigeon.tsv")
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
        classification_txt = step03.path("Classify_Filter/{prefix}.pigeon_classification.txt"),
        sorted_annotation = step03.path("Prepare/{prefix}.pigeon.sorted.gtf")
    output:
        filtered_gtf = step03.path("Classify_Filter/{prefix}.pigeon.sorted.filtered.gtf"),
        pigeon_summary = step03.path("Classify_Filter/{prefix}.pigeon_classification.filtered_lite_classification.txt")
    log: step03.log("Pigeon/filter_{prefix}.log")
    benchmark: step03.benchmark("Pigeon/filter_{prefix}.tsv")
    params:
        pigeon_output = lambda w, input: add_suffix_to_filename(input.sorted_annotation, ".filtered_lite")
    threads: 16
    conda: "../envs/pigeon.yaml"
    shell:
        "pigeon filter -j {threads} --log-level INFO "
        "{input.classification_txt} --isoforms {input.sorted_annotation} 2>&1 | tee {log}; "
        "mv {params.pigeon_output} {output.filtered_gtf}"

register_test_data_link(rules.pigeon_filter.output.filtered_gtf)
register_test_data_link(rules.pigeon_filter.output.pigeon_summary)

## Debug: benchmark loading time
_log(f"\t+ {step03.name} imported in {time.perf_counter() - _start_time:.3f}s")


# rule download_sqanti3:
#     output:
#         sq3_dir = directory(sq3_path),
#         sq3_qc = f"{sq3_path}/sqanti3_qc.py",
#         sq3_filter = f"{sq3_path}/sqanti3_filter.py",
#         sq3_rescue = f"{sq3_path}/sqanti3_rescue.py",
#         sq3_env = f"{sq3_path}/SQANTI3.conda_env.yml",
#         default_rules = f"{sq3_path}/src/utilities/filter/filter_default.json"
#     params:
#         version = sq3_version
#     log: log_path("download_sqanti3.log")
#     conda: "../envs/download_scripts.yaml"
#     script: "../scripts/03a-Download_sqanti3.sh"

# rule sqanti3_qc:
#     input: 
#         sq3_qc = f"{sq3_path}/sqanti3_qc.py",
#         isoform = lambda wc: rules.manual_filter.output if wc.sq3_ref != "_reference" else ref_annotation,
#         ref_gtf = ref_annotation,
#         ref_fasta = ref_genome,
#     output:
#         classification = sqanti3_path("Sqanti3_QC{sq3_ref}/{prefix}_classification.txt"),
#         corrected_gff = sqanti3_path("Sqanti3_QC{sq3_ref}/{prefix}_corrected.gtf"),
#         corrected_fasta = sqanti3_path("Sqanti3_QC{sq3_ref}/{prefix}_corrected.fasta"),
#     log: log_path("Sqanti3_QC{sq3_ref}/{prefix}_qc.log")
#     benchmark: benchmark_path("Sqanti3_QC{sq3_ref}/{prefix}_qc.tsv")
#     params:
#         out_dir = sqanti3_path("Sqanti3_QC{sq3_ref}")
#     wildcard_constraints: 
#         sq3_ref = ".{0}|_.+"
#     threads: 16
#     conda: "../envs/sqanti3.yaml" # ".snakemake/conda/sq3_env" ~ Bug in current snakemake version https://github.com/snakemake/snakemake/issues/3192
#     shell: 
#         "python {input.sq3_qc} --force_id_ignore --skipORF --report skip "
#         "-o {wildcards.prefix} -d {params.out_dir} -n {threads} "
#         "{input.isoform} {input.ref_gtf} {input.ref_fasta} 2>&1 | tee {log}"
#     # shell: 
#     #     "python {input.sq3_qc} --force_id_ignore --skipORF --report skip --isoforms {input.isoform} "
#     #     "--refGTF {input.ref_gtf} --refFasta {input.ref_fasta} "
#     #     "-o {wildcards.prefix} -d {params.out_dir} -n {threads} 2>&1 | tee {log}"

# rule sqanti3_filter:
#     input:
#         sq3_filter = f"{sq3_path}/sqanti3_filter.py",
#         classification = sqanti3_path("Sqanti3_QC/{prefix}_classification.txt"),
#         corrected_gtf = sqanti3_path("Sqanti3_QC/{prefix}_corrected.gtf")
#     output: 
#         filtered_gtf = sqanti3_path("Sqanti3_Filter/{prefix}.filtered.gtf"),
#         filtered_classification = sqanti3_path("Sqanti3_Filter/{prefix}_RulesFilter_result_classification.txt")
#     log: log_path("Sqanti3_Filter/{prefix}_filter.log")
#     benchmark: benchmark_path("Sqanti3_Filter/{prefix}_filter.tsv")
#     params:
#         out_dir = sqanti3_path("Sqanti3_Filter")
#     threads: 16
#     conda: "../envs/sqanti3.yaml"
#     shell: 
#         "python {input.sq3_filter} rules "
#         "--gtf {input.corrected_gtf} "
#         "-o {wildcards.prefix} -d {params.out_dir} "
#         "{input.classification} 2>&1 | tee {log}"
#     # shell: 
#     #     "python {input.sq3_filter} rules --sqanti_class {input.classification} "
#     #     "--filter_gtf {input.corrected_gtf} "
#     #     "-o {wildcards.prefix} -d {params.out_dir} -c {threads} 2>&1 | tee {log}"

# rule sqanti3_rescue:
#     input:
#         sq3_rescue = f"{sq3_path}/sqanti3_rescue.py",
#         filter_class = sqanti3_path("Sqanti3_Filter/{prefix}_RulesFilter_result_classification.txt"),
#         ref_gtf = ref_annotation,
#         ref_fasta = ref_genome,
#         rescue_isoforms = sqanti3_path("Sqanti3_QC/{prefix}_corrected.fasta"),
#         rescue_gtf = sqanti3_path("Sqanti3_Filter/{prefix}.filtered.gtf"),
#         ref_classif = sqanti3_path("Sqanti3_QC_reference/sq3.reference_classification.txt"),
#         sq3_rules = os.path.join(sq3_path, "src/utilities/filter/filter_default.json"),
#     output: 
#         rescued_gtf = sqanti3_path("Sqanti3_Rescue/{prefix}_rescued.gtf")
#     log: log_path("Sqanti3_Rescue/{prefix}_rescue.log")
#     benchmark: benchmark_path("Sqanti3_Rescue/{prefix}_rescue.tsv")
#     params:
#         out_dir = sqanti3_path("Sqanti3_Rescue")
#     threads: 16
#     conda: "../envs/sqanti3.yaml"
#     shell:
#         "python {input.sq3_rescue} rules "
#         "--refGTF {input.ref_gtf} --refGenome {input.ref_fasta} "
#         "--isoforms {input.rescue_isoforms} --gtf {input.rescue_gtf} "
#         "--refClassif {input.ref_classif} "
#         "--mode full -j {input.sq3_rules} "
#         "-o {wildcards.prefix} -d {params.out_dir} {input.filter_class} 2>&1 | tee {log}"
#     # shell:
#     #     "python {input.sq3_rescue} rules --filter_class {input.filter_class} "
#     #     "--refGTF {input.ref_gtf} --refFasta {input.ref_fasta} "
#     #     "--rescue_isoforms {input.rescue_isoforms} --rescue_gtf {input.rescue_gtf} "
#     #     "--refClassif {input.ref_classif} "
#     #     "--mode full -j {input.sq3_rules} "
#     #     "-o {wildcards.prefix} -d {params.out_dir} -c {threads} 2>&1 | tee {log}"
