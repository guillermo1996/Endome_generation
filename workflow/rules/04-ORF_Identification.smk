################################################################################
## ORF Identification
################################################################################

## Debug: benchmark loading time
_start_time = time.perf_counter()

## Variables
################################################################################
### Register the step: generates the hash suffix, creates the parameter rule and returns the output_path helpers.
step04 = register_step(
    name="04-ORF_Identification",
    params={
        "orfannotate_presets": build_tool_settings(config, "orfannotate_settings", "orfannotate_preset"),
        "orf_filter_presets": build_tool_settings(config, "orf_filter_settings"),
    },
    extra_params={
        "orf_filter_preset": config["orf_filter_preset"]
    }
)

### Resolved settings used directly in the rule bodies below
orfannotate_settings = resolve_preset(config, "orfannotate_settings", "orfannotate_preset")
orf_filter_settings = resolve_preset(config, "orf_filter_settings", "orf_filter_preset")

### ORFannotate Download path
orfannotate_version = orfannotate_settings["version"]
orfannotate_tar_url = f"https://github.com/egustavsson/ORFannotate/archive/refs/tags/{orfannotate_version}.tar.gz"
orfannotate_env_url = f"https://raw.githubusercontent.com/egustavsson/ORFannotate/refs/tags/{orfannotate_version}/ORFannotate.conda_env.yml"

## Functions
################################################################################
def orfannotate_input(wildcards):
    if orfannotate_settings.get("gffread_clean", False):
        return rules.gffread.output.gtf
    else:
        return rules.pigeon_filter.output.filtered_gtf

## Rules
################################################################################
rule download_ORF_annotate:
    message: "--- Downloading and Extracting ORF Annotate ---"
    output:
        orfanntoate_dir = directory(f"tools/ORFannotate-{orfannotate_version}"),
        script = f"tools/ORFannotate-{orfannotate_version}/ORFannotate.py"
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
        gtf = step03.path("Classify_Filter/{prefix}.sorted.filtered.gtf")
    output:
        gtf = step04.path("gffread/{prefix}.gtf")
    log: step04.log("ORF_Identification/gffread_{prefix}.log")
    benchmark: step04.benchmark("ORFannotate/gffread_{prefix}.tsv")
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
        summary_tsv = step04.path("ORFannotate/{prefix}/ORFannotate_summary.tsv"),
        gtf = step04.path("ORFannotate/{prefix}/ORFannotate_annotated.gtf"),
        protein_fa = step04.path("ORFannotate/{prefix}/protein.fa"),
        utr3_fa = step04.path("ORFannotate/{prefix}/utr3.fa"),
        utr5_fa = step04.path("ORFannotate/{prefix}/utr5.fa")
    log: step04.log("ORFannotate/{prefix}.log")
    benchmark: step04.benchmark("ORFannotate/{prefix}.tsv")
    params:
        out_dir = lambda w, output: Path(output.summary_tsv).parent
    threads: 1
    conda: orfannotate_env_url
    shell:
        "python {input.script} --gtf {input.gtf} --fa {input.ref_genome} --outdir {params.out_dir} 2>&1 | tee {log}"

# rule ORF_annotate_clean:
#     message: """--- Clean ORFannotate output with gffread ---"""
#     input:
#         gtf = step04.path("ORFannotate/{prefix}/ORFannotate_annotated.gtf")
#     output:
#         gtf = step04.path("ORFannotate/{prefix}/ORFannotate_annotated_clean.gtf"),
#     log: step04.log("ORF_Identification/gffread_clean_{prefix}.log")
#     benchmark: step04.benchmark("ORFannotate/gffread_clean_{prefix}.tsv")
#     conda:  "../envs/gffread.yaml"
#     shell:
#         "gffread -E {input.gtf} -T -o {output.gtf} 2>&1 | tee {log}"

rule ORF_categorization:
    message: """--- Transcript categorization ---"""
    input:
        gtf = step04.path("ORFannotate/{prefix}/ORFannotate_annotated.gtf"),
        ref_annotation = ref_annotation,
        protein_fa = step04.path("ORFannotate/{prefix}/protein.fa"),
        orf_summary = step04.path("ORFannotate/{prefix}/ORFannotate_summary.tsv"),
        pigeon_summary = step03.path("Classify_Filter/{prefix}.pigeon_classification.filtered_lite_classification.txt")
    output:
        gtf = step04.path("ORF_Category/{prefix}.orf_annotated.gtf"),
        isoform_summary = step04.path("ORF_Category/{prefix}.isoform_summary.tsv"),
    log: step04.log("ORF_Category/{prefix}.orf_annotated.log")
    benchmark: step04.benchmark("ORF_Category/{prefix}.orf_annotated.tsv")
    conda: "../envs/r.yaml"
    script: "../scripts/04a-ORF_Categorization.R"


# def select_source_dataset(wildcards):
#     parts = wildcards.prefix.split(".")

#     if parts[0] == "gencode":
#         return ref_annotation
#     else:
#         return rules.ORF_categorization.output.gtf

rule ORF_filtration:
    message: """--- Transcript filtration ---"""
    input:
        gtf = rules.ORF_categorization.output.gtf
    output:
        gtf_filter = step04.path("ORF_Filtration/{prefix}.{orf_filter}.orf_filter.gtf"),
    log: step04.log("ORF_Category/{prefix}.{orf_filter}_orf.log")
    benchmark: step04.benchmark("ORF_Category/{prefix}.{orf_filter}_orf.tsv")
    params:
        main_config = lambda wc: wc.orf_filter,
        valid_ref_gene_type = orf_filter_settings["valid_ref_gene_type"],
        valid_ref_tx_type = orf_filter_settings["valid_ref_tx_type"],
        valid_orfannotate_type = orf_filter_settings["valid_orfannotate_type"],
        in_ref_filter = orf_filter_settings["in_ref_filter"],
    conda: "../envs/r.yaml"
    script: "../scripts/04b-ORF_Filtration.R"

register_test_data_link(rules.ORF_annotate.output.summary_tsv)
register_test_data_link(rules.ORF_annotate.output.gtf)
register_test_data_link(rules.ORF_annotate.output.protein_fa)
register_test_data_link(rules.ORF_annotate.output.utr3_fa)
register_test_data_link(rules.ORF_annotate.output.utr5_fa)
register_test_data_link(rules.ORF_categorization.output.isoform_summary)
register_test_data_link(rules.ORF_categorization.output.gtf)
register_test_data_link(rules.ORF_filtration.output.gtf_filter)

## Debug: benchmark loading time
_log(f"\t+ {step04.name} imported in {time.perf_counter() - _start_time:.3f}s")
