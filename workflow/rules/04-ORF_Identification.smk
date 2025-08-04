################################################################################
## ORF Identification
################################################################################

## Variables
################################################################################
orfannotate_path = "tools/ORFannotate"
orf_path = lambda x: Path(results_output_path) / "04-ORF_Identification" / x

### Configurations

use_pigeon = config["use_Pigeon"]
use_sq3 = config["use_Sqanti3"]
sq3_rescue_run = config["sq3_rescue_run"]
gffread_clean = config["gffread_clean"]

## Only one of the previous can be true
use_sq3 = False if use_pigeon and use_sq3 else use_sq3

## Functions
################################################################################
# def gffread_input(wildcards):


#     prefix = wildcards.prefix

#     # Return the output of each rules with the prefix name as priority (this is mostly for debugging purposes)
#     if prefix == "sq3.filter":
#         return rules.sqanti3_filter.output.filtered_gtf
#     elif prefix == "sq3.rescue":
#         return rules.sqanti3_rescue.output.rescued_gtf
#     elif prefix == "pigeon.filter":
#         return rules.pigeon_filter.output.filtered_gtf

#     # If prefix does not fall in the previous categories, prioritize the configuration parameters
#     if use_pigeon:
#         return rules.pigeon_filter.output.filtered_gtf
#     if use_sq3:
#         return rules.sqanti3_rescue.output.rescued_gtf if sq3_rescue else rules.sqanti3_filter.output.filtered_gtf

#     # If nothing else matches, default to no artifact correction
#     return rules.manual_filter.output.gtf

def orfannotate_input(wildcards):
    if gffread_clean:
        return rules.pigeon_filter.output.filtered_gtf
    else:
        return rules.gffread.output.gtf

## Rules
################################################################################
rule gffread:
    message: """--- Cleaning gtf with gffread ----"""
    input: 
        gtf = pigeon_path("Classify_Filter/{prefix}.sorted.filtered.gtf")
    output: 
        gtf = orf_path("gffread/{prefix}.gtf")
    log: log_path("ORF_Identification/gffread_{prefix}.log")
    benchmark: benchmark_path("ORFannotate/gffread_{prefix}.tsv")
    conda: "../envs/gffread.yaml"
    shell:
        "gffread -E {input.gtf} -T -o {output.gtf} 2>&1 | tee {log}"

rule ORF_annotate:
    message: """--- ORF Identification with ORFannotate.py ----"""
    input:
        script = os.path.join(orfannotate_path, "ORFannotate.py"),
        gtf = orfannotate_input,
        ref_genome = ref_genome
    output: 
        summary_tsv = orf_path("ORFannotate/{prefix}/ORFannotate_summary.tsv"),
        gtf = orf_path("ORFannotate/{prefix}/ORFannotate_annotated.gtf"),
        protein_fa = orf_path("ORFannotate/{prefix}/protein.fa")
    log: log_path("ORFannotate/{prefix}.log")
    benchmark: benchmark_path("ORFannotate/{prefix}.tsv")
    params:
        out_dir = lambda w, output: Path(output.summary_tsv).parent
    threads: 1
    conda: "../envs/ORFannotate.yaml"
    shell:
        "python {input.script} --gtf {input.gtf} --fa {input.ref_genome} --outdir {params.out_dir} 2>&1 | tee {log}"

rule ORF_annotate_clean:
    message: """--- Clean ORFannotate output with gffread ----"""
    input:
        gtf = orf_path("ORFannotate/{prefix}/ORFannotate_annotated.gtf")
    output:
        gtf = orf_path("ORFannotate/{prefix}/ORFannotate_annotated_clean.gtf"),
    log: log_path("ORF_Identification/gffread_clean_{prefix}.log")
    benchmark: benchmark_path("ORFannotate/gffread_clean_{prefix}.tsv")
    conda: "../envs/gffread.yaml"
    shell:
        "gffread -E {input.gtf} -T -o {output.gtf} 2>&1 | tee {log}"

rule ORF_categorization:
    message: """--- Transcript categorization ----"""
    input:
        gtf = orf_path("ORFannotate/{prefix}/ORFannotate_annotated_clean.gtf"),
        ref_annotation = ref_annotation,
        protein_fa = orf_path("ORFannotate/{prefix}/protein.fa"),
        orf_summary = orf_path("ORFannotate/{prefix}/ORFannotate_summary.tsv"),
        pigeon_summary = pigeon_path("Classify_Filter/{prefix}.pigeon_classification.filtered_lite_classification.txt")
    output:
        gtf = orf_path("ORF_Category/{prefix}_orf.gtf"),
        isoform_summary = orf_path("ORF_Category/{prefix}_isoform_summary.tsv"),
    log: log_path("ORF_Category/{prefix}_orf.log")
    benchmark: benchmark_path("ORF_Category/{prefix}_orf.tsv")
    conda: "../envs/r.yaml"
    script: "../scripts/04a-ORF_Categorization.R"

rule ORF_filtration:
    message: """--- Transcript filtration ----"""
    input:
        gtf = rules.ORF_categorization.output.gtf
    output:
        gtf = orf_path("ORF_Filtration/{prefix}.orf_filter.gtf")
    log: log_path("ORF_Category/{prefix}_orf.log")
    benchmark: benchmark_path("ORF_Category/{prefix}_orf.tsv")
    params:
        valid_ref_gene_type = config["valid_ref_gene_type"],
        valid_ref_tx_type = config["valid_ref_tx_type"],
        valid_orfannotate_type = config["valid_orfannotate_type"],
        in_ref_filter = config["in_ref_filter"]
    conda: "../envs/r.yaml"
    script: "../scripts/04b-ORF_Filtration.R"