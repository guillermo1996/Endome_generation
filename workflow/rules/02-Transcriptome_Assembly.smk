################################################################################
## Transcriptome Assembly
################################################################################

## Variables
################################################################################
step02_name = "02-Transcriptome_Assembly"

### Configurations
stringtie_settings = config["stringtie_settings"][config["stringtie_profile"]]
gffcompare_settings = config["gffcompare_settings"][config["gffcompare_profile"]]

step02_params = {
    "stringtie_settings": stringtie_settings,
    "gffcompare_settings": gffcompare_settings,
    **({"toy_data": True} if config.get("use_toy_data", False) else {})
}

### Compute hash and set global parameters
global_params.update({step02_name: step02_params})
step02_hash = compute_hash(global_params)
global_params.update({f"{step02_name}_{step02_hash}": global_params.pop(step02_name)})

### Paths
transcriptome_path = lambda x: Path(results_path) / f"{step02_name}_{step02_hash}" / x
transcriptome_log_path = lambda x: Path(results_path) / f"{step02_name}_{step02_hash}" / log_path / x
transcriptome_benchmark_path = lambda x: Path(results_path) / f"{step02_name}_{step02_hash}" / benchmark_path / x

## Functions
################################################################################
create_save_params_rule(step02_name, transcriptome_path, step02_params, global_params)

def get_input_bam(wildcards):
    if config.get("use_toy_data", False):
        return alignment_path(f"Samtools_subsample/{wildcards.sample}_sorted.bam")

    if wildcards.dataset == "Ebbert":
        return alignment_path(f"Samtools_sort/{wildcards.sample}_sorted.bam")
    elif wildcards.dataset == "Wood":
        return Path(config["input_dir_wood"]) / wildcards.sample / "Mapping" / f"{wildcards.sample}_minimap.bam"

## Rules
################################################################################
rule stringtie_assembly:
    message: """--- StringTie Assembly - {wildcards.sample} ----"""
    input: 
        bam = get_input_bam,
    output: 
        gtf = transcriptome_path("StringTie_Sample_Assembly/{sample}.gtf")
    log: transcriptome_log_path("StringTie_Sample_Assembly/{sample}.log")
    benchmark: transcriptome_benchmark_path("StringTie_Sample_Assembly/{sample}.tsv")
    params:
        guided = f"-G {ref_annotation}" if stringtie_settings["assembly_guided"] else "",
        flags = stringtie_settings["assembly_flags"]
    threads: 8
    conda: "../envs/stringtie.yaml"
    shell:
        "stringtie {params.flags} -p {threads} {params.guided} -o {output} {input} 2>&1 | tee {log}"

rule write_gtf_list:
    message: """--- Write GTF list ----"""
    input: 
        gtf = lambda wc: expand(transcriptome_path("StringTie_Sample_Assembly/{sample}.gtf"), sample = generate_input_samples_df(wc.dataset, wc.group)["sample_id"].tolist(), dataset = wc.dataset, group = wc.group)
    output: 
        gtf_list = transcriptome_path("StringTie_Sample_Assembly/{dataset}.{group}.gtf_list.txt"),
    shell: 
        """
        printf '%s\n' {input} > {output.gtf_list}
        """

rule stringtie_merge:
    message: """--- StringTie Merge ----"""
    input: 
        gtf_list = rules.write_gtf_list.output.gtf_list
    output: 
        gtf = transcriptome_path("StringTie_Merge/{dataset}.{group}.merged.gtf"),
    log: transcriptome_log_path("StringTie_Merge/{dataset}.{group}.merged.log"),
    benchmark: transcriptome_benchmark_path("StringTie_Merge/{dataset}.{group}.merged.tsv")
    params:
        guided = f"-G {ref_annotation}" if stringtie_settings["merge_guided"] else "",
        flags = stringtie_settings["merge_flags"]
    threads: 16,
    conda: "../envs/stringtie.yaml",
    shell: "stringtie --merge {params.flags} {params.guided} -p {threads} -o {output} {input} 2>&1 | tee {log}"

rule gffcompare:
    message: """--- Running gffcompare ----"""
    input: 
        gtf = rules.stringtie_merge.output.gtf,
    output: 
        gtf = transcriptome_path("gffcompare/{dataset}.{group}.gffcompare.annotated.gtf"),
        stats = transcriptome_path("gffcompare/{dataset}.{group}.gffcompare.stats")
    log: transcriptome_log_path("gffcompare/{dataset}.{group}.gffcompare.annotated.log")
    benchmark: transcriptome_benchmark_path("gffcompare/{dataset}.{group}.gffcompare.annotated.tsv")
    params: 
        out_name = lambda wc: expand(transcriptome_path("gffcompare/{dataset}.{group}.gffcompare"), dataset = wc.dataset, group = wc.group),
        r_flag = f"-r {ref_annotation}",
        flags = gffcompare_settings["gffcompare_flags"]
    conda: "../envs/gffcompare.yaml",
    shell:
        """
        gffcompare {params.flags} {params.r_flag} -o {params.out_name} {input} 2>&1 | tee {log}

        if [ ! -f "{output.stats}" ] && [ -f "{params.out_name}" ]; then
            mv "{params.out_name}" "{output.stats}"
        fi
        """

rule manual_filter:
    message: """--- Chromosome and Strand filtration ----"""
    input: rules.gffcompare.output.gtf
    output: 
        gtf = transcriptome_path("gffcompare/{dataset}.{group}.gffcompare.annotated.clean.gtf")
    log: transcriptome_log_path("gffcompare/{dataset}.{group}.gffcompare.annotated.clean.log")
    benchmark: transcriptome_benchmark_path("gffcompare/{dataset}.{group}.gffcompare.annotated.clean.tsv")
    shell: "awk '$1 ~ /^(chr)?([1-9]|1[0-9]|2[0-2]|X|Y)$/ && $7 ~ /^[+-]$/' {input} > {output} 2> {log}"


# rule write_bam_list:
#     message: """--- Write BAM list ----"""
#     input: 
#         bam = lambda wc: expand(rules.samtools_subsample.output.bam, sample = generate_input_samples_df(wc.dataset, wc.group)["sample_id"].tolist(), dataset = wc.dataset, group = wc.group)
#     output: 
#         bam_list = transcriptome_path("IsoQuant_Assembly/{dataset}.{group}.bam_list.txt"),
#     shell: 
#         """
#         printf '%s\n' {input.bam} > {output.bam_list}
#         """

# rule isoquant_assembly:
#     message: """--- IsoQuant Assembly ----"""
#     input:
#         bam_list = rules.write_bam_list.output.bam_list
#     output:
#         "{dataset}.{group}.flag"
#         # gtf = transcriptome_path("IsoQuant/{dataset}.{group}.transcript_models.gtf")
#     log: transcriptome_logs_path("Isoquant/{dataset}.{group}.test.log")
#     benchmark: transcriptome_benchmark_path("Isoquant/{dataset}.{group}.test.tsv")
#     params:
#         ref_genome = ref_genome,
#         genedb = ref_annotation
#         # flags = "--data_type nanopore --complete_gemedb"
#         # out_dir = lambda w, output: Path(output.summary_tsv).parent
#     threads: 8,
#     conda: "../envs/isoquant.yaml",
#     shell: 
#         """
#         isoquant.py --reference {params.ref_genome} --genedb {params.genedb} --complete_genedb \\
#         --bam $(cat {input.bam_list}) \\
#         -o test \\
#         --threads 8 \\
#         --process_only_chr chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 \\
#         --high_memory \\
#         --force \\
#         --data_type nanopore

#         touch {output}
#         """
