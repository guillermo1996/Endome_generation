################################################################################
## Transcriptome Assembly
################################################################################

## Variables
################################################################################

### Configurations
ref_in_assembly = config["ref_in_assembly"]
long_read_processing_in_assembly = config["long_read_processing_in_assembly"]

ref_in_merge = config["ref_in_merge"]
long_read_processing_in_merge = config["long_read_processing_in_merge"]

k_flag = config["k_flag"]

## Functions
################################################################################
def get_input_bam(wildcards):
    if config["input_dataset"] == "Ebbert":
        return alignment_path(f"Samtools_sort.k{k_flag}/{wildcards.sample}_sorted.bam")
    elif config["input_dataset"] == "Wood":
        return Path(config["input_dir_wood"]) / wildcards.sample / "Mapping" / f"{wildcards.sample}_minimap.bam"

## Rules
################################################################################
rule stringtie_assembly:
    message: """--- StringTie Assembly - {wildcards.sample} ----"""
    input: 
        bam = get_input_bam
    output: 
        gtf = transcriptome_path("StringTie_Sample_Assembly/{sample}.gtf")
    log: log_path("StringTie_Sample_Assembly/{sample}.log")
    benchmark: benchmark_path("StringTie_Sample_Assembly/{sample}.tsv")
    params:
        G_flag = f"-G {ref_annotation}" if ref_in_assembly else "",
        L_flag = f"-L" if long_read_processing_in_assembly else "",
    threads: 8
    conda: "../envs/stringtie.yaml"
    shell:
        "stringtie --rf -p {threads} {params.L_flag} {params.G_flag} -o {output} {input} 2>&1 | tee {log}"

rule write_gtf_list:
    message: """--- Write GTF list ----"""
    input: 
        gtf = lambda wc: expand(transcriptome_path("StringTie_Sample_Assembly/{sample}.gtf"), sample = generate_input_samples(wc.dataset, wc.group)["sample_id"].tolist(), dataset = wc.dataset, group = wc.group)
    output: 
        gtf_list = transcriptome_path("StringTie_Sample_Assembly/{dataset}.{group}.gtf_list.txt"),
        build_config = transcriptome_path("StringTie_Sample_Assembly/{dataset}.{group}.build_config.txt")
    params: 
        minimap2_k_flag = config["k_flag"]
    shell: 
        """
        printf '%s\n' {input} > {output.gtf_list}
        echo "k_flag: {params.minimap2_k_flag}" > {output.build_config}
        """

rule stringtie_merge:
    message: """--- StringTie Merge ----"""
    input: 
        gtf_list = rules.write_gtf_list.output.gtf_list
    output: 
        gtf = transcriptome_path("StringTie_Merge/{dataset}.{group}.merged_genome.gtf"),
    log: log_path("StringTie_Merge/{dataset}.{group}.merged_genome.log"),
    benchmark: benchmark_path("StringTie_Merge/{dataset}.{group}.merged_genome.tsv")
    params:
        G_flag = f"-G {ref_annotation}" if ref_in_merge else "",
        L_flag = f"-L" if long_read_processing_in_merge else "",
    threads: 16,
    conda: "../envs/stringtie.yaml",
    shell: "stringtie --merge {params.L_flag} {params.G_flag} -p {threads} -o {output} {input} 2>&1 | tee {log}"

rule gffcompare:
    message: """--- Running gffcompare ----"""
    input: 
        gtf = rules.stringtie_merge.output.gtf,
    output: 
        gtf = transcriptome_path("gffcompare/{dataset}.{group}.gffcompare.annotated.gtf")
    log: log_path("gffcompare/{dataset}.{group}.gffcompare.annotated.log")
    benchmark: benchmark_path("gffcompare/{dataset}.{group}.gffcompare.annotated.tsv")
    params: 
        r_flag = f"-r {ref_annotation}",
        output_name = lambda wc: expand(transcriptome_path("gffcompare/{dataset}.{group}.gffcompare"), dataset = wc.dataset, group = wc.group),
    conda: "../envs/gffcompare.yaml",
    shell: "gffcompare -V {params.r_flag} -o {params.output_name} {input} 2>&1 | tee {log}"

rule manual_filter:
    message: """--- Chromosome and Strand filtration ----"""
    input: rules.gffcompare.output.gtf
    output: 
        gtf = transcriptome_path("gffcompare/{dataset}.{group}.gffcompare.annotated.clean.gtf")
    log: log_path("gffcompare/{dataset}.{group}.gffcompare.annotated.clean.log")
    benchmark: benchmark_path("gffcompare/{dataset}.{group}.gffcompare.annotated.clean.tsv")
    shell: "awk '$1 ~ /^(chr)?([1-9]|1[0-9]|2[0-2]|X|Y)$/ && $7 ~ /^[+-]$/' {input} > {output} 2> {log}"