################################################################################
## Alignment & Sorting
################################################################################

## Variables
################################################################################

### Configurations
minimap2_settings = config["minimap2_settings"][config["minimap2_profile"]]

step01_params = {
    **minimap2_settings
}

### Compute the step hash
step01_hash = compute_hash(step01_params)

### Paths
alignment_path = lambda x: Path(results_path) / f"01-Alignment_{step01_hash}" / x
alignment_logs_path = lambda x: log_path(x, step01_hash)
alignment_benchmark_path = lambda x: benchmark_path(x, step01_hash)

## Functions
################################################################################
create_save_params_rule(step_num=1, step_name="01-Alignment", step_dir=alignment_path(""), step_params=step01_params)

def get_input_fastq(wildcards):
    endome_samples_df = generate_input_samples_df(wildcards.dataset, wildcards.group)
    row = endome_samples_df[endome_samples_df["sample_id"] == wildcards.sample]

    return row.iloc[0]["path"]

## Rules
################################################################################
rule minimap2_align:
    message: """--- Runing Minimap2 k{params.k} - {wildcards.sample} ----"""
    input:
        genome = ref_genome,
        fastq = get_input_fastq,
        params_file = alignment_path("parameters.yaml")
    output: 
        sam = temporary(alignment_path("Minimap2/{sample}.sam"))
    log: alignment_logs_path("Minimap2/{sample}.log")
    benchmark: alignment_benchmark_path("Minimap2/{sample}.tsv")
    params:
        k = step01_params["minimap2_kmer"],
        flags = step01_params["minimap2_flags"]
    threads: 8
    conda: "../envs/minimap2.yaml",
    shell: "minimap2 {params.flags} -k {params.k} -t {threads} -o {output} {input.genome} {input.fastq} 2>&1 | tee {log}"

rule samtools_sort:
    message: """--- Runing Samtool Sort - {wildcards.sample} ----"""
    input: 
        sam = rules.minimap2_align.output.sam
    output: 
        bam = alignment_path("Samtools_sort/{sample}_sorted.bam"),
    log: alignment_logs_path("Samtools_sort/{sample}.log"),
    benchmark: alignment_benchmark_path("Samtools_sort/{sample}.tsv"),
    resources:
        mem="10G"
    threads: 6
    conda: "../envs/samtools.yaml"
    shell:
        "samtools view -b -u -@{threads} {input} | "
        "samtools sort -@{threads} -m {resources.mem} -T tmp_{wildcards.sample} -o {output.bam}"

rule samtools_subsample:
    message: """--- Subsampling BAM files - {wildcards.sample} ----"""
    input:
        bam = rules.samtools_sort.output.bam
    output:
        bam = alignment_path("Samtools_subsample/{sample}_sorted.bam"),
        bai = alignment_path("Samtools_subsample/{sample}_sorted.bam.bai")
    log: alignment_logs_path("Samtools_subsample/{sample}.log"),
    benchmark: alignment_benchmark_path("Samtools_subsample/{sample}.tsv"),
    params:
        seed = config.get("subsample_seed", 0),
        prop = config.get("subsample_prop", "0001")
    resources:
        mem="10G"
    threads: 6
    conda: "../envs/samtools.yaml"
    shell: 
        """
        samtools view -s {params.seed}.{params.prop} -b {input.bam} | \\
        samtools sort -@{threads} -m {resources.mem} -T tmp_{wildcards.sample} -o {output.bam}

        samtools index {output.bam}
        """
    
