################################################################################
## Alignment & Sorting
################################################################################

## Variables
################################################################################
step01_name = "01-Alignment"

### Settings
minimap2_settings = config["minimap2_settings"][config["minimap2_profile"]]
step01_params = {
    "minimap2_settings": minimap2_settings
}

### Compute hash and set global parameters
global_params.update({step01_name: step01_params})
step01_hash = compute_hash(global_params)
global_params.update({f"{step01_name}_{step01_hash}": global_params.pop(step01_name)})

### Paths
alignment_path = lambda x: Path(results_path) / f"{step01_name}_{step01_hash}" / x
alignment_log_path = lambda x: Path(results_path) / f"{step01_name}_{step01_hash}" / log_path / x
alignment_benchmark_path = lambda x: Path(results_path) / f"{step01_name}_{step01_hash}" / benchmark_path / x

## Functions
################################################################################
create_save_params_rule(step01_name, alignment_path, step01_params, global_params)

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
    output: 
        sam = temporary(alignment_path("Minimap2/{sample}.sam"))
    log: alignment_log_path("Minimap2/{sample}.log")
    benchmark: alignment_benchmark_path("Minimap2/{sample}.tsv")
    params:
        k = minimap2_settings["minimap2_kmer"],
        flags = minimap2_settings["minimap2_flags"]
    threads: 8
    conda: "../envs/minimap2.yaml",
    shell: "minimap2 {params.flags} -k {params.k} -t {threads} -o {output} {input.genome} {input.fastq} 2>&1 | tee {log}"

rule samtools_sort:
    message: """--- Runing Samtool Sort - {wildcards.sample} ----"""
    input: 
        sam = rules.minimap2_align.output.sam
    output: 
        bam = alignment_path("Samtools_sort/{sample}_sorted.bam"),
    log: alignment_log_path("Samtools_sort/{sample}.log"),
    benchmark: alignment_benchmark_path("Samtools_sort/{sample}.tsv"),
    resources:
        mem="10G"
    threads: 6
    conda: "../envs/samtools.yaml"
    shell:
        "samtools view -b -u -@{threads} {input} | "
        "samtools sort -@{threads} -m {resources.mem} -T tmp_{wildcards.sample} -o {output.bam}"

def subsamples_samples(wildcards):
    if wildcards.dataset == "Ebbert":
        return rules.samtools_sort.output.bam
    elif wildcards.dataset == "Wood":
        return Path(config["input_dir_wood"]) / wildcards.sample / "Mapping" / f"{wildcards.sample}_minimap.bam"

rule samtools_subsample:
    message: """--- Subsampling BAM files - {wildcards.sample} ----"""
    input:
        bam = subsamples_samples
    output:
        bam = alignment_path("Samtools_subsample/{sample}_sorted.bam"),
        bai = alignment_path("Samtools_subsample/{sample}_sorted.bam.bai")
    log: alignment_log_path("Samtools_subsample/{sample}.log"),
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
    
