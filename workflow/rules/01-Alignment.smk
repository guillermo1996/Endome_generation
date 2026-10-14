################################################################################
## Alignment & Sorting
################################################################################

## Variables
################################################################################
### Register the step: generates the hash suffix, creates the parameter rule and returns the output_path helpers.
step01 = register_step(
    name="01-Alignment",
    params={
        "minimap2_settings": build_tool_settings(config, "minimap2_settings", "minimap2_preset"),
    },
)

### Resolved settings used directly in the rule bodies below
minimap2_settings = resolve_preset(config, "minimap2_settings", "minimap2_preset")

## Functions
################################################################################
def get_input_fastq(wildcards):
    sample = get_input_sample(wildcards)

    if sample["file_type"] != "fastq":
        raise ValueError(f"{wildcards.sample} is a '{sample['file_type']}' input, so it has no alignment step")
    return sample["path"]

def get_toy_samples(wildcards):
    sample = get_input_sample(wildcards)
    return sample["path"] if sample["file_type"] == "bam" else rules.samtools_sort.output.bam

## Rules
################################################################################
rule minimap2_align:
    message: """--- Runing Minimap2 k{params.k} - {wildcards.sample} ----"""
    input:
        genome = ref_genome,
        fastq = get_input_fastq,
    output: 
        # sam = temporary(step01.path("Minimap2/{sample}.sam")),
        sam = step01.path("Minimap2/{sample}.sam"),
    log: step01.logs("Minimap2/{sample}.log")
    benchmark: step01.benchmark("Minimap2/{sample}.tsv")
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
        bam = step01.path("Samtools_sort/{sample}_sorted.bam"),
    log: step01.logs("Samtools_sort/{sample}.log"),
    benchmark: step01.benchmark("Samtools_sort/{sample}.tsv"),
    resources:
        mem_gb=10
    threads: 6
    conda: "../envs/samtools.yaml"
    shell:
        "samtools view -b -u -@{threads} {input} | "
        "samtools sort -@{threads} -m {resources.mem_gb}G -T tmp_{wildcards.sample} -o {output.bam}"

rule samtools_subsample:
    message: """--- Subsampling BAM files - {wildcards.sample} ----"""
    input:
        bam = get_toy_samples
    output:
        bam = step01.path("Samtools_subsample/{sample}_sorted.bam"),
        bai = step01.path("Samtools_subsample/{sample}_sorted.bam.bai")
    log: step01.logs("Samtools_subsample/{sample}.log"),
    benchmark: step01.benchmark("Samtools_subsample/{sample}.tsv"),
    params:
        seed = config.get("subsample_seed", 0),
        prop = config.get("subsample_prop", "0001")
    resources:
        mem_gb=20
    threads: 6
    conda: "../envs/samtools.yaml"
    shell:
        """
        samtools view -s {params.seed}.{params.prop} -b {input.bam} | \\
        samtools sort -@{threads} -m {resources.mem_gb}G -T tmp_{wildcards.sample} -o {output.bam}

        samtools index {output.bam}
        """