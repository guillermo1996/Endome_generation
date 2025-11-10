################################################################################
## Alignment & Sorting
################################################################################

## Variables
################################################################################

### Configurations

## Functions
################################################################################

# def get_fastq_files(wildcards):
#     file_name = f"{sample_prefix}{wildcards.sample}{sample_suffix}"
#     file_path = os.path.join(input_dir, file_name)
#     return(file_path)

def get_input_fastq(wildcards):
    endome_samples_df = generate_input_samples(wildcards.dataset, wildcards.group)
    row = endome_samples_df[endome_samples_df["sample_id"] == wildcards.sample]

    return row.iloc[0]["path"]


## Rules
################################################################################
rule minimap2_align:
    message: """--- Runing Minimap2 k{wildcards.k} - {wildcards.sample} ----"""
    input:
        genome = ref_genome,
        fastq = get_input_fastq
    output: 
        sam = temporary(alignment_path("Minimap2.k{k}/{sample}.sam"))
    log: log_path("Minimap2.k{k}/{sample}.log")
    benchmark: benchmark_path("Minimap2.k{k}/{sample}.tsv")
    params:
        extra = config["minimap2_flags"],
    threads: 8
    conda: "../envs/minimap2.yaml",
    shell: "minimap2 {params.extra} -k {wildcards.k} -t {threads} -o {output} {input.genome} {input.fastq} 2>&1 | tee {log}"

rule samtools_sort:
    message: """--- Runing Samtool Sort - {wildcards.sample} ----"""
    input: 
        sam = rules.minimap2_align.output.sam
    output: 
        bam = alignment_path("Samtools_sort.k{k}/{sample}_sorted.bam"),
    log: log_path("Samtools_sort.k{k}/{sample}.log"),
    benchmark: benchmark_path("Samtools_sort.k{k}/{sample}.tsv"),
    resources:
        mem="10G"
    threads: 6
    conda: "../envs/samtools.yaml"
    shell:
        "samtools view -b -u -@{threads} {input} | "
        "samtools sort -@{threads} -m {resources.mem} -T tmp_{wildcards.sample} -o {output}"