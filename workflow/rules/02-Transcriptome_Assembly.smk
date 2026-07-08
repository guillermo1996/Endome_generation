################################################################################
## Transcriptome Assembly
################################################################################

## Debug: benchmark loading time
_start_time = time.perf_counter()

## Variables
################################################################################
### Register the step: generates the hash suffix, creates the parameter rule and returns the output_path helpers.
step02 = register_step(
    name="02-Transcriptome_Assembly",
    params={
        "stringtie_presets": build_tool_settings(config, "stringtie_settings", "stringtie_preset"),
        "gffcompare_presets": build_tool_settings(config, "gffcompare_settings", "gffcompare_preset"),
    },
)

### Resolved settings used directly in the rule bodies below
stringtie_settings = resolve_preset(config, "stringtie_settings", "stringtie_preset")
gffcompare_settings = resolve_preset(config, "gffcompare_settings", "gffcompare_preset")

## Functions
################################################################################
def get_input_bam(wildcards):
    if config.get("use_toy_data", False):
        return step01.path(f"Samtools_subsample/{wildcards.sample}_sorted.bam")

    if wildcards.dataset == "Ebbert":
        return step01.path(f"Samtools_sort/{wildcards.sample}_sorted.bam")
    elif wildcards.dataset == "Wood":
        return Path(config["input_dir_wood"]) / wildcards.sample / "Mapping" / f"{wildcards.sample}_minimap.bam"


## Rules
################################################################################
rule stringtie_assembly:
    message: """--- StringTie Assembly - {wildcards.sample} ----"""
    input: 
        bam = get_input_bam,
    output: 
        gtf = step02.path("StringTie_Sample_Assembly/{sample}.gtf")
    log: step02.log("StringTie_Sample_Assembly/{sample}.log")
    benchmark: step02.benchmark("StringTie_Sample_Assembly/{sample}.tsv")
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
        gtf = lambda wc: expand(step02.path("StringTie_Sample_Assembly/{sample}.gtf"), sample = generate_input_samples_df(wc.dataset, wc.group)["sample_id"].tolist(), dataset = wc.dataset, group = wc.group)
    output: 
        gtf_list = step02.path("StringTie_Sample_Assembly/{dataset}.{group}.gtf_list.txt"),
    shell: 
        """
        printf '%s\n' {input} > {output.gtf_list}
        """

rule stringtie_merge:
    message: """--- StringTie Merge ----"""
    input: 
        gtf_list = rules.write_gtf_list.output.gtf_list
    output: 
        gtf = step02.path("StringTie_Merge/{dataset}.{group}.merged.gtf"),
    log: step02.log("StringTie_Merge/{dataset}.{group}.merged.log"),
    benchmark: step02.benchmark("StringTie_Merge/{dataset}.{group}.merged.tsv")
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
        gtf = step02.path("gffcompare/{dataset}.{group}.gffcompare.annotated.gtf"),
        stats = step02.path("gffcompare/{dataset}.{group}.gffcompare.stats")
    log: step02.log("gffcompare/{dataset}.{group}.gffcompare.annotated.log")
    benchmark: step02.benchmark("gffcompare/{dataset}.{group}.gffcompare.annotated.tsv")
    params: 
        out_name = lambda wc: expand(step02.path("gffcompare/{dataset}.{group}.gffcompare"), dataset = wc.dataset, group = wc.group),
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
        gtf = step02.path("gffcompare/{dataset}.{group}.gffcompare.annotated.clean.gtf")
    log: step02.log("gffcompare/{dataset}.{group}.gffcompare.annotated.clean.log")
    benchmark: step02.benchmark("gffcompare/{dataset}.{group}.gffcompare.annotated.clean.tsv")
    shell: "awk '$1 ~ /^(chr)?([1-9]|1[0-9]|2[0-2]|X|Y)$/ && $7 ~ /^[+-]$/' {input} > {output} 2> {log}"
