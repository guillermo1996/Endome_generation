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
        "isomatch_presets": build_tool_settings(config, "isomatch_settings", "isomatch_preset"),
        "gffcompare_presets": build_tool_settings(config, "gffcompare_settings", "gffcompare_preset"),
    },
    extra_params={
        "transcript_merge_method": config["transcript_merge_method"],
    },
)

### Resolved settings used directly in the rule bodies below
stringtie_settings = resolve_preset(config, "stringtie_settings", "stringtie_preset")
isomatch_settings = resolve_preset(config, "isomatch_settings", "isomatch_preset")
gffcompare_settings = resolve_preset(config, "gffcompare_settings", "gffcompare_preset")

### Which tool's merged GTF feeds into gffcompare
transcript_merge_method = config["transcript_merge_method"]

### ORFannotate Download path
isomatch_version = isomatch_settings["version"]
isomatch_tar_url = f"https://github.com/zhengxinchang/isomatch/releases/download/{isomatch_version}/isomatch-{isomatch_version}-linux-x86_64.tar.gz"

## Functions
################################################################################
def get_input_bam(wildcards):
    if config.get("use_toy_data", False):
        return step01.path(f"Samtools_subsample/{wildcards.sample}_sorted.bam")

    if wildcards.dataset == "Ebbert":
        return step01.path(f"Samtools_sort/{wildcards.sample}_sorted.bam")
    elif wildcards.dataset == "Wood":
        return Path(config["input_dir_wood"]) / wildcards.sample / "Mapping" / f"{wildcards.sample}_minimap.bam"

def get_filter_input(wildcards):
    if wildcards.merge_method == "stringtie":
        return rules.gffcompare.output.gtf
    elif wildcards.merge_method == "isomatch":
        return rules.isomatch_classify.output.annotated_gtf
    elif wildcards.merge_method == "ref":
        return ref_annotation

## Rules
################################################################################
rule download_isomatch:
    message: "--- Downloading and Extracting isomatch ---"
    output:
        isomatch_dir = directory(f"tools/isomatch-{isomatch_version}"),
        executable = f"tools/isomatch-{isomatch_version}/isomatch"
    params:
        url = isomatch_tar_url
    shell:
        """
        mkdir -p {output.isomatch_dir}
        curl -L {params.url} | tar -xz -C {output.isomatch_dir}
        """
### StringTie Rules
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
        gtf = lambda wc: expand(step02.path("StringTie_Sample_Assembly/{sample}.gtf"), sample = generate_input_samples_df(wc.dataset, wc.group)["sample_id"].tolist(), dataset = wc.dataset, group = wc.group)
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

### Isomatch rules
rule isomatch_index_reference:
    message: "--- Indexing reference annotation with isomatch ---"
    input:
        gtf = ref_annotation,
        executable = rules.download_isomatch.output.executable,
        ref_genome = ref_genome
    output:
        gtf = step02.path("isomatch_ref_index/ref_annotation.gtf"),
        isomx = step02.path("isomatch_ref_index/ref_annotation.gtf.isomx"),
        isoms = step02.path("isomatch_ref_index/ref_annotation.gtf.isoms"),
    log: step02.log("isomatch_ref_index/ref_annotation.log")
    benchmark: step02.benchmark("isomatch_ref_index/ref_annotation.tsv")
    params:
        index_dir = lambda wildcards, output: str(Path(output.gtf).parent)
    threads: 1,
    shell:
        """
        set -euo pipefail
        ln -sf $(readlink -f {input.gtf}) {output.gtf}
        {input.executable} index --ref-fa {input.ref_genome} {output.gtf} --out {output.isomx} 2>&1 | tee {log}
        rm -rf {params.index_dir}/.isomatch-index-*
        """

rule isomatch_index:
    message: """--- Indexing with isomatch ----"""
    input:
        gtf = rules.stringtie_assembly.output.gtf,
        executable = rules.download_isomatch.output.executable,
        ref_genome = ref_genome
    output:
        gtf = step02.path("isomatch_index/{sample}.gtf"),
        isomx = step02.path("isomatch_index/{sample}.gtf.isomx"),
        isoms = step02.path("isomatch_index/{sample}.gtf.isoms"),
        info = step02.path("isomatch_index/{sample}.gtf.isomx.info.json")
    log: step02.log("isomatch_index/{sample}.log")
    benchmark: step02.benchmark("isomatch_index/{sample}.tsv")
    threads: 1,
    shell: 
        """
        set -euo pipefail
        ln -sf $(readlink -f {input.gtf}) {output.gtf}
        {input.executable} index --ref-fa {input.ref_genome} {output.gtf} --out {output.isomx} 2>&1 | tee {log}
        """

rule isomatch_merge:
    message: """--- Merging with isomatch ----"""
    input:
        gtf = lambda wc: expand(step02.path("isomatch_index/{sample}.gtf"), sample = generate_input_samples_df(wc.dataset, wc.group)["sample_id"].tolist(), dataset = wc.dataset, group = wc.group),
        executable = rules.download_isomatch.output.executable,
        ref_genome = ref_genome
    output:
        gtf = step02.path("isomatch_merge/{dataset}.{group}.merged.gtf.gz"),
        info = step02.path("isomatch_merge/{dataset}.{group}.merged_info.json")
    log: step02.log("isomatch_merge/{dataset}.{group}.merged.log")
    benchmark: step02.benchmark("isomatch_merge/{dataset}.{group}.merged.tsv")
    params:
        out_prefix = str(step02.path("isomatch_merge/{dataset}.{group}")),
        index_dir = str(step02.path("isomatch_index")),
        flags = isomatch_settings["merge_flags"]
    threads: 1,
    shell:
        """
        set -euo pipefail
        {input.executable} merge -d 3 -a 3 {params.flags} --ref-fa {input.ref_genome} --out {params.out_prefix} {input.gtf} 2>&1 | tee {log}
        rm -rf {params.index_dir}/.isomatch-index-*
        """

rule isomatch_chop:
    message: """--- Chopping with isomatch ----"""
    input:
        gtf = rules.isomatch_merge.output.gtf,
        executable = rules.download_isomatch.output.executable
    output:
        gtf = step02.path("isomatch_merge/{dataset}.{group}.chopped.gtf.gz")
    log: step02.log("isomatch_merge/{dataset}.{group}.chopped.log")
    benchmark: step02.benchmark("isomatch_merge/{dataset}.{group}.chopped.tsv")
    params:
        out_prefix = str(step02.path("isomatch_merge/{dataset}.{group}"))
    threads: 1,
    shell:
        """
        set -euo pipefail
        {input.executable} tools chop --out {params.out_prefix} {input.gtf} 2>&1 | tee {log}
        """

rule isomatch_classify:
    message: """--- Classifying with isomatch ----"""
    input:
        gtf = rules.isomatch_chop.output.gtf,
        ref_gtf = rules.isomatch_index_reference.output.gtf,
        executable = rules.download_isomatch.output.executable,
        ref_genome = ref_genome
    output:
        classification = step02.path("isomatch_classify/{dataset}.{group}.classification.txt.gz"),
        annotated_gtf = step02.path("isomatch_classify/{dataset}.{group}.annotated.gtf.gz"),
        info = step02.path("isomatch_classify/{dataset}.{group}.classify_info.json")
    log: step02.log("isomatch_classify/{dataset}.{group}.log")
    benchmark: step02.benchmark("isomatch_classify/{dataset}.{group}.tsv")
    params:
        out_prefix = str(step02.path("isomatch_classify/{dataset}.{group}")),
        merge_dir = lambda wildcards, input: str(Path(input.gtf).parent)
    threads: 1,
    shell:
        """
        set -euo pipefail
        {input.executable} classify --ref-gtf {input.ref_gtf} --ref-fa {input.ref_genome} --out {params.out_prefix} {input.gtf} 2>&1 | tee {log}
        rm -rf {params.merge_dir}/.isomatch-index-*
        """

rule isomatch_fix:
    message: """--- Reconciling isomatch attributes with gffcompare schema ----"""
    input:
        gtf = rules.isomatch_classify.output.annotated_gtf
    output:
        gtf = step02.path("isomatch_classify/{dataset}.{group}.fixed.gtf")
    log: step02.log("isomatch_classify/{dataset}.{group}.fixed.log")
    benchmark: step02.benchmark("isomatch_classify/{dataset}.{group}.fixed.tsv")
    conda: "../envs/r.yaml"
    threads: 1
    script: "../scripts/02a-isomatch_fix.R"

### Common rules
rule manual_filter:
    message: """--- Chromosome and Strand filtration ----"""
    input: get_filter_input
    output: 
        gtf = step02.path("transcriptome_assembly/{dataset}.{group}.{merge_method}.annotated.clean.gtf")
    log: step02.log("transcriptome_assembly/{dataset}.{group}.{merge_method}.annotated.clean.log")
    benchmark: step02.benchmark("transcriptome_assembly/{dataset}.{group}.{merge_method}.annotated.clean.tsv")
    shell: "zcat -f {input} | awk '$1 ~ /^(chr)?([1-9]|1[0-9]|2[0-2]|X|Y)$/ && $7 ~ /^[+-]$/' > {output} 2> {log}"
    # shell: "awk '$1 ~ /^(chr)?([1-9]|1[0-9]|2[0-2]|X|Y)$/ && $7 ~ /^[+-]$/' {input} > {output} 2> {log}"
