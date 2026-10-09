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
        "stringtie_assembly_settings": build_tool_settings(config, "stringtie_assembly_settings", "stringtie_assembly_preset"),
        "merge_settings": build_tool_settings(config, "merge_settings"),
        "gffcompare_settings": build_tool_settings(config, "gffcompare_settings", "gffcompare_preset"),
    }
)

### Resolved settings used directly in the rule bodies below
stringtie_assembly_settings = resolve_preset(config, "stringtie_assembly_settings", "stringtie_assembly_preset")
merge_settings = config["merge_settings"]
gffcompare_settings = resolve_preset(config, "gffcompare_settings", "gffcompare_preset")

### isomatch Download path
isomatch_version = config.get("isomatch_version", "v0.6.0")
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
    if wildcards.merge_method == "ref":
        return ref_annotation
    elif wildcards.merge_method.startswith("st"):
        return rules.gffcompare.output.gtf
    elif wildcards.merge_method.startswith("iso"):
        return rules.isomatch_fix.output.gtf
    raise ValueError(f"Unknown merge method: {wildcards.merge_method}")

## Rules
################################################################################

### StringTie Rules
rule stringtie_assembly:
    message: """--- StringTie Assembly - {wildcards.sample} ----"""
    input: 
        bam = get_input_bam,
    output: 
        gtf = step02.path("StringTie_Sample_Assembly/{sample}.gtf")
    log: step02.logs("StringTie_Sample_Assembly/{sample}.log")
    benchmark: step02.benchmark("StringTie_Sample_Assembly/{sample}.tsv")
    params:
        guided = f"-G {ref_annotation}" if stringtie_assembly_settings["assembly_guided"] else "",
        flags = stringtie_assembly_settings["assembly_flags"]
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
        gtf = step02.path("StringTie_Merge/{dataset}.{group}.{merge_method}.merged.gtf"),
    log: step02.logs("StringTie_Merge/{dataset}.{group}.{merge_method}.merged.log"),
    benchmark: step02.benchmark("StringTie_Merge/{dataset}.{group}.{merge_method}.merged.tsv")
    params:
        guided = lambda wc: f"-G {ref_annotation}" if merge_settings[wc.merge_method].get("use_reference", False) else "",
        flags = lambda wc: merge_settings[wc.merge_method].get("merge_flags", "")
    threads: 16,
    conda: "../envs/stringtie.yaml",
    shell: "stringtie --merge {params.flags} {params.guided} -p {threads} -o {output} {input} 2>&1 | tee {log}"

rule gffcompare:
    message: """--- Running gffcompare ----"""
    input:
        gtf = rules.stringtie_merge.output.gtf,
    output: 
        gtf = step02.path("gffcompare/{dataset}.{group}.{merge_method}.gffcompare.annotated.gtf"),
        stats = step02.path("gffcompare/{dataset}.{group}.{merge_method}.gffcompare.stats")
    log: step02.logs("gffcompare/{dataset}.{group}.{merge_method}.gffcompare.annotated.log")
    benchmark: step02.benchmark("gffcompare/{dataset}.{group}.{merge_method}.gffcompare.annotated.tsv")
    params: 
        out_name = lambda wc: expand(step02.path("gffcompare/{dataset}.{group}.{merge_method}.gffcompare"), 
            dataset = wc.dataset, group = wc.group, merge_method = wc.merge_method),
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

### Common rules
rule manual_filter:
    message: """--- Chromosome and Strand filtration ----"""
    input: get_filter_input
    output: 
        gtf = step02.path("transcriptome_assembly/{dataset}.{group}.{merge_method}.annotated.clean.gtf")
    log: step02.logs("transcriptome_assembly/{dataset}.{group}.{merge_method}.annotated.clean.log")
    benchmark: step02.benchmark("transcriptome_assembly/{dataset}.{group}.{merge_method}.annotated.clean.tsv")
    shell:
        """
        {{
            echo "##gff-version 2"
            echo "##source ENDome manual_filter ({wildcards.merge_method})"
            echo "##input {input}"
            echo "##date $(date +%Y-%m-%d)"
            zcat -f {input} | awk '$1 ~ /^(chr)?([1-9]|1[0-9]|2[0-2]|X|Y)$/ && $7 ~ /^[+-]$/'
        }} > {output} 2> {log}
        """

### Isomatch rules
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

rule download_isomatch_guide:
    message: "--- Downloading isomatch guide evidence ---"
    output:
        tes = f"tools/isomatch-{isomatch_version}-guides/human.grch38.tes.bed",
        tss = f"tools/isomatch-{isomatch_version}-guides/human.grch38.tss.bed",
    params:
        tes_url = f"https://raw.githubusercontent.com/zhengxinchang/isomatch/{isomatch_version}/evidence/human.grch38.tes.bed",
        tss_url = f"https://raw.githubusercontent.com/zhengxinchang/isomatch/{isomatch_version}/evidence/human.grch38.tss.bed"
    shell:
        """
        curl -fsSL -o {output.tes} {params.tes_url}
        curl -fsSL -o {output.tss} {params.tss_url}
        """

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
    log: step02.logs("isomatch_ref_index/ref_annotation.log")
    benchmark: step02.benchmark("isomatch_ref_index/ref_annotation.tsv")
    params:
        index_dir = lambda wildcards, output: str(Path(output.gtf).parent)
    threads: 1, #rm -rf {params.index_dir}/.isomatch-index-*
    shell:
        """
        ln -sf $(readlink -f {input.gtf}) {output.gtf}
        {input.executable} index --ref-fa {input.ref_genome} {output.gtf} --out {output.isomx} 2>&1 | tee {log}
        """

rule isomatch_filter_sample:
    message: """--- Filtering sample assembly for isomatch - {wildcards.sample} ----"""
    input:
        gtf = rules.stringtie_assembly.output.gtf
    output:
        gtf = step02.path("isomatch_index/{merge_method}/{sample}.filter.gtf")
    log: step02.logs("isomatch_index/{merge_method}/{sample}.filter.log")
    benchmark: step02.benchmark("isomatch_index/{merge_method}/{sample}.filter.tsv")
    params:
        min_fpkm = lambda wc: merge_settings[wc.merge_method].get("min_fpkm", 0),
        min_tpm = lambda wc: merge_settings[wc.merge_method].get("min_tpm", 0)
    conda: "../envs/r.yaml"
    threads: 1
    script: "../scripts/02a-isomatch_filter_sample.R"

rule isomatch_index:
    message: """--- Indexing with isomatch ----"""
    input:
        gtf = rules.isomatch_filter_sample.output.gtf,
        executable = rules.download_isomatch.output.executable,
        ref_genome = ref_genome
    output:
        isomx = step02.path("isomatch_index/{merge_method}/{sample}.filter.gtf.isomx"),
        isoms = step02.path("isomatch_index/{merge_method}/{sample}.filter.gtf.isoms"),
        info = step02.path("isomatch_index/{merge_method}/{sample}.filter.gtf.isomx.info.json")
    log: step02.logs("isomatch_index/{merge_method}/{sample}.index.log")
    benchmark: step02.benchmark("isomatch_index/{merge_method}/{sample}.index.tsv")
    threads: 1,
    shell: 
        """
        {input.executable} index --ref-fa {input.ref_genome} {input.gtf} --out {output.isomx} 2>&1 | tee {log}
        """

rule isomatch_merge:
    message: """--- Merging with isomatch ----"""
    input:
        gtf = lambda wc: expand(step02.path("isomatch_index/{merge_method}/{sample}.filter.gtf"), 
            sample = generate_input_samples_df(wc.dataset, wc.group)["sample_id"].tolist(), 
            dataset = wc.dataset, group = wc.group, merge_method = wc.merge_method),
        index = lambda wc: expand(step02.path("isomatch_index/{merge_method}/{sample}.filter.gtf.isomx"), 
            sample = generate_input_samples_df(wc.dataset, wc.group)["sample_id"].tolist(), 
            dataset = wc.dataset, group = wc.group, merge_method = wc.merge_method), 
        ref_gtf = lambda wc: rules.isomatch_index_reference.output.gtf if merge_settings[wc.merge_method].get("use_reference", False) else [],
        guide_tss = lambda wc: rules.download_isomatch_guide.output.tss if merge_settings[wc.merge_method].get("use_bed_files", False) else [],
        guide_tes = lambda wc: rules.download_isomatch_guide.output.tes if merge_settings[wc.merge_method].get("use_bed_files", False) else [],
        executable = rules.download_isomatch.output.executable,
        ref_genome = ref_genome
    output:
        gtf = step02.path("isomatch_merge/{dataset}.{group}.{merge_method}.merged.gtf.gz"),
        info = step02.path("isomatch_merge/{dataset}.{group}.{merge_method}.merged_info.json"),
        track = step02.path("isomatch_merge/{dataset}.{group}.{merge_method}.track.tsv.gz"),
        present_absent = step02.path("isomatch_merge/{dataset}.{group}.{merge_method}.present_absent.tsv.gz")
    log: step02.logs("isomatch_merge/{dataset}.{group}.{merge_method}.merged.log")
    benchmark: step02.benchmark("isomatch_merge/{dataset}.{group}.{merge_method}.merged.tsv")
    params:
        out_prefix = str(step02.path("isomatch_merge/{dataset}.{group}.{merge_method}")),
        index_dir = str(step02.path("isomatch_index")),
        ref_index_dir = str(step02.path("isomatch_ref_index")),
        flags = lambda wc: merge_settings[wc.merge_method].get("merge_flags", ""),
        guide_flags = lambda wc, input: f"--guide-tss {input.guide_tss} --guide-tes {input.guide_tes}" if merge_settings[wc.merge_method].get("use_bed_files", False) else []
    threads: 1,  # rm -rf {params.index_dir}/.isomatch-index-* {params.ref_index_dir}/.isomatch-index-*
    shell:
        """
        {input.executable} merge -d 3 -a 3 {params.flags} {params.guide_flags} --ref-fa {input.ref_genome} --out {params.out_prefix} {input.gtf} {input.ref_gtf} 2>&1 | tee {log}
        """
       
rule isomatch_chop:
    message: """--- Chopping with isomatch ----"""
    input:
        gtf = rules.isomatch_merge.output.gtf,
        executable = rules.download_isomatch.output.executable
    output:
        gtf = step02.path("isomatch_merge/{dataset}.{group}.{merge_method}.chopped.gtf.gz")
    log: step02.logs("isomatch_merge/{dataset}.{group}.{merge_method}.chopped.log")
    benchmark: step02.benchmark("isomatch_merge/{dataset}.{group}.{merge_method}.chopped.tsv")
    params:
        out_prefix = str(step02.path("isomatch_merge/{dataset}.{group}.{merge_method}"))
    threads: 1,
    shell:
        """
        {input.executable} tools chop --out {params.out_prefix} {input.gtf} 2>&1 | tee {log}
        """

rule isomatch_classify:
    message: """--- Classifying with isomatch ----"""
    input:
        gtf = rules.isomatch_chop.output.gtf,
        ref_gtf = rules.isomatch_index_reference.output.gtf,
        guide_tss = lambda wc: rules.download_isomatch_guide.output.tss if merge_settings[wc.merge_method].get("use_bed_files", False) else [],
        guide_tes = lambda wc: rules.download_isomatch_guide.output.tes if merge_settings[wc.merge_method].get("use_bed_files", False) else [],
        executable = rules.download_isomatch.output.executable,
        ref_genome = ref_genome
    output:
        classification = step02.path("isomatch_classify/{dataset}.{group}.{merge_method}.classification.txt.gz"),
        annotated_gtf = step02.path("isomatch_classify/{dataset}.{group}.{merge_method}.annotated.gtf.gz"),
        info = step02.path("isomatch_classify/{dataset}.{group}.{merge_method}.classify_info.json")
    log: step02.logs("isomatch_classify/{dataset}.{group}.{merge_method}.log")
    benchmark: step02.benchmark("isomatch_classify/{dataset}.{group}.{merge_method}.tsv")
    params:
        out_prefix = str(step02.path("isomatch_classify/{dataset}.{group}.{merge_method}")),
        merge_dir = lambda wc, input: str(Path(input.gtf).parent),
        guide_flags = lambda wc, input: f"--guide-tss {input.guide_tss} --guide-tes {input.guide_tes}" if merge_settings[wc.merge_method].get("use_bed_files", False) else []
    threads: 1, #rm -rf {params.merge_dir}/.isomatch-index-*,
    shell:
        """
        {input.executable} classify {params.guide_flags} --ref-gtf {input.ref_gtf} --ref-fa {input.ref_genome} --out {params.out_prefix} {input.gtf} 2>&1 | tee {log}
        """
        
rule isomatch_fix:
    message: """--- Reconciling isomatch attributes with gffcompare schema ----"""
    input:
        gtf = rules.isomatch_classify.output.annotated_gtf,
        present_absent = rules.isomatch_merge.output.present_absent,
        track = rules.isomatch_merge.output.track
    output:
        gtf = step02.path("isomatch_classify/{dataset}.{group}.{merge_method}.fixed.gtf")
    log: step02.logs("isomatch_classify/{dataset}.{group}.{merge_method}.fixed.log")
    benchmark: step02.benchmark("isomatch_classify/{dataset}.{group}.{merge_method}.fixed.tsv")
    params:
        ref_source = lambda wc: Path(str(rules.isomatch_index_reference.output.gtf)).name if merge_settings[wc.merge_method].get("use_reference", False) else "",
        min_sample_cnt = lambda wc: merge_settings[wc.merge_method].get("min_sample_cnt", 0)
    conda: "../envs/r.yaml"
    threads: 1
    script: "../scripts/02b-isomatch_fix.R"

register_test_data_link(rules.isomatch_classify.output.annotated_gtf)
register_test_data_link(rules.isomatch_merge.output.present_absent)
register_test_data_link(rules.isomatch_merge.output.track)
