################################################################################
## PacBio Processing
################################################################################

## Variables
################################################################################
pacbio_settings = resolve_preset(config, "pacbio_settings", "pacbio_preset")

## Functions
################################################################################
def get_input_hifi_bam(wildcards):
    sample = get_input_sample(wildcards)

    if sample["file_type"] != "hifi_bam":
        raise ValueError(f"{wildcards.sample} is a '{sample['file_type']}' input, so it has no PacBio processing")
    return sample["path"]

def get_flnc_reads(wildcards):
    """Full-length reads of a sample: made by pacbio_refine (hifi_bam) or given in the sample sheet (flnc)"""
    sample = get_input_sample(wildcards)

    if sample["file_type"] == "hifi_bam":
        return rules.pacbio_refine.output.bam
    if sample["file_type"] == "flnc":
        return sample["path"]
    raise ValueError(f"{wildcards.sample} is a '{sample['file_type']}' input, so it has no PacBio alignment")

def get_lima_pair(primers: str) -> str:
    """Name of the primer pair in lima's output file: <5' primer>--<3' primer>."""
    names = [line[1:].split()[0] for line in open(primers) if line.startswith(">")]
    primers_5p = [name for name in names if name.endswith("_5p")]
    primers_3p = [name for name in names if name.endswith("_3p")]
    if len(primers_5p) != 1 or len(primers_3p) != 1:
        raise ValueError(f"{primers} must contain one 5' (*_5p) and one 3' (*_3p) primer, so lima writes a single output")
    return f"{primers_5p[0]}--{primers_3p[0]}"

lima_pair = get_lima_pair(pacbio_settings["primers"])

def get_pacbio_samples(wildcards) -> list:
    """Sample IDs of one dataset x group"""
    return generate_input_samples_df(wildcards.dataset, wildcards.group)["sample_id"].tolist()

def get_pb_setting(key: str, default):
    """Value of `key` in the pb* merge preset of the wildcards"""
    return lambda wildcards: config["merge_settings"][wildcards.merge_method].get(key, default)


## Rules
################################################################################
rule pacbio_hifi_filter:
    message: """--- Keeping the HiFi reads - {wildcards.sample} ----"""
    input:
        bam = get_input_hifi_bam
    output:
        bam = temp(step01.path("PacBio/HiFi/{sample}.hifi.bam"))
    log: step01.logs("PacBio/HiFi/{sample}.log")
    benchmark: step01.benchmark("PacBio/HiFi/{sample}.tsv")
    params:
        min_rq = pacbio_settings["min_rq"]
    threads: 4
    conda: "../envs/samtools.yaml"
    shell: "samtools view -@ {threads} -b -e '[rq]>={params.min_rq}' -o {output.bam} {input.bam} > {log} 2>&1"

rule pacbio_lima:
    message: """--- Removing the cDNA primers with lima - {wildcards.sample} ----"""
    input:
        bam = rules.pacbio_hifi_filter.output.bam,
        primers = pacbio_settings["primers"]
    output:
        bam = temp(step01.path(f"PacBio/lima/{{sample}}.fl.{lima_pair}.bam")),
        pbi = temp(step01.path(f"PacBio/lima/{{sample}}.fl.{lima_pair}.bam.pbi")),
        summary = step01.path("PacBio/lima/{sample}.fl.lima.summary"),
        counts = step01.path("PacBio/lima/{sample}.fl.lima.counts")
    log: step01.logs("PacBio/lima/{sample}.log")
    benchmark: step01.benchmark("PacBio/lima/{sample}.tsv")
    params:
        out_bam = lambda wc, output: output.summary.replace(".lima.summary", ".bam"),
        flags = pacbio_settings["lima_flags"]
    threads: 16
    conda: "../envs/isoseq.yaml"
    shell: "lima {params.flags} -j {threads} {input.bam} {input.primers} {params.out_bam} > {log} 2>&1"

rule pacbio_refine:
    message: """--- Removing poly(A) tails and concatemers - {wildcards.sample} ----"""
    input:
        bam = rules.pacbio_lima.output.bam,
        primers = pacbio_settings["primers"]
    output:
        bam = step01.path("PacBio/FLNC/{sample}.flnc.bam"),
        pbi = step01.path("PacBio/FLNC/{sample}.flnc.bam.pbi"),
        report = step01.path("PacBio/FLNC/{sample}.flnc.filter_summary.report.json")
    log: step01.logs("PacBio/FLNC/{sample}.log")
    benchmark: step01.benchmark("PacBio/FLNC/{sample}.tsv")
    params:
        flags = pacbio_settings["refine_flags"]
    threads: 16
    conda: "../envs/isoseq.yaml"
    shell: "isoseq refine {params.flags} -j {threads} {input.bam} {input.primers} {output.bam} > {log} 2>&1"

rule pbmm2_index:
    message: """--- Indexing the genome for pbmm2 ----"""
    input:
        genome = ref_genome
    output:
        mmi = step01.path("PacBio/pbmm2_index/genome.isoseq.mmi")
    log: step01.logs("PacBio/pbmm2_index/genome.log")
    benchmark: step01.benchmark("PacBio/pbmm2_index/genome.tsv")
    threads: 8
    conda: "../envs/isoseq.yaml"
    shell: "pbmm2 index --preset ISOSEQ -j {threads} {input.genome} {output.mmi} > {log} 2>&1"

rule pbmm2_align:
    message: """--- Aligning the FLNC reads with pbmm2 - {wildcards.sample} ----"""
    input:
        mmi = rules.pbmm2_index.output.mmi,
        reads = get_flnc_reads
    output:
        bam = step01.path("PacBio/pbmm2/{sample}.aligned.bam"),
        bai = step01.path("PacBio/pbmm2/{sample}.aligned.bam.bai")
    log: step01.logs("PacBio/pbmm2/{sample}.log")
    benchmark: step01.benchmark("PacBio/pbmm2/{sample}.tsv")
    params:
        flags = pacbio_settings["pbmm2_flags"]
    threads: 16
    conda: "../envs/isoseq.yaml"
    shell: "pbmm2 align {params.flags} -j {threads} {input.mmi} {input.reads} {output.bam} > {log} 2>&1"

rule pacbio_read_ids:
    message: """--- Listing the aligned reads - {wildcards.sample} ----"""
    input:
        bam = rules.pbmm2_align.output.bam
    output:
        ids = step01.path("PacBio/pbmm2/{sample}.read_ids.tsv.gz")
    log: step01.logs("PacBio/pbmm2/{sample}.read_ids.log")
    conda: "../envs/samtools.yaml"
    shell:
        """
        samtools view -F 0x904 {input.bam} 2> {log} | awk -v s={wildcards.sample} 'BEGIN {{OFS = "\\t"}} {{print $1, s}}' | gzip > {output.ids}
        """

rule pacbio_merge_alignments:
    message: """--- Pooling the PacBio alignments of {wildcards.dataset}.{wildcards.group} ----"""
    input:
        bam = lambda wc: expand(rules.pbmm2_align.output.bam, sample = get_pacbio_samples(wc), dataset = wc.dataset)
    output:
        bam = temp(step01.path("PacBio/collapse/{dataset}.{group}.{merge_method}.aligned.bam"))
    log: step01.logs("PacBio/collapse/{dataset}.{group}.{merge_method}.merge.log")
    benchmark: step01.benchmark("PacBio/collapse/{dataset}.{group}.{merge_method}.merge.tsv")
    threads: 8
    conda: "../envs/samtools.yaml"
    shell: "samtools merge -f -@ {threads} -o {output.bam} {input.bam} > {log} 2>&1"

rule pacbio_collapse:
    message: """--- Collapsing the PacBio reads of {wildcards.dataset}.{wildcards.group} ----"""
    input:
        bam = rules.pacbio_merge_alignments.output.bam
    output:
        gtf = step01.path("PacBio/collapse/{dataset}.{group}.{merge_method}.collapsed.gtf"),
        read_stat = step01.path("PacBio/collapse/{dataset}.{group}.{merge_method}.collapsed.read_stat.txt")
    log: step01.logs("PacBio/collapse/{dataset}.{group}.{merge_method}.collapse.log")
    benchmark: step01.benchmark("PacBio/collapse/{dataset}.{group}.{merge_method}.collapse.tsv")
    params:
        gff = lambda wc, output: output.gtf.replace(".gtf", ".gff"),
        flags = get_pb_setting("collapse_flags", "")
    threads: 16
    conda: "../envs/isoseq.yaml"
    shell:
        """
        isoseq collapse {params.flags} -j {threads} {input.bam} {params.gff} > {log} 2>&1
        mv {params.gff} {output.gtf}
        """

rule pacbio_pigeon_prepare_reference:
    message: """--- Preparing the reference annotation for pigeon ----"""
    input:
        gtf = ref_annotation
    output:
        link_gtf = step01.path("PacBio/pigeon_ref/ref.pigeon.gtf"),
        sorted_gtf = step01.path("PacBio/pigeon_ref/ref.pigeon.sorted.gtf"),
        pgi = step01.path("PacBio/pigeon_ref/ref.pigeon.sorted.gtf.pgi")
    log: step01.logs("PacBio/pigeon_ref/prepare.log")
    benchmark: step01.benchmark("PacBio/pigeon_ref/prepare.tsv")
    conda: "../envs/pigeon.yaml"
    shell:
        """
        ln -sf $(readlink -f {input.gtf}) {output.link_gtf}
        pigeon prepare {output.link_gtf} > {log} 2>&1
        """

rule pacbio_pigeon_prepare:
    message: """--- Preparing the collapsed transcripts for pigeon - {wildcards.dataset}.{wildcards.group} ----"""
    input:
        gtf = rules.pacbio_collapse.output.gtf
    output:
        sorted_gtf = step01.path("PacBio/collapse/{dataset}.{group}.{merge_method}.collapsed.sorted.gtf"),
        pgi = step01.path("PacBio/collapse/{dataset}.{group}.{merge_method}.collapsed.sorted.gtf.pgi")
    log: step01.logs("PacBio/pigeon/{dataset}.{group}.{merge_method}.prepare.log")
    benchmark: step01.benchmark("PacBio/pigeon/{dataset}.{group}.{merge_method}.prepare.tsv")
    conda: "../envs/pigeon.yaml"
    shell: "pigeon prepare {input.gtf} > {log} 2>&1"

rule pacbio_pigeon_classify:
    message: """--- Classifying the collapsed transcripts with pigeon - {wildcards.dataset}.{wildcards.group} ----"""
    input:
        sorted_gtf = rules.pacbio_pigeon_prepare.output.sorted_gtf,
        pgi = rules.pacbio_pigeon_prepare.output.pgi,
        reference = rules.pacbio_pigeon_prepare_reference.output.sorted_gtf,
        reference_pgi = rules.pacbio_pigeon_prepare_reference.output.pgi,
        genome = ref_genome
    output:
        classification = step01.path("PacBio/pigeon/{dataset}.{group}.{merge_method}.pigeon_classification.txt")
    log: step01.logs("PacBio/pigeon/{dataset}.{group}.{merge_method}.classify.log")
    benchmark: step01.benchmark("PacBio/pigeon/{dataset}.{group}.{merge_method}.classify.tsv")
    params:
        out_dir = lambda wc, output: Path(output.classification).parent,
        out_prefix = "{dataset}.{group}.{merge_method}.pigeon"
    threads: 16
    conda: "../envs/pigeon.yaml"
    shell:
        "pigeon classify --gene-id -d {params.out_dir} -o {params.out_prefix} -j {threads} "
        "{input.sorted_gtf} {input.reference} {input.genome} > {log} 2>&1"
        
rule pacbio_build_transcriptome:
    message: """--- Building the PacBio transcriptome - {wildcards.dataset}.{wildcards.group}.{wildcards.merge_method} ----"""
    input:
        gtf = rules.pacbio_collapse.output.gtf,
        classification = rules.pacbio_pigeon_classify.output.classification,
        read_stat = rules.pacbio_collapse.output.read_stat,
        read_ids = lambda wc: expand(rules.pacbio_read_ids.output.ids, sample = get_pacbio_samples(wc), dataset = wc.dataset),
        ref_gtf = ref_annotation
    output:
        gtf = step01.path("PacBio/transcriptome/{dataset}.{group}.{merge_method}.pacbio.gtf")
    log: step01.logs("PacBio/transcriptome/{dataset}.{group}.{merge_method}.log")
    benchmark: step01.benchmark("PacBio/transcriptome/{dataset}.{group}.{merge_method}.tsv")
    params:
        use_reference = get_pb_setting("use_reference", False),
        min_tpm = get_pb_setting("min_tpm", 1.0),
        min_sample_reads = get_pb_setting("min_sample_reads", 2),
        min_sample_cnt = get_pb_setting("min_sample_cnt", 2),
        end_match_bp = get_pb_setting("end_match_bp", 50),
    threads: 1
    conda: "../envs/r.yaml"
    script: "../scripts/S01b-PacBio_Build_Transcriptome.R"