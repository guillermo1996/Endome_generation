################################################################################
## 06 - UTR Quantification
################################################################################

## Debug: benchmark loading time
_start_time = time.perf_counter()

## Variables
################################################################################
### Register the step: hash suffix, parameters.yaml rule, path helpers.
step06 = register_step(
    name="06-UTR_Quantification",
    params={
        "scutrquant_presets": build_tool_settings(config, "scutrquant_settings", "scutrquant_preset"),
    },
)

### Resolved settings used directly in the rule bodies below
scutrquant_settings = resolve_preset(config, "scutrquant_settings", "scutrquant_preset")
scutrquant_dataset_name = scutrquant_settings["dataset_name"]
scutrquant_genome = gtxcutr_settings["genome"]

### output_type may be a scalar ("txs") or a list -> normalise to a list
_output_type = scutrquant_settings.get("output_type", ["txs"])
scutrquant_output_type = _output_type if isinstance(_output_type, list) else [_output_type]

### tmp_dir / bx_whitelist are optional -> fall back sensibly if not configured
scutrquant_tmp_dir = scutrquant_settings.get("tmp_dir") or config.get("scutrquant_tmp_dir") or "/tmp"
scutrquant_bx_whitelist = scutrquant_settings.get("bx_whitelist") or config.get("scutrquant_bx_whitelist") or None

### Which BUS file feeds counting (barcode-corrected or raw), set by correct_bus
_count_bus = "output.corrected.sorted.bus" if scutrquant_settings["correct_bus"] else "output.sorted.bus"

### Input samples
input_samples = pd.read_csv(config["input_scUTRquant_samples"], index_col="sample_id")

### scUTRquant module. The "endome" target (config/utrome_config.yaml) is a stub:
### every reused rule below has its I/O overridden, so scUTRquant never reads the
### target's file paths — the stub only satisfies scUTRquant's target validation.
scutrquant_config = {
    "dataset_name": scutrquant_dataset_name,
    "sample_file": config["input_scUTRquant_samples"],
    "sample_regex": ".*",
    "targets_config": "config/utrome_config.yaml",
    "target": "endome",
    "tech": scutrquant_settings["tech"],
    "strand": scutrquant_settings["strand"],
    "min_umis": scutrquant_settings["min_umis"],
    "correct_bus": scutrquant_settings["correct_bus"],
    "bx_whitelist": scutrquant_bx_whitelist,
    "cell_annots": None,
    "cell_annots_key": "cell_id",
    "exclude_unannotated_cells": False,
    "output_type": scutrquant_output_type,
    "output_format": ["sce"],
    "use_hdf5": False,
    "include_reports": True,
    "tmp_dir": scutrquant_tmp_dir,
}

module scUTRquant:
    snakefile: github("Mayrlab/scUTRquant", path="Snakefile", tag="v0.5.1")
    config: scutrquant_config
    skip_validation: True

## Functions (read the sample sheet; output/input paths are written inline)
################################################################################
def get_sequence_files(wildcards):
    return input_samples.files[wildcards.sample_id].split(";")

def get_file_type(wildcards):
    return "--bam" if input_samples.file_type[wildcards.sample_id] == "bam" else ""

## Rules
################################################################################
rule kallisto_index:
    message: "--- Building kallisto index for the ENDome ---"
    input:
        fa = step05.path("gtxcutr/{prefix}.{orf_filter}.gtxcutr.w{width}.{txEnd}.fa.gz")
    output:
        kdx = step06.path("kallisto_index/{prefix}.{orf_filter}.gtxcutr.w{width}.{txEnd}.kdx")
    log: step06.log("kallisto_index//{prefix}.{orf_filter}.gtxcutr.w{width}.{txEnd}.log")
    benchmark: step06.benchmark("kallisto_index/{prefix}.{orf_filter}.gtxcutr.w{width}.{txEnd}.tsv")
    params:
        tmp_dir = lambda wildcards, output: output.kdx + ".tmp"
    conda: "../envs/kallisto.yaml"
    threads: 4
    shell:
        "kallisto index --threads {threads} --tmp {params.tmp_dir} -i {output.kdx} {input.fa} 2>&1 | tee {log}"

use rule generate_tx_merge from scUTRquant as sq_generate_tx_merge with:
    input:
        tsv = step05.path("gtxcutr/{prefix}.{orf_filter}.gtxcutr.w{width}.{txEnd}.merge.tsv")
    output:
        step06.path("scUTRquant/{prefix}.{orf_filter}.w{width}.{txEnd}.tx_merge.tsv")

use rule generate_gene_merge from scUTRquant as sq_generate_gene_merge with:
    input:
        tsv = step05.path("gtxcutr/{prefix}.{orf_filter}.gtxcutr.w{width}.{txEnd}.merge.tsv")
    output:
        step06.path("scUTRquant/{prefix}.{orf_filter}.w{width}.{txEnd}.gene_merge.tsv")

use rule kallisto_bus from scUTRquant as sq_kallisto_bus with:
    input:
        idx = rules.kallisto_index.output.kdx,
        files = get_sequence_files
    output:
        bus = temp(step06.path("kallisto/{prefix}.{orf_filter}.w{width}.{txEnd}/{sample_id}/output.bus")),
        ec = step06.path("kallisto/{prefix}.{orf_filter}.w{width}.{txEnd}/{sample_id}/matrix.ec"),
        tx = step06.path("kallisto/{prefix}.{orf_filter}.w{width}.{txEnd}/{sample_id}/transcripts.txt")
    params:
        tech = scutrquant_settings["tech"],
        strand = scutrquant_settings["strand"],
        bam = get_file_type

use rule bustools_sort from scUTRquant as sq_bustools_sort with:
    input:
        step06.path("kallisto/{prefix}.{orf_filter}.w{width}.{txEnd}/{sample_id}/output.bus")
    output:
        step06.path("kallisto/{prefix}.{orf_filter}.w{width}.{txEnd}/{sample_id}/output.sorted.bus")
    params:
        tmpDir = lambda wildcards, output: output[0] + ".tmp"

use rule bustools_whitelist from scUTRquant as sq_bustools_whitelist with:
    input:
        step06.path("kallisto/{prefix}.{orf_filter}.w{width}.{txEnd}/{sample_id}/output.sorted.bus")
    output:
        step06.path("kallisto/{prefix}.{orf_filter}.w{width}.{txEnd}/{sample_id}/whitelist.txt")

use rule bustools_correct from scUTRquant as sq_bustools_correct with:
    input:
        bus = step06.path("kallisto/{prefix}.{orf_filter}.w{width}.{txEnd}/{sample_id}/output.sorted.bus"),
        bxs = scutrquant_bx_whitelist if scutrquant_bx_whitelist else step06.path("kallisto/{prefix}.{orf_filter}.w{width}.{txEnd}/{sample_id}/whitelist.txt")
    output:
        temp(step06.path("kallisto/{prefix}.{orf_filter}.w{width}.{txEnd}/{sample_id}/output.corrected.bus"))

use rule bustools_correct_sort from scUTRquant as sq_bustools_correct_sort with:
    input:
        step06.path("kallisto/{prefix}.{orf_filter}.w{width}.{txEnd}/{sample_id}/output.corrected.bus")
    output:
        temp(step06.path("kallisto/{prefix}.{orf_filter}.w{width}.{txEnd}/{sample_id}/output.corrected.sorted.bus"))
    params:
        tmpDir = lambda wildcards, output: output[0] + ".tmp"

use rule bustools_count_txs from scUTRquant as sq_bustools_count_txs with:
    input:
        bus = step06.path("kallisto/{prefix}.{orf_filter}.w{width}.{txEnd}/{sample_id}/" + _count_bus),
        txs = step06.path("kallisto/{prefix}.{orf_filter}.w{width}.{txEnd}/{sample_id}/transcripts.txt"),
        ec = step06.path("kallisto/{prefix}.{orf_filter}.w{width}.{txEnd}/{sample_id}/matrix.ec"),
        merge = step06.path("scUTRquant/{prefix}.{orf_filter}.w{width}.{txEnd}.tx_merge.tsv")
    output:
        mtx = step06.path("kallisto/{prefix}.{orf_filter}.w{width}.{txEnd}/{sample_id}/txs.mtx"),
        txs = step06.path("kallisto/{prefix}.{orf_filter}.w{width}.{txEnd}/{sample_id}/txs.genes.txt"),
        bxs = step06.path("kallisto/{prefix}.{orf_filter}.w{width}.{txEnd}/{sample_id}/txs.barcodes.txt")

use rule bustools_count_genes from scUTRquant as sq_bustools_count_genes with:
    input:
        bus = step06.path("kallisto/{prefix}.{orf_filter}.w{width}.{txEnd}/{sample_id}/" + _count_bus),
        txs = step06.path("kallisto/{prefix}.{orf_filter}.w{width}.{txEnd}/{sample_id}/transcripts.txt"),
        ec = step06.path("kallisto/{prefix}.{orf_filter}.w{width}.{txEnd}/{sample_id}/matrix.ec"),
        merge = step06.path("scUTRquant/{prefix}.{orf_filter}.w{width}.{txEnd}.gene_merge.tsv")
    output:
        mtx = step06.path("kallisto/{prefix}.{orf_filter}.w{width}.{txEnd}/{sample_id}/genes.mtx"),
        txs = step06.path("kallisto/{prefix}.{orf_filter}.w{width}.{txEnd}/{sample_id}/genes.genes.txt"),
        bxs = step06.path("kallisto/{prefix}.{orf_filter}.w{width}.{txEnd}/{sample_id}/genes.barcodes.txt")

use rule mtxs_to_sce_txs from scUTRquant as sq_mtxs_to_sce_txs with:
    input:
        bxs = expand(str(step06.path("kallisto/{prefix}.{orf_filter}.w{width}.{txEnd}/{sample_id}/txs.barcodes.txt")), sample_id=input_samples.index.values, allow_missing=True),
        txs = expand(str(step06.path("kallisto/{prefix}.{orf_filter}.w{width}.{txEnd}/{sample_id}/txs.genes.txt")), sample_id=input_samples.index.values, allow_missing=True),
        mtxs = expand(str(step06.path("kallisto/{prefix}.{orf_filter}.w{width}.{txEnd}/{sample_id}/txs.mtx")), sample_id=input_samples.index.values, allow_missing=True),
        gtf = step05.path("gtxcutr/{prefix}.{orf_filter}.gtxcutr.w{width}.{txEnd}.gtf"),
        tx_annots = [],
        cell_annots = []
    output:
        sce = step06.path("scUTRquant/{prefix}.{orf_filter}.w{width}.{txEnd}." + scutrquant_dataset_name + ".txs.Rds")
    params:
        genome = scutrquant_genome,
        sample_ids = list(input_samples.index.values),
        min_umis = scutrquant_settings["min_umis"],
        cell_annots_key = "cell_id",
        exclude_unannotated_cells = False,
        tmp_dir = scutrquant_tmp_dir,
        use_hdf5 = False

use rule mtxs_to_sce_genes from scUTRquant as sq_mtxs_to_sce_genes with:
    input:
        bxs = expand(str(step06.path("kallisto/{prefix}.{orf_filter}.w{width}.{txEnd}/{sample_id}/genes.barcodes.txt")), sample_id=input_samples.index.values, allow_missing=True),
        genes = expand(str(step06.path("kallisto/{prefix}.{orf_filter}.w{width}.{txEnd}/{sample_id}/genes.genes.txt")), sample_id=input_samples.index.values, allow_missing=True),
        mtxs = expand(str(step06.path("kallisto/{prefix}.{orf_filter}.w{width}.{txEnd}/{sample_id}/genes.mtx")), sample_id=input_samples.index.values, allow_missing=True),
        gtf = step05.path("gtxcutr/{prefix}.{orf_filter}.gtxcutr.w{width}.{txEnd}.gtf"),
        gene_annots = [],
        cell_annots = []
    output:
        sce = step06.path("scUTRquant/{prefix}.{orf_filter}.w{width}.{txEnd}." + scutrquant_dataset_name + ".genes.Rds")
    params:
        genome = scutrquant_genome,
        sample_ids = list(input_samples.index.values),
        min_umis = scutrquant_settings["min_umis"],
        cell_annots_key = "cell_id",
        exclude_unannotated_cells = False,
        tmp_dir = scutrquant_tmp_dir,
        use_hdf5 = False

use rule report_umis_per_cell from scUTRquant as sq_report_umis_per_cell with:
    input:
        bxs = step06.path("kallisto/{prefix}.{orf_filter}.w{width}.{txEnd}/{sample_id}/txs.barcodes.txt"),
        txs = step06.path("kallisto/{prefix}.{orf_filter}.w{width}.{txEnd}/{sample_id}/txs.genes.txt"),
        mtx = step06.path("kallisto/{prefix}.{orf_filter}.w{width}.{txEnd}/{sample_id}/txs.mtx")
    output:
        step06.path("qc/{prefix}.{orf_filter}.w{width}.{txEnd}/{sample_id}.umi_count.html")
    params:
        min_umis = scutrquant_settings["min_umis"],
        sq_version = "0.5.1"

## Debug: benchmark loading time
_log(f"\t+ {step06.name} imported in {time.perf_counter() - _start_time:.3f}s")