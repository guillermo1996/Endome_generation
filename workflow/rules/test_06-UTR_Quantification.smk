################################################################################
## 06 - UTR Quantification (scUTRquant on scRNA/snRNA-seq)
################################################################################
## Quantifies UTR-isoform usage of the generated ENDome on single-cell data by
## REUSING the scUTRquant pipeline (https://github.com/Mayrlab/scUTRquant).
##
## scUTRquant is imported as a Snakemake `module` straight from GitHub (no
## vendoring), and its rules are reused via `use rule X from scUTRquant as sq_X
## with:` so that all outputs land inside this step's (hashed) folder, keeping
## the register_step / parameters.yaml conventions of steps 01-05.
##
## The scUTRquant "target" transcriptome IS the ENDome: the gtxcutr GTF and
## merge table from step 05 map directly onto a scUTRquant target (the merge
## table header `tx_in / tx_out / gene_out` is exactly what scUTRquant's
## generate_{tx,gene}_merge consume). The kallisto index is built here from the
## gtxcutr FASTA. The target YAML (config/utrome_config.yaml) is GENERATED at
## parse time from the step-05 output paths, so it never needs hand-editing.

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

### Resolved settings used directly below
scutrquant_settings = resolve_preset(config, "scutrquant_settings", "scutrquant_preset")
scutrquant_dataset_name = config["scutrquant_dataset_name"]
scutrquant_tmp_dir = config["scutrquant_tmp_dir"]

### Single-cell sample sheet (sample_id, file_type, files)
samples = pd.read_csv(config["scutrquant_sample_file"], index_col="sample_id")

### Genome is inherited from the step-05 gtxcutr settings (avoid duplication).
scutrquant_genome = gtxcutr_settings["genome"]

### Which count matrices to produce. For DTU on UTR-isoform bins (DRIMSeq) only
### "txs" is needed; "genes" is optional. Default: txs only ("minimal pipeline").
scutrquant_output_type = config.get("scutrquant_output_type", ["txs"])

### Barcode whitelist. None -> derive it from the data with `bustools whitelist`.
### For 10x libraries the recommended approach is the official 10x whitelist
### (e.g. 3M-february-2018.txt for 10x v3); set scutrquant_bx_whitelist to its path.
scutrquant_bx_whitelist = config.get("scutrquant_bx_whitelist") or None

## Output-path templates (wildcards match step 05: prefix={dataset}.{group})
################################################################################
ENDOME = "{prefix}.{orf_filter}.w{width}.{txEnd}"
GTXCUTR = "{prefix}.{orf_filter}.gtxcutr.w{width}.{txEnd}"

def sq_path(rel):
    """Per-sample kallisto/bustools work area inside the step folder."""
    return str(step06.path("kallisto/" + ENDOME + "/{sample_id}/" + rel))

def utr_path(rel):
    """Per-ENDome UTR merge tables."""
    return str(step06.path("utrs/" + ENDOME + "/" + rel))

def sce_path(rel):
    """SingleCellExperiment outputs (aggregated over samples)."""
    return str(step06.path("sce/" + ENDOME + "/" + rel))

def qc_path(rel):
    """Per-sample QC reports."""
    return str(step06.path("qc/" + ENDOME + "/" + rel))

def step05_gtxcutr(ext):
    """A step-05 gtxcutr output file (gtf / fa.gz / merge.tsv)."""
    return str(step05.path("gtxcutr/" + GTXCUTR + "." + ext))

kallisto_index_out = str(step06.path("kallisto_index/" + GTXCUTR + ".kdx"))

## scUTRquant module configuration
################################################################################
## scUTRquant requires a "targets" YAML and a `target` name. Normally this maps
## a target to a prebuilt UTRome (gtf/kdx/merge_tsv). Here every reused rule has
## its inputs/outputs OVERRIDDEN to point at the step-05/step-06 ENDome files
## directly, so scUTRquant never reads the target paths. The target file
## (config/utrome_config.yaml) is therefore just a static stub whose only job is
## to satisfy scUTRquant's parse-time validation (`target` must be a known key).
scutrquant_targets_config = "config/utrome_config.yaml"
scutrquant_target = "endome"

scutrquant_config = {
    "dataset_name": scutrquant_dataset_name,
    "sample_file": config["scutrquant_sample_file"],
    "sample_regex": ".*",
    "targets_config": scutrquant_targets_config,
    "target": scutrquant_target,
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

## Functions (mirror scUTRquant input helpers, bound to our sample sheet)
################################################################################
def get_sequence_files(wildcards):
    return samples.files[wildcards.sample_id].split(";")

def get_file_type(wildcards):
    return "--bam" if samples.file_type[wildcards.sample_id] == "bam" else ""

## Rules
################################################################################
## Build the kallisto index from the gtxcutr FASTA. scUTRquant targets normally
## ship a prebuilt index; we build ours with the SAME kallisto build scUTRquant
## uses for `kallisto bus` (see envs/scutrquant-kallisto-bustools.yaml) so the
## index format is compatible.
rule kallisto_index:
    message: "--- Building kallisto index for the ENDome ---"
    input:
        fa = step05_gtxcutr("fa.gz")
    output:
        kdx = kallisto_index_out
    log: str(step06.log("kallisto_index/" + GTXCUTR + ".log"))
    benchmark: str(step06.benchmark("kallisto_index/" + GTXCUTR + ".tsv"))
    conda: "../envs/scutrquant-kallisto-bustools.yaml"
    threads: 4
    shell:
        "kallisto index -i {output.kdx} {input.fa} 2>&1 | tee {log}"

## --- Reused scUTRquant rules (I/O rewired into the step folder) -------------

## Merge tables (tx-level and gene-level) from the gtxcutr merge table.
use rule generate_tx_merge from scUTRquant as sq_generate_tx_merge with:
    input:
        tsv = step05_gtxcutr("merge.tsv")
    output:
        utr_path("tx_merge.tsv")

use rule generate_gene_merge from scUTRquant as sq_generate_gene_merge with:
    input:
        tsv = step05_gtxcutr("merge.tsv")
    output:
        utr_path("gene_merge.tsv")

## Pseudoalignment.
use rule kallisto_bus from scUTRquant as sq_kallisto_bus with:
    input:
        idx = kallisto_index_out,
        files = get_sequence_files
    output:
        bus = temp(sq_path("output.bus")),
        ec = sq_path("matrix.ec"),
        tx = sq_path("transcripts.txt")
    params:
        tech = scutrquant_settings["tech"],
        strand = scutrquant_settings["strand"],
        bam = get_file_type

## Barcode processing.
use rule bustools_sort from scUTRquant as sq_bustools_sort with:
    input:
        sq_path("output.bus")
    output:
        sq_path("output.sorted.bus")
    params:
        tmpDir = lambda wcs: scutrquant_tmp_dir + "/bs-sort-" + wcs.sample_id

## Derive a whitelist from the data. Only used (and only added to the DAG) when
## no fixed whitelist is configured; otherwise the official whitelist is used.
use rule bustools_whitelist from scUTRquant as sq_bustools_whitelist with:
    input:
        sq_path("output.sorted.bus")
    output:
        sq_path("whitelist.txt")

use rule bustools_correct from scUTRquant as sq_bustools_correct with:
    input:
        bus = sq_path("output.sorted.bus"),
        bxs = scutrquant_bx_whitelist if scutrquant_bx_whitelist else sq_path("whitelist.txt")
    output:
        temp(sq_path("output.corrected.bus"))

use rule bustools_correct_sort from scUTRquant as sq_bustools_correct_sort with:
    input:
        sq_path("output.corrected.bus")
    output:
        temp(sq_path("output.corrected.sorted.bus"))
    params:
        tmpDir = lambda wcs: scutrquant_tmp_dir + "/bs-corrsort-" + wcs.sample_id

## Which BUS file feeds counting (barcode-corrected or not).
_count_bus = (sq_path("output.corrected.sorted.bus")
              if scutrquant_settings["correct_bus"]
              else sq_path("output.sorted.bus"))

## Counting (UTR-isoform level and gene level).
use rule bustools_count_txs from scUTRquant as sq_bustools_count_txs with:
    input:
        bus = _count_bus,
        txs = sq_path("transcripts.txt"),
        ec = sq_path("matrix.ec"),
        merge = utr_path("tx_merge.tsv")
    output:
        mtx = sq_path("txs.mtx"),
        txs = sq_path("txs.genes.txt"),
        bxs = sq_path("txs.barcodes.txt")

use rule bustools_count_genes from scUTRquant as sq_bustools_count_genes with:
    input:
        bus = _count_bus,
        txs = sq_path("transcripts.txt"),
        ec = sq_path("matrix.ec"),
        merge = utr_path("gene_merge.tsv")
    output:
        mtx = sq_path("genes.mtx"),
        txs = sq_path("genes.genes.txt"),
        bxs = sq_path("genes.barcodes.txt")

## --- Aggregate per-sample matrices into SingleCellExperiment objects --------
use rule mtxs_to_sce_txs from scUTRquant as sq_mtxs_to_sce_txs with:
    input:
        bxs = expand(sq_path("txs.barcodes.txt"), sample_id=samples.index.values, allow_missing=True),
        txs = expand(sq_path("txs.genes.txt"), sample_id=samples.index.values, allow_missing=True),
        mtxs = expand(sq_path("txs.mtx"), sample_id=samples.index.values, allow_missing=True),
        gtf = step05_gtxcutr("gtf"),
        tx_annots = [],
        cell_annots = []
    output:
        sce = sce_path(scutrquant_dataset_name + ".txs.Rds")
    params:
        genome = scutrquant_genome,
        sample_ids = list(samples.index.values),
        min_umis = scutrquant_settings["min_umis"],
        cell_annots_key = "cell_id",
        exclude_unannotated_cells = False,
        tmp_dir = scutrquant_tmp_dir,
        use_hdf5 = False

use rule mtxs_to_sce_genes from scUTRquant as sq_mtxs_to_sce_genes with:
    input:
        bxs = expand(sq_path("genes.barcodes.txt"), sample_id=samples.index.values, allow_missing=True),
        genes = expand(sq_path("genes.genes.txt"), sample_id=samples.index.values, allow_missing=True),
        mtxs = expand(sq_path("genes.mtx"), sample_id=samples.index.values, allow_missing=True),
        gtf = step05_gtxcutr("gtf"),
        gene_annots = [],
        cell_annots = []
    output:
        sce = sce_path(scutrquant_dataset_name + ".genes.Rds")
    params:
        genome = scutrquant_genome,
        sample_ids = list(samples.index.values),
        min_umis = scutrquant_settings["min_umis"],
        cell_annots_key = "cell_id",
        exclude_unannotated_cells = False,
        tmp_dir = scutrquant_tmp_dir,
        use_hdf5 = False

## --- Per-sample QC report ---------------------------------------------------
use rule report_umis_per_cell from scUTRquant as sq_report_umis_per_cell with:
    input:
        bxs = sq_path("txs.barcodes.txt"),
        txs = sq_path("txs.genes.txt"),
        mtx = sq_path("txs.mtx")
    output:
        qc_path("{sample_id}.umi_count.html")
    params:
        min_umis = scutrquant_settings["min_umis"],
        sq_version = "0.5.1"

## Debug: benchmark loading time
_log(f"\t+ {step06.name} imported in {time.perf_counter() - _start_time:.3f}s")
