################################################################################
## 06 (ALT) - UTR Quantification — SINGLE-ENDome, minimal-code variant
################################################################################
## Alternative implementation of step 06 for evaluation against the combinatorial
## 06-UTR_Quantification.smk. Instead of rewiring every scUTRquant rule's I/O,
## this version quantifies ONE chosen ENDome and reuses scUTRquant WHOLESALE via
## `use rule * from ...`, relocating its native data/ layout under the step folder
## with the module `prefix:` directive.
##
## Trade-offs vs the combinatorial version:
##   + tiny: one `use rule *` + one kallisto_index rule + a 1-entry target YAML.
##   + reuses scUTRquant's own input functions / output names verbatim.
##   - NO combinatorial fan-out: it processes a single ENDome (the first value of
##     each combinatorial config list). Multiple variants would collide on the
##     single scUTRquant {target}.
##   - layout is scUTRquant-native (data/<type>/{target}/...), not the per-combo
##     filenames of the main step.
##
## Output root: {main_output_path}/06-UTR_quantification_ALT/{ENDOME_NAME}/
##
## Distinct module/rule names (scUTRquant_alt, kallisto_index_alt, sqalt_*) let
## this file coexist with 06-UTR_Quantification.smk. To evaluate, add
##   include: "rules/06-UTR_quantification_ALT.smk"
## to the Snakefile and run:  snakemake utr_quant_alt
##
## Why this reuses scUTRquant's kallisto|bustools verbatim: the index is built
## with the SAME kallisto (merv::kallisto=0.46.2sq, see the conda env) so its
## format matches scUTRquant's `kallisto bus`.

## Debug: benchmark loading time
_start_time = time.perf_counter()

## Variables — pick the single ENDome (first value of each combinatorial list)
################################################################################
def _first(value, default):
    if isinstance(value, list):
        return value[0] if value else default
    return value if value is not None else default

_alt_dataset = _first(config.get("input_dataset"), "Ebbert")
_alt_group = _first(config.get("input_group"), "control")
_alt_orf = _first(config.get("orf_filter_preset"), "pc")
_alt_width = _first(config.get("trunc_width"), 500)
_alt_txEnd = _first(config.get("trunc_site"), "3p")

ENDOME_NAME = f"{_alt_dataset}.{_alt_group}.{_alt_orf}.w{_alt_width}.{_alt_txEnd}"
GTXCUTR_STEM = f"{_alt_dataset}.{_alt_group}.{_alt_orf}.gtxcutr.w{_alt_width}.{_alt_txEnd}"

scutrquant_settings = resolve_preset(config, "scutrquant_settings", "scutrquant_preset")
scutrquant_sample_file = config["input_scUTRquant_samples"]
scutrquant_dataset_name = scutrquant_settings["dataset_name"]
# tmp_dir / bx_whitelist are optional; fall back sensibly if not configured.
scutrquant_tmp_dir = scutrquant_settings.get("tmp_dir") or config.get("scutrquant_tmp_dir") or "/tmp"
scutrquant_bx_whitelist = scutrquant_settings.get("bx_whitelist") or config.get("scutrquant_bx_whitelist") or None
# output_type may be a scalar ("txs") or a list; normalise to a list.
_alt_output_type = scutrquant_settings.get("output_type", ["txs"])
scutrquant_output_type = _alt_output_type if isinstance(_alt_output_type, list) else [_alt_output_type]
scutrquant_genome = gtxcutr_settings["genome"]

## Output root for this ALT step and the concrete step-05 inputs
################################################################################
alt_dir = Path(main_output_path) / "06-UTR_quantification_ALT" / ENDOME_NAME

def _step05_file(ext):
    """Concrete step-05 gtxcutr file for the chosen ENDome, as an ABSOLUTE path.
    Absolute so the module `prefix` skips it (prefix only rewrites relative paths)
    — otherwise these step-05 inputs would be wrongly relocated under alt_dir."""
    p = str(step05.path(f"gtxcutr/{GTXCUTR_STEM}.{ext}")).format(dataset=_alt_dataset, group=_alt_group)
    return os.path.abspath(p)

# Absolute, for the same reason: it is consumed as scUTRquant's target 'kdx'.
alt_kdx = os.path.abspath(str(alt_dir / "kallisto_index" / f"{ENDOME_NAME}.kdx"))

## Generate the (single-entry) scUTRquant target YAML at parse time
################################################################################
## scUTRquant reads these paths for real here (we do NOT override rule inputs).
## get_target_file is a *callable*, so Snakemake never applies the module prefix
## to it — these stay pointing at the step-05 / index files as written.
_alt_target = "endome"
_alt_targets_yaml = {
    _alt_target: {
        "path": "",
        "genome": scutrquant_genome,
        "gtf": _step05_file("gtf"),
        "kdx": alt_kdx,
        "merge_tsv": _step05_file("merge.tsv"),
        "tx_annots": None,
        "gene_annots": None,
        "tx_annots_csv": None,
        "gene_annots_csv": None,
    }
}
scutrquant_alt_targets_config = "config/utrome_config_ALT.yaml"
with open(scutrquant_alt_targets_config, "w") as _f:
    yaml.safe_dump(_alt_targets_yaml, _f, sort_keys=False, default_flow_style=False)

## scUTRquant module configuration
################################################################################
scutrquant_alt_config = {
    "dataset_name": scutrquant_dataset_name,
    "sample_file": scutrquant_sample_file,
    "sample_regex": ".*",
    "targets_config": scutrquant_alt_targets_config,
    "target": _alt_target,
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

## Import ALL scUTRquant rules, relocating their native data/ layout under alt_dir.
## prefix is skipped for absolute paths and callables (get_target_file inputs), so
## only scUTRquant's relative OUTPUTS (data/..., qc/...) are moved into alt_dir.
module scUTRquant_alt:
    snakefile: github("Mayrlab/scUTRquant", path="Snakefile", tag="v0.5.1")
    config: scutrquant_alt_config
    prefix: str(alt_dir)
    skip_validation: True

use rule * from scUTRquant_alt as sqalt_*

## The one rule scUTRquant lacks: build the index from the gtxcutr FASTA.
## Single ENDome -> fully concrete paths, no wildcards. Built with the same
## kallisto as scUTRquant's `kallisto bus` (env below) for index compatibility.
rule kallisto_index_alt:
    message: f"--- Building kallisto index for the ENDome ({ENDOME_NAME}) ---"
    input:
        fa = _step05_file("fa.gz")
    output:
        kdx = alt_kdx
    log: str(alt_dir / "kallisto_index" / f"{ENDOME_NAME}.log")
    benchmark: str(alt_dir / "kallisto_index" / f"{ENDOME_NAME}.benchmark.tsv")
    conda: "../envs/scutrquant-kallisto-bustools.yaml"
    shell:
        "kallisto index -i {output.kdx} {input.fa} 2>&1 | tee {log}"

## Convenience target: build the key outputs (txs SCE + tx->bin mapping).
rule utr_quant_alt:
    message: "--- scUTRquant (ALT) finished ---"
    input:
        sce = str(alt_dir / "data" / "sce" / _alt_target / f"{scutrquant_dataset_name}.txs.Rds"),
        tx_merge = str(alt_dir / "data" / "utrs" / _alt_target / "tx_merge.tsv")

## Debug: benchmark loading time
_log(f"\t+ 06-UTR_quantification_ALT imported in {time.perf_counter() - _start_time:.3f}s")
