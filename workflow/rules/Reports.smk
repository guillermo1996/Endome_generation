################################################################################
## Reports: rendered reports of the pipeline outputs
################################################################################

## Variables
################################################################################
### Axes of the report targets
endome_axes_pseudobulk = table_axes("dataset", "group", "endome_name", "pseudobulk_method")
endome_axes_dtu = table_axes("dataset", "group", "endome_name", "pseudobulk_method", "saturn_method")

### Transcriptomes (prefix) with ORFs: one information content report each
endome_axes_info = table_axes("dataset", "group", "prefix", where = "has_orfs")

### Transcriptomes with a step 02 assembly (pre-truncated GTFs have none)
endome_axes_isopedia = table_axes("dataset", "group", "merge_method", where = "has_orfs")

### Bin information scoring settings.
scoring_settings = config.get("scoring_settings", {}).get("default")

### Input function: `template` filled in for every ENDome (orf_filter x width x txEnd) of the prefix wildcard
def prefix_endomes(template):
    def _input(wildcards):
        rows = endome_table.loc[endome_table["has_orfs"] & (endome_table["prefix"] == wildcards.prefix), ["orf_filter", "width", "txEnd"]].drop_duplicates()
        if rows.empty:
            raise ValueError(f"'{wildcards.prefix}' has no ENDome with ORFs in pipeline_runs (pre-truncated GTFs have no reports)")
        n = len(rows)
        return expand(template, zip, dataset = [wildcards.dataset] * n,
                      prefix = [wildcards.prefix] * n, **rows.to_dict("list"))
    return _input

## Information content report (step 05): one per transcriptome, comparing its ENDomes
################################################################################
rule endome_info_report:
    message: "--- Rendering the information content report - {wildcards.prefix} ---"
    input:
        duckdb = prefix_endomes(rules.endome_build_db.output.duckdb),
        done = prefix_endomes(rules.endome_compute_scores.output.done)
    output:
        html = step05.path("reports/{prefix}.info_report.html")
    log: step05.logs("reports/{prefix}.info_report.log")
    benchmark: step05.benchmark("reports/{prefix}.info_report.tsv")
    params:
        merge_distance = txendcutr_settings["merge_distance"],
        protein_metric = scoring_settings["protein_metric"],
        length_scaling = scoring_settings["length_scaling"],
        member_weights = scoring_settings["member_weights"],
        component_weights = scoring_settings["component_weights"]
    conda: "../envs/r.yaml"
    script: "../reports/R05a-Information_Content_Report.Rmd"

## UTR usage reports (step 07)
################################################################################
rule pseudobulk_report:
    message: "--- Rendering the pseudobulk report of {wildcards.endome_name} ---"
    input:
        cells = expand(rules.process_scUTRquant.output.cells, pseudobulk_method = config.get("pseudobulk_preset", ["default"]), allow_missing = True),
        pseudobulk = expand(rules.pseudobulk_scUTRquant.output.pseudobulk, pseudobulk_method = config.get("pseudobulk_preset", ["default"]), allow_missing = True),
        covariates = sc_dataset["sample_table"],
    output:
        html = step07.path("reports/{endome_name}.pseudobulk_report.html")
    log: step07.logs("reports/{endome_name}.pseudobulk_report.log")
    benchmark: step07.benchmark("reports/{endome_name}.pseudobulk_report.tsv")
    params:
        preset_settings = pseudobulk_settings,
        satuRn = {**satuRn_settings["default"], **sc_dataset},
    conda: "../envs/r.yaml"
    script: "../reports/R07a-Pseudobulk_Report.Rmd"

rule satuRn_report:
    message: "--- Rendering the satuRn report of {wildcards.endome_name} ---"
    input:
        results = expand(rules.satuRn_DTU.output.results_rds,
            pseudobulk_method = config.get("pseudobulk_preset", ["default"]),
            saturn_method = config.get("satuRn_preset", ["default"]), allow_missing = True),
        diagnostics = expand(rules.satuRn_DTU.output.diagnostics_rds,
            pseudobulk_method = config.get("pseudobulk_preset", ["default"]),
            saturn_method = config.get("satuRn_preset", ["default"]), allow_missing = True),
        pseudobulk = expand(rules.pseudobulk_scUTRquant.output.pseudobulk,
            pseudobulk_method = config.get("pseudobulk_preset", ["default"]), allow_missing = True),
        n_cells = expand(rules.pseudobulk_scUTRquant.output.n_cells,
            pseudobulk_method = config.get("pseudobulk_preset", ["default"]), allow_missing = True),
        sample_table = sc_dataset["sample_table"],
    output:
        html = step07.path("reports/{endome_name}.satuRn_report.html")
    log: step07.logs("reports/{endome_name}.satuRn_report.log")
    benchmark: step07.benchmark("reports/{endome_name}.satuRn_report.tsv")
    params:
        pseudobulk_settings = pseudobulk_settings,
        satuRn_settings = satuRn_settings,
        comparison = {key: sc_dataset.get(key) for key in ["stratify_by", "celltypes", "contrasts", "covariates", "exclude_samples"]},
        genes_of_interest = ["APP", "SNCA"],
        top_genes = 6,
    conda: "../envs/r.yaml"
    script: "../reports/R07b-satuRn_Report.Rmd"

## Isopedia validation report (step 02)
################################################################################
rule isopedia_report:
    message: "--- Rendering the Isopedia report - {wildcards.dataset}.{wildcards.group}.{wildcards.merge_method} ---"
    input:
        samples = rules.isopedia_summary.output.samples,
        transcripts = rules.isopedia_summary.output.transcripts,
        groups = rules.isopedia_summary.output.groups,
        gtf = rules.manual_filter.output.gtf
    output:
        html = step02.path("isopedia/reports/{dataset}.{group}.{merge_method}.isopedia_report.html")
    log: step02.logs("isopedia/{dataset}.{group}.{merge_method}.isopedia_report.log")
    benchmark: step02.benchmark("isopedia/{dataset}.{group}.{merge_method}.isopedia_report.tsv")
    params:
        isopedia_version = isopedia_version,
        isopedia_index = isopedia_index_dir,
        isopedia_settings = isopedia_settings
    conda: "../envs/r.yaml"
    script: "../reports/R02a-Isopedia_Report.Rmd"

## Targets
################################################################################
rule isopedia:
    input:
        transcripts = expand(rules.isopedia_summary.output.transcripts, zip, **endome_axes_isopedia),
        report = expand(rules.isopedia_report.output.html, zip, **endome_axes_isopedia)

rule info_report:
    input:
        html = expand(rules.endome_info_report.output.html, zip, **endome_axes_info)

rule pseudobulk:
    input:
        pseudobulk = expand(rules.pseudobulk_scUTRquant.output.pseudobulk, zip, **endome_axes_pseudobulk),
        report = expand(rules.pseudobulk_report.output.html, zip, **endome_axes_sce)

rule satuRn:
    input:
        results = expand(rules.satuRn_DTU.output.results, zip, **endome_axes_dtu),
        report = expand(rules.satuRn_report.output.html, zip, **endome_axes_sce)
