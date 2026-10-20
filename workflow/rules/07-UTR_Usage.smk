################################################################################
## 07 - UTR Usage
################################################################################
## Differential usage of the ENDome bins. Step by step:
##   1. process_scUTRquant    (S07a) -> annotated nuclei and their QC
##   2. pseudobulk_scUTRquant (S07b) -> processed SCE and pseudobulk (sample x cell type)
##   3. satuRn_DTU            (S07c) -> satuRn + stageR per stratum x cell type
##   4. pseudobulk_report / satuRn_report (R07a / R07b)
## The cell annotation and sample tables are dataset-specific and are prepared
## outside the pipeline

## Variables
################################################################################
### Input tables and what to compare (config "sc_dataset")
sc_dataset = config["sc_dataset"]

### Register the step: hash suffix, parameters.yaml rule, path helpers.
step07 = register_step(
    name="07-UTR_Usage",
    params={
        "pseudobulk_settings": build_tool_settings(config, "pseudobulk_settings"),
        "satuRn_settings": build_tool_settings(config, "satuRn_settings"),
    },
)

### Resolved settings (the step hash is derived from these)
satuRn_settings = config["satuRn_settings"]
pseudobulk_settings = config["pseudobulk_settings"]

## Rules
################################################################################
rule process_scUTRquant:
    message: "--- Annotating and filtering the nuclei of {wildcards.endome_name} ---"
    input:
        sce = step06.path(f"data/sce/{{endome_name}}/{scUTRquant_dataset_name}.txs.Rds"),
        cell_annotations = sc_dataset["cell_annotations"],
    output:
        cells = step07.path("cells/{endome_name}.{pseudobulk_method}.cells.tsv"),
        qc = step07.path("cells/{endome_name}.{pseudobulk_method}.cells_qc.tsv"),
    log: step07.logs("cells/{endome_name}.{pseudobulk_method}.cells.log")
    benchmark: step07.benchmark("cells/{endome_name}.{pseudobulk_method}.cells.tsv")
    threads: 1
    params:
        exclude_samples = sc_dataset.get("exclude_samples"),
        effective_txs_mode = lambda wc: pseudobulk_settings[wc.pseudobulk_method].get("effective_txs_mode", "fixed"),
        min_effective_txs = lambda wc: pseudobulk_settings[wc.pseudobulk_method].get("min_effective_txs"),
        mad_below_median = lambda wc: pseudobulk_settings[wc.pseudobulk_method].get("mad_below_median"),
    conda: "../envs/r.yaml"
    script: "../scripts/S07a-Process_scUTRquant.R"

rule pseudobulk_scUTRquant:
    message: "--- Building the processed SCE and pseudobulk of {wildcards.endome_name} ---"
    input:
        sce = rules.process_scUTRquant.input.sce,
        cells = rules.process_scUTRquant.output.cells,
        gtf = step05.path("txendcutr/{endome_name}.txendcutr.gtf"),
        merge_tsv = step05.path("txendcutr/{endome_name}.txendcutr.merge.tsv"),
    output:
        sce = step07.path("sce/{endome_name}.{pseudobulk_method}.processed.txs.rds"),
        pseudobulk = step07.path("pseudobulk/{endome_name}.{pseudobulk_method}.pseudobulk.rds"),
        n_cells = step07.path("pseudobulk/{endome_name}.{pseudobulk_method}.pseudobulk_n_cells.tsv"),
    log: step07.logs("pseudobulk/{endome_name}.{pseudobulk_method}.pseudobulk.log")
    benchmark: step07.benchmark("pseudobulk/{endome_name}.{pseudobulk_method}.pseudobulk.tsv")
    threads: 4
    conda: "../envs/r.yaml"
    script: "../scripts/S07b-Pseudobulk.R"

rule satuRn_DTU:
    message: "--- Testing differential bin usage of {wildcards.endome_name} with satuRn ---"
    input:
        pseudobulk = rules.pseudobulk_scUTRquant.output.pseudobulk,
        n_cells = rules.pseudobulk_scUTRquant.output.n_cells,
        sample_table = sc_dataset["sample_table"],
    output:
        results = step07.path("DTU/{endome_name}.{pseudobulk_method}.{saturn_method}.satuRn.condition.tsv"),
        diagnostics = step07.path("DTU/{endome_name}.{pseudobulk_method}.{saturn_method}.satuRn.condition.diagnostics.tsv"),
        results_rds = step07.path("DTU/{endome_name}.{pseudobulk_method}.{saturn_method}.satuRn.condition.rds"),
        diagnostics_rds = step07.path("DTU/{endome_name}.{pseudobulk_method}.{saturn_method}.satuRn.condition.diagnostics.rds"),
        diagplots = directory(step07.path("DTU/{endome_name}.{pseudobulk_method}.{saturn_method}.satuRn.condition.diagplots")),
    log: step07.logs("DTU/{endome_name}.{pseudobulk_method}.{saturn_method}.satuRn.condition.log")
    benchmark: step07.benchmark("DTU/{endome_name}.{pseudobulk_method}.{saturn_method}.satuRn.condition.tsv")
    threads: 16
    params:
        exclude_samples = sc_dataset.get("exclude_samples"),
        stratify_by = sc_dataset.get("stratify_by"),
        celltypes = sc_dataset.get("celltypes"),
        contrasts = sc_dataset["contrasts"],
        covariates = sc_dataset.get("covariates"),
        min_cells = lambda w: satuRn_settings[w.saturn_method]["min_cells"],
        min_genes = lambda w: satuRn_settings[w.saturn_method]["min_genes"],
        min_samples_per_group = lambda w: satuRn_settings[w.saturn_method]["min_samples_per_group"],
        dmfilter = lambda w: satuRn_settings[w.saturn_method]["dmfilter"],
        p_value = lambda w: satuRn_settings[w.saturn_method]["p_value"],
        alpha = lambda w: satuRn_settings[w.saturn_method]["alpha"]
    conda: "../envs/r.yaml"
    script: "../scripts/S07c-satuRn_DTU.R"

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
    script: "../scripts/R07a-Pseudobulk_Report.Rmd"

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
    script: "../scripts/R07b-satuRn_Report.Rmd"


register_test_data_link(rules.process_scUTRquant.output.cells)
register_test_data_link(rules.process_scUTRquant.output.qc)
register_test_data_link(rules.pseudobulk_scUTRquant.output.pseudobulk)
register_test_data_link(rules.pseudobulk_scUTRquant.output.n_cells)
register_test_data_link(rules.satuRn_DTU.output.results)
register_test_data_link(rules.satuRn_DTU.output.diagnostics)
register_test_data_link(rules.satuRn_DTU.output.results_rds)
register_test_data_link(rules.satuRn_DTU.output.diagnostics_rds)