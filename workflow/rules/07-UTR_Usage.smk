################################################################################
## 06b - UTR usage analysis
################################################################################

## Debug: benchmark loading time
_start_time = time.perf_counter()

## Variables
################################################################################
### Register the step: hash suffix, parameters.yaml rule, path helpers.
step07 = register_step(
    name="07-UTR_Usage",
    params={
        "pseudobulk_settings": build_tool_settings(config, "pseudobulk_settings"),
        "satuRn_settings": build_tool_settings(config, "satuRn_settings"),
    },
)

### Resolved settings (the step hash is derived from these)
# pseudobulk_settings = resolve_preset(config, "pseudobulk_settings", "pseudobulk_preset")
# satuRn_settings = resolve_preset(config, "satuRn_settings", "satuRn_preset")

satuRn_settings = config["satuRn_settings"]
pseudobulk_settings = config["pseudobulk_settings"]

annotations = satuRn_settings["default"]["annotations"]
sample_metadata = satuRn_settings["default"]["sample_metadata"]
sample_covariates = satuRn_settings["default"]["sample_covariates"]
exclude_samples = satuRn_settings["default"]["exclude_samples"]

## Functions
################################################################################

## Rules
################################################################################
rule process_scUTRquant:
    message: "--- Processing scUTRquant output ---"
    input:
        sce = step06.path(f"data/sce/{{endome_name}}/{scUTRquant_dataset_name}.txs.Rds"),
        annotations = annotations,
        sample_metadata = sample_metadata,
    output:
        cells = step07.path("cells/{endome_name}.{pseudobulk_method}.cells.tsv"),
        qc = step07.path("cells/{endome_name}.{pseudobulk_method}.cells_qc.tsv"),
    log: step07.logs("cells/{endome_name}.{pseudobulk_method}.cells.log")
    benchmark: step07.benchmark("cells/{endome_name}.{pseudobulk_method}.cells.tsv")
    threads: 1
    params:
        exclude_samples = exclude_samples,
        effective_txs_mode = lambda wc: pseudobulk_settings[wc.pseudobulk_method].get("effective_txs_mode", "fixed"),
        min_effective_txs = lambda wc: pseudobulk_settings[wc.pseudobulk_method].get("min_effective_txs", 300),
        mad_below_median = lambda wc: pseudobulk_settings[wc.pseudobulk_method].get("mad_below_median", 3),
        cell_type_level = lambda wc: pseudobulk_settings[wc.pseudobulk_method]["cell_type_level"]
    conda: "../envs/r.yaml"
    script: "../scripts/07a-Process_scUTRquant.R"

rule pseudobulk_scUTRquant:
    message: "--- Building the processed SCE and pseudobulk ---"
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
    params:
        cell_type_level = lambda wc: pseudobulk_settings[wc.pseudobulk_method]["cell_type_level"]
    conda: "../envs/r.yaml"
    script: "../scripts/07b-Pseudobulk.R"

rule satuRn_DTU:
    message: "--- Testing differential bin usage between conditions with satuRn ---"
    input:
        pseudobulk = rules.pseudobulk_scUTRquant.output.pseudobulk,
        n_cells = rules.pseudobulk_scUTRquant.output.n_cells,
        sample_covariates = sample_covariates,
    output:
        results = step07.path("DTU/{endome_name}.{pseudobulk_method}.{saturn_method}.satuRn.condition.tsv.gz"),
        diagnostics = step07.path("DTU/{endome_name}.{pseudobulk_method}.{saturn_method}.satuRn.condition.diagnostics.tsv"),
        diagplots = directory(step07.path("DTU/{endome_name}.{pseudobulk_method}.{saturn_method}.satuRn.condition.diagplots")),
    log: step07.logs("DTU/{endome_name}.{pseudobulk_method}.{saturn_method}.satuRn.condition.log")
    benchmark: step07.benchmark("DTU/{endome_name}.{pseudobulk_method}.{saturn_method}.satuRn.condition.tsv")
    threads: 16
    params:
        exclude_samples = lambda w: satuRn_settings[w.saturn_method]["exclude_samples"],
        regions = lambda w: satuRn_settings[w.saturn_method]["regions"],
        celltypes = lambda w: satuRn_settings[w.saturn_method]["celltypes"],
        contrasts = lambda w: satuRn_settings[w.saturn_method]["contrasts"],
        covariates = lambda w: satuRn_settings[w.saturn_method]["covariates"],
        min_cells = lambda w: satuRn_settings[w.saturn_method]["min_cells"],
        min_genes = lambda w: satuRn_settings[w.saturn_method]["min_genes"],
        min_samples_per_group = lambda w: satuRn_settings[w.saturn_method]["min_samples_per_group"],
        dmfilter = lambda w: satuRn_settings[w.saturn_method]["dmfilter"],
        p_value = lambda w: satuRn_settings[w.saturn_method]["p_value"],
        alpha = lambda w: satuRn_settings[w.saturn_method]["alpha"]
    conda: "../envs/r.yaml"
    script: "../scripts/07c-satuRn_DTU.R"

rule pseudobulk_report:
    message: "--- Rendering the pseudobulk report ---"
    input:
        cells = expand(rules.process_scUTRquant.output.cells, pseudobulk_method = config.get("pseudobulk_preset", ["default"]), allow_missing = True),
        pseudobulk = expand(rules.pseudobulk_scUTRquant.output.pseudobulk, pseudobulk_method = config.get("pseudobulk_preset", ["default"]), allow_missing = True),
        covariates = sample_covariates,
    output:
        html = step07.path("reports/{endome_name}.pseudobulk_report.html")
    log: step07.logs("reports/{endome_name}.pseudobulk_report.log")
    benchmark: step07.benchmark("reports/{endome_name}.pseudobulk_report.tsv")
    params:
        preset_settings = pseudobulk_settings,
        satuRn = satuRn_settings["default"],
    conda: "../envs/r.yaml"
    script: "../scripts/07d-Pseudobulk_Report.Rmd"

rule satuRn_report:
    message: "--- Rendering the satuRn report ---"
    input:
        results = expand(rules.satuRn_DTU.output.results, 
            pseudobulk_method = config.get("pseudobulk_preset", ["default"]),
            saturn_method = config.get("satuRn_preset", ["default"]), allow_missing = True),
        diagnostics = expand(rules.satuRn_DTU.output.diagnostics, 
            pseudobulk_method = config.get("pseudobulk_preset", ["default"]),
            saturn_method = config.get("satuRn_preset", ["default"]), allow_missing = True),
        pseudobulk = expand(rules.pseudobulk_scUTRquant.output.pseudobulk, 
            pseudobulk_method = config.get("pseudobulk_preset", ["default"]), allow_missing = True),
        sample_covariates = sample_covariates,
    output:
        html = step07.path("reports/{endome_name}.satuRn_report.html")
    log: step07.logs("reports/{endome_name}.satuRn_report.log")
    benchmark: step07.benchmark("reports/{endome_name}.satuRn_report.tsv")
    params:
        pseudobulk_settings = pseudobulk_settings,
        satuRn_settings = satuRn_settings,
        genes_of_interest = ["APP", "SNCA"],
        top_genes = 6,
    conda: "../envs/r.yaml"
    script: "../scripts/07e-satuRn_Report.Rmd"


register_test_data_link(rules.process_scUTRquant.output.cells)
register_test_data_link(rules.process_scUTRquant.output.qc)
register_test_data_link(rules.pseudobulk_scUTRquant.output.pseudobulk)
register_test_data_link(rules.pseudobulk_scUTRquant.output.n_cells)

# rule satuRn_DTU:
#     message: "--- Testing differential bin usage between conditions with satuRn ---"
#     input:
#         pseudobulk = rules.pseudobulk_scUTRquant.output.pseudobulk,
#         n_cells = rules.pseudobulk_scUTRquant.output.n_cells,
#         covariates = satuRn_settings["covariates_table"],
#     output:
#         results = step06.path("satuRn/{endome_name}.satuRn.condition.tsv.gz"),
#         diagnostics = step06.path("satuRn/{endome_name}.satuRn.condition.diagnostics.tsv"),
#         diagplots = directory(step06.path("satuRn/{endome_name}.satuRn.condition.diagplots")),
#     log: step06.logs("satuRn/{endome_name}.satuRn.condition.log")
#     benchmark: step06.benchmark("satuRn/{endome_name}.satuRn.condition.tsv")
#     threads: 16
#     params:
#         regions = satuRn_settings["regions"],
#         celltypes = satuRn_settings["celltypes"],
#         exclude_samples = satuRn_settings["exclude_samples"],
#         min_genes = satuRn_settings["min_genes"],
#         min_cells = satuRn_settings["min_cells"],
#         min_samples_per_group = satuRn_settings["min_samples_per_group"],
#         covariates_sample_col = satuRn_settings["covariates_sample_col"],
#         group_col = satuRn_settings["group_col"],
#         contrasts = satuRn_settings["contrasts"],
#         fit_mode = satuRn_settings["fit_mode"],
#         covariates = satuRn_settings["covariates"],
#         scale_covariates = satuRn_settings["scale_covariates"],
#         dmfilter = satuRn_settings["dmfilter"],
#         p_fallback = satuRn_settings["p_fallback"],
#         alpha = satuRn_settings["alpha"],
#     conda: "../envs/r.yaml"
#     script: "../scripts/06e-satuRn_DTU.R"

# register_test_data_link(rules.process_scUTRquant.output.cells)
# register_test_data_link(rules.pseudobulk_scUTRquant.output.sce)
# register_test_data_link(rules.pseudobulk_scUTRquant.output.pseudobulk)
# register_test_data_link(rules.pseudobulk_scUTRquant.output.n_cells)