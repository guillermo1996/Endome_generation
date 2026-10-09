################################################################################
## 06 - UTR Quantification
################################################################################
## This step prepares everything scUTRquant needs to quantify the ENDome, then
## scUTRquant is run STANDALONE against the generated config, e.g.:
##   snakemake -s tools/scUTRquant/Snakefile \
##             --configfile <step06>/scUTRquant/<endome>.scUTRquant_config.yaml
##
## It produces, per ENDome ({prefix}.{orf_filter}.w{width}.{txEnd}):
##   1. kallisto_index            -> the kallisto index (.kdx)
##   2. generate_endome_target    -> a scUTRquant "targets" YAML entry
##   3. generate_scutrquant_config-> the user's template config + injected
##                                   `target` and `targets_config`
from snakemake.common.configfile import load_configfile

## Debug: benchmark loading time
_start_time = time.perf_counter()

## Variables
################################################################################
### Register the step: hash suffix, parameters.yaml rule, path helpers.
step06 = register_step(
    name="06-UTR_Quantification",
    params={
        "scUTRquant_settings": build_tool_settings(config, "scUTRquant_settings", "scUTRquant_preset"),
    },
)

### Resolved settings (the step hash is derived from these)
scUTRquant_settings = resolve_preset(config, "scUTRquant_settings", "scUTRquant_preset")

### scUTRquant Download path
scUTRquant_version = scUTRquant_settings["version"]
scUTRquant_tar_url = f"https://github.com/Mayrlab/scUTRquant/archive/refs/tags/{scUTRquant_version}.tar.gz"

### scUTRquant template config and samples (used to name the run_scUTRquant outputs)
scUTRquant_template = load_configfile(config["scUTRquant_config_template"])
scUTRquant_dataset_name = scUTRquant_template["dataset_name"]
scUTRquant_output_types = scUTRquant_template["output_type"]

## Functions
################################################################################
def endome_target_files(wildcards):
    return expand(step06.path("ENDomes/{dataset}.{group}.{merge_method}.{orf_filter}.w{width}.{txEnd}.ENDome_targets.yml"),
        dataset = wildcards.dataset, group = wildcards.group, merge_method = merge_method,
        orf_filter = pc_filter, width = trunc_width, txEnd = trunc_site)

## Rules
################################################################################
rule download_scUTRquant:
    message: "--- Downloading and Extracting scUTRquant ---"
    output:
        scUTRquant_dir = directory(f"tools/scUTRquant-{scUTRquant_version}"),
        snakefile = f"tools/scUTRquant-{scUTRquant_version}/Snakefile"
    params:
        url = scUTRquant_tar_url
    shell:
        """
        mkdir -p {output.scUTRquant_dir}
        curl -L {params.url} | tar -xz -C {output.scUTRquant_dir} --strip-components=1
        """

rule generate_endome_target:
    message: "--- Writing scUTRquant target entry for {wildcards.endome_name} ---"
    input:
        gtf = step05.path("txendcutr/{endome_name}.txendcutr.gtf"),
        merge = step05.path("txendcutr/{endome_name}.txendcutr.merge.tsv"),
        kdx = step05.path("kallisto_index/{endome_name}.kdx"),
    output:
        targets = step06.path("Targets/{endome_name}.ENDome_targets.yml")
    run:
        name = f"{wildcards.endome_name}"
        base = os.path.abspath(str(step05.path()).format(dataset=wildcards.dataset, group=wildcards.group))
        workdir = os.path.abspath(str(step06.path()).format(dataset=wildcards.dataset, group=wildcards.group))

        entry = {
            name: {
                "path": os.path.relpath(base, workdir) + "/",
                "genome": scUTRquant_settings["genome"],
                "gtf": os.path.relpath(input.gtf, base),
                "kdx": os.path.relpath(input.kdx, base),
                "merge_tsv": os.path.relpath(input.merge, base),
                "tx_annots": None,
                "tx_annots_csv": None,
                "gene_annots": None,
                "gene_annots_csv": None,
            }
        }
        with open(output.targets, "w") as f:
            yaml.safe_dump(entry, f, sort_keys=False, default_flow_style=False)
            
rule generate_scutrquant_config:
    message: "--- Writing resolved scUTRquant config ---"
    input:
        template = config["scUTRquant_config_template"],
        targets = rules.generate_endome_target.output.targets, 
        samples = config["input_scUTRquant_samples"],
    output:
        config_file = step06.path("Configs/{endome_name}.scUTRquant_config.yaml")
    run:
        # Relative paths are resolved by scUTRquant from its working directory (run_scUTRquant's --directory)
        base = os.path.abspath(str(step06.path()).format(dataset=wildcards.dataset, group=wildcards.group))
        cfg = load_configfile(input.template)
        with open(input.targets) as f:
            catalog = yaml.safe_load(f)
        cfg["target"] = wildcards.endome_name
        cfg["targets_config"] = os.path.relpath(input.targets, base)
        cfg["sample_file"] = os.path.relpath(input.samples, base)
        for k in ("bx_whitelist", "cell_annots"):
            if cfg.get(k):
                cfg[k] = os.path.relpath(os.path.expanduser(cfg[k]), base)
        cfg["tmp_dir"] = os.path.abspath(os.path.expanduser(cfg["tmp_dir"]))
        with open(output.config_file, "w") as f:
            yaml.safe_dump(cfg, f, sort_keys=False, default_flow_style=False)

rule run_scUTRquant:
    input:
        config_file = rules.generate_scutrquant_config.output.config_file,
        snakefile = rules.download_scUTRquant.output.snakefile,
    output:
        sce = expand(step06.path("data/sce/{endome_name}/{dataset_name}.{output_type}.Rds"),
            dataset_name=scUTRquant_dataset_name, output_type=scUTRquant_output_types, allow_missing=True), 
    log: step06.logs("scUTRquant/{endome_name}.run_scUTRquant.log")
    threads: 16
    params:
        workdir = lambda w: str(step06.path()).format(dataset=w.dataset, group=w.group),
        configfile = lambda w, input: os.path.abspath(input.config_file),
        conda_prefix = os.path.abspath(".snakemake/conda"),   # reuse outer envs
    shell:
        """
        snakemake -s {input.snakefile} \
            --configfile {params.configfile} \
            --config target={wildcards.endome_name} \
            --directory {params.workdir} \
            --nolock \
            --use-conda --conda-prefix {params.conda_prefix} \
            --cores {threads} 2>&1 | tee {log}
        """ 

rule pseudoalignment_rates:
    message: "--- Plotting scUTRquant pseudoalignment rates for {wildcards.prefix} ---"
    input:
        config_file = rules.generate_scutrquant_config.output.config_file,
        scUTRquant = rules.run_scUTRquant.output,   # run_info.json files live in scUTRquant's data/kallisto
    output:
        tsv = step06.path("pseudoalignment/{prefix}.pseudoalignment.tsv"),
        png = step06.path("pseudoalignment/{prefix}.pseudoalignment.png"),
    log: step06.logs("pseudoalignment/{prefix}.pseudoalignment.log")
    params:
        min_samples_boxplot = 3,
    conda: "../envs/r.yaml"
    script: "../scripts/06b-Pseudoalignment_Rates.R"

register_test_data_link(rules.run_scUTRquant.output.sce)

# rule generate_endome_target:
#     message: "--- Writing scUTRquant target entry for {wildcards.prefix}.{wildcards.orf_filter}.w{wildcards.width}.{wildcards.txEnd} ---"
#     input:
#         gtf = step05.path("gtxcutr/{prefix}.{orf_filter}.gtxcutr.w{width}.{txEnd}.gtf"),
#         merge = step05.path("gtxcutr/{prefix}.{orf_filter}.gtxcutr.w{width}.{txEnd}.merge.tsv"),
#         kdx = step05.path("kallisto_index/{prefix}.{orf_filter}.w{width}.{txEnd}.kdx"),
#     output:
#         targets = step06.path("scUTRquant/ENDomes/{prefix}.{orf_filter}.w{width}.{txEnd}.ENDome_targets.yml")
#     run:
#         name = f"{wildcards.prefix}.{wildcards.orf_filter}.w{wildcards.width}.{wildcards.txEnd}"
#         base = str(step05.path()).format(dataset=wildcards.dataset, group=wildcards.group)
#         entry = {
#             name: {
#                 "path": base + "/",
#                 "genome": gtxcutr_settings["genome"],
#                 "gtf": os.path.relpath(input.gtf, base),
#                 "kdx": os.path.relpath(input.kdx, base),
#                 "merge_tsv": os.path.relpath(input.merge, base),
#                 "tx_annots": None,
#                 "tx_annots_csv": None,
#                 "gene_annots": None,
#                 "gene_annots_csv": None,
#             }
#         }
#         with open(output.targets, "w") as f:
#             yaml.safe_dump(entry, f, sort_keys=False, default_flow_style=False)
# rule generate_scutrquant_config:
#     message: "--- Writing resolved scUTRquant config for {wildcards.dataset}.{wildcards.group} ---"
#     input:
#         template = config["scUTRquant_config_template"],
#         targets = step06.path("scUTRquant/ENDome_targets.yml"),
#     output:
#         config_file = step06.path("scUTRquant/scUTRquant_config.yaml")
#     run:
#         cfg = load_configfile(input.template)
#         with open(input.targets) as f:
#             catalog = yaml.safe_load(f)
#         cfg["target"] = list(catalog.keys())
#         cfg["targets_config"] = input.targets
#         with open(output.config_file, "w") as f:
#             yaml.safe_dump(cfg, f, sort_keys=False, default_flow_style=False)
# rule run_scUTRquant:
#     message: "--- Running scUTRquant (nested) for {wildcards.dataset}.{wildcards.group} ---"
#     input:
#         config_file = step06.path("scUTRquant/scUTRquant_config.yaml")
#     output:
#         flag = touch(step06.path("scUTRquant/scUTRquant.done"))
#     log: step06.logs("scUTRquant/run_scUTRquant.log")
#     threads: 30
#     shell:
#         "snakemake -s tools/scUTRquant/Snakefile "
#         "--configfile {input.config_file} "
#         "--use-conda --nolock --cores {threads} "
#         "2>&1 | tee {log}"

# rule generate_endome_target2:
#     message: "--- Writing scUTRquant target entry for {wildcards.prefix}.{wildcards.orf_filter}.w{wildcards.width}.{wildcards.txEnd} ---"
#     input:
#         gtf = step05.path("gtxcutr/{prefix}.{orf_filter}.gtxcutr.w{width}.{txEnd}.gtf"),
#         merge = step05.path("gtxcutr/{prefix}.{orf_filter}.gtxcutr.w{width}.{txEnd}.merge.tsv"),
#         kdx = step05.path("kallisto_index/{prefix}.{orf_filter}.w{width}.{txEnd}.kdx"),
#     output:
#         targets = step06.path("scUTRquant/ENDomes/{prefix}.{orf_filter}.w{width}.{txEnd}.ENDome_targets.yml")
#     run:
#         name = f"{wildcards.prefix}.{wildcards.orf_filter}.w{wildcards.width}.{wildcards.txEnd}"
#         base = os.path.abspath(str(step05.path()).format(dataset=wildcards.dataset, group=wildcards.group))
        
#         entry = {
#             name: {
#                 "path": base + "/",
#                 "genome": scUTRquant_settings["genome"],
#                 "gtf": os.path.relpath(input.gtf, base),
#                 "kdx": os.path.relpath(input.kdx, base),
#                 "merge_tsv": os.path.relpath(input.merge, base),
#                 "tx_annots": None,
#                 "tx_annots_csv": None,
#                 "gene_annots": None,
#                 "gene_annots_csv": None,
#             }
#         }
#         with open(output.targets, "w") as f:
#             yaml.safe_dump(entry, f, sort_keys=False, default_flow_style=False)





# This approach works, but it generates a single config file, which implies that
# any change in a specific endome can turn into a global change


# rule merge_endome_targets:
#     message: "--- Merging scUTRquant target catalog for {wildcards.dataset}.{wildcards.group} ---"
#     input: endome_target_files
#     output:
#         catalog = step06.path("scUTRquant/ENDome_targets.yml")
#     run:
#         catalog = {}
#         for t in input:
#             with open(t) as f:
#                 catalog.update(yaml.safe_load(f))
#         with open(output.catalog, "w") as f:
#             yaml.safe_dump(catalog, f, sort_keys=False, default_flow_style=False)
# rule generate_scutrquant_config:
#     message: "--- Writing resolved scUTRquant config for {wildcards.dataset}.{wildcards.group} ---"
#     input:
#         template = config["scUTRquant_config_template"],
#         targets = rules.merge_endome_targets.output.catalog,
#         samples = config["input_scUTRquant_samples"],
#     output:
#         config_file = step06.path("scUTRquant/scUTRquant_config.yaml")
#     run:
#         cfg = load_configfile(input.template)
#         with open(input.targets) as f:
#             catalog = yaml.safe_load(f)
#         cfg["target"] = list(catalog.keys())
#         cfg["targets_config"] = os.path.abspath(input.targets)
#         cfg["sample_file"] = os.path.abspath(input.samples)
#         for k in ("bx_whitelist", "cell_annots"):
#             if cfg.get(k):
#                 cfg[k] = os.path.abspath(os.path.expanduser(cfg[k]))
#         cfg["tmp_dir"] = os.path.abspath(os.path.expanduser(cfg["tmp_dir"]))
#         with open(output.config_file, "w") as f:
#             yaml.safe_dump(cfg, f, sort_keys=False, default_flow_style=False)

# rule run_scUTRquant:
#     message: "--- Running scUTRquant (nested) for {wildcards.dataset}.{wildcards.group} ---"
#     input:
#         config_file = rules.generate_scutrquant_config.output.config_file,
#         snakefile = rules.download_scUTRquant.output.snakefile,
#     output:
#         flag = step06.path("scUTRquant/scUTRquant.done"),
#     log: step06.logs("scUTRquant/run_scUTRquant.log")
#     threads: 64
#     params:
#         workdir = lambda w, output: os.path.dirname(output.flag),
#         configfile = lambda w, input: os.path.abspath(input.config_file),
#         conda_prefix = os.path.abspath(".snakemake/conda"),   # reuse outer envs
#     shell:
#         """
#         snakemake -s {input.snakefile} \
#             --configfile {params.configfile} \
#             --directory {params.workdir} \
#             --use-conda --conda-prefix {params.conda_prefix} \
#             --cores {threads} 2>&1 | tee {log}
#         touch {output.flag}
#         """