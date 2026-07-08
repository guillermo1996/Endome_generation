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
        "scUTRquant_presets": build_tool_settings(config, "scUTRquant_settings", "scUTRquant_preset"),
    },
)

### Resolved settings (the step hash is derived from these)
scUTRquant_settings = resolve_preset(config, "scUTRquant_settings", "scUTRquant_preset")

### scUTRquant Download path
scUTRquant_version = scUTRquant_settings["version"]
scUTRquant_tar_url = f"https://github.com/Mayrlab/scUTRquant/archive/refs/tags/{scUTRquant_version}.tar.gz"

## Functions
################################################################################
def endome_target_files(wildcards):
    return expand(step06.path("scUTRquant/ENDomes/{dataset}.{group}.{orf_filter}.w{width}.{txEnd}.ENDome_targets.yml"),
                  dataset=wildcards.dataset, group=wildcards.group,
                  orf_filter=pc_filter, width=trunc_width, txEnd=trunc_site)

## Rules
################################################################################
rule download_scUTRquant:
    message: "--- Downloading and Extracting scUTRquant ---"
    output:
        scUTRquant_dir = directory(f"tools/scUTRquant-{orfannotate_version}"),
        snakefile = f"tools/scUTRquant-{orfannotate_version}/Snakefile"
    params:
        url = scUTRquant_tar_url
    shell:
        """
        mkdir -p {output.scUTRquant_dir}
        curl -L {params.url} | tar -xz -C {output.scUTRquant_dir} --strip-components=1
        """

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
#     log: step06.log("scUTRquant/run_scUTRquant.log")
#     threads: 30
#     shell:
#         "snakemake -s tools/scUTRquant/Snakefile "
#         "--configfile {input.config_file} "
#         "--use-conda --nolock --cores {threads} "
#         "2>&1 | tee {log}"

rule generate_endome_target:
    message: "--- Writing scUTRquant target entry for {wildcards.prefix}.{wildcards.orf_filter}.w{wildcards.width}.{wildcards.txEnd} ---"
    input:
        gtf = step05.path("gtxcutr/{prefix}.{orf_filter}.gtxcutr.w{width}.{txEnd}.gtf"),
        merge = step05.path("gtxcutr/{prefix}.{orf_filter}.gtxcutr.w{width}.{txEnd}.merge.tsv"),
        kdx = step05.path("kallisto_index/{prefix}.{orf_filter}.w{width}.{txEnd}.kdx"),
    output:
        targets = step06.path("scUTRquant/ENDomes/{prefix}.{orf_filter}.w{width}.{txEnd}.ENDome_targets.yml")
    run:
        name = f"{wildcards.prefix}.{wildcards.orf_filter}.w{wildcards.width}.{wildcards.txEnd}"
        base = os.path.abspath(str(step05.path()).format(dataset=wildcards.dataset, group=wildcards.group))
        
        entry = {
            name: {
                "path": base + "/",
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

rule merge_endome_targets:
    message: "--- Merging scUTRquant target catalog for {wildcards.dataset}.{wildcards.group} ---"
    input: endome_target_files
    output:
        catalog = step06.path("scUTRquant/ENDome_targets.yml")
    run:
        catalog = {}
        for t in input:
            with open(t) as f:
                catalog.update(yaml.safe_load(f))
        with open(output.catalog, "w") as f:
            yaml.safe_dump(catalog, f, sort_keys=False, default_flow_style=False)

rule generate_scutrquant_config:
    message: "--- Writing resolved scUTRquant config for {wildcards.dataset}.{wildcards.group} ---"
    input:
        template = config["scUTRquant_config_template"],
        targets = step06.path("scUTRquant/ENDome_targets.yml"),
    output:
        config_file = step06.path("scUTRquant/scUTRquant_config.yaml")
    run:
        cfg = load_configfile(input.template)
        with open(input.targets) as f:
            catalog = yaml.safe_load(f)
        cfg["target"] = list(catalog.keys())
        cfg["targets_config"] = os.path.abspath(input.targets)
        for k in ("sample_file", "bx_whitelist"):
            if cfg.get(k):
                cfg[k] = os.path.abspath(os.path.expanduser(cfg[k]))
        cfg["tmp_dir"] = os.path.abspath(os.path.expanduser(cfg["tmp_dir"]))   # also fixes "~/tmp"
        with open(output.config_file, "w") as f:
            yaml.safe_dump(cfg, f, sort_keys=False, default_flow_style=False)

rule run_scUTRquant:
    message: "--- Running scUTRquant (nested) for {wildcards.dataset}.{wildcards.group} ---"
    input:
        config_file = step06.path("scUTRquant/scUTRquant_config.yaml"),
        snakefile = rules.download_scUTRquant.output.snakefile,
    output:
        flag = step06.path("scUTRquant/scUTRquant.done"),
    log: step06.log("scUTRquant/run_scUTRquant.log")
    threads: 32
    params:
        workdir = lambda w, output: os.path.dirname(output.flag),
        configfile = lambda w, input: os.path.abspath(input.config_file),
        conda_prefix = os.path.abspath(".snakemake/conda"),   # reuse outer envs
    shell:
        """
        snakemake -s {input.snakefile} \
            --configfile {params.configfile} \
            --directory {params.workdir} \
            --use-conda --conda-prefix {params.conda_prefix} \
            --cores {threads} 2>&1 | tee {log}
        touch {output.flag}
        """