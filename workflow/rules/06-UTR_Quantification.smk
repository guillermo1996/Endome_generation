################################################################################
## 06 - UTR Quantification
################################################################################
## Quantifies every ENDome on the single-cell dataset (config "sc_dataset") with
## scUTRquant
from snakemake.common.configfile import load_configfile

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
sc_dataset = config.get("sc_dataset")

### scUTRquant download path
scUTRquant_version = scUTRquant_settings["version"]
scUTRquant_tar_url = f"https://github.com/Mayrlab/scUTRquant/archive/refs/tags/{scUTRquant_version}.tar.gz"

### scUTRquant names its outputs after these keys of the template config
scUTRquant_template = load_configfile(sc_dataset["scUTRquant_config_template"])
scUTRquant_dataset_name = scUTRquant_template["dataset_name"]
scUTRquant_output_types = scUTRquant_template["output_type"]

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

        # ENDome patches to scUTRquant:
        ## 1) target-specific tmp dir for bustool_sort
        sed -i 's|"/bs-utrome-sort" + wcs.sample_id|"/bs-utrome-sort" + wcs.target + "-" + wcs.sample_id|' {output.snakefile}
        ## 2) mark output.sorted.bus as temp() to save disk space
        sed -i '/^rule bustools_sort:/,/^rule /s|"\\(data/kallisto/.*/output\\.sorted\\.bus\\)"|temp("\\1")|' {output.snakefile}
        ## Print if successful
        [ "$(grep -c 'wcs.target + "-" + wcs.sample_id' {output.snakefile})" -eq 2 ] || {{ echo "scUTRquant tmp_dir patch failed" >&2; exit 1; }}
        grep -q 'temp("data/kallisto/.*/output.sorted.bus")' {output.snakefile} || {{ echo "scUTRquant temp() patch failed" >&2; exit 1; }}
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
        # scUTRquant resolves "path" from its working directory (the step06
        # folder) and the file paths from "path" (the step05 folder)
        base = os.path.abspath(str(step05.path()).format(dataset=wildcards.dataset, group=wildcards.group))
        workdir = os.path.abspath(str(step06.path()).format(dataset=wildcards.dataset, group=wildcards.group))
        entry = {
            wildcards.endome_name: {
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
    message: "--- Writing the scUTRquant config for {wildcards.endome_name} ---"
    input:
        template = sc_dataset["scUTRquant_config_template"],
        targets = rules.generate_endome_target.output.targets,
        samples = sc_dataset["sample_sheet"],
    output:
        config_file = step06.path("Configs/{endome_name}.scUTRquant_config.yaml")
    run:
        # Paths relative to scUTRquant's working directory (the step06 folder)
        base = os.path.abspath(str(step06.path()).format(dataset=wildcards.dataset, group=wildcards.group))
        cfg = load_configfile(input.template)
        cfg["target"] = wildcards.endome_name
        cfg["targets_config"] = os.path.relpath(input.targets, base)
        cfg["sample_file"] = os.path.relpath(input.samples, base)
        for k in ("bx_whitelist", "cell_annots"):
            if cfg.get(k):
                cfg[k] = os.path.relpath(os.path.expanduser(cfg[k]), base)
        cfg["tmp_dir"] = os.path.abspath(os.path.expanduser(cfg["tmp_dir"]))
        with open(output.config_file, "w") as f:
            yaml.safe_dump(cfg, f, sort_keys=False, default_flow_style=False)

rule prepare_scUTRquant_envs:
    message: "--- Creating scUTRquant conda environments ---"
    input:
        # Any ENDome config works: the envs are the same for every target
        config_file = lambda w: ancient(endome_targets[0]),
        snakefile = rules.download_scUTRquant.output.snakefile
    output:
        flag = touch(f".snakemake/conda/scUTRquant-{scUTRquant_version}.envs.done")
    log: f".snakemake/conda/scUTRquant-{scUTRquant_version}.envs.log"
    threads: 1
    params:
        workdir = lambda w, input: os.path.dirname(os.path.dirname(os.path.abspath(input.config_file))),
        configfile = lambda w, input: os.path.abspath(input.config_file),
        conda_prefix = os.path.abspath(".snakemake/conda"),
    shell:
        """
        snakemake -s {input.snakefile} \
            --configfile {params.configfile} \
            --directory {params.workdir} \
            --nolock \
            --use-conda --conda-prefix {params.conda_prefix} \
            --conda-create-envs-only \
            --cores {threads} 2>&1 | tee {log}
        """

rule run_scUTRquant:
    message: "--- Running scUTRquant for {wildcards.endome_name} ---"
    input:
        config_file = rules.generate_scutrquant_config.output.config_file,
        snakefile = rules.download_scUTRquant.output.snakefile,
        envs = ancient(rules.prepare_scUTRquant_envs.output.flag),
    output:
        sce = expand(step06.path("data/sce/{endome_name}/{dataset_name}.{output_type}.Rds"),
            dataset_name=scUTRquant_dataset_name, output_type=scUTRquant_output_types, allow_missing=True),
    log: step06.logs("scUTRquant/{endome_name}.run_scUTRquant.log")
    threads: 36
    params:
        workdir = lambda w: str(step06.path()).format(dataset=w.dataset, group=w.group),
        configfile = lambda w, input: os.path.abspath(input.config_file),
        conda_prefix = os.path.abspath(".snakemake/conda"),
        # kallisto bus spends most of its time loading the index (single-threaded),
        # so it is faster to run several libraries at once with fewer threads
        # each than one library with all of them (benchmarks/scUTRquant_threads)
        kallisto_bus_threads = lambda w, threads: max(1, min(10, threads-1))
    shell:
        """
        snakemake -s {input.snakefile} \
            --configfile {params.configfile} \
            --config target={wildcards.endome_name} \
            --directory {params.workdir} \
            --nolock --rerun-incomplete \
            --use-conda --conda-prefix {params.conda_prefix} \
            --set-threads kallisto_bus={params.kallisto_bus_threads} bustools_sort=1 bustools_correct_sort=1 \
            --cores {threads} 2>&1 | tee {log}
        """

rule collect_run_info:
    message: "--- Collecting kallisto pseudoalignment rates for {wildcards.endome_name} ---"
    input:
        sce = rules.run_scUTRquant.output.sce,
    output:
        tsv = step06.path("pseudoalignment/{endome_name}.run_info.tsv")
    params:
        kallisto_dir = lambda w: str(step06.path("data/kallisto/{endome_name}")).format(dataset=w.dataset, group=w.group, endome_name=w.endome_name),
    run:
        import json, csv
        fields = ["sample_id", "n_processed", "n_pseudoaligned", "n_unique", "p_pseudoaligned", "p_unique", "kallisto_version"]
        with open(output.tsv, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=fields, delimiter="\t", extrasaction="ignore")
            writer.writeheader()
            for run_info in sorted(Path(params.kallisto_dir).glob("*/run_info.json")):
                with open(run_info) as j:
                    writer.writerow({"sample_id": run_info.parent.name, **json.load(j)})

rule pseudoalignment_rates:
    message: "--- Plotting scUTRquant pseudoalignment rates for {wildcards.dataset}.{wildcards.group}.{wildcards.merge_method} ---"
    input:
        run_info = expand(step06.path("pseudoalignment/{dataset}.{group}.{merge_method}.{orf_filter}.w{width}.{txEnd}.run_info.tsv"),
            orf_filter = config["orf_filter_preset"], width = config["trunc_width"], txEnd = config["trunc_site"], allow_missing = True)
    output:
        tsv = step06.path("pseudoalignment/{dataset}.{group}.{merge_method}.pseudoalignment.tsv"),
        p_pseudoaligned = step06.path("pseudoalignment/{dataset}.{group}.{merge_method}.pseudoalignment.png"),
        p_unique = step06.path("pseudoalignment/{dataset}.{group}.{merge_method}.unique.png"),
    log: step06.logs("pseudoalignment/{dataset}.{group}.{merge_method}.pseudoalignment.log")
    params:
        min_samples_boxplot = 3,
        dataset_name = scUTRquant_dataset_name,
        regions = ["ACG", "IPL"]    # split the plots by these strings of the library names (null = no split)
    conda: "../envs/r.yaml"
    script: "../scripts/S06a-Pseudoalignment_Rates.R"


register_test_data_link(rules.generate_scutrquant_config.output.config_file)
register_test_data_link(rules.collect_run_info.output.tsv)
register_test_data_link(rules.run_scUTRquant.output.sce)
