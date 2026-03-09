################################################################################
## 06 - Evaluation
################################################################################

## Variables
################################################################################
step06_name = "06-Evaluation"

### Configurations
evaluation_samples = pd.read_csv(config['evaluation_sample_file'], index_col='sample_id')
kb_settings = config["kb_settings"][config["kb_profile"]]

step06_params = {
    "kb_settings": kb_settings,
    **({"toy_data": True} if config.get("use_toy_data", False) else {})
}

### Compute hash and set global parameters
global_params.update({step06_name: step06_params})
step06_hash = compute_hash(global_params)
global_params.update({f"{step06_name}_{step06_hash}": global_params.pop(step06_name)})

### Paths
evaluation_path = lambda x: Path(results_path) / f"{step06_name}-{step06_hash}" / x
evaluation_log_path = lambda x: Path(results_path) / f"{step06_name}-{step06_hash}" / log_path / x
evaluation_benchmark_path = lambda x: Path(results_path) / f"{step06_name}-{step06_hash}" / benchmark_path / x

## Functions
################################################################################
create_save_params_rule(step06_name, evaluation_path, step06_params, global_params)

def get_input_fa(wildcards):
    if wildcards.prefix == "gencode.none" and wildcards.txEnd == "0p" and wildcards.orf_filter == "all":
        return "/home/drihome/MRGuillermoPerez/RytenLab-Research/Resources/GENCODE/gencode.v48.transcripts.fa.gz"
    elif wildcards.prefix == "gencode.none" and wildcards.txEnd == "0p" and wildcards.orf_filter == "pc":
        return "/home/drihome/MRGuillermoPerez/RytenLab-Research/Resources/GENCODE/gencode.v48.pc_transcripts.fa.gz"
    else:
        return rules.gtxcutr_truncation.output.fa

def get_sequence_files(wildcards):
    return evaluation_samples.files[wildcards.sample_id].split(';')

## Rules
################################################################################
rule kallisto_index:
    input:
        fa = get_input_fa
    output:
        kdx = evaluation_path("kallisto_index/{prefix}.{orf_filter}.gtxcutr.w{width}.{txEnd}.kdx")
    log: evaluation_log_path("kallisto_index/{prefix}.{orf_filter}.gtxcutr.w{width}.{txEnd}.log")
    benchmark: evaluation_benchmark_path("kallisto_index/{prefix}.{orf_filter}.gtxcutr.w{width}.{txEnd}.tsv")
    params:
        tmp_dir = "tmp_kallisto_index_{prefix}_{orf_filter}_{width}_{txEnd}"
    conda: "../envs/kallisto.yaml"
    threads: 8
    shell:
        "kallisto index --threads {threads} --tmp {params.tmp_dir} -i {output.kdx} {input.fa} 2>&1 | tee {log}"
        # "kallisto index -i {output.kdx} {input.fa} 2>&1 | tee {log}"

rule kallisto_bus:
    input:
        kdx = rules.kallisto_index.output.kdx,
        files = get_sequence_files
    output:
        run_info = evaluation_path("kallisto_bus/{prefix}.{orf_filter}.w{width}.{txEnd}/{sample_id}/run_info.json"),
        bus = evaluation_path("kallisto_bus/{prefix}.{orf_filter}.w{width}.{txEnd}/{sample_id}/output.bus"),
        ec = evaluation_path("kallisto_bus/{prefix}.{orf_filter}.w{width}.{txEnd}/{sample_id}/matrix.ec"),
        tx = evaluation_path("kallisto_bus/{prefix}.{orf_filter}.w{width}.{txEnd}/{sample_id}/transcripts.txt"),
    params:
        tech = kb_settings["tech"],
        strand = kb_settings["strand"],
        out_dir = lambda w, output: Path(output.run_info).parent
    log: evaluation_log_path("kallisto_bus/{prefix}.{orf_filter}.w{width}.{txEnd}/{sample_id}.log")
    benchmark: evaluation_benchmark_path("kallisto_bus/{prefix}.{orf_filter}.w{width}.{txEnd}/{sample_id}.tsv")
    conda: "../envs/kallisto.yaml"
    threads: 12
    shell:
        """
        kallisto bus -t {threads} -i {input.kdx} -x {params.tech} {params.strand} -o {params.out_dir} --verbose {input.files} 2>&1 | tee {log}
        """

# module scUTRquant:
#     snakefile:
#         "../../tools/scUTRquant-0.5.0/Snakefile"
#         # github("Mayrlab/scUTRquant", path="Snakefile", tag="v0.5.0")
#     skip_validation: True
#     config: config["scUTRquant_config"]

# # use rule * from scUTRquant as scUTRquant_*

# use rule bustools_sort from scUTRquant as scUTR_bustools_sort with:
#     input: rules.kallisto_bus.output.bus
#     # input: evaluation_path("kallisto_bus/{prefix}.{orf_filter}.w{width}.{txEnd}/{sample_id}/output.bus"),
#     output: evaluation_path("kallisto_bus/{prefix}.{orf_filter}.w{width}.{txEnd}/{sample_id}/output.sorted.bus"),
    
# use rule bustools_whitelist from scUTRquant as scUTR_bustools_whitelist with:
#     input: rules.scUTR_bustools_bus.output
#     # input: evaluation_path("kallisto_bus/{prefix}.{orf_filter}.w{width}.{txEnd}/{sample_id}/output.sorted.bus"),
#     output: evaluation_path("kallisto_bus/{prefix}.{orf_filter}.w{width}.{txEnd}/{sample_id}/whitelist.txt"),

# def get_whitelist(wildcards):
#     if not config['bx_whitelist']:
#         return rules.scUTR_bustools_whitelist.output
#         # return evaluation_path("kallisto_bus/{prefix}.{orf_filter}.w{width}.{txEnd}/{sample_id}/whitelist.txt"),
#     else:
#         return config['bx_whitelist']

# use rule bustools_correct from scUTRquant as scUTR_bustools_correct with:
#     input: 
#         bus = rules.scUTR_bustools_sort.output
#         # bus = evaluation_path("kallisto_bus/{prefix}.{orf_filter}.w{width}.{txEnd}/{sample_id}/output.sorted.bus"),
#         bxs = get_whitelist
#     output: temp(evaluation_path("kallisto_bus/{prefix}.{orf_filter}.w{width}.{txEnd}/{sample_id}/output.corrected.bus"))

# use rule bustools_correct_sort from scUTRquant as scUTR_bustools_correct_sort with:
#     input: rules.scUTR_bustools_correct.output
#     # input: evaluation_path("kallisto_bus/{prefix}.{orf_filter}.w{width}.{txEnd}/{sample_id}/output.corrected.bus"),
#     output: temp(evaluation_path("kallisto_bus/{prefix}.{orf_filter}.w{width}.{txEnd}/{sample_id}/output.corrected.sorted.bus"))

# use rule bustools_count_txs from scUTRquant as scUTR_bustools_count_txs with:
#     input:
#         bus = get_input_busfile,
#         txs = rules.kallisto_bus.output.tx,
#         ec = rules.kallisto_bus.output.ec,
#         merge = 
#     input: rules.scUTR_bustools_correct.output
#     # input: evaluation_path("kallisto_bus/{prefix}.{orf_filter}.w{width}.{txEnd}/{sample_id}/output.corrected.bus"),
#     output: temp(evaluation_path("kallisto_bus/{prefix}.{orf_filter}.w{width}.{txEnd}/{sample_id}/output.corrected.sorted.bus"))


# use rule generate_tx_merge from scUTRquant as scUTR_generate_tx_merge with:
#     output: 

# rule generate_tx_merge:
#     input:
#         tsv=get_target_file('merge_tsv')
#     output:
#         "data/utrs/{target}/tx_merge.tsv"
#     shell:
#         """
#         tail -n+2 {input.tsv} | cut -f1,2 > {output}
#         """