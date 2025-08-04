################################################################################
## 06 - Evaluation
################################################################################

## Variables
################################################################################
evaluation_path = lambda x: Path(results_output_path) / "06-Evaluation" / x

evaluation_samples = pd.read_csv(config['evaluation_sample_file'], index_col='sample_id')

## Functions
################################################################################
def get_sequence_files(wildcards):
    return evaluation_samples.files[wildcards.sample_id].split(';')

## Rules
################################################################################
rule kallisto_index:
    input:
        fa = rules.gtxcutr_truncation.output.fa
    output:
        kdx = evaluation_path("kallisto_index/{prefix}.gtxcutr.w{width}.{txEnd}.kdx")
    log: log_path("kallisto_index/{prefix}.gtxcutr.w{width}.{txEnd}.log")
    benchmark: benchmark_path("kallisto_index/{prefix}.gtxcutr.w{width}.{txEnd}.tsv")
    params:
        tmp_dir = "tmp_kallisto_index_{prefix}_{width}_{txEnd}"
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
        run_info = evaluation_path("06-Evaluation/kallisto/{prefix}.w{width}.{txEnd}/{sample_id}/run_info.json"),
    params:
        tech = config["tech"],
        strand = config["strand"],
        out_dir = lambda w, output: Path(output.run_info).parent
    log: log_path("kallisto/{prefix}.w{width}.{txEnd}/{sample_id}.log")
    benchmark: benchmark_path("kallisto/{prefix}.w{width}.{txEnd}/{sample_id}.tsv")
    conda: "../envs/kallisto.yaml"
    threads: 12
    shell:
        """
        kallisto bus -t {threads} -i {input.kdx} -x {params.tech} {params.strand} -o {params.out_dir} --verbose {input.files} 2>&1 | tee {log}
        """