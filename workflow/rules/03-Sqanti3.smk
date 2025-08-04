################################################################################
## Sqanti3 Artifact Removal
################################################################################

## Variables
################################################################################
sqanti3_path = lambda x: os.path.join(results_path, "03-Sqanti3", x)

### Configurations
sq3_path = config["sq3_path"]
sq3_version = config["sq3_version"]

## Functions
################################################################################

## Rules
################################################################################
rule download_sqanti3:
    output:
        sq3_dir = directory(sq3_path),
        sq3_qc = f"{sq3_path}/sqanti3_qc.py",
        sq3_filter = f"{sq3_path}/sqanti3_filter.py",
        sq3_rescue = f"{sq3_path}/sqanti3_rescue.py",
        sq3_env = f"{sq3_path}/SQANTI3.conda_env.yml",
        default_rules = f"{sq3_path}/src/utilities/filter/filter_default.json"
    params:
        version = sq3_version
    log: log_path("download_sqanti3.log")
    conda: "../envs/download_scripts.yaml"
    script: "../scripts/03a-Download_sqanti3.sh"

rule sqanti3_qc:
    input: 
        sq3_qc = f"{sq3_path}/sqanti3_qc.py",
        isoform = lambda wc: rules.manual_filter.output if wc.sq3_ref != "_reference" else ref_annotation,
        ref_gtf = ref_annotation,
        ref_fasta = ref_genome,
    output:
        classification = sqanti3_path("Sqanti3_QC{sq3_ref}/{prefix}_classification.txt"),
        corrected_gff = sqanti3_path("Sqanti3_QC{sq3_ref}/{prefix}_corrected.gtf"),
        corrected_fasta = sqanti3_path("Sqanti3_QC{sq3_ref}/{prefix}_corrected.fasta"),
    log: log_path("Sqanti3_QC{sq3_ref}/{prefix}_qc.log")
    benchmark: benchmark_path("Sqanti3_QC{sq3_ref}/{prefix}_qc.tsv")
    params:
        out_dir = sqanti3_path("Sqanti3_QC{sq3_ref}")
    wildcard_constraints: 
        sq3_ref = ".{0}|_.+"
    threads: 16
    conda: "../envs/sqanti3.yaml" # ".snakemake/conda/sq3_env" ~ Bug in current snakemake version https://github.com/snakemake/snakemake/issues/3192
    shell: 
        "python {input.sq3_qc} --force_id_ignore --skipORF --report skip "
        "-o {wildcards.prefix} -d {params.out_dir} -n {threads} "
        "{input.isoform} {input.ref_gtf} {input.ref_fasta} 2>&1 | tee {log}"
    # shell: 
    #     "python {input.sq3_qc} --force_id_ignore --skipORF --report skip --isoforms {input.isoform} "
    #     "--refGTF {input.ref_gtf} --refFasta {input.ref_fasta} "
    #     "-o {wildcards.prefix} -d {params.out_dir} -n {threads} 2>&1 | tee {log}"

rule sqanti3_filter:
    input:
        sq3_filter = f"{sq3_path}/sqanti3_filter.py",
        classification = sqanti3_path("Sqanti3_QC/{prefix}_classification.txt"),
        corrected_gtf = sqanti3_path("Sqanti3_QC/{prefix}_corrected.gtf")
    output: 
        filtered_gtf = sqanti3_path("Sqanti3_Filter/{prefix}.filtered.gtf"),
        filtered_classification = sqanti3_path("Sqanti3_Filter/{prefix}_RulesFilter_result_classification.txt")
    log: log_path("Sqanti3_Filter/{prefix}_filter.log")
    benchmark: benchmark_path("Sqanti3_Filter/{prefix}_filter.tsv")
    params:
        out_dir = sqanti3_path("Sqanti3_Filter")
    threads: 16
    conda: "../envs/sqanti3.yaml"
    shell: 
        "python {input.sq3_filter} rules "
        "--gtf {input.corrected_gtf} "
        "-o {wildcards.prefix} -d {params.out_dir} "
        "{input.classification} 2>&1 | tee {log}"
    # shell: 
    #     "python {input.sq3_filter} rules --sqanti_class {input.classification} "
    #     "--filter_gtf {input.corrected_gtf} "
    #     "-o {wildcards.prefix} -d {params.out_dir} -c {threads} 2>&1 | tee {log}"

rule sqanti3_rescue:
    input:
        sq3_rescue = f"{sq3_path}/sqanti3_rescue.py",
        filter_class = sqanti3_path("Sqanti3_Filter/{prefix}_RulesFilter_result_classification.txt"),
        ref_gtf = ref_annotation,
        ref_fasta = ref_genome,
        rescue_isoforms = sqanti3_path("Sqanti3_QC/{prefix}_corrected.fasta"),
        rescue_gtf = sqanti3_path("Sqanti3_Filter/{prefix}.filtered.gtf"),
        ref_classif = sqanti3_path("Sqanti3_QC_reference/sq3.reference_classification.txt"),
        sq3_rules = os.path.join(sq3_path, "src/utilities/filter/filter_default.json"),
    output: 
        rescued_gtf = sqanti3_path("Sqanti3_Rescue/{prefix}_rescued.gtf")
    log: log_path("Sqanti3_Rescue/{prefix}_rescue.log")
    benchmark: benchmark_path("Sqanti3_Rescue/{prefix}_rescue.tsv")
    params:
        out_dir = sqanti3_path("Sqanti3_Rescue")
    threads: 16
    conda: "../envs/sqanti3.yaml"
    shell:
        "python {input.sq3_rescue} rules "
        "--refGTF {input.ref_gtf} --refGenome {input.ref_fasta} "
        "--isoforms {input.rescue_isoforms} --gtf {input.rescue_gtf} "
        "--refClassif {input.ref_classif} "
        "--mode full -j {input.sq3_rules} "
        "-o {wildcards.prefix} -d {params.out_dir} {input.filter_class} 2>&1 | tee {log}"
    # shell:
    #     "python {input.sq3_rescue} rules --filter_class {input.filter_class} "
    #     "--refGTF {input.ref_gtf} --refFasta {input.ref_fasta} "
    #     "--rescue_isoforms {input.rescue_isoforms} --rescue_gtf {input.rescue_gtf} "
    #     "--refClassif {input.ref_classif} "
    #     "--mode full -j {input.sq3_rules} "
    #     "-o {wildcards.prefix} -d {params.out_dir} -c {threads} 2>&1 | tee {log}"
