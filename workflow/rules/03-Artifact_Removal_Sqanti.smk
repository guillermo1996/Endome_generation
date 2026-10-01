################################################################################
## Artifact Removal - SQANTI3
################################################################################
## Alternative to 03-Artifact_Removal.smk (pigeon) that uses SQANTI3 instead:
## QC -> rules filter -> rescue. Include this file INSTEAD of the pigeon one,
## never both: they register the same step (03-Artifact_Removal) and share the
## `step03` name.
##
## Step 04 reads the pigeon outputs (rules.pigeon_filter.output.* and the
## Classify_Filter/ paths), so it has to point to the SQANTI3 outputs below
## before this file can replace the pigeon one.
##
## Config settings to add (config/config.yaml, section 03):
##
##   sqanti3_preset: default
##
##   sqanti3_settings:
##     default:
##       version: "v5.3.6"                                     # SQANTI3 release, downloaded into tools/SQANTI3-{version}/
##       qc_flags: "--force_id_ignore --skipORF --report skip" # extra flags for sqanti3_qc.py
##       filter_flags: ""                                      # extra flags for sqanti3_filter.py rules
##       rescue_mode: "full"                                   # automatic | full
##
## The command-line flags below follow the SQANTI3 v5.3 interface. Check them
## against `python sqanti3_qc.py --help` (and filter/rescue) for the version used.
################################################################################

## Debug: benchmark loading time
_start_time = time.perf_counter()

## Variables
################################################################################
### Register the step: generates the hash suffix, creates the parameter rule and returns the output_path helpers.
step03 = register_step(
    name="03-Artifact_Removal",
    params={
        "sqanti3_settings": build_tool_settings(config, "sqanti3_settings", "sqanti3_preset"),
    },
)

### Resolved settings used directly in the rule bodies below
sqanti3_settings = resolve_preset(config, "sqanti3_settings", "sqanti3_preset")

### SQANTI3 Download path
sqanti3_version = sqanti3_settings["version"]
sqanti3_dir = f"tools/SQANTI3-{sqanti3_version}"
sqanti3_tar_url = f"https://github.com/ConesaLab/SQANTI3/archive/refs/tags/{sqanti3_version}.tar.gz"

## Functions
################################################################################
def sqanti3_input(wildcards):
    parts = wildcards.prefix.split(".")

    if wildcards.prefix == "ref":
        return ref_annotation
    elif wildcards.prefix == "gencode":
        return expand(rules.manual_filter.output.gtf, dataset = wildcards.prefix, group = "none", merge_method = "ref")[0]
    else:
        return expand(rules.manual_filter.output.gtf, dataset = parts[0], group = parts[1], merge_method = parts[2])[0]

## Rules
################################################################################
rule download_sqanti3:
    message: "--- Downloading and Extracting SQANTI3 ---"
    output:
        sqanti3_dir = directory(sqanti3_dir),
        qc = f"{sqanti3_dir}/sqanti3_qc.py",
        filter = f"{sqanti3_dir}/sqanti3_filter.py",
        rescue = f"{sqanti3_dir}/sqanti3_rescue.py",
        filter_rules = f"{sqanti3_dir}/src/utilities/filter/filter_default.json"
    params:
        url = sqanti3_tar_url
    shell:
        """
        mkdir -p {output.sqanti3_dir}
        curl -L {params.url} | tar -xz -C {output.sqanti3_dir} --strip-components=1
        """

### The reference is classified once per step folder; rescue needs it.
rule sqanti3_qc_reference:
    message: "--- SQANTI3 QC - Reference annotation ---"
    input:
        qc = rules.download_sqanti3.output.qc,
        ref_gtf = ref_annotation,
        ref_fasta = ref_genome
    output:
        classification = step03.path("Sqanti3_QC_reference/reference_classification.txt")
    log: step03.logs("Sqanti3_QC_reference/reference_qc.log")
    benchmark: step03.benchmark("Sqanti3_QC_reference/reference_qc.tsv")
    params:
        out_dir = lambda w, output: str(Path(output.classification).parent),
        flags = sqanti3_settings.get("qc_flags", "")
    threads: 16
    conda: "../envs/sqanti3.yaml"
    shell:
        "python {input.qc} --isoforms {input.ref_gtf} --refGTF {input.ref_gtf} --refFasta {input.ref_fasta} "
        "{params.flags} -o reference -d {params.out_dir} -t {threads} 2>&1 | tee {log}"

rule sqanti3_qc:
    message: "--- SQANTI3 QC - {wildcards.prefix} ---"
    input:
        qc = rules.download_sqanti3.output.qc,
        isoforms = sqanti3_input,
        ref_gtf = ref_annotation,
        ref_fasta = ref_genome
    output:
        classification = step03.path("Sqanti3_QC/{prefix}_classification.txt"),
        corrected_gtf = step03.path("Sqanti3_QC/{prefix}_corrected.gtf"),
        corrected_fasta = step03.path("Sqanti3_QC/{prefix}_corrected.fasta")
    log: step03.logs("Sqanti3_QC/{prefix}_qc.log")
    benchmark: step03.benchmark("Sqanti3_QC/{prefix}_qc.tsv")
    params:
        out_dir = lambda w, output: str(Path(output.classification).parent),
        flags = sqanti3_settings.get("qc_flags", "")
    threads: 16
    conda: "../envs/sqanti3.yaml"
    shell:
        "python {input.qc} --isoforms {input.isoforms} --refGTF {input.ref_gtf} --refFasta {input.ref_fasta} "
        "{params.flags} -o {wildcards.prefix} -d {params.out_dir} -t {threads} 2>&1 | tee {log}"

rule sqanti3_filter:
    message: "--- SQANTI3 Filter - {wildcards.prefix} ---"
    input:
        filter = rules.download_sqanti3.output.filter,
        classification = rules.sqanti3_qc.output.classification,
        corrected_gtf = rules.sqanti3_qc.output.corrected_gtf
    output:
        filtered_gtf = step03.path("Sqanti3_Filter/{prefix}.filtered.gtf"),
        filtered_classification = step03.path("Sqanti3_Filter/{prefix}_RulesFilter_result_classification.txt")
    log: step03.logs("Sqanti3_Filter/{prefix}_filter.log")
    benchmark: step03.benchmark("Sqanti3_Filter/{prefix}_filter.tsv")
    params:
        out_dir = lambda w, output: str(Path(output.filtered_gtf).parent),
        flags = sqanti3_settings.get("filter_flags", "")
    threads: 16
    conda: "../envs/sqanti3.yaml"
    shell:
        "python {input.filter} rules --sqanti_class {input.classification} --filter_gtf {input.corrected_gtf} "
        "{params.flags} -o {wildcards.prefix} -d {params.out_dir} -c {threads} 2>&1 | tee {log}"

rule sqanti3_rescue:
    message: "--- SQANTI3 Rescue - {wildcards.prefix} ---"
    input:
        rescue = rules.download_sqanti3.output.rescue,
        filter_rules = rules.download_sqanti3.output.filter_rules,
        filter_class = rules.sqanti3_filter.output.filtered_classification,
        rescue_gtf = rules.sqanti3_filter.output.filtered_gtf,
        rescue_isoforms = rules.sqanti3_qc.output.corrected_fasta,
        ref_classif = rules.sqanti3_qc_reference.output.classification,
        ref_gtf = ref_annotation,
        ref_fasta = ref_genome
    output:
        rescued_gtf = step03.path("Sqanti3_Rescue/{prefix}_rescued.gtf")
    log: step03.logs("Sqanti3_Rescue/{prefix}_rescue.log")
    benchmark: step03.benchmark("Sqanti3_Rescue/{prefix}_rescue.tsv")
    params:
        out_dir = lambda w, output: str(Path(output.rescued_gtf).parent),
        mode = sqanti3_settings.get("rescue_mode", "full")
    threads: 16
    conda: "../envs/sqanti3.yaml"
    shell:
        "python {input.rescue} rules --filter_class {input.filter_class} "
        "--refGTF {input.ref_gtf} --refFasta {input.ref_fasta} "
        "--rescue_isoforms {input.rescue_isoforms} --rescue_gtf {input.rescue_gtf} "
        "--refClassif {input.ref_classif} --mode {params.mode} -j {input.filter_rules} "
        "-o {wildcards.prefix} -d {params.out_dir} -c {threads} 2>&1 | tee {log}"

## Debug: benchmark loading time
_log(f"\t+ {step03.name} imported in {time.perf_counter() - _start_time:.3f}s")
