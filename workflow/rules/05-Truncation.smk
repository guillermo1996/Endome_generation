################################################################################
## 05 - ENDome truncation
################################################################################

## Debug: benchmark loading time
_start_time = time.perf_counter()

## Variables
################################################################################
### Register the step: generates the hash suffix, creates the parameter rule and returns the output_path helpers.
step05 = register_step(
    name="05-Truncation",
    params={
        "txendcutr_settings": build_tool_settings(config, "txendcutr_settings", "txendcutr_preset"),
        "mmseqs2_settings": build_tool_settings(config, "mmseqs2_settings", "mmseqs2_preset"),
    },
    extra_params={
        "scoring_settings": config["scoring_settings"]
    }
)

### Resolved settings used directly in the rule bodies below
txendcutr_settings = resolve_preset(config, "txendcutr_settings", "txendcutr_preset")
mmseq2_settings = resolve_preset(config, "mmseqs2_settings", "mmseqs2_preset")

### Bin information scoring settings.
scoring_settings = config.get("scoring_settings", {}).get("default")

## Truncation Rules
################################################################################
rule txendcutr_truncation:
    message: """--- Transcriptome Truncation ---"""
    input:
        gtf = rules.ORF_filtration.output.gtf_filter
    output:
        gtf = step05.path("txendcutr/{prefix}.{orf_filter}.w{width}.{txEnd}.txendcutr.gtf"),
        fa = step05.path("txendcutr/{prefix}.{orf_filter}.w{width}.{txEnd}.txendcutr.fa.gz"),
        transcript_overlap = step05.path("txendcutr/{prefix}.{orf_filter}.w{width}.{txEnd}.txendcutr.overlaps.tsv"),
        merge_table = step05.path("txendcutr/{prefix}.{orf_filter}.w{width}.{txEnd}.txendcutr.merge.tsv")
    log: step05.logs("txendcutr/{prefix}.{orf_filter}.w{width}.{txEnd}.txendcutr.log")
    benchmark: step05.benchmark("txendcutr/{prefix}.{orf_filter}.w{width}.{txEnd}.txendcutr.tsv")
    params:
        mergeDist = txendcutr_settings["merge_distance"],
        genome = txendcutr_settings["genome"],
        script = lambda wc: "05a-txendcutr.R" if wc.width.isdigit() else "05a-txendcutr_UTR.R"
    conda: "../envs/r.yaml"
    threads: 12
    script: "../scripts/{params.script}" # Modified `txendcutr.R` script to include my version of the package

## Bin Information Content Rules
################################################################################
rule endome_export_proteins:
    message: """--- Exporting the ORF protein set for MMseqs2 ---"""
    input:
        protein_fa = rules.ORF_annotate.output.protein_fa,
        gtf_filter = rules.ORF_filtration.output.gtf_filter
    output:
        faa = step05.path("MMseqs2/{prefix}.{orf_filter}.proteins.faa"),
        protein_map = step05.path("MMseqs2/{prefix}.{orf_filter}.protein_map.tsv"),
    log: step05.logs("MMseqs2/{prefix}.{orf_filter}.export_proteins.log")
    benchmark: step05.benchmark("MMseqs2/{prefix}.{orf_filter}.export_proteins.tsv")
    conda: "../envs/r.yaml"
    script: "../scripts/05b-Export_Proteins.R"

rule endome_mmseqs_search:
    message: """--- All-vs-all MMseqs2 alignment of the ORF proteins ---"""
    input:
        faa = rules.endome_export_proteins.output.faa
    output:
        pairs = step05.path("MMseqs2/{prefix}.{orf_filter}.pairs.tsv")
    log: step05.logs("MMseqs2/{prefix}.{orf_filter}.mmseqs2.log")
    benchmark: step05.benchmark("MMseqs2/{prefix}.{orf_filter}.mmseqs2.tsv")
    params:
        sensitivity = mmseq2_settings["sensitivity"],
        evalue = mmseq2_settings["evalue"],
        max_seqs = mmseq2_settings["max_seqs"],
        max_accept = mmseq2_settings["max_accept"],
        max_rejected = mmseq2_settings["max_rejected"],
        alignment_mode = mmseq2_settings["alignment_mode"],
        coverage = mmseq2_settings["coverage"],
        cov_mode = mmseq2_settings["cov_mode"]
    conda: "../envs/mmseqs2.yaml"
    threads: 64
    shadow: "minimal"
    shell:
        """
        mmseqs createdb {input.faa} protDB 2>&1 | tee {log}

        mmseqs search protDB protDB alnDB tmp \
            -s {params.sensitivity} \
            -e {params.evalue} \
            --min-ungapped-score 0 \
            --max-seqs {params.max_seqs} \
            --max-accept {params.max_accept} \
            --max-rejected {params.max_rejected} \
            --alignment-mode {params.alignment_mode} \
            -a 1 \
            --add-self-matches 1 \
            -c {params.coverage} \
            --cov-mode {params.cov_mode} \
            --threads {threads} 2>&1 | tee -a {log}

        mmseqs convertalis protDB protDB alnDB {output.pairs} \
            --format-mode 4 \
            --format-output "query,target,fident,nident,alnlen,mismatch,qlen,tlen,qcov,tcov,raw,bits,evalue" \
            --threads {threads} 2>&1 | tee -a {log}
        """

rule endome_build_db:
    message: """--- Building the ENDome DuckDB ---"""
    input:
        gtf_filter = rules.ORF_filtration.output.gtf_filter,
        isoform_summary = rules.ORF_categorization.output.isoform_summary,
        protein_fa = rules.ORF_annotate.output.protein_fa,
        utr5_fa = rules.ORF_annotate.output.utr5_fa,
        utr3_fa = rules.ORF_annotate.output.utr3_fa,
        merge_table = rules.txendcutr_truncation.output.merge_table,
        overlap_table = rules.txendcutr_truncation.output.transcript_overlap,
        ref_genome = ref_genome
    output:
        duckdb = step05.path("DuckDB/{prefix}.{orf_filter}.w{width}.{txEnd}.duckdb")
    log: step05.logs("DuckDB/{prefix}.{orf_filter}.w{width}.{txEnd}.build_db.log")
    benchmark: step05.benchmark("DuckDB/{prefix}.{orf_filter}.w{width}.{txEnd}.build_db.tsv")
    conda: "../envs/r.yaml"
    threads: 4
    script: "../scripts/05c-Build_ENDome_DB.R"

rule endome_compute_scores:
    message: """--- Computing the Bin Information Content scores ---"""
    input:
        duckdb = rules.endome_build_db.output.duckdb,
        pairs = step05.path("MMseqs2/{prefix}.{orf_filter}.pairs.tsv"),
        protein_map = step05.path("MMseqs2/{prefix}.{orf_filter}.protein_map.tsv"),
    output:
        done = touch(step05.path("DuckDB/{prefix}.{orf_filter}.w{width}.{txEnd}.scores.done"))
    log: step05.logs("characterization/{prefix}.{orf_filter}.w{width}.{txEnd}.scores.log")
    benchmark: step05.benchmark("characterization/{prefix}.{orf_filter}.w{width}.{txEnd}.scores.tsv")
    params:
        protein_metric = scoring_settings["protein_metric"],
        length_scaling = scoring_settings["length_scaling"],
        member_weights = scoring_settings["member_weights"],
        component_weights = scoring_settings["component_weights"]
    conda: "../envs/r.yaml"
    threads: 4
    script: "../scripts/05d-Compute_Scores.R"

rule endome_report:
    message: """--- Rendering the ENDome characterization report ---"""
    input:
        duckdb = rules.endome_build_db.output.duckdb,
        done = rules.endome_compute_scores.output.done
    output:
        html = step05.path("characterization/reports/{prefix}.{orf_filter}.w{width}.{txEnd}.html")
    log: step05.logs("characterization/{prefix}.{orf_filter}.w{width}.{txEnd}.report.log")
    params:
        merge_distance = txendcutr_settings["merge_distance"],
        protein_metric = scoring_settings["protein_metric"],
        length_scaling = scoring_settings["length_scaling"],
        member_weights = scoring_settings["member_weights"],
        component_weights = scoring_settings["component_weights"]
    conda: "../envs/r.yaml"
    script: "../scripts/05e-ENDome_Report.Rmd"


## Pseudoalignment RUles
################################################################################
rule kallisto_index:
    message: "--- Building kallisto index for the ENDome ---"
    input:
        fa = step05.path("txendcutr/{prefix}.{orf_filter}.w{width}.{txEnd}.txendcutr.fa.gz")
    output:
        kdx = step05.path("kallisto_index/{prefix}.{orf_filter}.w{width}.{txEnd}.kdx")
    log: step05.logs("kallisto_index/{prefix}.{orf_filter}.w{width}.{txEnd}.log")
    benchmark: step05.benchmark("kallisto_index/{prefix}.{orf_filter}.w{width}.{txEnd}.tsv")
    conda: "../envs/scutrquant-kallisto-bustools.yaml"
    shell:
        "kallisto index -i {output.kdx} {input.fa} 2>&1 | tee {log}"

# To Do: Check Tama Collapse approach to solve the cascade binning issue in gtxcutr: https://github.com/GenomeRIK/tama/wiki/Tama-Collapse

register_test_data_link(rules.txendcutr_truncation.output.gtf)
register_test_data_link(rules.txendcutr_truncation.output.transcript_overlap)
register_test_data_link(rules.txendcutr_truncation.output.merge_table)
register_test_data_link(rules.endome_export_proteins.output.faa)
register_test_data_link(rules.endome_export_proteins.output.protein_map)
register_test_data_link(rules.endome_mmseqs_search.output.pairs)
register_test_data_link(rules.endome_build_db.output.duckdb)

## Debug: benchmark loading time
_log(f"\t+ {step05.name} imported in {time.perf_counter() - _start_time:.3f}s")
