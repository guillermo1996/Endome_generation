import pandas as pd
import os
import re

from snakemake.utils import min_version
from snakemake.utils import validate
from sys import stderr, stdout
from datetime import datetime
from pathlib import Path

# Set the minimum snakemake version
min_version("8.30.0")

wildcard_constraints:
    dataset = "Wood|Ebbert",
    group = "control|case|control_case",

################################################################################
## Load reference configuration
################################################################################
ref_genome = config["ref_genome"]
ref_annotation = config["ref_annotation"]

################################################################################
## Define samples to generate the ENDome
################################################################################
def generate_input_samples(dataset, group):
    endome_samples = []
    input_groups = ["control", "case"] if group == "control_case" else [group]

    if dataset == "Ebbert":
        for input_group in input_groups:
            samples = config.get(f"{input_group}_samples_ebbert", [])
            for sample_id in samples:
                sample_file = f"{sample_id}.fastq.gz"
                sample_path = Path(config["input_dir_ebbert"]) / sample_file

                endome_samples.append({
                    "sample_id": sample_id,
                    "group": input_group,
                    "path": sample_path,
                })

    elif dataset == "Wood":
        for input_group in input_groups:
            samples = config.get(f"{input_group}_samples_wood", [])
            for sample_id in samples:
                sample_file = f"{sample_id}.fastq"
                sample_path = Path(config["input_dir_wood"]) / sample_id / "Pychopper" / f"{sample_id}_full_length_reads.fastq"

                endome_samples.append({
                    "sample_id": sample_id,
                    "group": input_group,
                    "path": sample_path,
                })

    endome_samples_df = pd.DataFrame(endome_samples)

    return(endome_samples_df)

evaluation_samples = pd.read_csv(config['evaluation_sample_file'], index_col='sample_id')

################################################################################
## Define the Output paths
################################################################################
main_output_path = config["main_output_path"]
project_output_path = f"{{dataset}}.{{group}}.k{config["k_flag"]}"
results_output_path = Path(main_output_path) / project_output_path

log_path = lambda x: Path(results_output_path) / config["log_path"] / x
benchmark_path = lambda x: Path(results_output_path) / config["benchmark_path"] / x


rule print_input_samples:
    output: "test.log"
    run:
        print(generate_input_samples("Ebbert", "control_case"))

rule generate_input_df:
    output: 
        tsv = (Path(results_output_path) / "{dataset}.{group}.endome_samples_df.tsv").as_posix()
    run:
        endome_samples_df = generate_input_samples(wildcards.dataset, wildcards.group)
        endome_samples_df.to_csv(output.tsv, sep="\t", index=False)