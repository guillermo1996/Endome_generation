import pandas as pd
import os
import re
import json
import hashlib
import yaml

from sys import stderr, stdout
from datetime import datetime
from pathlib import Path

################################################################################
## Load reference configuration
################################################################################
ref_genome = config["ref_genome"]
ref_annotation = config["ref_annotation"]

################################################################################
## Define the Output paths
################################################################################
main_output_path = config["main_output_path"] # "results"
project_output_path = f"{{dataset}}.{{group}}"
results_path = Path(main_output_path) / project_output_path

log_path = lambda x, hash = "0000": Path(results_path) / config["log_path"] / f"{str(Path(x).parent)}_{hash}" / Path(x).name
benchmark_path = lambda x, hash = "0000": Path(results_path) / config["benchmark_path"] / f"{str(Path(x).parent)}_{hash}" / Path(x).name

################################################################################
## Helper Functions
################################################################################
def generate_input_samples_df(dataset, group):
    endome_samples = []
    
    #####################################
    ## Remove after testing
    # if config.get("use_toy_data", False):
    #     return pd.DataFrame({"sample_id": ["HG00154.chrom11.ILLUMINA.bwa.GBR.low_coverage.20120522", "HG00154.chrom20.ILLUMINA.bwa.GBR.low_coverage.20120522"]})
    #####################################

    # Experimental groups to use (if control_case is provided, then use both)
    input_groups = ["control", "case"] if group == "control_case" else [group]
    
    # Logic for the Ebbert dataset
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

    # Logic for the Wood dataset
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

    # Return the data.frame
    endome_samples_df = pd.DataFrame(endome_samples)
    return(endome_samples_df)

def create_save_params_rule(step_num, step_name, step_dir, step_params):
    # Dynamically create a rule to save parameters for a specific step
    rule:
        name: f"save_step{step_num}_params"
        message: f"--- Saving parameters to parameters.yaml ----"
        output:
            params_file = f"{step_dir}/parameters.yaml"
        run:
            cumulative_params = get_cumulative_params(step_num)
            previous_hash = compute_hash(get_cumulative_params(step_num-1))

            param_data = {
                "step_dir": step_name,
                "run_id": compute_hash(cumulative_params),
                "step_parameters": step_params,
                "all_parameters": cumulative_params
                # "hash_input": json.dumps(cumulative_params, sort_keys=True)
            }

            if previous_hash:
                param_data["inherited_from_run_id"] = previous_hash
            
            with open(output.params_file, "w") as f:
                yaml.dump(param_data, f, default_flow_style=False, sort_keys=False)

def compute_hash(params_dict):
    """Generate a short hash from parameters"""
    param_str = json.dumps(params_dict, sort_keys=True)
    return hashlib.md5(param_str.encode()).hexdigest()[:4]

def get_cumulative_params(step_num):
    cumulative = {}
    for i in range(1, step_num + 1):
        varname = f"step{i:02d}_params"
        if varname in globals():
            cumulative.update(globals()[varname])
    return cumulative
