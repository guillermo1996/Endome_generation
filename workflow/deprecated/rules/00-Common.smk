import pandas as pd
import os
import re
import json
import hashlib
import yaml
import copy

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
log_path = config["log_path"]
benchmark_path = config["benchmark_path"]

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
yaml.Dumper.ignore_aliases = lambda *args: True

def create_save_params_rule(step_name, step_dir, step_params, global_params):
    def previous_params(global_params, step_params):
        new_global_params = global_params.copy()
        for k in step_params.keys():
            new_global_params.pop(k, None)

    def filter_steps(d, step_num):
        out = {}
        for key, val in d.items():
            prefix = key.split("-")[0]        # e.g. "01"
            if prefix.isdigit() and int(prefix) <= step_num:
                out[key] = val
        return out

    # Extract step number from prefix "XX-Name"
    match = re.match(r"^(\d{2})-", step_name)
    if not match:
        raise ValueError(f"step_name must start with 'NN-' (00-99): received {step_name}")
    step_num = int(match.group(1))
    
    # Resolve the directory by calling the function )
    step_dir = step_dir("")
    previous_global_dict = previous_params(global_params, step_params)

    # Dynamically create a rule to save parameters for a specific step
    rule:
        name: f"save_step{step_num}_params"
        message: f"--- Saving parameters to parameters.yaml ---"
        output:
            params_file = f"{step_dir}/parameters.yaml"
        run:
            previous_hash = compute_hash(previous_global_dict)

            param_data = {
                "step_dir": step_name,
                # "run_id": compute_hash(global_params),
                # "inherited_from_run_id": previous_hash if previous_global_dict else "",
                # "step_parameters": step_params,
                "all_parameters": filter_steps(global_params, step_num)
            }
            with open(output.params_file, "w") as f:
                yaml.dump(param_data, f, default_flow_style=False, sort_keys=False)

def compute_hash(params_dict, prev_step_hash = ""):
    """Generate a short hash from parameters"""
    if not params_dict:
        return ""
        
    # Ignore step naming (only care about the configurations)
    values_only = list(params_dict.values())
    values_only = [json.dumps(v, sort_keys=True) for v in values_only]
    values_only.sort()

    merged = "[" + ",".join(values_only) + "]"
    current_hash = hashlib.md5(merged.encode()).hexdigest()[:4]
    # if prev_step_hash != "":
    #     current_hash = f"{prev_step_hash}.{current_hash}"
    return current_hash

# def pretty_print_dict(d):
#     #take empty string
#     pretty_dict = ''  
    
#     #get items for dict
#     for k, v in d.items():
#         pretty_dict += f'{k}: \n'
#         for value in v:
#             pretty_dict += f'    {value}: {v[value]}\n'
#     #return result
#     return pretty_dict

def dict_to_readable(d, indent=2):
    lines = []
    for i, (key, value) in enumerate(d.items()):
        if indent == 0 and i > 0:  # Add blank line before top-level keys (except first)
            lines.append("")
        
        if isinstance(value, dict):
            lines.append("  " * indent + f"{key}:")
            lines.append(dict_to_readable(value, indent + 1))
        else:
            lines.append("  " * indent + f"{key}: {value}")
    return "\n".join(lines)


################################################################################
## Global variables
################################################################################
global_params = {}