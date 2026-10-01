################################################################################
## 00 - Common: shared settings, helpers and step registration
################################################################################

## Debug: benchmark loading time
_start_time = time.perf_counter()

import pandas as pd
import re
import json
import hashlib
import yaml
import copy

from sys import stderr
from pathlib import Path

################################################################################
## Load reference configuration
################################################################################
ref_genome = config["ref_genome"]
ref_annotation = config["ref_annotation"]

################################################################################
## Define the Output paths
################################################################################
main_output_path = config["main_output_path"]
project_output_path = "{dataset}.{group}"

results_path = Path(main_output_path) / project_output_path
log_path = config["log_path"]
benchmark_path = config["benchmark_path"]

################################################################################
## Global variables
################################################################################
### Settings of every registered step, filled in by register_step()
global_params = {}
global_extra_params = {}

################################################################################
## Helper Functions - For Analysis
################################################################################
def generate_input_samples_df(dataset: str, group: str) -> pd.DataFrame:
    """Build the sample sheet for one dataset × group combination.

    Args:
        dataset (str): Dataset name, e.g. ``"Ebbert"`` or ``"Wood"``.
        group (str): Sample group: ``"control"``, ``"case"`` or
            ``"control_case"`` (which expands to both).

    Returns:
        pd.DataFrame: One row per sample with columns ``sample_id``, ``group``
        and ``path``. Empty if the dataset/group has no configured samples.
    """
    endome_samples = []

    def add_to_endome(samples: list, sample_id: str, sample_group: str, sample_path: Path) -> None:
        samples.append({
            "sample_id": sample_id,
            "group": sample_group,
            "path": sample_path
        })

    # Experimental groups to use
    group_list = ["control", "case"] if group == "control_case" else [group]

    # Logic for the Ebbert dataset
    if dataset == "Ebbert":
        for input_group in group_list:
            ebbert_samples = config.get(f"{input_group}_samples_ebbert", [])

            for sample_id in ebbert_samples:
                sample_file = f"{sample_id}.fastq.gz"
                sample_path = Path(config["input_dir_ebbert"]) / sample_file

                add_to_endome(endome_samples, sample_id, input_group, sample_path)
    elif dataset == "Wood":
        for input_group in group_list:
            wood_samples = config.get(f"{input_group}_samples_wood", [])

            for sample_id in wood_samples:
                sample_file = f"{sample_id}.fastq"
                sample_path = Path(config["input_dir_wood"]) / sample_id / "Pychopper" / f"{sample_id}_full_length_reads.fastq"

                add_to_endome(endome_samples, sample_id, input_group, sample_path)
    
    # Return the data.frame
    endome_samples_df = pd.DataFrame(endome_samples)
    return endome_samples_df

################################################################################
## Helper Functions - For Snakemake
################################################################################
def compute_hash(params_dict: dict) -> str:
    """Build the 5-character hash that is appended to a step's folder name.

    ``params_dict`` holds the settings of the current step and of every step
    registered before it, so a change in any earlier step also changes the
    hash of all later steps.

    Only setting values are hashed. The folder names of the steps and the
    names of the selected presets are left out, so renaming a preset without
    changing its values keeps the same hash.

    Args:
        params_dict (dict): Settings per step, e.g.
            ``{"01-Alignment": {"minimap2_settings": {...}}, ...}``.

    Returns:
        str: The first 5 characters of the MD5 digest, or "" when
        ``params_dict`` is empty.
    """
    if not params_dict:
        return ""

    # Drop the step names and pool the tool blocks together.
    tool_params = {}
    for step_body in params_dict.values():
        tool_params.update(step_body)

    # Remove the preset names, so only the values are hashed.
    clean_params = strip_preset_names(tool_params)

    # Sort the keys, so the order of the settings in the config doesn't matter.
    merged = json.dumps(clean_params, sort_keys=True)

    return hashlib.md5(merged.encode()).hexdigest()[:5]

def strip_preset_names(value: object) -> object:
    """Remove every ``"preset"`` key (the preset names), so only the setting
    values are hashed.

    Args:
        value (object): Any JSON-like value (dict, list or scalar).

    Returns:
        object: The same structure without the ``"preset"`` keys.
    """
    if isinstance(value, dict):
        return {k: strip_preset_names(v) for k, v in value.items() if k != "preset"}
    if isinstance(value, list):
        return [strip_preset_names(v) for v in value]
    return value

def resolve_preset(config: dict, presets_key: str, preset_key: str) -> dict:
    """Look up and validate the active preset for a tool.

    Validates that both the selector key and the named preset exist, raising a
    clear error otherwise.

    Args:
        config (dict): The Snakemake ``config`` mapping.
        presets_key (str): Key of the presets block, e.g. ``"minimap2_settings"``.
        preset_key (str): Key of the active-preset selector, e.g.
            ``"minimap2_preset"``.

    Returns:
        dict: The resolved settings dict for the selected preset.

    Raises:
        KeyError: If the selector or presets block is missing from ``config``.
        ValueError: If the selected preset name is absent from the presets block.
    """
    if preset_key not in config:
        raise KeyError(f"Missing preset selector '{preset_key}' in config.")
    if presets_key not in config:
        raise KeyError(f"Missing presets block '{presets_key}' in config.")

    preset_name = config[preset_key]
    available = config[presets_key]

    if preset_name not in available:
        raise ValueError(
            f"Preset '{preset_name}' (selected via '{preset_key}') not found in "
            f"'{presets_key}'. Available presets: {sorted(available)}"
        )

    return available[preset_name]

def build_tool_settings(config: dict, presets_key: str, preset_key: str = None) -> dict:
    """Collect a tool's settings for the step's hash.

    With ``preset_key``, only the preset selected in the config is recorded.
    Use this when the preset is not part of the output filenames.

    Without ``preset_key``, every preset in the block is recorded. Use this when
    the preset is part of the filenames (e.g. ``orf_filter``).

    Args:
        config (dict): The Snakemake ``config`` mapping.
        presets_key (str): Key of the presets block, e.g. ``"minimap2_settings"``.
        preset_key (str): Key of the active-preset selector, e.g.
            ``"minimap2_preset"``. When None, all presets in ``presets_key`` are
            recorded. Defaults to None.

    Returns:
        dict: ``{"preset": <name(s)>, "params": <settings>}``

    Raises:
        KeyError: If the presets block is missing from ``config``.
        ValueError: If the selected preset name is absent from the presets
        block.
    """
    if preset_key is None:
        if presets_key not in config:
            raise KeyError(f"Missing presets block '{presets_key}' in config.")
        available = config[presets_key]
        return {
            "preset": sorted(available),  # human-readable
            "params": available,          # hashed values: every preset
        }

    return {
        "preset": config[preset_key],                          # human-readable
        "params": resolve_preset(config, presets_key, preset_key),  # hashed values
    }

class Step:
    """Output-path helper bound to a single pipeline step.

    Bundles a step's folder name with helpers that build output, log and
    benchmark paths underneath it.

    Attributes:
        name (str): Final folder name, including the hash suffix when hashing is
            enabled (e.g. ``"01-Alignment-1d7a6"``).
        hash (str): The parameter hash, or "" when hashing is disabled.
    """

    def __init__(
        self,
        name: str,
        results_path: Path,
        log_path: str,
        benchmark_path: str,
        hash: str = "",
    ) -> None:
        """Initialise the step path helper.

        Args:
            name (str): Final step folder name.
            results_path (Path): Base results path (may contain ``{wildcards}``).
            log_path (str): Log sub-folder name, relative to the step folder.
            benchmark_path (str): Benchmark sub-folder name, relative to the step.
            hash (str): The parameter hash, or "". Defaults to "".
        """
        self.name = name
        self.hash = hash
        self._results_path = results_path
        self._log_path = log_path
        self._benchmark_path = benchmark_path

    def path(self, relpath: str = "") -> Path:
        """Build a path to an output file/dir inside the step folder.

        Returns:
            Path: ``{results_path}/{name}/{relpath}``.
        """
        return Path(self._results_path) / self.name / relpath

    def logs(self, relpath: str = "") -> Path:
        """Build a path inside the step's log sub-folder.

        Returns:
            Path: ``{results_path}/{name}/{log_path}/{relpath}``.
        """
        return Path(self._results_path) / self.name / self._log_path / relpath

    def benchmark(self, relpath: str = "") -> Path:
        """Build a path inside the step's benchmark sub-folder.

        Returns:
            Path: ``{results_path}/{name}/{benchmark_path}/{relpath}``.
        """
        return Path(self._results_path) / self.name / self._benchmark_path / relpath

def register_step(name: str, params: dict, extra_params: dict = None, save_params: bool = True) -> Step:
    """Set up a step's output folder. Every rule file calls this once, at the top.

    It does three things:

    1. Names the folder. With ``use_hash: True``, the folder gets a hash built
       from the settings of this step and of all the steps registered before it
       (e.g. ``01-Alignment-befcf``). With ``use_hash: False``, it is just the
       step name.
    2. Keeps a record of the settings, which is written to ``parameters.yaml``
       inside the folder.
    3. Returns a ``Step`` object, used in the rules to build the output, log and
       benchmark paths inside that folder.

    Args:
        name (str): Step name, starting with its number (e.g. ``"01-Alignment"``).
        params (dict): ``{presets_key: build_tool_settings(...)}`` for each tool
            of the step. Used for the hash.
        extra_params (dict): Values only recorded in ``parameters.yaml``.
            Defaults to None.
        save_params (bool): Whether to create the rule that writes
            ``parameters.yaml``. Defaults to True.

    Returns:
        Step: Path helper for the step's folder.
    """
    global_params[name] = params

    step_hash = ""
    final_name = name
    if config.get("use_hash", False):
        step_hash = compute_hash(global_params)
        final_name = f"{name}-{step_hash}"
        # Re-key under the hashed name so downstream dumps reflect the folder.
        global_params[final_name] = global_params.pop(name)

    if extra_params:
        global_extra_params[final_name] = extra_params

    step = Step(final_name, results_path, log_path, benchmark_path, hash=step_hash)

    if save_params:
        create_save_params_rule(step)

    return step

def create_save_params_rule(step: Step) -> None:
    """Create the rule that writes ``parameters.yaml`` into a step's folder.

    The file lists the settings of this step and of all the steps before it.
    It is rewritten when those settings change, or on every run with
    ``always_save_params: True``.

    Args:
        step (Step): The step to write the file for.
    """
    # Step number, taken from the "NN-" prefix of the name
    match = re.match(r"^(\d{2})-", step.name)
    if not match:
        raise ValueError(f"step name must start with 'NN-' (00-99): received {step.name}")
    step_num = int(match.group(1))

    # Copy the settings now, before the later steps are registered.
    params_snapshot = copy.deepcopy(global_params)

    # Add the extra_params, which are recorded but not hashed.
    for step_name, extra in global_extra_params.items():
        params_snapshot.setdefault(step_name, {})
        params_snapshot[step_name] = {**params_snapshot[step_name], **extra}

    param_data = {
        "step_dir": step.name,
        "all_parameters": params_snapshot,
    }
    params_yaml = yaml.dump(param_data, default_flow_style=False, sort_keys=False)

    # Snakemake reruns the rule when its params change: a timestamp changes on
    # every run, the digest only when the settings do.
    always_save = config.get("always_save_params", False)
    params_digest = hashlib.md5(params_yaml.encode()).hexdigest()

    # One rule per step, e.g. "save_step01_params"
    rule:
        name: f"save_step{step_num}_params"
        message: "--- Saving parameters to parameters.yaml ---"
        output:
            params_file = str(step.path("parameters.yaml"))
        params:
            _trigger = (lambda wildcards: time.time()) if always_save else params_digest
        run:
            with open(output.params_file, "w") as f:
                f.write(params_yaml)

################################################################################
## Temporary rules
################################################################################
import itertools

test_datasets = config.get("input_dataset", ["Ebbert"])
test_groups = config.get("input_group", ["control"])
test_merge_methods = config.get("merge_preset", ["st_ref"])
_test_wildcard_values = {
    "orf_filter": config.get("orf_filter_preset", ["pc"]),
    "width": config.get("trunc_width", [500]),
    "txEnd": config.get("trunc_site", ["3p"]),
}
_test_data_links = {}

def register_test_data_link(rule_output) -> None:
    """Register a rule output to be linked into data/test_data/ by the `test_data` rule.

    Every dataset x group x merge method gets its own folder
    (data/test_data/{dataset}.{group}.{merge_method}/) and the files keep a generic
    "test" prefix, so a script's interactive block only needs to change the folder.
    Every combination of orf_filter, width and txEnd is linked.
    """
    template = str(rule_output)
    ## Step 02 outputs name the combination with three wildcards instead of {prefix}
    name_template = Path(template).name.replace("{dataset}.{group}.{merge_method}", "{prefix}")

    for dataset, group, merge_method in itertools.product(test_datasets, test_groups, test_merge_methods):
        prefix = f"{dataset}.{group}.{merge_method}"
        for values in itertools.product(*_test_wildcard_values.values()):
            wildcards = dict(zip(_test_wildcard_values.keys(), values))
            source = expand(template, dataset=dataset, group=group, merge_method=merge_method, prefix=prefix, **wildcards)[0]
            dest = Path("data/test_data") / prefix / expand(name_template, prefix="test", **wildcards)[0]
            _test_data_links[str(dest)] = source

rule test_data:
    message: "--- Refreshing test-data snapshot in data/test_data/ ---"
    run:
        for target, source in _test_data_links.items():
            src = Path(source)
            if not src.exists():
                print(f"[test_data] skipping {target}: source not built ({src})", file=stderr)
                continue
            dst = Path(target)
            dst.parent.mkdir(parents=True, exist_ok=True)
            if dst.is_symlink() or dst.exists():
                dst.unlink()
            # Relative link, so it works from any mount of the repo (/home/drihome/... or /home/...)
            dst.symlink_to(os.path.relpath(src, dst.parent))
            print(f"[test_data] linked {target} -> {src}", file=stderr)


## Debug: benchmark loading time
_log(f"\t+ 00-Common imported in {time.perf_counter() - _start_time:.3f}s")