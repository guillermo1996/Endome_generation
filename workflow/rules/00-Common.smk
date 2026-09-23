## Debug: benchmark loading time
_start_time = time.perf_counter()

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
main_output_path = config["main_output_path"]
project_output_path = f"{{dataset}}.{{group}}"

results_path = Path(main_output_path) / project_output_path
log_path = config["log_path"]
benchmark_path = config["benchmark_path"]

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
def compute_hash(params_dict: dict, prev_step_hash: str = "") -> str:
    """Generate a short, deterministic hash from a step-keyed parameter dict.

    The hash is computed over the *structured* parameter dict so that tool
    namespaces are preserved (two tools that happen to share a leaf-key name
    no longer collide). Only the preset *name* is stripped — renaming a preset
    while keeping the same values yields the same hash. Step-name keys are
    dropped (and the tool-namespaced bodies merged) so that the hash depends on
    parameter values, not on the step folder names (which themselves carry the
    hash).

    Args:
        params_dict (dict): Provenance dict keyed by step name, e.g.
            ``{"01-Alignment": {"minimap2_presets": {...}}, ...}``.
        prev_step_hash (str): Reserved for future explicit hash chaining;
            currently unused. Defaults to "".

    Returns:
        str: A 5-character hex digest, or "" when ``params_dict`` is empty.
    """
    if not params_dict:
        return ""

    # Drop the step-name keys, merging the tool-namespaced bodies. Tool keys
    # (e.g. "minimap2_presets") are unique across steps, so merging is safe and
    # keeps each tool's params under its own namespace.
    tool_params = {}
    for step_body in params_dict.values():
        tool_params.update(step_body)

    # Strip preset *names* but keep the nested structure intact.
    clean_params = strip_preset_names(tool_params)

    # Sort keys to ensure order-independence, then hash.
    merged = json.dumps(clean_params, sort_keys=True)

    return hashlib.md5(merged.encode()).hexdigest()[:5]

def strip_preset_names(value: object) -> object:
    """Recursively drop any "preset" metadata key, preserving nested structure.

    Keeps the hash sensitive to *which parameters* are set and *under which
    tool* (tool namespaces and the "params" wrapper are retained), but
    insensitive to the human-readable preset name.

    Args:
        value (object): Any JSON-like value (dict, list or scalar).

    Returns:
        object: The same structure with every ``"preset"`` key removed.
    """
    if isinstance(value, dict):
        return {k: strip_preset_names(v) for k, v in value.items() if k != "preset"}
    if isinstance(value, list):
        return [strip_preset_names(v) for v in value]
    return value

def resolve_preset(config: dict, presets_key: str, preset_key: str) -> dict:
    """Look up and validate the active preset for a tool.

    Validates that both the selector key and the named preset exist, raising a
    clear error otherwise. This is the single entry point for reading a tool's
    settings, so a typo fails fast at workflow-load time.

    Args:
        config (dict): The Snakemake ``config`` mapping.
        presets_key (str): Key of the presets block, e.g. ``"minimap2_presets"``.
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
    """Wrap a tool's preset(s) with their name(s) for provenance dumps.

    When ``preset_key`` is given, only the selected preset is recorded. When it
    is omitted, the whole presets block is recorded instead - use this for tools
    whose preset varies as a filename wildcard (e.g. ``orf_filter``), where every
    preset may be materialised within the same run and so all of them must be
    part of the step's provenance and hash.

    Args:
        config (dict): The Snakemake ``config`` mapping.
        presets_key (str): Key of the presets block, e.g. ``"minimap2_presets"``.
        preset_key (str): Key of the active-preset selector, e.g.
            ``"minimap2_preset"``. When None, all presets in ``presets_key`` are
            recorded. Defaults to None.

    Returns:
        dict: ``{"preset": <name(s)>, "params": <settings>}``. The ``preset``
        entry is for human-readable dumps only; ``params`` drives the hash. With
        no ``preset_key``, ``preset`` is the sorted list of preset names and
        ``params`` is the full ``{name: settings}`` block.

    Raises:
        KeyError: If the presets block (or the selector, when given) is missing
            from ``config``.
        ValueError: If the selected preset name is absent from the presets block.
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

    Bundles a step's (optionally hashed) folder name with helpers that build
    output, log and benchmark paths underneath it. One instance is returned by
    :func:`register_step` per rule module, replacing the per-module
    ``*_path`` / ``*_log_path`` / ``*_benchmark_path`` lambdas.

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
            name (str): Final (possibly hashed) step folder name.
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

        Args:
            relpath (str): Path relative to the step folder. Defaults to "".

        Returns:
            Path: ``{results_path}/{name}/{relpath}``.
        """
        return Path(self._results_path) / self.name / relpath

    def log(self, relpath: str = "") -> Path:
        """Build a path inside the step's log sub-folder.

        Args:
            relpath (str): Path relative to the log folder. Defaults to "".

        Returns:
            Path: ``{results_path}/{name}/{log_path}/{relpath}``.
        """
        return Path(self._results_path) / self.name / self._log_path / relpath

    def benchmark(self, relpath: str = "") -> Path:
        """Build a path inside the step's benchmark sub-folder.

        Args:
            relpath (str): Path relative to the benchmark folder. Defaults to "".

        Returns:
            Path: ``{results_path}/{name}/{benchmark_path}/{relpath}``.
        """
        return Path(self._results_path) / self.name / self._benchmark_path / relpath

def register_step(name: str, params: dict, extra_params: dict = None, save_params: bool = True) -> Step:
    """Register a pipeline step and return its path helper.

    Performs the boilerplate shared by every rule module, so a module only has
    to declare its name and presets:

    1. Records ``params`` in the global provenance dict ``global_params``.
    2. When ``config["use_hash"]`` is truthy, computes a parameter hash over the
       accumulated ``global_params`` and appends it to the folder name.
    3. Records ``extra_params`` (if given) in ``global_extra_params`` under the
       same (possibly hashed) step name. These are dumped into the step's
       ``parameters.yaml`` for provenance exactly like ``params``, but are never
       fed into ``compute_hash`` - use this for things that should be visible in
       the audit trail but shouldn't force a folder change (e.g. a selector like
       ``transcript_merge_method`` that's meant to vary as a filename wildcard
       instead of a hash-scoped setting).
    4. Optionally creates the rule that writes the step's ``parameters.yaml``.

    Args:
        name (str): Base step name of the form ``"NN-Title"`` (e.g.
            ``"01-Alignment"``); the ``NN-`` prefix is required.
        params (dict): Mapping of ``{presets_key: build_tool_settings(...)}`` for
            every tool configured in this step. Drives the hash.
        extra_params (dict): Additional provenance-only values to record in
            ``parameters.yaml``. Never contributes to the hash, for this step or
            any step registered afterwards. Defaults to None.
        save_params (bool): When True, also create the ``parameters.yaml`` dump
            rule for the step. Defaults to True.

    Returns:
        Step: Path helper bound to the step's (possibly hashed) folder name.
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
    """Create the rule that dumps a step's resolved parameters to YAML.

    The generated rule writes ``parameters.yaml`` into the step folder,
    capturing the provenance of every step registered up to and including this
    one (later steps are excluded, even though they'll still be registered
    into ``global_params`` before any rule actually runs).

    How often the rule fires is controlled by ``config["always_save_params"]``:

    * True - the rule's ``params`` carry a fresh timestamp, so the "params"
      rerun-trigger fires and every ``snakemake`` invocation rewrites all the
      ``parameters.yaml`` files.
    * False (default) - the ``params`` carry a digest of the YAML that would be
      written, so the rule only reruns when the file is missing or when the
      recorded parameters actually changed.

    Args:
        step (Step): The step whose ``parameters.yaml`` should be produced.

    Returns:
        None.

    Raises:
        ValueError: If ``step.name`` does not start with an ``NN-`` prefix.
    """
    # Extract step number from prefix "NN-Name"
    match = re.match(r"^(\d{2})-", step.name)
    if not match:
        raise ValueError(f"step name must start with 'NN-' (00-99): received {step.name}")
    step_num = int(match.group(1))

    # global_params keeps growing as later steps register, and the `run:`
    # block below only executes once the whole workflow has been parsed - so
    # without snapshotting here, every step's dump would end up containing
    # every other step's params too.
    params_snapshot = copy.deepcopy(global_params)

    # Merge in provenance-only extras (see register_step's extra_params) after
    # the snapshot is taken - they were never part of global_params, so they
    # never reached compute_hash for this step or any step registered so far.
    for step_name, extra in global_extra_params.items():
        params_snapshot.setdefault(step_name, {})
        params_snapshot[step_name] = {**params_snapshot[step_name], **extra}

    param_data = {
        "step_dir": step.name,
        "all_parameters": params_snapshot,
    }
    params_yaml = yaml.dump(param_data, default_flow_style=False, sort_keys=False)

    # The value fed to the "params" rerun-trigger decides how often the dump is
    # refreshed: a timestamp changes on every invocation (always rewrite), a
    # digest of the payload only changes when the parameters themselves do
    # (rewrite only when needed).
    always_save = config.get("always_save_params", False)
    params_digest = hashlib.md5(params_yaml.encode()).hexdigest()

    # Dynamically create a rule to save parameters for a specific step
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
## Global variables
################################################################################
global_params = {}
global_extra_params = {}

################################################################################
## Temporary rules
################################################################################
test_datasets = config.get("input_dataset", ["Ebbert"])
test_groups = config.get("input_group", ["control"])
test_merge_method = config.get("transcript_merge_method", ["stringtie"])[0]
_test_wildcard_defaults = {
    "orf_filter": config.get("orf_filter_preset", "pc"),
    "width": config.get("trunc_width", [500])[0],
    "txEnd": config.get("trunc_site", ["3p"])[0],
}
_test_data_links = []

def register_test_data_link(rule_output) -> None:
    template = str(rule_output)

    # The file keeps the generic "test" prefix, so the dataset x group combination
    # is carried by a sub-folder instead (data/test_data/{dataset}.{group}/).
    dest_name = expand(Path(template).name, prefix="test", **_test_wildcard_defaults)[0]

    for dataset in test_datasets:
        for group in test_groups:
            prefix = f"{dataset}.{group}.{test_merge_method}"
            source = expand(template, dataset=dataset, group=group, prefix=prefix, **_test_wildcard_defaults)[0]
            dest = Path("data/test_data") / f"{dataset}.{group}" / dest_name
            _test_data_links.append((source, str(dest)))

rule test_data:
    message: "--- Refreshing test-data snapshot in data/test_data/ ---"
    run:
        for source, target in _test_data_links:
            src = Path(source)
            if not src.exists():
                print(f"[test_data] skipping {target}: source not built ({src})", file=stderr)
                continue
            dst = Path(target)
            dst.parent.mkdir(parents=True, exist_ok=True)
            if dst.is_symlink() or dst.exists():
                dst.unlink()
            dst.symlink_to(src.resolve())
            print(f"[test_data] linked {target} -> {src}", file=stderr)

## Debug: benchmark loading time
_log(f"\t+ 00-Common imported in {time.perf_counter() - _start_time:.3f}s")