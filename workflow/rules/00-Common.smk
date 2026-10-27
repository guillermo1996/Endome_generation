################################################################################
## 00 - Common: shared settings, helpers and step registration
################################################################################

import pandas as pd
import re
import json
import hashlib
import yaml
import copy
import itertools

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
registered_steps = []

################################################################################
## Load the pipeline inputs
################################################################################
### Config keys of pipeline_inputs and their wildcard in the file names
wildcard_names = {
    "merge_preset": "merge_method",
    "orf_filter_preset": "orf_filter",
    "trunc_site": "txEnd",
    "trunc_width": "width",
    "pseudobulk_preset": "pseudobulk_method",
    "satuRn_preset": "saturn_method",
}

def as_list(value) -> list:
    return value if isinstance(value, list) else [value]

def build_endome_table(defaults: dict, inputs: list) -> pd.DataFrame:
    """One row per ENDome to build. Every key becomes a column, and keys given as a list are combined."""
    rows = []
    for entry in inputs:
        entry = {**defaults, **entry}
        if entry["pre-assembled_gtf"] or entry["pre-truncated_gtf"]:
            entry["merge_preset"] = "ext"
        if entry["pre-truncated_gtf"]:
            entry["orf_filter_preset"] = "ext"

        values = {key: as_list(value) for key, value in entry.items()}
        rows += [dict(zip(values, combo)) for combo in itertools.product(*values.values())]

    table = pd.DataFrame(rows, dtype = object).rename(columns = wildcard_names)
    table["width"] = table["width"].astype(str)
    table["prefix"] = table["dataset"] + "." + table["group"] + "." + table["merge_method"]
    table["endome_name"] = table["prefix"] + "." + table["orf_filter"] + ".w" + table["width"] + "." + table["txEnd"]
    table["has_orfs"] = table["pre-truncated_gtf"].isna()

    ### First step each entry runs: pre-truncated GTFs start at 05, pre-assembled GTFs at 02 (manual_filter)
    table["first_step"] = 1
    table.loc[table["pre-assembled_gtf"].notna(), "first_step"] = 2
    table.loc[table["pre-truncated_gtf"].notna(), "first_step"] = 5

    return table.drop_duplicates().reset_index(drop = True)

def table_axes(*cols, where = None) -> dict:
    """Unique rows of `cols` as lists, for expand(template, zip, **table_axes(...)).
    With `where`, only the rows where that True/False column is True."""
    table = endome_table if where is None else endome_table[endome_table[where]]
    sub = table[list(cols)].drop_duplicates()
    return {col: sub[col].tolist() for col in cols}

def get_input_setting(prefix: str, key: str):
    """Value of `key` in pipeline_inputs for one transcriptome (dataset.group.merge_method)"""
    values = endome_table.loc[endome_table["prefix"] == prefix, key]
    if values.empty:
        raise ValueError(f"'{prefix}' is not in pipeline_inputs: check the dataset, group and entry point of the requested file")
    return values.iloc[0]

################################################################################
## Long-read sample sheet
################################################################################
def load_endome_samples(endome_table: pd.DataFrame, sample_sheet: str = None) -> pd.DataFrame:
    """Load and check the long-read sample sheet.

    Only the entries of ``pipeline_inputs`` without a ``gtf`` start from reads.
    Without them the sheet is not read, so ``input_endome_samples`` can be
    missing or null.

    Args:
        endome_table (pd.DataFrame): Output of ``build_endome_table()``.
        sample_sheet (str): Path of the sample sheet (``input_endome_samples``).

    Returns:
        pd.DataFrame: The sample sheet, or an empty table with its columns when
        no entry starts from reads.

    Raises:
        ValueError: If the sheet is needed but missing, or its content is invalid.
    """
    required_cols = ["dataset", "group", "sample_id", "file_type", "path"]

    ### Only the entries without a pre-assembled or pre-truncated GTF start from reads
    from_reads = endome_table["pre-assembled_gtf"].isna() & endome_table["pre-truncated_gtf"].isna()
    read_inputs = endome_table.loc[from_reads, ["dataset", "group"]].drop_duplicates()

    if read_inputs.empty:
        return pd.DataFrame(columns = required_cols)

    if not sample_sheet or not Path(sample_sheet).is_file():
        raise ValueError(f"Entries without `pre-assembled_gtf` start from reads, so `input_endome_samples` must point to "
                         f"the long-read sample sheet (got: {sample_sheet})")
    samples = pd.read_csv(sample_sheet, sep = "\t", dtype = str)

    ### Ensure all data is provided
    missing_cols = set(required_cols) - set(samples.columns)
    if missing_cols:
        raise ValueError(f"Missing columns in {sample_sheet}: {sorted(missing_cols)}")

    ### Ensure data is valid (only for the datasets that are run)
    selected = samples[samples["dataset"].isin(read_inputs["dataset"])]
    if not selected["file_type"].isin(["fastq", "bam"]).all():
        raise ValueError(f"file_type must be `fastq` or `bam` for every sample of the selected datasets")
    if selected["path"].isna().any():
        raise ValueError(f"Samples without a path: {selected[selected['path'].isna()]['sample_id'].tolist()}")
    if selected.duplicated(["dataset", "sample_id"]).any():
        raise ValueError(f"sample_id must be unique within a dataset")

    ### Every dataset x group that starts from reads has samples
    for dataset, group in read_inputs.itertuples(index = False):
        group_list = ["control", "case"] if group == "control_case" else [group]
        if not ((selected["dataset"] == dataset) & selected["group"].isin(group_list)).any():
            raise ValueError(f"No samples for {dataset}.{group} in {sample_sheet}. "
                             f"Check `group` in pipeline_inputs, or add a `pre-assembled_gtf` if it is an external transcriptome")
    return samples

################################################################################
## Helper Functions - For Analysis
################################################################################
def generate_input_samples_df(dataset: str, group: str) -> pd.DataFrame:
    """Samples of one dataset x group combination from the sample sheet.
    
    Args:
        dataset (str): Dataset name, as in the ``dataset`` columns.
        group (str): Sample group

    Returns:
        pd.DataFrame: One row per sample with columns ``sample_id``, ``group``, ``file_type`` and ``path``.
    """
    group_list = ["control", "case"] if group == "control_case" else [group]
    selected = endome_samples[(endome_samples["dataset"] == dataset) & (endome_samples["group"].isin(group_list))]
    return selected[["sample_id", "group", "file_type", "path"]].reset_index(drop = True)

def get_input_sample(wildcards) -> pd.Series:
    samples_df = generate_input_samples_df(wildcards.dataset, wildcards.group)
    return samples_df[samples_df["sample_id"] == wildcards.sample].iloc[0]

def generate_param_files():
    """``parameters.yaml`` of every registered step, for the datasets that run it"""
    targets = []
    for step in registered_steps:
        rows = endome_table.loc[endome_table["first_step"] <= int(step.name[:2]), ["dataset", "group"]].drop_duplicates()
        targets += expand(step.path("parameters.yaml"), zip, dataset = rows["dataset"].tolist(), group = rows["group"].tolist())
    return targets

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
        registered_steps.append(step)

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
    "pseudobulk_method": config.get("pseudobulk_preset", ["default"]),
    "saturn_method": config.get("satuRn_preset", ["default"]),
}
_test_data_links = {}

def register_test_data_link(rule_output) -> None:
    """Register a rule output to be linked into data/test_data/ by the `test_data` rule.

    Every dataset x group x merge method gets its own folder
    (data/test_data/{dataset}.{group}.{merge_method}/) and the files keep a generic
    "test" prefix, so a script's interactive block only needs to change the folder.
    Every combination of the wildcards in `_test_wildcard_values` that appear in
    the template (orf_filter, width, txEnd, pseudobulk_method, saturn_method) is linked.

    `rule_output` can be a single template or a list of them (e.g. an output
    entry built with expand(..., allow_missing=True)). Templates with an
    {endome_name} wildcard are supported: it stands for
    {prefix}.{orf_filter}.w{width}.{txEnd}. When it only appears in a folder
    (e.g. scUTRquant's data/sce/{endome_name}/hardy_snRNAseq.txs.Rds), it is
    prepended to the file name so each ENDome gets its own link
    (test.{orf_filter}.w{width}.{txEnd}.hardy_snRNAseq.txs.Rds).
    """
    if isinstance(rule_output, (list, tuple)):
        for single_output in rule_output:
            register_test_data_link(single_output)
        return

    endome_template = "{prefix}.{orf_filter}.w{width}.{txEnd}"
    template = str(rule_output)
    name_template = Path(template).name
    if "{endome_name}" in template and "{endome_name}" not in name_template:
        name_template = "{endome_name}." + name_template
    template = template.replace("{endome_name}", endome_template)
    name_template = name_template.replace("{endome_name}", endome_template)
    ## Step 02 outputs name the combination with three wildcards instead of {prefix}
    name_template = name_template.replace("{dataset}.{group}.{merge_method}", "{prefix}")

    ## Only the wildcards used by this template are expanded
    used_wildcards = {k: v for k, v in _test_wildcard_values.items() if "{" + k + "}" in template}

    for dataset, group, merge_method in itertools.product(test_datasets, test_groups, test_merge_methods):
        prefix = f"{dataset}.{group}.{merge_method}"
        for values in itertools.product(*used_wildcards.values()):
            wildcards = dict(zip(used_wildcards.keys(), values))
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