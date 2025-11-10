#!/usr/bin/env bash

# Redirect stderr and stdout to snakemake log file
exec 2> "${snakemake_log[0]}"
exec 1> "${snakemake_log[0]}"

# Exit on any command failure
set -o pipefail

# Accept sq3_version as first argument, default to 5.3.6 if not provided
sq3_version="${snakemake_params[version]}"
echo "Downloading SQANTI3 version ${sq3_version}..."

# Create tools directory and enter it
mkdir -p tools
pushd tools

# Download, extract and remove SQANTI3 archive
wget -q -O sq3_tmp.zip "https://github.com/ConesaLab/SQANTI3/archive/refs/tags/v${sq3_version}.zip"
unzip -q -o sq3_tmp.zip
rm sq3_tmp.zip

# Remove existing Sqanti3 directory if it exists
rm -rf Sqanti3

# Rename downloaded directory to standard name
mv "SQANTI3-${sq3_version}" Sqanti3

# Remove unnecessary directories to save space
rm -rf Sqanti3/test
rm -rf Sqanti3/example
rm -rf Sqanti3/data

# Copy conda environment file to workflow/envs
cp Sqanti3/SQANTI3.conda_env.yml ../workflow/envs/sqanti3.yaml

# Return to original directory
popd
echo "SQANTI3 v${sq3_version} successfully installed in tools/Sqanti3"
