#!/usr/bin/env bash
set -euo pipefail

# Snakemake test runner for multiple assembly configurations
# Usage: ./snakescript.sh

TESTS=(
    "st_ungui_gui"
    "st_ungui_gui_strict"
    "st_gui_gui"
    "st_ungui_ungui"
)

echo "=========================================="
echo "Running Snakemake Assembly Tests"
echo "Tests: ${#TESTS[@]}"
echo "=========================================="
echo

for test in "${TESTS[@]}"; do
    echo ">>> Running test: ${test}"
    echo "-------------------------------------------"
    
    if snakemake \
        --config assembly_test="${test}"; then
        echo "✓ Test '${test}' completed successfully"
    else
        echo "✗ Test '${test}' FAILED" >&2
        exit 1
    fi
    
    echo
done

echo "=========================================="
echo "All tests completed successfully!"
echo "=========================================="
