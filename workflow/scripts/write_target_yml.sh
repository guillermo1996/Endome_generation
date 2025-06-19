#!/usr/bin/env bash

set -euo pipefail

TARGET_PATH=$1
TARGET_GTF=

path: "/home/MRGuillermoPerez/RytenLab-Research/38-Endome_generation/results/gtxcutr/"
  genome: "hg38"
  gtf: "sq3.annotated.gtxcutr.w500.3p.gtf"
  kdx: "sq3.annotated.gtxcutr.w500.3p.kdx"
  merge_tsv: "sq3.annotated.gtxcutr.w500.3p.merge.tsv"

TARGET_FILE=$1       # config/targets.yml
PREFIX=$2            # sq3.annotated
GTX_PATH=$3          # results/<build>/gtxcutr
WIDTH=$4             # 500
TXEND=$5             # 3p
GENOME="hg38"

GTF_FILE="${PREFIX}.gtxcutr.w${WIDTH}.${TXEND}.gtf"
KDX_FILE="${PREFIX}.gtxcutr.w${WIDTH}.${TXEND}.kdx"
MERGE_FILE="${PREFIX}.gtxcutr.w${WIDTH}.${TXEND}.merge.tsv"

TAG="endome_${TXEND}_w${WIDTH}"

# mkdir -p "$(dirname "$TARGET_FILE")"

{
    echo "${TAG}:"
    echo "  path: \"$GTX_PATH\""
    echo "  genome: \"$GENOME\""
    echo "  gtf: \"$GTF_FILE\""
    echo "  kdx: \"$KDX_FILE\""
    echo "  merge_tsv: \"$MERGE_FILE\""
    echo "  tx_annots: null"
    echo "  tx_annots_csv: null"
    echo "  gene_annots: null"
    echo "  gene_annots_csv: null"
    echo "  download_script: null"
} >> "$TARGET_FILE"
