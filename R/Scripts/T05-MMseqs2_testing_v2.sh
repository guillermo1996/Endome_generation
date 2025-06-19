#!/bin/bash

mkdir -p mmseqs_tmp
mkdir -p mmseqs_results

echo "transcript_id,num_pairs,mean_identity" > transcript_similarity_scores.csv

for i in $(seq 1 15); do
    fasta="ORFs_fasta/bin_${i}.fa"
    if [[ -f "$fasta" ]]; then
        base=$(basename "$fasta" .fa)

        # Create database
        mmseqs createdb "$fasta" mmseqs_tmp/${base}_db

        # All-vs-all search
        mmseqs search mmseqs_tmp/${base}_db mmseqs_tmp/${base}_db mmseqs_tmp/${base}_res mmseqs_tmp/tmp_${base} --max-seqs 3000

        # Convert to readable format
        mmseqs convertalis mmseqs_tmp/${base}_db mmseqs_tmp/${base}_db mmseqs_tmp/${base}_res mmseqs_results/${base}_results.tsv --format-output "query,target,pident"

        # Filter and calculate mean identity
        awk -v tx=$base '
            $1 != $2 { sum += $3; count++ }
            END {
                mean = (count > 0) ? sum / count : 1;
                print tx "," count "," mean
            }' mmseqs_results/${base}_results.tsv >> transcript_similarity_scores.csv
    else
        echo "Skipping missing file: $fasta"
    fi
done
