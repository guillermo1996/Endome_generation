#!/bin/bash

rm -rf test/*

mmseqs createdb ORFs_fasta/bin_1.fa test/bin.db
mmseqs search test/bin.db test/bin.db test/bin_res test/tmp
mmseqs convertalis test/bin.db test/bin.db test/bin_res test/results.tsv --format-output "query,target,pident,raw,evalue"
