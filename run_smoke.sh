#!/usr/bin/env bash
set -euo pipefail

Rscript 02_methods/R/core_pipeline/00_io.R
Rscript 02_methods/R/core_pipeline/04_make_core_table.R

test -f 04_results/logs/00_io_summary.txt
test -f 04_results/tables/core_instability_signal_table.tsv

echo "[OK] smoke run complete"
