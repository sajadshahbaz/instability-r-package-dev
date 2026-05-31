#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")/.."

Rscript NM2_instability_signal/02_methods/R/00_io.R
Rscript NM2_instability_signal/02_methods/R/01_pool_vs_strat.R
Rscript NM2_instability_signal/02_methods/R/02_plots_pool_vs_strat.R
Rscript NM2_instability_signal/02_methods/R/03_tail_divergence.R
Rscript NM2_instability_signal/02_methods/R/04_make_core_table.R

md5sum NM2_instability_signal/04_results/tables/*.tsv \
  > NM2_instability_signal/04_results/logs/md5_tables.txt

echo "[OK] smoke run complete"
