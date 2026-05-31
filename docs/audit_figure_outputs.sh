#!/usr/bin/env bash
set -euo pipefail

echo "=== Figure script audit ==="
echo "Repository: $(pwd)"
echo

for s in 02_methods/R/figure_generation/*.R; do
  echo "---- $s ----"
  echo "Declared output paths:"
  grep -nE "ggsave|png\\(|pdf\\(|tiff\\(|svg\\(|write.*figure|OUT_|outdir|out_dir|file.path" "$s" || true
  echo
done

echo "=== Existing generated image files in 04_results ==="
find 04_results -type f \( -iname "*.png" -o -iname "*.pdf" -o -iname "*.tif" -o -iname "*.tiff" -o -iname "*.svg" \) -print | sort

echo
echo "=== Existing manuscript figures ==="
find 03_figures/main -type f -print | sort
