#!/usr/bin/env bash
set -euo pipefail

ROOT="NM2_instability_signal"

# Refuse to overwrite unless user explicitly deletes it
if [[ -e "$ROOT" ]]; then
  echo "[ERROR] '$ROOT' already exists. Delete it or choose a different location."
  exit 1
fi

# Core structure
mkdir -p \
  "$ROOT/00_config" \
  "$ROOT/01_data/raw" \
  "$ROOT/01_data/derived" \
  "$ROOT/02_methods/R" \
  "$ROOT/02_methods/bash" \
  "$ROOT/03_runs" \
  "$ROOT/04_results/figures" \
  "$ROOT/04_results/tables" \
  "$ROOT/04_results/logs" \
  "$ROOT/08_theory" \
  "$ROOT/tmp"

# Keep empty dirs tracked (Git-friendly)
touch \
  "$ROOT/01_data/raw/.gitkeep" \
  "$ROOT/01_data/derived/.gitkeep" \
  "$ROOT/03_runs/.gitkeep" \
  "$ROOT/04_results/figures/.gitkeep" \
  "$ROOT/04_results/tables/.gitkeep" \
  "$ROOT/04_results/logs/.gitkeep" \
  "$ROOT/tmp/.gitkeep"

# Minimal config (edit paths later)
cat > "$ROOT/00_config/config.yaml" <<'YAML'
# NM2_instability_signal configuration
# Fill in absolute paths or keep data inside 01_data/.

paths:
  vst: "01_data/raw/vst_counts_filtered.tsv"
  counts: "01_data/raw/counts.tsv"
  meta: "01_data/raw/metadata.csv"

analysis:
  condition_col: "condition"
  subgroup_col: "baseline_block"
  group_col: "group_type"      # e.g., treated/control if you use it
  sample_id_col: "sample_id"   # change if your metadata uses a different name

params:
  seed: 1
  boot: 500
  topn: 150

labels:
  uv:
    short: "short"
    long:  "long"
  des:
    short: "fast"
    long:  "slow"
YAML

# Real-data checklist (reviewer-proof discipline)
cat > "$ROOT/REALDATA_CHECKLIST.md" <<'MD'
# REALDATA_CHECKLIST (NM2_instability_signal)

This checklist is a guardrail against self-inflicted chaos.

## Data integrity
- [ ] `vst_counts_filtered.tsv` is present and readable.
- [ ] `counts.tsv` is present and readable (optional for some analyses).
- [ ] `metadata.csv` is present and readable.
- [ ] Sample IDs match between metadata and VST (intersection > 0, no duplicates).
- [ ] Condition labels are normalized (lowercase, trimmed).

## Subgroup sanity
- [ ] Subgroup column exists (e.g., `baseline_block`).
- [ ] For each condition, subgroup levels exist (e.g., UV: short/long; DES: fast/slow).
- [ ] Per-subgroup sample size is reported (n per subgroup, treated/control if applicable).
- [ ] Any subgroup with n < 3 is flagged and handled (weights / fallback / excluded).

## Reproducibility
- [ ] Seed is fixed and recorded.
- [ ] Bootstrap count is fixed and recorded.
- [ ] Output folders include timestamp or are explicitly overwritten.
- [ ] Every run writes:
  - [ ] `run_summary.txt`
  - [ ] key tables (`instability_scores.tsv`, `pool_vs_strat_summary.tsv`)
  - [ ] key figures (`I1_vs_variance.png`, etc.)

## Artifacts
- [ ] UV run produces all tables + figures.
- [ ] DES run produces all tables + figures.
- [ ] Figures include labels + axes + legends and are publication-ready.
MD

# Stub theory docs
cat > "$ROOT/08_theory/definitions.md" <<'MD'
# Definitions (Draft)

This file will contain formal definitions of stability and instability functionals,
and how they differ from variance-based dispersion.
MD

cat > "$ROOT/08_theory/propositions.md" <<'MD'
# Propositions (Draft)

This file will contain short, testable propositions about instability as signal,
with empirical tests defined in 04_results outputs.
MD

# Convenience runner stub (we’ll fill it later)
cat > "$ROOT/run_all.sh" <<'SH2'
#!/usr/bin/env bash
set -euo pipefail
echo "Stub: run_all.sh will orchestrate UV/DES runs once methods scripts exist."
SH2
chmod +x "$ROOT/run_all.sh"

echo "[OK] Created project tree at: $ROOT"
echo "Next: copy your files into $ROOT/01_data/raw/ and update 00_config/config.yaml."
