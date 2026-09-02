# Original UV/DES equivalence oracle

This directory contains the pre-package oracle for the unchanged original
UV/DES Instability implementation. It is deliberately limited to the formal
v1 inputs (preprocessed VST expression and metadata) and the four authoritative
UV/DES core scripts.

The oracle does not implement scientific calculations. It runs the original
scripts unchanged in a newly created `/tmp` workspace and compares every
scientific TSV with the frozen repository references.

## Authoritative scripts

The runner invokes, in order:

1. `02_methods/R/core_pipeline/01_pool_vs_strat_rebuild_uv_des.R`
2. `02_methods/R/core_pipeline/02_plots_pool_vs_strat.R`
3. `02_methods/R/core_pipeline/03_tail_divergence.R`
4. `02_methods/R/core_pipeline/04_make_core_table.R`

Their hashes, the formal input hashes, and the frozen scientific-reference
hashes are locked in `uv_des_manifest.sha256`.

## Safety model

The runner:

- creates its workspace with `mktemp` below `/tmp`;
- rejects broad, repository-local, or pre-existing output workspaces;
- copies only `config.yaml`, metadata, and the VST matrix;
- sets both the working directory and `NM2_ROOT` to the workspace;
- invokes authoritative scripts by absolute repository path;
- captures Git HEAD, exact `git status --short`, and protected hashes before
  execution;
- repeats all repository protection checks on exit, including failed exits;
- preserves the workspace for diagnosis;
- never writes generated results to repository `04_results`.

The repository working tree need not be clean. Its post-run short status must
be byte-for-byte identical to its pre-run short status.

## Comparisons

Layer 1 requires byte-identical SHA-256 values for all scoped scientific TSVs.

Layer 2 parses every table and requires exact equality of dimensions, schema,
column and row order, types, identifiers, labels, logical values, special-value
masks, and numeric values. It also checks Top-150 order, Top-K membership/order,
ranks, K order, and final condition order. No numerical tolerance is present.

Logs and plots are checked only for presence and non-zero size.

## Requirements

- Bash, Git, `sha256sum`, `mktemp`, and `Rscript`
- R packages `data.table` and `ggplot2`
- locale `en_US.UTF-8`
- no package installation or environment mutation is performed

## Controlled execution

Only run after explicit authorization:

```bash
bash oracle/run_original_uv_des_oracle.sh
```

Success and failure both print the exact preserved workspace. Detailed console,
environment, checksum, and parsed-comparison records remain in that workspace.
Any upstream failure prevents subsequent scientific stages.
