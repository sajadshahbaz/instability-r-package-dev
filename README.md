# instability

`instability` provides the validated ultraviolet (UV) and desiccation (DES)
Instability computation for already-preprocessed variance-stabilized expression
data and associated metadata.

## Installation

Install from a local source checkout:

    R CMD INSTALL /path/to/instability-r-package-dev

Or install a built source tarball:

    R CMD INSTALL instability_0.1.0.tar.gz

## Public API

The package exports one function:

    compute_uv_des_instability(vst, metadata)

`vst` must be a `data.frame` or `data.table` containing already-preprocessed
variance-stabilized expression data. Its first column contains feature
identifiers, and the remaining named columns contain samples.

`metadata` must be a `data.frame` or `data.table` containing these columns:

- `sample_id`
- `condition`
- `baseline_block`
- `group_type`

## Minimal usage

    vst <- data.table::fread("path/to/vst.tsv")
    metadata <- data.table::fread("path/to/metadata.csv")

    result <- instability::compute_uv_des_instability(
      vst = vst,
      metadata = metadata
    )

## Return value

The function returns an ordinary named list containing these eleven
`data.table` objects, in order:

1. `uv_pool_vs_strat_all`
2. `uv_pool_vs_strat_top150`
3. `uv_pool_vs_strat_summary`
4. `des_pool_vs_strat_all`
5. `des_pool_vs_strat_top150`
6. `des_pool_vs_strat_summary`
7. `pool_vs_strat_metrics`
8. `tail_jaccard_metrics`
9. `rank_displacement_pool_vs_strat`
10. `tail_headline_summary`
11. `core_instability_signal_table`

## Scope limitations

Version 0.1.0 does not provide:

- upstream raw-count filtering;
- normalization;
- VST generation;
- Gamma analysis;
- homogeneous-null controls;
- figure-generation workflows;
- annotation or NCBI workflows;
- synthetic benchmark workflows;
- archived or historical manuscript workflows.

The UV/DES scientific semantics are intentionally preserved from the validated
source implementation.

## Provenance

- Authoritative scientific source repository:
  https://github.com/sajadshahbaz/NM2_instability_signal
- Authoritative scientific source commit:
  `6a261c2aa665eb457a94aa1cead2cc8b39be1987`
- Baseline tag: `v1.1.1`
- Package development repository:
  https://github.com/sajadshahbaz/instability-r-package-dev
- Scientific-source repository provenance DOI:
  https://doi.org/10.5281/zenodo.21110170

The Zenodo DOI above identifies provenance for the scientific-source
repository. It is not a DOI for package version 0.1.0.

## Validation

The validated package implementation reproduced all eleven frozen UV/DES
scientific reference tables byte-for-byte and passed exact parsed and semantic
comparison. This validation does not extend to conditions or workflows outside
the package scope described above.
