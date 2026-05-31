# Repository Structure Policy

## Canonical Repository Root

```text
NM2_instability_signal/
```

All scripts, figures, tables, logs, and manuscript-facing outputs should be generated and stored within the canonical repository structure.

## Archive Material

Legacy scripts, exploratory analyses, intermediate outputs, and development-stage artifacts may be retained in archive directories for transparency and historical reference.

Archive materials are not part of the minimal reproducible workflow and should not be used to regenerate manuscript results unless explicitly stated.

## Repository Rules

* All manuscript-facing analyses must be generated from the canonical repository workflow.
* Final figures, tables, and supplementary outputs should be written to their designated repository locations.
* Historical and exploratory files should remain separated from the reproducible analysis pipeline.
* Any changes affecting manuscript results should be documented and traceable through repository scripts and logs.

