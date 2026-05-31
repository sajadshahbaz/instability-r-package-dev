# Raw Reference Genome Data

## Description

This directory contains information and reference resources associated with the genome assembly used for annotation and functional characterization analyses in the NM2 instability-signal study.

Reference assembly:

```text
Ramazzottius varieornatus
Assembly accession: GCA_001949185.1
Source: NCBI GenBank
```

The original assembly, annotations, and associated metadata are publicly available through NCBI.

---

## Data Source

Assembly accession:

```text
GCA_001949185.1
```

Organism:

```text
Ramazzottius varieornatus
```

Assembly provider:

```text
NCBI GenBank
```

Original assembly page:

```text
https://www.ncbi.nlm.nih.gov/assembly/GCA_001949185.1
```

---

## Purpose Within This Repository

The NM2 repository does not reproduce genome assembly, genome annotation, or RNA-seq alignment from raw sequencing reads.

Instead, processed expression matrices and curated metadata are used for the instability analyses described in the manuscript.

Reference annotation resources derived from assembly GCA_001949185.1 are used for:

* gene mapping
* annotation matching
* functional characterization
* Figure 5 annotation analyses
* supplementary annotation resources

---

## Files Included in This Repository

The following reference files are included because they are required for reproducing annotation-related analyses and associated repository outputs:

```text
cds_from_genomic.fna
protein.faa
genomic.gff
genomic.gtf
```

These files are sufficient for reproducing the Figure 5 annotation workflow and associated outputs included in this repository.

---

## Files Omitted From GitHub

The following large files were used during exploratory annotation development but are not required for reproducing the final NM2 analyses, figures, supplementary outputs, or repository audit results:

```text
GCA_001949185.1_Rvar_4.0_genomic.fna
genomic.gbff
```

These files can be obtained directly from the original NCBI assembly distribution.

They may also be preserved in the archival Zenodo release associated with this repository.

Removal of these files does not affect reproduction of:

* manuscript figures
* supplementary figures
* supplementary tables
* annotation outputs
* repository audit results
* instability-versus-variance analyses

---

## Reproducibility

Users may reproduce repository outputs by:

1. Using the processed datasets distributed with this repository.
2. Using the included reference annotation resources listed above.
3. Downloading the complete assembly package directly from NCBI if additional genome-level analyses are desired.
4. Using the archival Zenodo release accompanying this repository.

Repository-wide reproducibility instructions are available in:

```text
README.md
REPOSITORY_MAP.md
docs/repository_structure.md
docs/reproducibility.md
```

---

## Citation

Hashimoto T, Horikawa DD, Saito Y, Kuwahara H, Kozuka-Hata H, Shin-I T, et al.

Genome sequencing and analysis of the tardigrade Ramazzottius varieornatus.

Nature Communications. 2016;7:12808.

---

## Notes

This repository is intended to reproduce the NM2 instability-signal analyses and manuscript outputs.

The files distributed here are not intended to replace the complete NCBI assembly package and should be considered a minimal reference resource necessary for reproducing repository analyses.

