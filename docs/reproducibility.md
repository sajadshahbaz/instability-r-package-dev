

## Manuscript reproducibility checkpoint

The submitted manuscript workflow is reproduced with:

```bash
Rscript 02_methods/R/run_nm2_pipeline.R
```

The master runner executes each analysis script in a separate `Rscript --vanilla` process to avoid shared-session side effects from open graphics devices, file connections, or global objects.

The final manuscript reproducibility run completed successfully with exit code 0 and is logged in:

```text
04_results/logs/final_manuscript_reproducibility_run.log
```

A frozen copy of this log is stored in:

```text
docs/checkpoints/final_manuscript_reproducibility_run.log
```

Figure 5e requires supplied precomputed annotation resources:

```text
04_results/figure5e_global_functional_structure/annotation_from_scratch/figure5e_pfam.domtblout
04_results/figure5e_global_functional_structure/annotation_from_scratch/figure5e_interproscan.tsv
```

These files are included as annotation resources. The manuscript pipeline does not rerun Pfam/HMMER or InterProScan.

The required NCBI-derived files used by the manuscript scripts are:

```text
01_data/raw/ncbi_dataset/ncbi_dataset/data/GCA_001949185.1/cds_from_genomic.fna
01_data/raw/ncbi_dataset/ncbi_dataset/data/GCA_001949185.1/protein.faa
```

The full genomic FASTA and GenBank files are not required by the manuscript reproduction pipeline.
