# NM2 synthetic validation design (locked)

## Purpose
Demonstrate when pooled summaries are valid and when they become misleading under internal heterogeneity.

## Scenarios
1. Coherent homogeneous shift
2. Opposing subgroup shift
3. Temporal phase switch

## Dataset size
- 1000 genes total
- 700 null genes
- 100 coherent-shift genes
- 100 opposing-subgroup genes
- 100 phase-switch genes

## Replicates
- 3 replicates per subgroup or timepoint

## Noise levels
- low
- moderate

## Comparators
- pooled limma
- stratified or time-aware limma
- instability
- variance proxy

## Expected behavior

### Coherent shift
pooled limma: strong signal  
instability: low

### Opposing subgroup shift
pooled limma: attenuated or misleading  
stratified limma: subgroup signal recovered  
instability: high

### Temporal phase switch
pooled limma: blurred average  
time-aware limma: phase signal recovered  
instability: high
