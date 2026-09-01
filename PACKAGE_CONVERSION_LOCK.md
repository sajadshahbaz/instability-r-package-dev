# Instability R Package Conversion Lock

## Purpose

This repository is an isolated development copy created solely for converting
the published Instability implementation into a standard R package.

The package conversion must preserve the scientific behavior of the original
implementation.

## Authoritative Source

Original repository:
https://github.com/sajadshahbaz/NM2_instability_signal

Source branch:
main

Source commit:
6a261c2aa665eb457a94aa1cead2cc8b39be1987

Source tag present at baseline:
v1.1.1

Development repository:
https://github.com/sajadshahbaz/instability-r-package-dev

## Baseline Validation

Baseline smoke test:

    bash run_smoke.sh

Result:

    PASS

Metadata samples: 111
VST samples: 111
Sample overlap: 111
Metadata-only samples: 0
VST-only samples: 0

Baseline core scientific output:

    04_results/tables/core_instability_signal_table.tsv

Baseline SHA-256:

    6991e7024cefac76d15fabbc68e3884adffc3e86267bfce1aa1f4e49c2274d92

## Scientific Behavior Lock

The package conversion MUST NOT silently change:

- algorithms
- mathematical equations
- scientific logic
- thresholds
- default parameter values
- preprocessing semantics
- input semantics
- output semantics
- ranking behavior
- numerical behavior
- interpretation of Instability
- scientific claims

Any proposed change that could alter scientific behavior must be reported
separately and must not be implemented without explicit approval.

## Permitted Engineering Work

The conversion may introduce or restructure engineering components required
for a conventional R package, including:

- DESCRIPTION
- NAMESPACE
- R/
- roxygen2 documentation
- man/
- tests/
- testthat infrastructure
- examples
- vignettes
- dependency declarations
- package build/check configuration
- .Rbuildignore
- package installation workflows

Refactoring is permitted only when scientific behavior is preserved.

## Known Portability Issues at Baseline

The source repository contains historical absolute filesystem paths in some
logs, precomputed outputs, and scripts.

These paths must not be silently interpreted as authorization to change
scientific behavior.

Packaging-related portability corrections must be distinguished from
scientific modifications.

## Acceptance Requirement

The packaged implementation must be independently compared against the
baseline implementation using representative inputs.

Acceptance requires equivalent scientific outputs under explicitly documented
comparison criteria and tolerances.

Package installation, tests, examples, documentation, and R CMD check must
also be independently verified.

## AI Development Boundary

An AI coding system may assist with software engineering and package
construction in this development repository.

It is not authorized to redefine, improve, optimize, simplify, or otherwise
modify the scientific method without explicit human approval.

The original NM2_instability_signal repository is outside the AI development
boundary and must remain unchanged.
