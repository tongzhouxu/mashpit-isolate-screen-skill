# QC policy 1.0.0

The executable policy is `config/qc-policy.json`. Changes require a policy-version increment and validation against representative isolates.

## Evidence basis

NCBI Pathogen Detection documents a 20× minimum read coverage, a 10,000 bp minimum assembly N50, and a maximum of 500 contigs for most organisms, with an exception for Escherichia/Shigella. It also applies species-specific genome-length and contamination checks. See [NCBI Pathogen Detection data processing](https://www.ncbi.nlm.nih.gov/pathogens/docs/data_processing/).

The organism length ranges for Salmonella, Listeria, Campylobacter, and E. coli use the ranges reported in the 2025 inter-laboratory foodborne-pathogen workflow comparison. Salmonella's range is also consistent with the 4.3–5.3 Mb distribution reported by SISTR. Cronobacter uses a deliberately broad 4–5 Mb screening envelope and requires empirical validation before production use. See [inter-laboratory study](https://doi.org/10.3389/fmicb.2025.1629731) and [SISTR assessment](https://doi.org/10.1371/journal.pone.0192573).

## State rules

`FAIL` stops Mashpit. It includes malformed input, N50 below 10 kb, excessive organism-specific contigs, length outside the organism range, more than 5% ambiguous bases, or estimated read coverage below 20×.

`WARN` proceeds with a visible warning. It includes N50 from 10–20 kb, more than 300 contigs while still below the fail ceiling, 1–5% ambiguous bases, or post-fastp Q30 fraction below 0.80.

`PASS` means no implemented fail or warning rule fired. It does not guarantee purity, correctness, or fitness for epidemiological inference.

The broad pre-classification 1–10 Mb length envelope only catches grossly invalid bacterial-isolate inputs. Final length assessment happens after user-directed or local-MLST routing.

Contamination is not claimed to be measured in MVP 1.0. MLST statuses `MIXED` and `BAD` stop automatic routing, but this is not a validated mixed-isolate contamination assay.
