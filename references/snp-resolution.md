# SNP resolution 1.0.0

Optional, opt-in refinement of a Mashpit candidate using ska2 pairwise SNP distances. Governed by `config/snp-resolution-policy.json` (selection thresholds) and the `snp_resolution` block of `config/workflow.json` (fixed ska2 command profile). Enable with `--snp-resolve`.

## Why this exists

Mashpit's MinHash similarity is a coarse, sketch-resolution screen. A single unambiguous top cluster at Mash resolution is not proof of SNP-level closeness, and two near-tied clusters cannot be told apart by Mash alone. SNP resolution runs ska2 (split k-mer analysis) between the query assembly and representative genomes from the relevant cluster(s) to get an actual pairwise SNP count.

## When it runs

Whenever `--snp-resolve` is set and Mashpit returned a candidate (`mashpit_result.best_candidate` is present) — regardless of whether that candidate was unambiguous. This is deliberate: an unambiguous Mash hit can still be many SNPs away from its representatives, so the confirmation step is not limited to the ambiguous case.

## Network dependency

Unlike the rest of the screen, this step is **not** fully local. A Mashpit database only retains sourmash signatures, not the representative assemblies (they are sketched and discarded during `mashpit build`), so resolving SNPs requires re-downloading the relevant representative genomes from NCBI via the pinned `datasets` CLI. This is why the step is opt-in rather than automatic: default screens keep sensitive query data fully local, and `--snp-resolve` is an explicit choice to reach out to NCBI for public reference genomes (the query sequence itself is never uploaded).

## Target selection

1. Read the raw `<sample>_cluster_candidates.csv` rows (not just the summarized `mashpit_result`). The relevant cluster set is the top-ranked cluster plus any cluster Mashpit flagged `near_top`.
2. For each relevant cluster, read `<sample>_representative_matches.csv`, sort its representatives by `similarity_score`, and keep the top `max_representatives_per_cluster`.
3. Round-robin across relevant clusters (so each keeps at least one representative) and cap the total at `max_total_genomes`.

## Genome retrieval and SNP distance

Representative genomes are downloaded with `datasets download genome accession --include genome --dehydrated` followed by `datasets rehydrate`, retried up to `download_attempts` times for genomes that fail. The query assembly plus every successfully downloaded representative are built into one ska2 split-kmer file (`ska build -f <name-path list> -k <kmer_size>`, pinned k-mer size in `workflow.json`), then compared with `ska distance`, which reports the number of SNPs differing between every pair.

## Interpretation

The parser ranks representatives by SNP distance to the query and reports the nearest one, its cluster, and the raw SNP count — no invented "close/distant" cutoff, matching [mashpit-interpretation.md](mashpit-interpretation.md)'s policy of reporting structure over invented labels. It also states whether the SNP-nearest genome's cluster agrees with Mashpit's top candidate; a disagreement is surfaced as a warning, not silently resolved one way or the other.

## Limitations

- ska2 pairwise SNP distance is itself a screening refinement, not a validated outbreak-confirmation pipeline (unlike, e.g., CFSAN SNP Pipeline or an accredited cgMLST scheme). Always recommend a validated, organism-appropriate high-resolution comparison for actual outbreak confirmation, exactly as for the Mash result.
- Only representative genomes that Mashpit already selected during database build are compared — this step does not search NCBI independently, and its resolution is bounded by however sparse or dense the underlying Mashpit database's representative set is for that cluster.
- A failed or partial download (fewer than two usable genomes, including the query) skips SNP resolution with a warning; it does not fail the underlying Mash screen, which remains authoritative on its own.
