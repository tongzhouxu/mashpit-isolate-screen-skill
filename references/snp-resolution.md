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
2. For each relevant cluster, read `<sample>_representative_matches.csv`, sort its representatives by `similarity_score`, and keep the top `max_representatives_per_cluster` (default 100).
3. Round-robin across relevant clusters (so each keeps at least one representative) and cap the total at `max_total_genomes` (default 100). This is what keeps a cluster with hundreds or thousands of representatives (real example: a Listeria cluster with 714 representatives locally, 2,378 isolates in NCBI's full record) from turning into a 1000+-genome ska2 run — confirmed against that exact cluster: 100 representatives selected, not 714, for 5,050 pairwise comparisons instead of ~255,000.

## Genome retrieval and SNP distance

Representative genomes are downloaded with `datasets download genome accession --include genome --dehydrated` followed by `datasets rehydrate` — one batched call for every selected accession, not one request per genome — retried up to `download_attempts` times for genomes that fail. The query assembly plus every successfully downloaded representative are built into one ska2 split-kmer file (`ska build -f <name-path list> -k <kmer_size>`, pinned k-mer size in `workflow.json`), then compared with `ska distance`, which reports the number of SNPs differing between every pair — the *entire* pairwise matrix (every representative against every other, not just against the query).

Measured at the current 100/100 caps against the 714-representative Listeria cluster above: ~42s total (versus ~14s at the old 5/20 caps), ~299MB downloaded (versus ~15MB). The download itself is cheap regardless of count (one batched request); the added time is mostly `ska build`/`ska distance` doing proportionally more real work.

## Interpretation

`interpret_snp_resolution.py` builds four views of the same pairwise data:

- **`ranked`**: every compared representative, sorted by SNP distance to the query, with its cluster.
- **`cluster_summary`**: per relevant cluster, the count of genomes compared and the min/median/max SNP distance to the query — the "which cluster is really closest" table.
- **`confidence`**: the nearest cluster's minimum distance versus the next-nearest cluster's minimum distance, and their ratio. No strength label ("high/medium/low confidence") is invented — only the raw distances and ratio, consistent with [mashpit-interpretation.md](mashpit-interpretation.md)'s policy against unvalidated cutoffs. When only one cluster was compared, there is no second distance to compare against and the statement says so.
- **`newick_tree`**: a Neighbor-Joining tree (pure-Python Saitou-Nei implementation in `build_snp_tree.py`, no new pinned tool) over every genome compared — query plus every downloaded representative — built from the full pairwise matrix, not just the query's row. Verified against a hand-computed additive 4-taxon case, where NJ is guaranteed to recover the exact tree.

It also states whether the SNP-nearest genome's cluster agrees with Mashpit's top candidate; a disagreement is surfaced as a warning, not silently resolved one way or the other.

`generate_report.py` renders all of this into `report.md`, written alongside `result.json` on every run (not just `--snp-resolve` ones): organism determination, the Mashpit candidate table, the ska2 cluster/genome tables and confidence statement, and the Newick tree, in plain language for a reader who doesn't want to open JSON.

## Limitations

- ska2 pairwise SNP distance is itself a screening refinement, not a validated outbreak-confirmation pipeline (unlike, e.g., CFSAN SNP Pipeline or an accredited cgMLST scheme). Always recommend a validated, organism-appropriate high-resolution comparison for actual outbreak confirmation, exactly as for the Mash result.
- Only representative genomes that Mashpit already selected during database build are compared — this step does not search NCBI independently, and its resolution is bounded by however sparse or dense the underlying Mashpit database's representative set is for that cluster.
- A failed or partial download (fewer than two usable genomes, including the query) skips SNP resolution with a warning; it does not fail the underlying Mash screen, which remains authoritative on its own.
