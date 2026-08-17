# Mashpit result interpretation

The parser sorts candidate rows by numeric score and applies only structural labels:

| Rule | Label |
|---|---|
| no candidates | No candidate |
| Mashpit `near_top` marks another candidate cluster | Ambiguous top candidates |
| one top candidate | Top-ranked candidate |

Mashpit derives `near_top` from its fixed query profile. The Skill reports the raw score and does not label it strong, moderate, or weak because no validated biological cutoff has been supplied. The query profile's `--threshold 0.85` controls Mashpit's local tree construction, not which rows appear in `cluster_candidates.csv` — Mashpit's `--number`/`--tie-tolerance-hashes` group and return the top hits regardless of how low their absolute score is. A query from a genuinely unrelated organism can therefore still produce a `near_top`-flagged "candidate" at a near-zero score (observed directly: a Salmonella genome queried against the Cronobacter database returned a `Top-ranked`/`Ambiguous` result at score 0.007, with Mashpit's own log noting "Top query similarity is smaller than the threshold"). The parser surfaces this as a separate `below_threshold` flag (best score under the same `--threshold`) with its own warning — independent of, and in addition to, the `near_top` ambiguity check — instead of leaving that Mashpit-native signal only in the log.

Report the best candidate cluster, score, alternatives, ambiguity, QC, database name/version, and warnings. Always state that Mashpit similarity is a screening result and does not establish an epidemiological or outbreak relationship. Recommend a validated organism-appropriate high-resolution SNP workflow and review of epidemiological metadata for confirmation.

If the candidate file is absent, malformed, ambiguous beyond the fixed rule, or inconsistent with its schema, do not reconstruct an answer from terminal logs.
