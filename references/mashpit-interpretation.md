# Mashpit result interpretation

The parser sorts candidate rows by numeric score and applies only structural labels:

| Rule | Label |
|---|---|
| no candidates | No candidate |
| Mashpit `near_top` marks another candidate cluster | Ambiguous top candidates |
| one top candidate | Top-ranked candidate |

Mashpit derives `near_top` from its fixed query profile. The Skill reports the raw score and does not label it strong, moderate, or weak because no validated biological cutoff has been supplied. The query profile's `--threshold 0.85` controls Mashpit's local tree construction and is not presented as an outbreak-relatedness cutoff.

Report the best candidate cluster, score, alternatives, ambiguity, QC, database name/version, and warnings. Always state that Mashpit similarity is a screening result and does not establish an epidemiological or outbreak relationship. Recommend a validated organism-appropriate high-resolution SNP workflow and review of epidemiological metadata for confirmation.

If the candidate file is absent, malformed, ambiguous beyond the fixed rule, or inconsistent with its schema, do not reconstruct an answer from terminal logs.
