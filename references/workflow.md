# Fixed workflow

`scripts/screen_isolate.py` owns routing; subordinate scripts own computation.

Assembly input:

```text
validate FASTA → broad assembly QC → user organism or local MLST scheme detection
→ select one supported database → organism-specific QC → run Mashpit once
→ structured Mashpit parser → provenance
```

Paired Illumina input adds:

```text
validate complete R1/R2 FASTQ → fastp → SKESA → QUAST → assembly-input path
```

Commands and parameters come only from `config/workflow.json`. Mashpit is pinned to upstream commit `538d3421302fe6dd129780605b8ff5dedbf4c046c`. The fixed profile is:

```text
mashpit query ASSEMBLY DATABASE --number 200 --threshold 0.85 --tie-tolerance-hashes 2
```

Mashpit writes its selected-database query within the retained `mashpit/` directory. The parser consumes the pinned commit's native `<sample>_representative_matches.csv` and `<sample>_cluster_candidates.csv` files. Missing, multiple, or malformed result files are errors.

The assembly component returns:

```json
{
  "status": "PASS",
  "assembly_path": "assembly.fasta",
  "read_qc": {},
  "commands": [],
  "workflow_version": "1.1.0"
}
```

Any future replacement assembly Skill must satisfy this contract and must not weaken input validation or QC stop rules.

`--organism` directly selects a supported database when the organism is already known. Without it, local `mlst` 2.35.0 auto-detects a bundled PubMLST scheme. The classifier's complete CSV, stderr log, normalized routing JSON, and command are retained.
