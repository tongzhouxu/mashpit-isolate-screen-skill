# Organism and database routing

If the user supplies `--organism`, use that supported key directly. This is an explicit user assertion; the workflow does not independently verify it.

If the option is omitted, run the pinned local classifier:

```text
mlst --full --csv --quiet ASSEMBLY
```

`mlst` auto-detects a scheme from the PubMLST schemes bundled with version 2.35.0. `config/mlst-routing.json` maps known schemes to the five supported Mashpit database groups. The assembly remains local; the workflow does not use the PubMLST REST service or update schemes during a run.

Use MLST's native status rather than inventing a numeric threshold:

| MLST status | Routing action |
|---|---|
| `PERFECT`, `NOVEL`, `OK` | proceed |
| `MISSING` | proceed with a prominent warning |
| `NONE`, `MIXED`, `BAD` | stop as uncertain |
| a valid but unmapped scheme | stop as unsupported |

After routing, query only the selected Mashpit database. A missing selected database is a setup error; never query or substitute an unrelated database.
