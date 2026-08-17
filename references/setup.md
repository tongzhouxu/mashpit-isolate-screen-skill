# Setup

Run analyses in the pinned container built from `container/Dockerfile`. Mashpit is installed from upstream commit `538d3421302fe6dd129780605b8ff5dedbf4c046c`, not from the older published PyPI artifact. Runtime package installation is not part of an analysis.

Start Docker Desktop, then build from the repository's parent directory:

```bash
docker build --tag mashpit-isolate-screen:local \
  --file mashpit-isolate-screen/container/Dockerfile mashpit-isolate-screen
```

Run an assembly with a known organism:

```bash
docker run --rm \
  --volume "/absolute/path/to/data:/data:ro" \
  --volume "$HOME/.mashpit/databases:/databases:ro" \
  --volume "/absolute/path/to/results:/results" \
  mashpit-isolate-screen:local \
  /data/sample.fasta --organism salmonella \
  --database-root /databases --output /results/sample-screen
```

Omit `--organism` to use local MLST routing. The output directory must not already exist.

Set `MASHPIT_DATABASE_ROOT` or pass `--database-root`. The directory must contain:

```text
databases/
├── salmonella/
├── ecoli_shigella/
├── listeria/
├── campylobacter/
└── cronobacter/
```

Each organism directory must contain the Mashpit `<name>.db` and `<name>.sig` files plus `database.json`:

```json
{
  "name": "salmonella",
  "version": "PDG-build-identifier",
  "build_date": "2026-08-01",
  "checksum": "sha256-of-primary-artifact",
  "checksum_file": "salmonella.db",
  "source": "NCBI Pathogen Detection via Mashpit build"
}
```

The wrapper checks the name, required metadata, required database files, and checksum target. It never downloads, builds, updates, or substitutes a database during a screen.

Only the selected organism database is required for a run. The directory names are placeholders until the user supplies the actual local Mashpit databases and metadata. Local `mlst` 2.35.0 and its bundled PubMLST schemes are installed in the image; no scheme download or sequence upload occurs at runtime.

The release process must replace `SET_AT_RELEASE` in `config/workflow.json` with the published container digest. Treat that placeholder as a setup failure for regulated or production use.
