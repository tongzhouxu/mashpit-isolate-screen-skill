# Setup

Run analyses in the pinned container built from `container/Dockerfile`. Mashpit is installed from upstream commit `538d3421302fe6dd129780605b8ff5dedbf4c046c`, not from the older published PyPI artifact. Runtime package installation is not part of an analysis.

Pull the published image:

```bash
docker pull ghcr.io/tongzhouxu/mashpit-isolate-screen-skill:latest
docker tag ghcr.io/tongzhouxu/mashpit-isolate-screen-skill:latest mashpit-isolate-screen:local
```

or build it locally from the repository root:

```bash
docker build --platform linux/amd64 --tag mashpit-isolate-screen:local --file container/Dockerfile .
```

`--platform linux/amd64` is required on Apple Silicon (`quast=5.3.0` has no `linux/arm64` build for the pinned Python 3.11); it runs fine under emulation there.

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

Pre-built, checksummed databases for all five are published at [Releases](../../../releases/tag/databases-v1) — download and verify before extracting:

```bash
mkdir -p ~/.mashpit/databases && cd ~/.mashpit/databases
for org in salmonella ecoli_shigella listeria campylobacter cronobacter; do
  curl -LO "https://github.com/tongzhouxu/mashpit-isolate-screen-skill/releases/download/databases-v1/${org}.tar.gz"
done
curl -LO https://github.com/tongzhouxu/mashpit-isolate-screen-skill/releases/download/databases-v1/checksums.sha256.txt
shasum -a 256 -c checksums.sha256.txt
for f in *.tar.gz; do tar -xzf "$f"; done
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

Only the selected organism database is required for a run — each is independent, so download just the one(s) you need. Local `mlst` 2.35.0 and its bundled PubMLST schemes are installed in the image; no scheme download or sequence upload occurs at runtime.

The release process must replace `SET_AT_RELEASE` in `config/workflow.json` with the published container digest. Treat that placeholder as a setup failure for regulated or production use.
