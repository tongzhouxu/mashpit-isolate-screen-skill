# Mashpit Isolate Screen

A local, deterministic Codex Skill for screening bacterial isolate assemblies or paired Illumina reads with Mashpit.

## Routing behavior

- With `--organism`: select that organism's local Mashpit database directly.
- Without `--organism`: run local `mlst` 2.35.0 against its bundled PubMLST schemes, map the selected scheme to a supported group, then query only that Mashpit database.
- Supported keys: `salmonella`, `ecoli`, `listeria`, `campylobacter`, and `cronobacter`.

No assembly is uploaded. Mashpit is installed from upstream commit `538d3421302fe6dd129780605b8ff5dedbf4c046c`.

## Build

Start Docker Desktop, then from the directory containing this repository run:

```bash
docker build \
  --tag mashpit-isolate-screen:local \
  --file mashpit-isolate-screen/container/Dockerfile \
  mashpit-isolate-screen
```

## Run with a known organism

```bash
docker run --rm \
  --volume "/absolute/path/to/data:/data:ro" \
  --volume "$HOME/.mashpit/databases:/databases:ro" \
  --volume "/absolute/path/to/results:/results" \
  mashpit-isolate-screen:local \
  /data/sample.fasta \
  --organism salmonella \
  --database-root /databases \
  --output /results/sample-screen
```

Omit `--organism salmonella` to invoke local MLST auto-detection. For paired reads, replace the assembly argument with both paths, for example `/data/sample_R1.fastq.gz /data/sample_R2.fastq.gz`.

The output directory must not already exist. Review `result.json` first and retain `provenance.json` plus the logs.

## Test without biological tools or databases

```bash
PYTHONPYCACHEPREFIX=/tmp/mashpit_pycache \
python3 -m unittest discover -s mashpit-isolate-screen/tests -v
```

These unit tests mock external bioinformatics execution. A real end-to-end test additionally requires one supported assembly and its corresponding Mashpit `.db`, `.sig`, and `database.json` files.
