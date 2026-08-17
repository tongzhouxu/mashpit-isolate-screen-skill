# Mashpit Isolate Screen Skill

A local, deterministic **agent skill** (works with Claude Code, OpenAI Codex, or any tool-calling agent that reads `SKILL.md`) for screening a bacterial isolate — raw paired Illumina reads or an assembly — against [Mashpit](https://github.com/tongzhouxu/mashpit), a MinHash-sketch database of NCBI Pathogen Detection SNP clusters. Optionally confirms a Mashpit candidate at true SNP resolution with [ska2](https://github.com/bacpop/ska.rust).

Every command, threshold, and parameter is fixed in version-controlled config (`config/*.json`) — the invoking LLM only runs one script and reports the result; it never constructs bioinformatics commands or invents a cutoff. See [SKILL.md](SKILL.md) for the full agent-facing contract.

Supports `salmonella`, `ecoli_shigella`, `listeria`, `campylobacter`, and `cronobacter`. Nothing is uploaded during a screen — the query stays on your machine. Mashpit is installed from a pinned upstream commit (`538d3421302fe6dd129780605b8ff5dedbf4c046c`), not the older published PyPI release.

## Setup

You need two things: the container image, and a Mashpit database for at least one organism.

**1. Get the image** — either pull the pre-built one:

```bash
docker pull --platform linux/amd64 ghcr.io/tongzhouxu/mashpit-isolate-screen-skill:latest
docker tag ghcr.io/tongzhouxu/mashpit-isolate-screen-skill:latest mashpit-isolate-screen:local
```

or build it yourself:

```bash
docker build --platform linux/amd64 --tag mashpit-isolate-screen:local --file container/Dockerfile .
```

`--platform linux/amd64` is required everywhere here, not just on Apple Silicon: the image is only published/buildable for `linux/amd64` because `quast=5.3.0` has no native `linux/arm64` build compatible with the pinned Python 3.11. It runs fine under emulation on Apple Silicon; omitting the flag there pulls/builds nothing since Docker defaults to your host's native architecture.

**2. Get a database** — download the pre-built ones from [Releases](../../releases/tag/databases-v1):

```bash
mkdir -p ~/.mashpit/databases && cd ~/.mashpit/databases
for org in salmonella ecoli_shigella listeria campylobacter cronobacter; do
  curl -LO "https://github.com/tongzhouxu/mashpit-isolate-screen-skill/releases/download/databases-v1/${org}.tar.gz"
done
curl -LO https://github.com/tongzhouxu/mashpit-isolate-screen-skill/releases/download/databases-v1/checksums.sha256.txt
shasum -a 256 -c checksums.sha256.txt   # verify before extracting
for f in *.tar.gz; do tar -xzf "$f"; done
```

You only need the organism(s) you actually plan to screen against — each is independent. `screen_isolate.py` additionally verifies a per-file checksum recorded in each organism's `database.json` on every run, so a corrupted or tampered database is caught automatically, not just at download time.

Database creation and updating are out of scope for this skill; see [mashpit](https://github.com/tongzhouxu/mashpit) itself for that.

## Run a screen

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

- Input is one assembly (`.fa`/`.fasta`/`.fna`), or two paired-end FASTQ files — replace the assembly argument with both paths (e.g. `/data/sample_R1.fastq.gz /data/sample_R2.fastq.gz`) to run the fixed fastp → SKESA → QUAST assembly workflow first.
- Omit `--organism` to auto-detect it with local `mlst` against its bundled PubMLST schemes instead of asserting it.
- The output directory must not already exist — nothing gets silently overwritten.

Add `--snp-resolve` to also download the relevant representative genomes from NCBI and compute exact pairwise SNP distances with ska2 once a Mashpit candidate is found — a Neighbor-Joining tree, a per-cluster distance summary, and a confidence comparison between the nearest and next-nearest cluster. This is the only step that reaches out to the network (for public reference genomes; the query itself is never uploaded), so it's opt-in. See [references/snp-resolution.md](references/snp-resolution.md).

### Reading the result

Three files land in the output directory:

- **`report.md`** — a plain-language summary for a non-technical reader: organism determination, Mashpit's candidate clusters and scores plus its own Mash-based tree, and (with `--snp-resolve`) the ska2 SNP tables, confidence statement, and SNP tree — both trees rendered as PNGs with the query highlighted.
- **`result.json`** — the structured, authoritative result. Use `status`, `stop_reason`, and `user_summary`.
- **`provenance.json`** — checksums, pinned tool versions, and every command actually run, for reproducibility.

## Use it without writing any code

If you're driving this through an agentic coding tool rather than the CLI directly, no code is required — you just ask in plain language and the agent invokes the pipeline for you, following [SKILL.md](SKILL.md)'s instructions:

- **Claude Code**: drop this repo (or a copy of it) into a skills directory Claude Code scans (project-level `.claude/skills/` or your user-level skills directory) and ask it to screen a genome in chat.
- **OpenAI Codex** (or any Codex-style agent): [agents/openai.yaml](agents/openai.yaml) already declares this skill's default invocation prompt.

What you still can't skip, regardless of which agent triggers it: Docker installed and a database downloaded locally (see Setup above) on whatever machine actually runs the agent's shell commands. A plain chat interface with no shell/Docker access (e.g. the base claude.ai chat, without a code-execution or agentic-coding capability) can't run this — it needs an agent that can actually execute local commands.

## Test without biological tools or databases

```bash
PYTHONPYCACHEPREFIX=/tmp/mashpit_pycache python3 -m unittest discover -s tests -v
```

These unit tests mock external bioinformatics execution and don't need Docker, Mashpit, or a real database. A real end-to-end run additionally needs the container and at least one downloaded database, per Setup above.

## Reference docs

- [SKILL.md](SKILL.md) — the agent-facing contract: when to run what, how to report results
- [references/setup.md](references/setup.md) — container/database contract in full
- [references/workflow.md](references/workflow.md) — the fixed pipeline, stage by stage
- [references/database-routing.md](references/database-routing.md) — organism routing rules
- [references/qc-policy.md](references/qc-policy.md) — QC thresholds and their evidence basis
- [references/mashpit-interpretation.md](references/mashpit-interpretation.md) — how a Mashpit result is labeled
- [references/snp-resolution.md](references/snp-resolution.md) — the optional ska2 SNP-resolution step
- [references/limitations.md](references/limitations.md) — scope and scientific limitations
- [CHANGELOG.md](CHANGELOG.md) — release history
- [CITATION.cff](CITATION.cff) — how to cite this skill and Mashpit itself
