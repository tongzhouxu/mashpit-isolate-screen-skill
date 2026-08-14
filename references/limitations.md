# Limitations

MVP 1.0 supports one bacterial isolate represented by paired-end Illumina FASTQ/FASTQ.GZ or an existing FASTA/FNA assembly. It does not support ONT, PacBio, hybrid, interleaved, single-end, or metagenomic input.

The workflow does not confirm outbreak membership, transmission, source attribution, clinical significance, or treatment decisions. Mashpit MinHash similarity can prioritize cluster representatives but does not replace validated SNP/cgMLST analysis and epidemiological investigation.

Local MLST auto-detection and the QC policy require validation with positive controls, unsupported near neighbors, contaminated isolates, and site-specific sequencing distributions. MLST scheme selection is not a general taxonomic classifier and can misassign close relatives or poor/mixed assemblies. A user-supplied organism bypasses this check. Cronobacter QC bounds especially require empirical calibration. The workflow retains warnings rather than masking these limitations.

The bundled `ecoli` MLST scheme is shared by Escherichia and Shigella in the upstream scheme metadata, so automatic routing cannot reliably distinguish those genera. Treat an automatically selected E. coli database as a routing choice, not definitive species identification.

Database creation and updating are outside scope. Results are only reproducible when the exact Mashpit databases, checksums, pinned source commit, container digest, and provenance are retained.

Optional `--snp-resolve` (see [references/snp-resolution.md](snp-resolution.md)) refines a Mashpit candidate with ska2 pairwise SNP distances against representative genomes re-downloaded from NCBI. It is a screening refinement bounded by whatever representatives the underlying Mashpit database happened to select, not a validated outbreak-confirmation pipeline, and it introduces a network dependency the rest of this workflow does not otherwise have.
