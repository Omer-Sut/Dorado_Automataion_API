# Real Integration Test Fixtures

These fixtures contain sanitized but real human Oxford Nanopore sequencing data. The lab authorized
their public redistribution on 2026-08-05 after sample/run identifiers were replaced with anonymous
test labels. Replacing metadata and identifiers does not change the retained genomic sequences.

The FASTQ and BAM are from unrelated runs and test different workflow stages. They must not be
compared with each other or treated as a matched sample:

- `fastq/barcode19/test_reads_barcode19.fastq.gz` tests NanoTel and NanoTel filtration.
- `bam/barcode07/test_aligned_barcode07.bam` tests an already-aligned BAM, alignment-summary
  generation, mapping, and BAM filtering. Its `.bai` is regenerated during preparation.

The preparation script replaces read and run identifiers and removes timestamps, user and device
names, filesystem paths, POD5 filenames, GPU details, and command lines. It preserves FASTQ
sequences and qualities and the BAM's sequence, quality, alignment, reference, and selected
alignment tags. See `manifest.json` for checksums and the exact policy.

These files do not test basecalling, demultiplexing, or remapping. The BAM has no `MM` or `ML`
modified-base tags, so it tests graceful handling of unavailable methylation data rather than
methylation calculations.

Regenerate the fixtures in an environment containing Python and samtools:

```text
python tests/fixtures/prepare_integration_fixtures.py tests/test_case.zip
```

The script refuses any source archive whose SHA256 differs from the approved archive. The original
ZIP and unsanitized extracted files must remain untracked.

By default, integration-test outputs use a temporary directory and are deleted. To retain a run for
manual inspection, set `REAL_INTEGRATION_OUTPUT_DIR` to a new, nonexistent directory. Persistent
outputs under `tests/integration_output/` are ignored by Git.

`expected_results.json` contains the team-approved scientific and workflow regression baseline.
It intentionally excludes generated timestamps, absolute paths, complete logs, and whole-output
checksums. Scientific-baseline approval and redistribution authorization remain separately recorded.
