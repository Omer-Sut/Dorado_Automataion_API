# Dorado Compatibility And Public Workflow Scope

This document records the supported user-facing workflow boundary and the external Dorado
requirements. Dorado is not installed by the Python package or the Bioconda recipe.

## Public GUI Workflows

The distributable application supports the three workflows exposed through the GUI:

1. POD5 input with Basecalling.
2. POD5 input with Basecalling followed by NanoTel analysis.
3. Existing FASTQ or BAM input with NanoTel analysis.

The backend contains additional standalone workflow methods and CLI subcommands for development and
internal composition. They are not part of the promised public user interface and should not be
advertised as separately supported workflows.

## Dorado Boundary

Dorado is required only for the two public workflows that start from POD5. Existing FASTQ/BAM input
does not invoke Dorado. The GUI-reachable Dorado commands are:

- `dorado basecaller`, which converts POD5 signal into BAM and can optionally call modified bases or
  align against the selected chromosome-end reference.
- `dorado demux`, which classifies the resulting reads and writes barcode-specific BAM files.

The workflow does not invoke `dorado aligner`. Mapping outside basecalling is performed with
minimap2 and samtools.

## Pinned Model And Dorado Version

The team has used one basecalling model and will keep it fixed for the initial package:

- `dna_r10.4.1_e8.2_400bps_sup@v5.2.0`

`v5.2.0` is the model version, not necessarily the Dorado application version. The exact output of
`dorado --version` from the lab's working installation has not yet been recorded in this repository.

- `PACKAGING TODO [DORADO-EXECUTABLE-VERSION]` - Record the exact tested Dorado executable version
  before release and use that version in installation and compatibility instructions.
- If the team changes either the Dorado executable version or the basecalling model, update this
  document, the default configuration, user installation instructions, compatibility checks, and
  POD5 integration-test record together.
- Do not describe a newer Dorado release as supported merely because it is available. Support starts
  after the GUI-reachable basecalling and demultiplexing workflow has been tested with it.

## Barcoding Kit

The current demultiplexing configuration supplies `--kit-name SQK-NBD114-24` to Dorado. The kit name
tells Dorado which barcode sequences and arrangement to search for when assigning reads to
`barcode01` through `barcode24`. It is not a software dependency or a file that this package ships.

Using the wrong kit name can reduce classification or place reads in `unclassified`, because Dorado
would search for barcodes that do not match the library preparation. The initial release therefore
needs one of these explicit decisions:

- Support only libraries prepared with `SQK-NBD114-24` and state that limitation to users; or
- Support additional kits and add a reviewed user-facing kit selection mechanism.

The current GUI has no kit selector, so `SQK-NBD114-24` remains the effective fixed kit until the
team decides otherwise.

## Methylation Decision

The GUI currently offers no modified bases, `5mCpG`, and `5mCpG + 5hmCpG`. Whether both methylation
modes belong in the supported public workflow remains `AWAITING DECISION`.

Before approving either mode, verify that its Dorado modified-base code is compatible with the
pinned SUP v5.2.0 model and test the complete GUI path using POD5 input. The existing sanitized BAM
fixture has no `MM`/`ML` tags and cannot validate modified-base calling or Modkit pileup results.
Bioconda's `ont-modkit` remains a runtime dependency so the downstream tool is available if the
methylation workflow is approved.

## Release Checks

- Record `dorado --version` from the installation used by the lab.
- Confirm whether `SQK-NBD114-24` is the only supported library-preparation kit.
- Resolve the supported GUI methylation modes.
- Run a small authorized POD5 fixture through basecalling and demultiplexing with the pinned Dorado
  executable and model.
- Verify the generated BAM, barcode outputs, summary, optional alignment, and modified-base tags.
- Ensure public documentation and entry points expose only the three approved GUI workflows.
