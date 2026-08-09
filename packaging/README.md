# Packaging Worklist

This directory tracks the work required to distribute the project as an installable Python
package and through Bioconda. The current files are a draft skeleton and must not be published.

Every unresolved field is marked with `PACKAGING TODO`. Items marked `RELEASE BLOCKER` must be
resolved before creating a public release or submitting a recipe to Bioconda.

## How to update a placeholder

1. Record the team's decision in this checklist.
2. Replace the temporary value or commented field in every affected file.
3. Remove the corresponding `PACKAGING TODO` only after the replacement is tested.

## Decision checklist

| ID | Status | Missing information | Where to fill it in |
| --- | --- | --- | --- |
| `PACKAGE-NAME` | Waiting for team | Final public package and command names | `pyproject.toml`, `packaging/bioconda/meta.yaml`, recipe directory name |
| `VERSION` | Waiting for team | Authoritative version source and tagging policy | `pyproject.toml`, package code, `meta.yaml`, Git release tag |
| `LICENSE` | Release blocker | Approved project license and SPDX identifier | Root `LICENSE`, `pyproject.toml`, `meta.yaml` |
| `AUTHORS` | Waiting for team | Author attribution | `pyproject.toml` |
| `MAINTAINERS` | Waiting for team | Release and Bioconda maintainers | `pyproject.toml`, `meta.yaml` |
| `PROJECT-URLS` | Waiting for public locations | Repository, documentation, issues, and publication links | `pyproject.toml`, `meta.yaml` |
| `PYTHON-DEPENDENCIES` | Initial audit complete | Python 3.10+ and direct `PySide6` dependency recorded; retest after entry-point work | `pyproject.toml`, `meta.yaml` |
| `R-DEPENDENCIES` | Declared and package-tested | NanoTel and `r_analysis` dependencies are recorded; missing packages now produce an actionable error instead of a runtime download | `R_DEPENDENCIES.md`, `meta.yaml` |
| `NANOTEL-DEPENDENCIES` | Conda package test passed | Direct dependencies and the approved summary/filtration baseline passed from the built package | `R_DEPENDENCIES.md`, `BIOCONDA_BUILD_REHEARSAL.md`, `meta.yaml` |
| `ENTRY-POINTS` | GUI launcher implemented | `telomere-analyzer` launches the installed GUI and exposes only help/version; final command name remains unresolved | `dorado_gui/main.py`, `pyproject.toml`, recipe tests |
| `PACKAGE-DATA` | Approved core resources included | Config, NanoTel, runtime R scripts, and icons are packaged; references remain partly unresolved | `RESOURCE_INVENTORY.md`, `pyproject.toml`, installed-resource tests |
| `SOURCE-URL` | Requires release | URL of the final tagged source archive | `meta.yaml` |
| `SOURCE-SHA256` | Requires release | SHA256 of that exact source archive | `meta.yaml` |
| `BIOCONDA-TESTS` | Local conda-build passed | Build, clean install, GUI smoke, resources, tools, R namespaces, and real FASTQ/BAM tests passed; official Bioconda CI remains | `BIOCONDA_BUILD_REHEARSAL.md`, `meta.yaml` |
| `DORADO` | External; compatibility audit in progress | SUP v5.2.0 model is pinned; record the tested Dorado executable version and kit policy | `DORADO_COMPATIBILITY.md`, runtime checks |
| `DORADO-MODEL` | External and pinned | Users obtain `dna_r10.4.1_e8.2_400bps_sup@v5.2.0` separately | `DORADO_COMPATIBILITY.md`, user documentation |
| `EXTERNAL-RESOURCE-CONFIG` | Code/design work required | Shared per-user settings, GUI selectors, precedence, validation, and cross-platform tests | `RESOURCE_INVENTORY.md`, config manager, GUI settings |
| `REFERENCES` | Partial decision | Human is CC0 pending derivation details; mouse awaits provenance; zebrafish is preferred to ship but temporarily excluded pending source/license and file review | `RESOURCE_INVENTORY.md`, package-data settings, notices |
| `NANOTEL` | Approved and integration-tested | Lab-owned, packaged, and verified against the approved real-data regression baseline | `RESOURCE_INVENTORY.md`, `R_DEPENDENCIES.md`, recipe |
| `RDS-MODEL` | Excluded | Awaiting a decision on deployment of the NanoTel statistical analysis | `RESOURCE_INVENTORY.md`, NanoTel runtime behavior |
| `MODKIT` | Approved Conda dependency | Install Bioconda's `ont-modkit`; methylation pileup still needs a tagged-BAM integration test | `meta.yaml`, integration tests |
| `WINDOWS` | Deferred separate distribution | Bioconda is the primary release; a native Windows executable will be packaged and validated later through a separate process | User documentation and future Windows release work |

## Current package boundary

The current skeleton packages Python source plus the approved default configuration, NanoTel script,
runtime R-analysis scripts, and GUI icons. References, trained models, executables, test data,
development artifacts, and generated output remain excluded unless the resource inventory approves
them explicitly.

Confirmed Conda runtime dependencies currently recorded in the draft recipe are `minimap2`,
`samtools`, `ont-modkit`, `pyside6`, `r-base`, and the packages loaded by the core R-analysis
scripts. Dorado and its pinned SUP v5.2.0 model are external; the tested executable version and
installation instructions remain unresolved.

The root `requirements.txt` is a snapshot of the development environment and includes transitive or
development tools that the application does not import directly. `pyproject.toml` is the authoritative
list of direct Python runtime dependencies for the distributable package.

Detailed decisions, provenance, sizes, and checksums for non-Python files are maintained in
[`RESOURCE_INVENTORY.md`](RESOURCE_INVENTORY.md).
NanoTel and `r_analysis` dependency classifications and remaining execution tests are maintained in
[`R_DEPENDENCIES.md`](R_DEPENDENCIES.md).
The successful local Conda build and installed-package evidence are recorded in
[`BIOCONDA_BUILD_REHEARSAL.md`](BIOCONDA_BUILD_REHEARSAL.md).

## Publication gate

Do not publish `package-name-pending`, version `0.0.0`, or a recipe containing unresolved release
blockers. Before release, build and install the package in a clean environment and inspect the
archives to confirm that only approved files are present.
