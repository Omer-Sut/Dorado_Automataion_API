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
| `R-DEPENDENCIES` | Source audit in progress | NanoTel is classified; active and possibly unused `r_analysis` dependencies are recorded for team classification | `R_DEPENDENCIES.md`, `meta.yaml` |
| `NANOTEL-DEPENDENCIES` | Clean-environment test passed | Direct dependencies and the approved summary/filtration baseline passed in an isolated Linux environment | `R_DEPENDENCIES.md`, `meta.yaml` |
| `ENTRY-POINTS` | CLI placeholder complete | `telomere-analyzer` is temporary; GUI command and final names remain unresolved | GUI entry module, `pyproject.toml`, recipe tests |
| `PACKAGE-DATA` | Approved core resources included | Config, NanoTel, runtime R scripts, and icons are packaged; references remain partly unresolved | `RESOURCE_INVENTORY.md`, `pyproject.toml`, installed-resource tests |
| `SOURCE-URL` | Requires release | URL of the final tagged source archive | `meta.yaml` |
| `SOURCE-SHA256` | Requires release | SHA256 of that exact source archive | `meta.yaml` |
| `BIOCONDA-TESTS` | Local integration passed | Real FASTQ/BAM tests passed with recipe dependencies; built-package CI and GUI coverage remain | `tests/test_real_integration.py`, `meta.yaml` |
| `DORADO` | Team/license decision | Supported Dorado versions and external installation instructions | README and runtime dependency checks |
| `DORADO-MODEL` | Redistribution decision | Whether ONT model files may ship; default is external | Package-data settings and user documentation |
| `EXTERNAL-RESOURCE-CONFIG` | Code/design work required | Shared per-user settings, GUI selectors, precedence, validation, and cross-platform tests | `RESOURCE_INVENTORY.md`, config manager, GUI settings |
| `REFERENCES` | Partial decision | Human is CC0 pending derivation details; mouse awaits provenance; zebrafish is preferred to ship but temporarily excluded pending source/license and file review | `RESOURCE_INVENTORY.md`, package-data settings, notices |
| `NANOTEL` | Approved and integration-tested | Lab-owned, packaged, and verified against the approved real-data regression baseline | `RESOURCE_INVENTORY.md`, `R_DEPENDENCIES.md`, recipe |
| `RDS-MODEL` | Excluded | Awaiting a decision on deployment of the NanoTel statistical analysis | `RESOURCE_INVENTORY.md`, NanoTel runtime behavior |
| `MODKIT` | Technical decision | Whether it is required, optional, or unsupported | `meta.yaml`, README, runtime dependency checks |
| `WINDOWS` | Support decision | Native executable, WSL2, best effort, or unsupported | User documentation and separate release process |

## Current package boundary

The current skeleton packages Python source plus the approved default configuration, NanoTel script,
runtime R-analysis scripts, and GUI icons. References, trained models, executables, test data,
development artifacts, and generated output remain excluded unless the resource inventory approves
them explicitly.

Confirmed Conda runtime dependencies currently recorded in the draft recipe are `minimap2`,
`samtools`, `pyside6`, `r-base`, and the packages loaded by the core R-analysis scripts. Dorado is
documented as external while its distribution status is unresolved.

The root `requirements.txt` is a snapshot of the development environment and includes transitive or
development tools that the application does not import directly. `pyproject.toml` is the authoritative
list of direct Python runtime dependencies for the distributable package.

Detailed decisions, provenance, sizes, and checksums for non-Python files are maintained in
[`RESOURCE_INVENTORY.md`](RESOURCE_INVENTORY.md).
NanoTel and `r_analysis` dependency classifications and remaining execution tests are maintained in
[`R_DEPENDENCIES.md`](R_DEPENDENCIES.md).

## Publication gate

Do not publish `package-name-pending`, version `0.0.0`, or a recipe containing unresolved release
blockers. Before release, build and install the package in a clean environment and inspect the
archives to confirm that only approved files are present.
