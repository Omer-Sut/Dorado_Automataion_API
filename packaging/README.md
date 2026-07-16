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
| `R-DEPENDENCIES` | Core audit complete | Core R-analysis dependencies recorded; retest when R resources are packaged | `meta.yaml`, README dependency documentation |
| `NANOTEL-DEPENDENCIES` | Waiting for NanoTel decision | Additional CRAN and Bioconductor dependencies used by `NanoTel.R` | `meta.yaml`, README dependency documentation |
| `ENTRY-POINTS` | CLI placeholder complete | `telomere-analyzer` is temporary; GUI command and final names remain unresolved | GUI entry module, `pyproject.toml`, recipe tests |
| `PACKAGE-DATA` | Audit required | Approved non-Python resources to install | `pyproject.toml`, installed-resource path code |
| `SOURCE-URL` | Requires release | URL of the final tagged source archive | `meta.yaml` |
| `SOURCE-SHA256` | Requires release | SHA256 of that exact source archive | `meta.yaml` |
| `BIOCONDA-TESTS` | Code/test work required | Fast tests of installed commands and resources | `meta.yaml` or a future `run_test.sh` |
| `DORADO` | Team/license decision | Supported Dorado versions and external installation instructions | README and runtime dependency checks |
| `DORADO-MODEL` | Redistribution decision | Whether ONT model files may ship; default is external | Package-data settings and user documentation |
| `REFERENCES` | Provenance review | Source and redistribution terms for every FASTA | Package-data settings and third-party notices |
| `NANOTEL` | Provenance review | Ownership and redistribution terms for `NanoTel.R` | Package-data settings and third-party notices |
| `RDS-MODEL` | Provenance review | Ownership, creation details, and redistribution permission | Package-data settings and documentation |
| `MODKIT` | Technical decision | Whether it is required, optional, or unsupported | `meta.yaml`, README, runtime dependency checks |
| `WINDOWS` | Support decision | Native executable, WSL2, best effort, or unsupported | User documentation and separate release process |

## Current package boundary

The first skeleton packages Python source only. It intentionally excludes R scripts, JSON
configuration, icons, references, trained models, executables, test data, development artifacts,
and generated output. Exclusion at this stage does not mean those files can never be included.

Confirmed Conda runtime dependencies currently recorded in the draft recipe are `minimap2`,
`samtools`, `pyside6`, `r-base`, and the packages loaded by the core R-analysis scripts. Dorado is
documented as external while its distribution status is unresolved.

The root `requirements.txt` is a snapshot of the development environment and includes transitive or
development tools that the application does not import directly. `pyproject.toml` is the authoritative
list of direct Python runtime dependencies for the distributable package.

## Publication gate

Do not publish `package-name-pending`, version `0.0.0`, or a recipe containing unresolved release
blockers. Before release, build and install the package in a clean environment and inspect the
archives to confirm that only approved files are present.
