# Draft Bioconda Recipe

This directory is the team's upstream working copy of the future Bioconda recipe. It is intentionally
incomplete and cannot be submitted while `PACKAGING TODO` release blockers remain.

## Files

- `meta.yaml` describes the package source, dependencies, tests, license, and maintainers.
- `build.sh` installs the upstream Python package into the Conda build environment.

The recipe does not compile minimap2, samtools, or Modkit. Bioconda installs their existing Conda
packages (`minimap2`, `samtools`, and `ont-modkit`) as runtime dependencies. Dorado and its
SUP v5.2.0 neural-network model remain external in this draft. Their tested compatibility and the
three supported GUI workflows are tracked in `../DORADO_COMPATIBILITY.md`.

## Relationship to Bioconda

Bioconda recipes are published from the separate `bioconda-recipes` repository. Once this upstream
package is release-ready:

1. Create a tagged source release in the project repository.
2. Add its archive URL and SHA256 to `meta.yaml`.
3. Resolve every release blocker listed in `../README.md`.
4. Copy the reviewed recipe to `bioconda-recipes/recipes/<final-package-name>/`.
5. Run Bioconda lint/build tests and open a pull request for maintainer review.

After acceptance, the published recipe in `bioconda-recipes` is authoritative. This local copy remains
useful for team review but must be kept synchronized deliberately if retained.

Bioconda is the primary supported distribution route. Native Windows support is deferred and will
use a separately built and tested executable rather than being treated as part of this recipe.

## Resource policy for this draft

The Python package includes the approved configuration, NanoTel script, runtime R scripts, and GUI
icons. Their installed paths are covered by package and recipe tests. References and trained models
remain excluded until their inventory blockers are resolved. Do not add the Dorado executable or the
315 MB ONT model without explicit permission. Current decisions and checksums are tracked in
`../RESOURCE_INVENTORY.md`.

Sanitized real-data fixtures are copied into the Bioconda test environment through
`test.source_files`. They are not installed in the Python wheel, source distribution, or runtime
environment. Conda-build does preserve them under the downloadable Conda artifact's `info/test`
metadata, where they contributed approximately 18.7 MB in the local rehearsal. The lab authorized
public redistribution on 2026-08-05 after sample/run identifiers were replaced with anonymous test
labels. See `../BIOCONDA_BUILD_REHEARSAL.md` for the complete local build record.
