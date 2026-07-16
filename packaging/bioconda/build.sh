#!/bin/bash
set -euxo pipefail

# Bioconda runs this script after downloading and unpacking the source archive.
# $PYTHON is the interpreter inside the Conda build environment. Installing with it places
# the package under $PREFIX, the isolated installation root created for this Conda package.
# --no-deps prevents pip from downloading dependencies because Conda installs everything
# declared in meta.yaml. --no-build-isolation keeps the build reproducible inside that environment.

# PACKAGING TODO [SOURCE-URL] - RELEASE BLOCKER
# This command will work only after meta.yaml points to a source archive containing pyproject.toml.

"${PYTHON}" -m pip install . --no-deps --no-build-isolation
