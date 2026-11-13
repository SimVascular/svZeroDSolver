#!/bin/bash
# Build the svMultiPhysics Apptainer image on a compute node (Apptainer is not
# installed on the login nodes). The digest pins the image used locally.
#SBATCH --job-name=svmp-image
#SBATCH --partition=devel
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=build_image_%j.log
set -euo pipefail
IMAGE=docker://simvascular/solver@sha256:d3ba3c7b0dfd8bf6e295e9b3ce5bac0ae5ed862d8a795f28f72d89d43cd0d625
export APPTAINER_CACHEDIR=${APPTAINER_CACHEDIR:-$PWD/.apptainer}
mkdir -p "$APPTAINER_CACHEDIR"
apptainer build --force svmp-solver.sif "$IMAGE"
apptainer exec svmp-solver.sif ls -l /build-trilinos/svMultiPhysics-build/bin/svmultiphysics
apptainer exec svmp-solver.sif mpirun --version | head -1
