#!/bin/bash

# NOTE: This script was written for the Stanford Sherlock HPC cluster.
# You may need to modify the following for your own cluster:
#   - --partition (currently: amarsden)
#   - --mail-user
#   - module load commands and version numbers

# Name of your job
#SBATCH --job-name=tuning_job

# Name of partition
#SBATCH --partition=amarsden

# Specify the name of the output file. The %j specifies the job ID
#SBATCH --output=tuning_job.o%j

# Specify a name of the error file. The %j specifies the job ID
#SBATCH --error=tuning_job.e%j

# The walltime you require for your simulation
#SBATCH --time=03:00:00

# Job priority. Leave as normal for now.
#SBATCH --qos=normal

# Number of nodes you are requesting for your job. You can have 24 processors per node, so plan accordingly
#SBATCH --nodes=1

# Amount of memory you require per node. The default is 4000 MB (or 4 GB) per node
#SBATCH --mem=8000

# Number of processors per node (for parallel differential_evolution)
#SBATCH --ntasks-per-node=24

# Send an email to this address when your job starts and finishes
#SBATCH --mail-user=your_email@your_institution.edu
#SBATCH --mail-type=begin
#SBATCH --mail-type=fail
#SBATCH --mail-type=end

# Load Modules
module purge
module load python/3.14.2
module load py-scipy/1.16.3_py314
module load py-pandas/2.3.3_py314
module load viz
module load py-matplotlib/3.10.8_py314
module load cmake

# Relative paths (submit from this directory)
VENV_DIR="../../venv"
REQUIREMENTS="../../requirements.txt"
SVZEROD_ROOT="../../../.."   # svZeroDSolver repo root (for pip install -e pysvzerod)

# Check if virtual environment exists and has Python; if not, create and install
if [[ ! -f "$VENV_DIR/bin/python" ]]; then
    echo "Virtual environment not found at $VENV_DIR. Creating and installing dependencies..."
    python3 -m venv "$VENV_DIR"
    "$VENV_DIR/bin/pip" install --upgrade pip
    # Install pip packages from requirements
    "$VENV_DIR/bin/pip" install -r "$REQUIREMENTS"
    # Install svZeroDSolver from source
    "$VENV_DIR/bin/pip" install -e "$SVZEROD_ROOT"
    echo "Virtual environment ready."
else
    echo "Using existing virtual environment at $VENV_DIR"
fi

# Run main.py with the venv's Python
source "$VENV_DIR/bin/activate"
python -u main.py

# Submit job with:
# sbatch ./tuning_job.sh
# Check status with:
# squeue -u $USER
