#!/bin/bash

# Add ~/bin to PATH if not already there
export PATH="$HOME/bin:$PATH"

# Source Conda setup from your local install in ~/bin/miniconda3
if [ -f "$HOME/bin/miniconda3/etc/profile.d/conda.sh" ]; then
    source "$HOME/bin/miniconda3/etc/profile.d/conda.sh"
    conda activate nfcore_env
else
    echo "ERROR: Conda not found in ~/bin/miniconda3. Aborting."
    return 1 2>/dev/null || exit 1
fi

# Optional: check if nf-core and nextflow work
echo "Environment 'nfcore_env' activated."
which nf-core
which nextflow