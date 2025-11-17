#!/bin/bash
# Launcher script for COLOSS GUI with conda environment

# Get the directory where this script is located
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# Source conda
CONDA_BASE=$(conda info --base 2>/dev/null)
if [ -z "$CONDA_BASE" ]; then
    # Try common conda locations
    if [ -f "$HOME/miniconda3/etc/profile.d/conda.sh" ]; then
        source "$HOME/miniconda3/etc/profile.d/conda.sh"
    elif [ -f "$HOME/anaconda3/etc/profile.d/conda.sh" ]; then
        source "$HOME/anaconda3/etc/profile.d/conda.sh"
    else
        echo "Error: conda not found!"
        echo "Please install conda or run compile.sh to set up the GUI"
        exit 1
    fi
else
    source "$CONDA_BASE/etc/profile.d/conda.sh"
fi

# Activate conda environment
conda activate coloss_gui

if [ $? -ne 0 ]; then
    echo "Error: Failed to activate coloss_gui environment"
    echo "Please run: ./compile.sh and choose option 2 to set up the GUI"
    exit 1
fi

# Change to the GUI directory
cd "$SCRIPT_DIR/coloss_gui"

# Run the GUI
echo "Starting COLOSS GUI..."
python main.py
