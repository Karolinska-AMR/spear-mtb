#!/bin/bash

set -e  # Exit on any error

SRC="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

echo "=================================================="
echo "SPEAR-MTB Setup Script"
echo "=================================================="
echo "Working directory: $SRC"
echo

echo "Step 1: Setting up conda environment..."
echo "- Deactivating current conda environment"
conda deactivate 2>/dev/null || true

echo "- Removing existing 'spear-mtb' environment (if exists)"
conda remove -y -n spear-mtb --all 2>/dev/null || true

echo "- Creating new conda environment from environment.yml"
if [ ! -f "./environment.yml" ]; then
    echo "ERROR: environment.yml not found in current directory!" >&2
    exit 1
fi
conda env create -f ./environment.yml

echo "- Activating spear-mtb environment"
source activate spear-mtb
echo "✓ Conda environment setup complete"
echo

echo "Step 2: Downloading and extracting reference databases..."
cd "$SRC"

# Check if assets.tar.gz already exists
if [ -f "assets.tar.gz" ]; then
    echo "- Found existing assets.tar.gz in current directory"
    read -p "Do you want to re-download? (y/N): " -n 1 -r
    echo
    if [[ $REPLY =~ ^[Yy]$ ]]; then
        echo "- Removing existing assets.tar.gz"
        rm "assets.tar.gz"
    else
        echo "- Using existing assets.tar.gz"
    fi
fi

# Download if not present
if [ ! -f "assets.tar.gz" ]; then
    echo "- Downloading assets from figshare..."
    wget -O "assets.tar.gz" https://figshare.com/ndownloader/files/55983806
    echo "✓ Download complete"
else
    echo "- Using existing assets.tar.gz"
fi

# Extract the tar.gz file
echo "- Extracting assets.tar.gz..."
echo "  This will create the required folder structure"
tar -xvzf assets.tar.gz

echo "✓ Extraction complete"

# Ask if user wants to remove the tar.gz file
echo
read -p "Remove assets.tar.gz file? (y/N): " -n 1 -r
echo
if [[ $REPLY =~ ^[Yy]$ ]]; then
    rm "assets.tar.gz"
    echo "- Removed assets.tar.gz"
else
    echo "- Keeping assets.tar.gz"
fi

echo
echo "=================================================="
echo "✓ Setup completed successfully!"
echo "=================================================="
echo "Next steps:"
echo "1. Activate the environment: conda activate spear-mtb"
echo "2. Verify the extracted folders are in place"
echo "3. Run the SPEAR-MTB workflow"
echo "=================================================="