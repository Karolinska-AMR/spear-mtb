#!/bin/bash

# Script to run Nextflow pipeline with specified parameters
# Usage: ./spear-mtb.sh --input_dir <path> --output_dir <path> [options]

set -e  # Exit on any error

# Function to display usage
usage() {
    cat << EOF
Usage: $0 --input_dir <path> --output_dir <path> [OPTIONS]

Required parameters:
    --input_dir <path>      Input directory path
    --output_dir <path>     Output directory path

Optional parameters:
    --working_dir <path>    Working directory (default: ./[source_dir]/.tmp)
    --config_file <path>    Nextflow config file
    --assets_dir <path>     Assets directory
    --profile <name>        Nextflow profile to use
    --archive_dir <path>    Archive directory
    --resume                Resume previous run 
    --skip_cryptic          Skip cryptic workflow 
    --ticket <value>        Ticket ID 
    -h, --help              Show this help message

Examples:
    $0 --input_dir /data/input --output_dir /data/output
    $0 --input_dir /data/input --output_dir /data/output --profile slurm --resume
EOF
}

# Function to generate random alphanumeric string
generate_ticket() {
    cat /dev/urandom | tr -dc 'a-zA-Z0-9' | fold -w 8 | head -n 1
}

# Function to convert to absolute path
get_absolute_path() {
    local path="$1"
    if [[ -z "$path" ]]; then
        echo ""
        return
    fi
    
    # If path is already absolute, return as is
    if [[ "$path" = /* ]]; then
        echo "$path"
    else
        # Convert relative path to absolute
        echo "$(cd "$(dirname "$path")" 2>/dev/null && pwd)/$(basename "$path")"
    fi
}

echo " _______  _______  _______  _______  _______         _______ _________ ______  "
echo "(  ____ \(  ____ )(  ____ \(  ___  )(  ____ )       (       )\__   __/(  ___ \ "
echo "| (    \/| (    )|| (    \/| (   ) || (    )|       | () () |   ) (   | (   ) )"
echo "| (_____ | (____)|| (__    | (___) || (____)| _____ | || || |   | |   | (__/ / "
echo "(_____  )|  _____)|  __)   |  ___  ||     __)(_____)| |(_)| |   | |   |  __ (  "
echo "      ) || (      | (      | (   ) || (\ (          | |   | |   | |   | (  \ \ "
echo "/\____) || )      | (____/\| )   ( || ) \ \__       | )   ( |   | |   | )___) )"
echo "\_______)|/       (_______/|/     \||/   \__/       |/     \|   )_(   |/ \___/ "
echo ""                                                                               


SRC="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
PSRC=$(cd "$SRC/.." && pwd)

# Initialize variables
INPUT_DIR=""
OUTPUT_DIR=""
WORKING_DIR="$PSRC/.tmp/"
CONFIG_FILE="$SRC/nextflow.config"
ASSETS_DIR="$SRC/assets" 
PROFILE="slurm"
ARCHIVE_DIR=""
RESUME=false
SKIP_CRYPTIC=false
TICKET=""

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --input_dir)
            INPUT_DIR="$2"
            shift 2
            ;;
        --output_dir)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        --working_dir)
            WORKING_DIR="$2"
            shift 2
            ;;
        --config_file)
            CONFIG_FILE="$2"
            shift 2
            ;;
        --assets_dir)
            ASSETS_DIR="$2"
            shift 2
            ;;
        --profile)
            PROFILE="$2"
            shift 2
            ;;
        --archive_dir)
            ARCHIVE_DIR="$2"
            shift 2
            ;;
        --resume)
            RESUME=true
            shift
            ;;
        --skip_cryptic)
            SKIP_CRYPTIC=true
            shift
            ;;
        --ticket)
            TICKET="$2"
            shift 2
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            echo "Error: Unknown parameter $1"
            usage
            exit 1
            ;;
    esac
done

# Check required parameters
if [[ -z "$INPUT_DIR" ]]; then
    echo "Error: --input_dir is required"
    usage
    exit 1
fi

if [[ -z "$OUTPUT_DIR" ]]; then
    echo "Error: --output_dir is required"
    usage
    exit 1
fi

# Generate default ticket if not provided
if [[ -z "$TICKET" ]]; then
    TICKET=$(generate_ticket)
fi
# Convert paths to absolute paths
echo "Converting paths to absolute paths..."
INPUT_DIR=$(get_absolute_path "$INPUT_DIR")
OUTPUT_DIR=$(get_absolute_path "$OUTPUT_DIR")
WORKING_DIR=$(get_absolute_path "$WORKING_DIR")

if [[ -n "$CONFIG_FILE" ]]; then
    CONFIG_FILE=$(get_absolute_path "$CONFIG_FILE")
fi

if [[ -n "$ASSETS_DIR" ]]; then
    ASSETS_DIR=$(get_absolute_path "$ASSETS_DIR")
fi

if [[ -n "$ARCHIVE_DIR" ]]; then
    ARCHIVE_DIR=$(get_absolute_path "$ARCHIVE_DIR")
fi

# Validate directories exist (for input_dir and config_file if specified)
if [[ ! -d "$INPUT_DIR" ]]; then
    echo "Error: Input directory '$INPUT_DIR' does not exist"
    exit 1
fi

if [[ -n "$CONFIG_FILE" && ! -f "$CONFIG_FILE" ]]; then
    echo "Error: Config file '$CONFIG_FILE' does not exist"
    exit 1
fi

# Create output and working directories if they don't exist
mkdir -p "$OUTPUT_DIR"
mkdir -p "$WORKING_DIR"

if [[ -n "$ARCHIVE_DIR" ]]; then
    mkdir -p "$ARCHIVE_DIR"
fi

# Print all parameters
echo "============================================="
echo "Nextflow Pipeline Configuration"
echo "============================================="
echo "Input Directory:    $INPUT_DIR"
echo "Output Directory:   $OUTPUT_DIR"
echo "Working Directory:  $WORKING_DIR"
echo "Config File:        ${CONFIG_FILE:-"Not specified"}"
echo "Assets Directory:   ${ASSETS_DIR:-"Not specified"}"
echo "Profile:            ${PROFILE:-"Not specified"}"
echo "Archive Directory:  ${ARCHIVE_DIR:-"Not specified"}"
echo "Resume:             $RESUME"
echo "Skip Cryptic:       $SKIP_CRYPTIC"
echo "Ticket:             $TICKET"
echo "============================================="

# Ask user for confirmation
read -p "Do you want to continue with these parameters? (y/N): " -n 1 -r
echo
if [[ ! $REPLY =~ ^[Yy]$ ]]; then
    echo "Pipeline execution cancelled."
    exit 0
fi

echo "Proceeding with pipeline execution..."


# Activate the environment
if ! source activate spear-mtb; then
    echo "Error: Failed to activate conda environment 'spear-mtb'"
    echo "Please ensure the environment exists by running: conda env list"
    exit 1
fi


# Build nextflow command
NEXTFLOW_CMD="nice -5 nextflow run \"$SRC/main.nf\" -process.log"

# Add parameters to the command
NEXTFLOW_CMD="$NEXTFLOW_CMD --input_dir '$INPUT_DIR'"
NEXTFLOW_CMD="$NEXTFLOW_CMD --out_dir '$OUTPUT_DIR'"
NEXTFLOW_CMD="$NEXTFLOW_CMD -work-dir '$WORKING_DIR'"
NEXTFLOW_CMD="$NEXTFLOW_CMD --ticket '$TICKET'"

# Add optional parameters
if [[ -n "$CONFIG_FILE" ]]; then
    NEXTFLOW_CMD="$NEXTFLOW_CMD -c '$CONFIG_FILE'"
fi

if [[ -n "$ASSETS_DIR" ]]; then
    NEXTFLOW_CMD="$NEXTFLOW_CMD --assets_dir '$ASSETS_DIR'"
fi

if [[ -n "$PROFILE" ]]; then
    NEXTFLOW_CMD="$NEXTFLOW_CMD -profile '$PROFILE'"
fi

if [[ -n "$ARCHIVE_DIR" ]]; then
    NEXTFLOW_CMD="$NEXTFLOW_CMD --archive_input true --archive_dir '$ARCHIVE_DIR'"
fi

if [[ "$RESUME" == true ]]; then
    NEXTFLOW_CMD="$NEXTFLOW_CMD -resume"
fi

if [[ "$SKIP_CRYPTIC" == true ]]; then
    NEXTFLOW_CMD="$NEXTFLOW_CMD --skip_cryptic"
fi

# Display the command that will be executed
echo "============================================="
echo "Executing Nextflow command:"
echo "$NEXTFLOW_CMD"
echo "============================================="

# Execute the nextflow pipeline
echo "Starting Nextflow pipeline..."

cd "$WORKING_DIR" || { echo "Error: Could not change to working directory '$WORKING_DIR'"; exit 1; }
eval $NEXTFLOW_CMD

echo "Pipeline execution completed!"