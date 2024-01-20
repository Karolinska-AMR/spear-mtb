#!/bin/bash


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

# default values
tmp_dir="$PSRC/.tmp/"
config_file="$SRC/nextflow.config"
assets_dir="$SRC/assets" 
out_dir="$PSRC/out"
profile='slurm'

# Set a default prefix for the nextflow trace report
default_prefix="$out_dir/$(date +'%y%m%d-%H%M%S')"

# Parse user-defined parameters
while getopts "t:c:a:p:o:f:" opt; do
  case $opt in
    t)
      tmp_dir=$(readlink -f "$OPTARG")
      ;;
    c)
      config_file=$(readlink -f "$OPTARG")
      ;;
    a)
      assets_dir=$(readlink -f "$OPTARG")
      ;;
    p)
      prefix="$OPTARG"
      ;;
    o)
      out_dir=$(readlink -f "$OPTARG")
      ;;
    f)
      profile="$OPTARG"
      ;;  
    \?)
      echo "Invalid option: -$OPTARG" >&2
      exit 1
      ;;
  esac
done

work_dir="$tmp_dir/work"
if [ ! -d "$work_dir" ]; then
    mkdir -p "$work_dir"
fi

shift "$((OPTIND - 1))"

# Check if the input is provided
if [ "$#" -ne 1 ]; then
  echo "Usage: ./spear-mtb.sh [-t tmp_dir] [-c config_file] [-a assets_dir] [-p prefix] [-f profile] [-o out_dir] input_directory"
  exit 1
fi

input_dir=$(readlink -f "$1")

if [ -z "$prefix" ]; then
  prefix="$default_prefix"
fi

trace_file="${prefix}_trace.txt"

echo "Using the following parameters:"
echo "  -Input directory: $input_dir"
echo "  -Output directory: $out_dir"
echo "  -Config file: $config_file"
echo "  -Assets directory: $assets_dir"
echo "  -Trace file: $trace_file"
echo "  -Temporary directory: $tmp_dir"
echo "  -Profile: $profile"
echo "----------------------"

echo "Running SPEAR-MTB..."

cd "$tmp_dir"
source activate spear-mtb
nextflow run "$SRC/main.nf" -profile "$profile" -c "$config_file" -w "$work_dir" --assets_dir "$assets_dir" --out_dir "$out_dir" --input_dir "$input_dir" -with-trace "$trace_file" -resume

echo ""
echo "    ______ _         _        __               __"
echo "   / ____/(_)____   (_)_____ / /_   ___   ____/ /"
echo "  / /_   / // __ \ / // ___// __ \ / _ \ / __  / "
echo " / __/  / // / / // /(__  )/ / / //  __// /_/ /  "
echo "/_/    /_//_/ /_//_//____//_/ /_/ \___/ \__,_/   "
                                                 


