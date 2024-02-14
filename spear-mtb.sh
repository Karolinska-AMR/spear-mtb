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
work_dir="$PSRC/.tmp/"
config_file="$SRC/nextflow.config"
assets_dir="$SRC/assets" 
out_dir="$PSRC/out"
profile='slurm'
archive_dir="$out_dir/raw_seqs"


# Parse user-defined parameters
while getopts "w:c:a:t:o:f:v:" opt; do
  case $opt in
    w)
      work_dir=$(readlink -f "$OPTARG")
      ;;
    c)
      config_file=$(readlink -f "$OPTARG")
      ;;
    a)
      assets_dir=$(readlink -f "$OPTARG")
      ;;
    t)
      ticket="$OPTARG"
      ;;
    o)
      out_dir=$(readlink -f "$OPTARG")
      ;;
    f)
      profile="$OPTARG"
      ;;  
    v)
      archive_dir=$(readlink -f "$OPTARG")
    ;;
    \?)
      echo "Invalid option: -$OPTARG" >&2
      exit 1
      ;;
  esac
done

work_dir="$work_dir/$ticket"
if [ ! -d "$work_dir" ]; then
    mkdir -p "$work_dir"
fi

trace_dir="$out_dir/traces"
if [ ! -d "$trace_dir" ]; then
    mkdir -p "$trace_dir"
fi

shift "$((OPTIND - 1))"

# Check if the input is provided
if [ "$#" -ne 1 ]; then
  echo "Usage: ./spear-mtb.sh [-w work_dir] [-c config_file] [-a assets_dir] [-t ticket] [-f profile] [-o out_dir] [-v archive_dir] input_directory"
  exit 1
fi

input_dir=$(readlink -f "$1")

# Set a default ticket for the nextflow trace report
default_ticket="$out_dir/$(date +'%y%m%d-%H%M%S')"

if [ -z "$ticket" ]; then
  ticket="$default_ticket"
fi

trace_file="$trace_dir/${ticket}.trace"
log_file="$trace_dir/${ticket}.log"

echo "Using the following parameters:"
echo "  -Input directory: $input_dir"
echo "  -Output directory: $out_dir"
echo "  -Archive [PE-reads] directory: $archive_dir"
echo "  -Config file: $config_file"
echo "  -Assets directory: $assets_dir"
echo "  -Trace file: $trace_file"
echo "  -Work directory: $work_dir"
echo "  -Profile: $profile"
echo "----------------------"

echo "Running SPEAR-MTB..."

cd "$work_dir"

nice -5 nextflow run "$SRC/main.nf" -profile "$profile" -c "$config_file" -w "$work_dir" --assets_dir "$assets_dir" --out_dir "$out_dir" --archive "$archive_dir" --input_dir "$input_dir" --ticket "$ticket" -with-trace "$trace_file" > "$log_file"

echo ""
echo "    ______ _         _        __               __"
echo "   / ____/(_)____   (_)_____ / /_   ___   ____/ /"
echo "  / /_   / // __ \ / // ___// __ \ / _ \ / __  / "
echo " / __/  / // / / // /(__  )/ / / //  __// /_/ /  "
echo "/_/    /_//_/ /_//_//____//_/ /_/ \___/ \__,_/   "
                                                 


