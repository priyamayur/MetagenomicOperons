#!/bin/bash

SCRIPT=$(realpath -s "$0")
SCRIPTPATH=$(dirname "$SCRIPT")

show_usage() {
    cat <<'EOF'
POEM - Core Operon Visualization

Usage:
  ./poem -s visualize-core-operon -o <poem_results_dir>
  ./pipeline/visualize_core_operon.sh -coo <core_operon_dir>

Arguments:
  -o <directory>    Root directory for POEM results (via ./poem)
  -coo <directory>  Core-operon output directory (direct script invocation)
  -h, --help        Show this help message

Required inputs:
  <core_operon_dir>/leiden_communities.pkl
  <core_operon_dir>/core_operons.csv
  cog-20.def.tab in the configured COG reference directory

Example:
  ./poem -s visualize-core-operon -o example_data/output_example
EOF
}

core_operon_path=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        -coo)
            if [[ $# -lt 2 || "$2" == -* || -z "$2" ]]; then
                echo "Error: Missing value for $1" >&2
                exit 1
            fi
            core_operon_path="$2"
            shift 2
            ;;
        -h|--help)
            show_usage
            exit 0
            ;;
        *)
            echo "Error: Unknown argument '$1'" >&2
            show_usage >&2
            exit 1
            ;;
    esac
done

if [[ -z "$core_operon_path" ]]; then
    echo "Error: Core-operon directory is required (-coo)." >&2
    exit 1
fi

source "$SCRIPTPATH/config_pipeline" || exit 1

leiden_graph="$core_operon_path/leiden_communities.pkl"
core_file="$core_operon_path/core_operons.csv"
cog_def_file="$cog_path/cog-20.def.tab"

for file in "$leiden_graph" "$core_file" "$cog_def_file"; do
    if [[ ! -s "$file" ]]; then
        printf 'Error: Required file not found or empty: %s\n' "$file" >&2
        exit 1
    fi
done

printf '\nPOEM - Starting core-operon visualization\nCore-operon directory: %s\n\n' "$core_operon_path"

python "$script_path/visualize_core_operon.py" \
    -g "$leiden_graph" \
    -cores "$core_file" \
    -cdef "$cog_def_file" \
    -o "$core_operon_path" || {
        echo "Error: Core-operon visualization failed." >&2
        exit 1
    }
