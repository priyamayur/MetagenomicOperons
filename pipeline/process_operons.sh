#!/bin/bash

show_usage() {
    cat <<'EOF'
POEM - Operon Processing

Usage:
  process_operons.sh -fpo <function_prediction_dir> -opo <operon_prediction_dir> -cdo <cdhit_dir>

Arguments:
  -fpo <directory>  Function prediction results
  -opo <directory>  Operon predictions and processing outputs
  -cdo <directory>  CD-HIT results
  -h, --help        Show this help message

POEM runs this step after predicting operons.
EOF
}

function_prediction_path=""
operon_output_path=""
cdhit_output_path=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        -h|--help)
            show_usage
            exit 0
            ;;
        -fpo|-opo|-cdo)
            if [[ $# -lt 2 || -z "$2" || "$2" == -* ]]; then
                printf 'Error: Missing directory for %s.\n' "$1" >&2
                exit 1
            fi

            case "$1" in
                -fpo) function_prediction_path="$2" ;;
                -opo) operon_output_path="$2" ;;
                -cdo) cdhit_output_path="$2" ;;
            esac
            shift 2
            ;;
        *)
            printf 'Error: Unknown argument %s.\n' "$1" >&2
            exit 1
            ;;
    esac
done

if [[ -z "$function_prediction_path" || -z "$operon_output_path" || -z "$cdhit_output_path" ]]; then
    printf 'Error: Operon processing requires -fpo, -opo, and -cdo.\n' >&2
    exit 1
fi

SCRIPT=$(realpath -s "$0") || exit 1
SCRIPTPATH=$(dirname "$SCRIPT") || exit 1
source "$SCRIPTPATH/config_pipeline" || exit 1

cog_reference="${cog_path}/cog-20.cog.csv"
diamond_matches="${function_prediction_path}/diamond_cog_reference.tsv"
operon_pairs="${operon_output_path}/gene_pair_names_to_operon.csv"
gene_clusters="${cdhit_output_path}/gene_to_clusters"
gene_to_cog="${operon_output_path}/gene_to_cogId.pkl"
processing_log="${operon_output_path}/process_operons.log"

for input_file in "$cog_reference" "$diamond_matches" "$operon_pairs" "$gene_clusters"; do
    if [[ ! -f "$input_file" || ! -r "$input_file" ]]; then
        printf 'Error: Required file is missing or unreadable: %s\n' "$input_file" >&2
        exit 1
    fi
done

: > "$processing_log" || exit 1
printf 'Operon processing log: %s\n' "$processing_log"

python "$script_path/process_gene_cog_relation.py" \
    -cp "$cog_reference" \
    -ca "$diamond_matches" \
    -o "$operon_output_path" >>"$processing_log" 2>&1 || {
        status=$?
        printf 'Error: Building the gene-to-COG mapping failed.\n' >&2
        tail -n 15 -- "$processing_log" >&2
        exit "$status"
    }

if [[ ! -s "$gene_to_cog" ]]; then
    printf 'Error: Gene-to-COG mapping is missing or empty: %s\n' "$gene_to_cog" >&2
    exit 1
fi

python "$script_path/annotate_operon_genes.py" \
    -og "$operon_pairs" \
    -gc "$gene_clusters" \
    -gcid "$gene_to_cog" \
    -o "$operon_output_path" >>"$processing_log" 2>&1 || {
        status=$?
        printf 'Error: Annotating predicted operons failed.\n' >&2
        tail -n 15 -- "$processing_log" >&2
        exit "$status"
    }