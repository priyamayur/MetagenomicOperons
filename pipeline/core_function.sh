
#!/bin/bash

# ==================================================
# 1. Setup
# ==================================================

SCRIPT=$(realpath -s "$0")
SCRIPTPATH=$(dirname "$SCRIPT")


# ==================================================
# 2. Help message
# ==================================================

show_usage() {
    cat <<'EOF'
POEM - Core Operon Identification

Usage:
  ./poem -s core-function -o <poem_results_dir>

Arguments:
  -o <directory>  Root directory for POEM pipeline results
  -h, --help      Show this help message

Input:
  POEM uses the processed operon predictions located at:
    <poem_results_dir>/operon_prediction_output/

Output:
  <poem_results_dir>/core_operon_output/

Example:
  ./poem -s core-function -o poem_stage_test
EOF
}


# ==================================================
# 3. Initialize variables
# ==================================================

operon_output_path=""
output_folder=""


# ==================================================
# 4. Parse command-line arguments
# ==================================================

while [[ $# -gt 0 ]]; do

    case "$1" in

        -opo|-o)
            if [ $# -lt 2 ] || [[ "$2" == -* ]]; then
                echo "Error: Missing value for $1"
                exit 1
            fi

            case "$1" in
                -opo) operon_output_path="$2" ;;
                -o) output_folder="$2" ;;
            esac

            shift 2
            ;;

        -h|--help)
            show_usage
            exit 0
            ;;

        *)
            echo "Error: Unknown argument '$1'"
            show_usage
            exit 1
            ;;

    esac

done


# ==================================================
# 5. Load configuration
# ==================================================

source "$SCRIPTPATH/config_pipeline" || exit 1


# ==================================================
# 6. Validate inputs
# ==================================================

if [ -z "$operon_output_path" ]; then
    echo "Error: Operon prediction directory is required (-opo)."
    exit 1
fi

if [ -z "$output_folder" ]; then
    echo "Error: Output directory is required (-o)."
    exit 1
fi

cog_pairs="${operon_output_path}/cog_pairs_df"

if [ ! -f "$cog_pairs" ]; then
    echo "Error: Processed COG pairs file not found:"
    echo "  $cog_pairs"
    exit 1
fi


# ==================================================
# 7. Prepare output directory
# ==================================================

# Remove previous core-operon results
if [ -d "$output_folder" ]; then
    echo "Removing previous core-operon results..."
    rm -rf -- "$output_folder" || exit 1
fi

mkdir -p "$output_folder" || exit 1

echo ""
echo "POEM - Starting core-operon identification"
echo "Input: $cog_pairs"
echo "Output directory: $output_folder"
echo ""


# ==================================================
# 8. Identify core operons
# ==================================================

core_log="${output_folder}/core_function.log"

echo "Identifying core operons..."

python "$script_path/find_community.py" \
    -cp "$cog_pairs" \
    -o "$output_folder" >"$core_log" 2>&1 || {
        status=$?
        echo "Error: Core-operon identification failed."
        tail -n 15 -- "$core_log" >&2
        exit "$status"
    }


# ==================================================
# 9. Validate outputs
# ==================================================

if [ ! -s "$output_folder/leiden_communities.pkl" ]; then
    echo "Error: Leiden communities file was not generated."
    exit 1
fi

if [ ! -s "$output_folder/core_operons.csv" ]; then
    echo "Error: Core-operon results file was not generated."
    exit 1
fi


# ==================================================
# 10. Completion
# ==================================================

echo ""
echo "Core-operon identification completed successfully."
echo "Output directory: $output_folder"
