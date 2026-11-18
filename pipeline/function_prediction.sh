
#!/bin/bash

# ==================================================
# 1. Setup
# ==================================================

SCRIPT=$(realpath -s "$0")
SCRIPTPATH=$(dirname "$SCRIPT")

source "$SCRIPTPATH/config_pipeline" || exit 1


# ==================================================
# 2. Help message
# ==================================================


show_usage() {
    cat <<'EOF'
POEM - Function Prediction

Usage:
  ./poem -s function-prediction -o <poem_results_dir> [-gf <proteins.faa>]

Arguments:
  -o <directory>  Root directory for POEM pipeline results
  -gf <file>      Protein sequences in FASTA format (optional)
  -h, --help      Show this help message

If you omit -gf, POEM uses the proteins from gene prediction under
the same results directory:
  <poem_results_dir>/gene_prediction_output/protein.translations.faa

Output:
  <poem_results_dir>/function_prediction_output/diamond_cog_reference.tsv

Example:
  ./poem -s function-prediction -o poem_stage_test
EOF
}


# ==================================================
# 3. Initialize variables
# ==================================================

gene_file=""
output_folder=""


# ==================================================
# 4. Parse command-line arguments
# ==================================================

while [[ $# -gt 0 ]]; do

    case "$1" in

        -gf|-o)
            if [ $# -lt 2 ] || [[ "$2" == -* ]]; then
                echo "Error: Missing value for $1"
                exit 1
            fi

            case "$1" in
                -gf) gene_file="$2" ;;
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
# 5. Validate inputs
# ==================================================

if [ -z "$gene_file" ]; then
    echo "Error: Predicted protein FASTA file is required (-gf)."
    exit 1
fi

if [ ! -s "$gene_file" ]; then
    echo "Error: Protein FASTA file not found or empty:"
    echo "  $gene_file"
    exit 1
fi

if [ -z "$output_folder" ]; then
    echo "Error: Output directory is required (-o)."
    exit 1
fi

# Reference COG sequences
cog_file="${cog_path}/cog-20.fa"

if [ ! -s "$cog_file" ]; then
    echo "Error: COG reference FASTA file not found or empty:"
    echo "  $cog_file"
    exit 1
fi


# ==================================================
# 6. Prepare output directory
# ==================================================

# Remove previous function prediction results
if [ -d "$output_folder" ]; then
    echo "Removing previous function prediction results..."
    rm -rf -- "$output_folder" || exit 1
fi

mkdir -p "$output_folder" || exit 1

echo ""
echo "POEM - Starting COG function prediction"
echo "Input proteins: $gene_file"
echo "COG reference: $cog_file"
echo "Output directory: $output_folder"
echo ""


# ==================================================
# 7. Build DIAMOND database
# ==================================================

echo "[1/2] Building the DIAMOND COG reference database."

diamond makedb \
    --in "$cog_file" \
    -d "$output_folder/cog_db" \
    --quiet || {
        echo "Error: DIAMOND database creation failed." >&2
        exit 1
    }


# ==================================================
# 8. Run DIAMOND BLASTP
# ==================================================

echo "[2/2] Running DIAMOND BLASTP."

diamond blastp \
    -d "$output_folder/cog_db" \
    -q "$gene_file" \
    -o "$output_folder/diamond_cog_reference.tsv" \
    --very-sensitive \
    --quiet || {
        echo "Error: DIAMOND BLASTP failed." >&2
        exit 1
    }


# ==================================================
# 9. Validate output
# ==================================================

diamond_output="${output_folder}/diamond_cog_reference.tsv"

if [ ! -f "$diamond_output" ]; then
    echo "Error: DIAMOND output file was not generated."
    exit 1
fi

echo ""
echo "Function prediction completed successfully."
echo "DIAMOND results: $diamond_output"
