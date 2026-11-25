
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
POEM - Operon Prediction and Processing

Usage:
  ./poem -s predict-operons -o <poem_results_dir> [--model <name>]

Arguments:
  -o <directory>  Root directory for POEM pipeline results
  --model <name>  Operon prediction classifier (optional).
                  Default: all_species (pooled classifier).
                  Select a species-specific classifier using its name from the model manifest file.

  -h, --help      Show this help message

Use the same results directory as the earlier stages. POEM predicts
operons from the assembled contigs and gene coordinates, then processes
the predictions using the DIAMOND and CD-HIT results.

Output:
  <poem_results_dir>/operon_prediction_output/

Example:
  ./poem -s predict-operons -o poem_stage_test
EOF
}


# ==================================================
# 3. Initialize variables
# ==================================================

gene_coordinates=""
assembled_file=""
output_folder=""
operon_model="all_species"


# ==================================================
# 4. Parse command-line arguments
# ==================================================

while [[ $# -gt 0 ]]; do

    case "$1" in

        -gc|-af|-o|--model)
            if [ $# -lt 2 ] || [[ "$2" == -* ]]; then
                echo "Error: Missing value for $1"
                exit 1
            fi

            case "$1" in
                -gc) gene_coordinates="$2" ;;
                -af) assembled_file="$2" ;;
                -o) output_folder="$2" ;;
                --model) operon_model="$2" ;;
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

if [ -z "$gene_coordinates" ]; then
    echo "Error: Gene coordinates file is required (-gc)."
    exit 1
fi

if [ ! -s "$gene_coordinates" ]; then
    echo "Error: Gene coordinates file not found or empty:"
    echo "  $gene_coordinates"
    exit 1
fi

if [ -z "$assembled_file" ]; then
    echo "Error: Assembled contigs file is required (-af)."
    exit 1
fi

if [ ! -s "$assembled_file" ]; then
    echo "Error: Assembled contigs FASTA file not found or empty:"
    echo "  $assembled_file"
    exit 1
fi

if [ -z "$output_folder" ]; then
    echo "Error: Output directory is required (-o)."
    exit 1
fi

# ==================================================
# 6. Prepare output directory
# ==================================================

# Preserve the existing output directory if it exists.
# Prediction and processing use this same directory.

mkdir -p "$output_folder" || exit 1

echo ""
echo "POEM - Starting operon prediction"
echo "Gene coordinates: $gene_coordinates"
echo "Assembled contigs: $assembled_file"
echo "Prediction model: $operon_model"
echo "Output directory: $output_folder"
echo ""


# ==================================================
# 7. Predict operons
# ==================================================

python "$script_path/predict_operons.py" \
    -g "$gene_coordinates" \
    --model-manifest "${operon_model_path}/model_manifest.json" \
    --model-name "$operon_model" \
    -af "$assembled_file" \
    -o "$output_folder" || {
        status=$?
        echo "Error: Operon prediction failed." >&2
        exit "$status"
    }


# ==================================================
# 8. Completion
# ==================================================

echo ""
echo "Operon prediction completed successfully."
