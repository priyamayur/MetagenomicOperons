
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
POEM - CD-HIT Protein Clustering

Usage:
  ./poem -s cdhit -o <poem_results_dir> [-gf <proteins.faa>]

Arguments:
  -o <directory>  Root directory for POEM pipeline results
  -gf <file>      Protein sequences in FASTA format (optional)
  -h, --help      Show this help message

If you omit -gf, POEM uses the proteins from gene prediction under
the same results directory:
  <poem_results_dir>/gene_prediction_output/protein.translations.faa

The sequence identity threshold is set by cdhit_threshold in config_pipeline.

CD-HIT outputs:
  <poem_results_dir>/cdhit_output/clustered_proteins.fasta
  <poem_results_dir>/cdhit_output/clustered_proteins.fasta.clstr

Example:
  ./poem -s cdhit -o poem_stage_test
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
    echo "Error: Protein FASTA file is required (-gf)."
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

if [ -z "$cdhit_threshold" ]; then
    echo "Error: cdhit_threshold is not defined in config_pipeline."
    exit 1
fi


# ==================================================
# 6. Prepare output directory
# ==================================================

# Remove previous CD-HIT results
if [ -d "$output_folder" ]; then
    echo "Removing previous CD-HIT results..."
    rm -rf -- "$output_folder" || exit 1
fi

mkdir -p "$output_folder" || exit 1

cdhit_output="${output_folder}/clustered_proteins.fasta"

echo ""
echo "POEM - Starting CD-HIT protein clustering"
echo "Input proteins: $gene_file"
echo "Sequence identity threshold: $cdhit_threshold"
echo "Output directory: $output_folder"
echo ""


# ==================================================
# 7. Run CD-HIT
# ==================================================

cdhit_log="${output_folder}/cdhit.log"
printf '[1/2] Running CD-HIT. Log: %s\n' "$cdhit_log"

cd-hit \
    -i "$gene_file" \
    -o "$cdhit_output" \
    -c "$cdhit_threshold" >"$cdhit_log" 2>&1 || {
        printf 'Error: CD-HIT clustering failed.\n' >&2
        tail -n 15 -- "$cdhit_log" >&2
        exit 1
    }

# ==================================================
# 8. Post-process CD-HIT output
# ==================================================

printf '[2/2] Post-processing CD-HIT clusters...\n'

cluster_file="${cdhit_output}.clstr"

if [ ! -s "$cluster_file" ]; then
    printf 'Error: CD-HIT cluster file not found or empty.\n' >&2
    exit 1
fi

python "$script_path/post_process_cdhit.py" \
    -i "$cluster_file" \
    -cf "$gene_file" \
    -o "$output_folder" || {
        printf 'Error: CD-HIT post-processing failed.\n' >&2
        exit 1
    }


# ==================================================
# 9. Completion
# ==================================================

echo ""

echo "Output directory: $output_folder"
