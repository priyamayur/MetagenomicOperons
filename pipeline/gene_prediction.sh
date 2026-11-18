
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
POEM - Gene Prediction

Usage:
  poem -s gene-prediction -o <poem_results_dir> [-af <contigs.fa>]

Arguments:
  -o <directory>  Root directory for POEM pipeline results.
  -af <file>      Assembled nucleotide sequences in FASTA format.
  -h, --help      Show this help message.

If you omit -af, use the same -o directory you used for POEM assembly.
POEM reads the assembled contigs from:
  <poem_results_dir>/assembly_output/poem_contigs.fa

If you supply -af, you can use a new POEM results directory.
In both cases, POEM writes the gene prediction results to:
  <poem_results_dir>/gene_prediction_output/

Examples:
  poem -s gene-prediction -o results/
  poem -s gene-prediction -af my_contigs.fa -o gene_check/
EOF
}


# ==================================================
# 3. Initialize variables
# ==================================================

assembled_contigs=""
output_folder=""


# ==================================================
# 4. Parse command-line arguments
# ==================================================

while [[ $# -gt 0 ]]; do

    case "$1" in

        -af|-o)
            if [ $# -lt 2 ] || [[ "$2" == -* ]]; then
                echo "Error: Missing value for $1"
                exit 1
            fi

            case "$1" in
                -af) assembled_contigs="$2" ;;
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

if [ -z "$assembled_contigs" ]; then
    echo "Error: Assembled contigs FASTA file is required (-af)."
    exit 1
fi

if [ ! -s "$assembled_contigs" ]; then
    echo "Error: Assembled contigs file not found or empty:"
    echo "  $assembled_contigs"
    exit 1
fi

if [ -z "$output_folder" ]; then
    echo "Error: Output directory is required (-o)."
    exit 1
fi


# ==================================================
# 6. Prepare output directory
# ==================================================

# Remove previous gene prediction results
if [ -d "$output_folder" ]; then
    echo "Removing previous gene prediction results..."
    rm -rf -- "$output_folder" || exit 1
fi

mkdir -p "$output_folder" || exit 1

echo ""
echo "POEM - Starting gene prediction"
echo "Input contigs: $assembled_contigs"
echo "Output directory: $output_folder"
echo ""


# ==================================================
# 7. Run Prodigal
# ==================================================

printf '[1/2] Running Prodigal gene prediction.\n'

prodigal \
    -i "$assembled_contigs" \
    -f gff \
    -o "$output_folder/gene.coords.gff" \
    -a "$output_folder/protein.translations.faa" \
    -d "$output_folder/dna.translations.faa" \
    -p meta \
    -q || {
        printf 'Error: Prodigal gene prediction failed.\n' >&2
        exit 1
    }

# ==================================================
# 8. Post-process Prodigal output
# ==================================================

printf '[2/2] Processing GFF gene coordinates.\n'

python "$script_path/process_gene_gff.py" \
    -i "$output_folder/gene.coords.gff" \
    -o "$output_folder" || {
        echo "Error: GFF post-processing failed."
        exit 1
    }

# ==================================================
# 9. Validate outputs
# ==================================================

if [ ! -s "$output_folder/protein.translations.faa" ]; then
    echo "Error: Predicted protein FASTA is missing or empty."
    exit 1
fi

if [ ! -s "$output_folder/gene_coordinates.csv" ]; then
    echo "Error: Gene coordinates CSV is missing or empty."
    exit 1
fi

echo ""
echo "Gene prediction completed successfully."
echo "Output directory: $output_folder"
