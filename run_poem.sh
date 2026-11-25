
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
POEM - Metagenomic Operon Prediction

Usage:
  poem -1 <R1.fastq> -2 <R2.fastq> -o <poem_results_dir> [-a assembler] [--model <name>]
  poem -af <contigs.fa> -o <poem_results_dir> [--model <name>]
  poem -af <contigs.fa> -gff <genes.gff> -pf <proteins.faa> -o <poem_results_dir> [--model <name>]
  poem --test 

Input options:
  -1 <file>       Read 1 paired-end FASTQ file
  -2 <file>       Read 2 paired-end FASTQ file
  -af <file>      Assembled contigs in FASTA format
  -gff <file>     Existing CDS annotations in GFF3 format
  -pf <file>      Existing protein translations in FASTA format

Required arguments:
  -o <directory>  Root directory for POEM pipeline results
                  (not required with --test)

Optional arguments:
  -a <assembler>  spades, megahit, or idba (default: spades)
                  Only applicable to paired-end reads
  --model <name>  Operon prediction classifier
                  Default: all_species (pooled classifier).
                  Species-specific classifiers can be selected
                  using their names in model_manifest.json.
  --test          Run POEM using the example dataset
  -h, --help      Show this help message

Pipeline:
  1. Assembly (if FASTQ files are provided)
  2. Gene prediction (skipped when -gff and -pf are supplied)
  3. Function prediction (DIAMOND BLASTP)
  4. Protein clustering (CD-HIT)
  5. Operon prediction and processing
  6. Core-operon identification

Examples:
  poem -1 R1.fq -2 R2.fq -o results/
  poem -1 R1.fq -2 R2.fq -a megahit -o results/
  poem -af contigs.fa -o results/
  poem -af contigs.fa -gff genes.gff -pf proteins.faa -o results/
  poem --test
EOF
}


# ==================================================
# 3. Initialize variables
# ==================================================

input_path1=""
input_path2=""
assembled_input=""
gff_input=""
protein_input=""
output_dir=""

assembler="spades"
assembler_specified=false
test_mode=false
operon_model="all_species"


# ==================================================
# 4. Parse command-line arguments
# ==================================================

if [ $# -eq 0 ]; then
    echo "Error: No arguments provided."
    show_usage
    exit 1
fi

while [[ $# -gt 0 ]]; do

    case "$1" in

        -1|-2|-af|-gff|-pf|-o|-a|--model)
            if [ $# -lt 2 ] || [[ "$2" == -* ]]; then
                echo "Error: Missing value for $1"
                exit 1
            fi

            case "$1" in
                -1)  input_path1="$2" ;;
                -2)  input_path2="$2" ;;
                -af) assembled_input="$2" ;;
                -gff) gff_input="$2" ;;
                -pf) protein_input="$2" ;;
                -o)  output_dir="$2" ;;
                --model) operon_model="$2" ;;
                -a)
                    assembler="$2"
                    assembler_specified=true
                    ;;
            esac

            shift 2
            ;;

        --test)
            test_mode=true
            shift
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
# 5. Determine input mode
# ==================================================

if [ "$test_mode" = true ]; then

    if [ -n "$input_path1" ] || [ -n "$input_path2" ] ||
       [ -n "$assembled_input" ] || [ -n "$gff_input" ] ||
       [ -n "$protein_input" ] || [ -n "$output_dir" ] ||
       [ "$assembler_specified" = true ]; then
        echo "Error: --test cannot be combined with other input options."
        exit 1
    fi

    echo "Running POEM in test mode..."

    source "$SCRIPTPATH/example_config" || exit 1

    input_mode="reads"

elif [ -n "$assembled_input" ]; then

    if [ -n "$input_path1" ] || [ -n "$input_path2" ]; then
        echo "Error: Provide either paired-end reads or assembled contigs, not both."
        exit 1
    fi

    if [ "$assembler_specified" = true ]; then
        echo "Error: -a is only valid when using paired-end reads."
        exit 1
    fi

    if [ -n "$gff_input" ] || [ -n "$protein_input" ]; then
        if [ -z "$gff_input" ] || [ -z "$protein_input" ]; then
            echo "Error: Option 3 requires BOTH -gff and -pf (with -af)." >&2
            exit 1
        fi
        input_mode="annotations"
    else
        input_mode="contigs"
    fi

else

    if [ -n "$gff_input" ] || [ -n "$protein_input" ]; then
        echo "Error: -gff and -pf require -af <contigs.fa>." >&2
        exit 1
    fi
    input_mode="reads"

fi


# ==================================================
# 6. Validate inputs
# ==================================================

if [ -z "$output_dir" ]; then
    echo "Error: POEM results directory is required (-o)."
    exit 1
fi

if [ "$input_mode" = "reads" ]; then

    if [ -z "$input_path1" ] || [ -z "$input_path2" ]; then
        echo "Error: Both paired-end reads (-1 and -2) are required."
        exit 1
    fi

    if [ ! -s "$input_path1" ]; then
        echo "Error: Read 1 file not found or empty: $input_path1"
        exit 1
    fi

    if [ ! -s "$input_path2" ]; then
        echo "Error: Read 2 file not found or empty: $input_path2"
        exit 1
    fi

    case "$assembler" in
        spades|megahit|idba) ;;
        *)
            echo "Error: Invalid assembler '$assembler'"
            echo "Supported assemblers: spades, megahit, idba"
            exit 1
            ;;
    esac

else

    if [ ! -s "$assembled_input" ]; then
        echo "Error: Assembled contigs file not found or empty:"
        echo "  $assembled_input"
        exit 1
    fi

    if [ "$input_mode" = "annotations" ]; then
        for input_file in "$gff_input" "$protein_input"; do
            if [ ! -s "$input_file" ]; then
                echo "Error: Annotation input file not found or empty: $input_file" >&2
                exit 1
            fi
        done
    fi

fi


# ==================================================
# 7. Prepare POEM results directory
# ==================================================

mkdir -p "$output_dir" || exit 1


# ==================================================
# 8. Display configuration
# ==================================================

echo ""
echo "=================================================="
echo "POEM - Metagenomic Operon Prediction"
echo "=================================================="
echo "Input mode: $input_mode"
echo "Results directory: $output_dir"
echo "Operon prediction model: $operon_model"

if [ "$input_mode" = "reads" ]; then
    echo "Read 1: $input_path1"
    echo "Read 2: $input_path2"
    echo "Assembler: $assembler"
else
    echo "Assembled contigs: $assembled_input"
    if [ "$input_mode" = "annotations" ]; then
        echo "Gene annotations: $gff_input"
        echo "Protein translations: $protein_input"
    fi
fi

echo ""


# ==================================================
# 9. Helper to execute POEM stages
# ==================================================

run_stage() {
    local stage="$1"
    local status

    shift

    "$SCRIPTPATH/run_poem_stage.sh" \
        -s "$stage" \
        -o "$output_dir" \
        "$@"

    status=$?

    if [ "$status" -ne 0 ]; then
        echo "Error: POEM stage '$stage' failed (exit code $status)."
        exit "$status"
    fi
}


# ==================================================
# 10. Assembly / Import assembled contigs
# ==================================================

if [ "$input_mode" = "reads" ]; then

    echo "------------------------ Step 1: Starting assembly using $assembler ------------------------"

    run_stage assembly \
        -1 "$input_path1" \
        -2 "$input_path2" \
        -a "$assembler"

    echo "------------------------ Step 1: Assembly completed successfully ------------------------"

elif [ "$input_mode" = "annotations" ]; then

    echo "------------------------ Step 1-2: Importing existing gene predictions ------------------------"

    assembly_folder="${output_dir}/assembly_output"
    gene_folder="${output_dir}/gene_prediction_output"
    mkdir -p "$assembly_folder" "$gene_folder" || exit 1

    # Protect external input if it resides inside the output folder.
    python "$SCRIPTPATH/scripts/normalize_fasta_headers.py" \
        -af "$assembled_input" \
        -o "$assembly_folder/poem_contigs.fa" \
        -map "$assembly_folder/contig_id_mapping.tsv" || exit 1

    python "$SCRIPTPATH/scripts/process_gene_gff.py" \
        -i "$gff_input" \
        -o "$gene_folder" \
        -af "$assembly_folder/poem_contigs.fa" \
        -map "$assembly_folder/contig_id_mapping.tsv" \
        -pf "$protein_input" || exit 1

    echo "Imported contigs, gene coordinates and protein sequences."

else

    echo "------------------------ Step 1: Using user-provided assembled contigs ------------------------"

    assembly_folder="${output_dir}/assembly_output"
    mkdir -p "$assembly_folder" || exit 1

    target_file="$assembly_folder/poem_contigs.fa"
    python "$SCRIPTPATH/scripts/normalize_fasta_headers.py" \
        -af "$assembled_input" \
        -o "$target_file" \
        -map "$assembly_folder/contig_id_mapping.tsv" || exit 1

    echo "Normalized contigs saved to: $target_file"

fi


# ==================================================
# 11. Gene prediction
# ==================================================

if [ "$input_mode" != "annotations" ]; then
    echo "------------------------ Step 2: Starting gene prediction ------------------------"
    run_stage gene-prediction
    echo "------------------------ Step 2: Gene prediction completed successfully ------------------------"
else
    echo "Step 2: Skipping gene prediction (user-supplied annotations)."
fi


# ==================================================
# 12. Function prediction
# ==================================================

echo "------------------------ Step 3: Starting function prediction ------------------------"

run_stage function-prediction

echo "------------------------ Step 3: Function prediction completed successfully ------------------------"


# ==================================================
# 13. CD-HIT protein clustering
# ==================================================

echo "------------------------ Step 4: Starting CD-HIT clustering ------------------------"

run_stage cdhit

echo "------------------------ Step 4: CD-HIT clustering completed successfully ------------------------"


# ==================================================
# 14. Operon prediction and processing
# ==================================================

echo "------------------------ Step 5: Starting operon prediction and processing ------------------------"

run_stage predict-operons --model "$operon_model"

echo "------------------------ Step 5: Operon prediction and processing completed successfully ------------------------"


# ==================================================
# 15. Core-operon identification
# ==================================================

echo "------------------------ Step 6: Starting core-operon identification ------------------------"

run_stage core-function

echo "------------------------ Step 6: Core-operon identification completed successfully ------------------------"

# ==================================================
# 16. Core-operon visualization
# ==================================================

echo "------------------------ Step 7: Starting core-operon visualization ------------------------"

if "$SCRIPTPATH/run_poem_stage.sh" \
    -s visualize-core-operon \
    -o "$output_dir"
then
    echo "Step 7: Core-operon visualization completed successfully."
else
    echo "Warning: Core-operon visualization failed."
    echo "All preceding POEM pipeline stages completed successfully."
fi

# ==================================================
# 17. Completion
# ==================================================

echo ""
echo "=================================================="
echo "POEM pipeline completed successfully."
echo "=================================================="
echo "Results directory: $output_dir"
echo ""
echo "To visualize the predicted core operons, run:"
echo "  ./poem -s visualize-core-operon -o \"$output_dir\""
