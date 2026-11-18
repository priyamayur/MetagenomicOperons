#!/bin/bash

# ==================================================
# 1. Setup
# ==================================================

SCRIPT=$(realpath "$0")
SCRIPTPATH=$(dirname "$SCRIPT")

# ==================================================
# 2. General help
# ==================================================

show_usage() {
    cat <<EOF
POEM - Metagenomic Operon Prediction

Usage:
  poem -s <stage> -o <output_dir> [options]
  poem -s <stage> -h

Available stages:
  assembly                Assemble paired-end reads
  gene-prediction         Predict genes from assembled contigs
  function-prediction     Assign COG functions to genes
  cdhit                   Cluster predicted protein sequences
  predict-operons         Predict operon boundaries
  core-function           Identify core operons
  visualize-core-operon   Visualize Leiden communities and core operons

Common arguments:
  -s, --stage <stage>      Pipeline stage to execute
  -o <directory>           Root directory for POEM pipeline results
  -h, --help               Show help message

For stage-specific help:
  poem -s assembly -h
  poem -s gene-prediction -h
  poem -s core-function -h

Examples:
  poem -s assembly -1 R1.fq -2 R2.fq -o results/
  poem -s gene-prediction -o results/
  poem -s predict-operons -o results/
EOF
}

# ==================================================
# 3. Map stage names to shell scripts
# ==================================================

get_stage_script() {
    case "$1" in
        assembly)              echo "assembly.sh" ;;
        gene-prediction)       echo "gene_prediction.sh" ;;
        function-prediction)   echo "function_prediction.sh" ;;
        cdhit)                 echo "cdhit.sh" ;;
        predict-operons)       echo "predict_operons.sh" ;;
        process-operons)       echo "process_operons.sh" ;;
        core-function)         echo "core_function.sh" ;;
        visualize-core-operon) echo "visualize_core_operon.sh" ;;
        *) return 1 ;;
    esac
}

# ==================================================
# 4. Initialize variables
# ==================================================

stage_value=""
output_folder=""
assembly_method="spades"

read_file_1=""
read_file_2=""
assembled_file=""
gene_file=""
gene_coordinates=""

function_prediction_output=""
operon_prediction_output=""
cdhit_prediction_output=""
core_operon_output=""

help_requested=false

# ==================================================
# 5. Parse command-line arguments
# ==================================================

if [ $# -eq 0 ]; then
    echo "Error: No arguments provided."
    show_usage
    exit 1
fi

while [[ $# -gt 0 ]]; do
    case "$1" in

        -h|--help)
            help_requested=true
            shift
            ;;

        -s|--stage|-a|-1|-2|-o|-af|-gf|-gc|-fpo|-opo|-cdo|-coo)

            if [ $# -lt 2 ] || [[ "$2" == -* ]]; then
              if [[ "$1" == "-s" || "$1" == "--stage" ]] &&
                [[ "$2" == "-h" || "$2" == "--help" ]]; then
                  show_usage
                  exit 0
              fi

              echo "Error: Missing value for $1"
              exit 1
            fi

            case "$1" in
                -s|--stage) stage_value="$2" ;;
                -a) assembly_method="$2" ;;
                -1) read_file_1="$2" ;;
                -2) read_file_2="$2" ;;
                -o) output_folder="$2" ;;
                -af) assembled_file="$2" ;;
                -gf) gene_file="$2" ;;
                -gc) gene_coordinates="$2" ;;
                -fpo) function_prediction_output="$2" ;;
                -opo) operon_prediction_output="$2" ;;
                -cdo) cdhit_prediction_output="$2" ;;
                -coo) core_operon_output="$2" ;;
            esac

            shift 2
            ;;

        *)
            echo "Error: Unknown argument $1"
            show_usage
            exit 1
            ;;
    esac
done

# ==================================================
# 6. Validate stage name
# ==================================================

if [ -z "$stage_value" ]; then
    if [ "$help_requested" = true ]; then
        show_usage
        exit 0
    fi

    echo "Error: Missing pipeline stage (-s)."
    show_usage
    exit 1
fi

stage_script=$(get_stage_script "$stage_value") || {
    echo "Error: Unknown pipeline stage '$stage_value'"
    show_usage
    exit 1
}

# ==================================================
# 7. Load configuration
# ==================================================

source "$SCRIPTPATH/config" || exit 1

if [ -z "$pipeline_dir" ]; then
    echo "Error: pipeline_dir is not defined in config."
    exit 1
fi

# ==================================================
# 8. Display stage-specific help
# ==================================================

if [ "$help_requested" = true ]; then
    "$pipeline_dir/$stage_script" -h
    exit $?
fi

# ==================================================
# 9. Validate output directory
# ==================================================

if [ -z "$output_folder" ]; then
    echo "Error: Missing output directory (-o)."
    echo "For help: poem -s $stage_value -h"
    exit 1
fi

# ==================================================
# 10. Activate Conda environment
# ==================================================

eval "$(conda shell.bash hook)" || exit 1
conda activate con_env || exit 1


# ==================================================
# 11. Define default output paths
# ==================================================

assembly_folder="${output_folder}/assembly_output"
gene_prediction_folder="${output_folder}/gene_prediction_output"
function_prediction_folder="${output_folder}/function_prediction_output"
cdhit_folder="${output_folder}/cdhit_output"
operon_prediction_folder="${output_folder}/operon_prediction_output"
core_operon_folder="${output_folder}/core_operon_output"

mkdir -p "$output_folder" || exit 1

echo ""
echo "POEM - Running stage: $stage_value"
echo "Output directory: $output_folder"
echo ""


# ==================================================
# 12. Execute selected pipeline stage
# ==================================================

case "$stage_value" in

    # ----------------------------------------------
    # Assembly
    # ----------------------------------------------

    assembly)

        if [ -z "$read_file_1" ] || [ -z "$read_file_2" ]; then
            echo "Error: Assembly requires -1 and -2."
            echo "For help: poem -s assembly -h"
            exit 1
        fi

        case "$assembly_method" in
            spades|megahit|idba) ;;
            *)
                echo "Error: Invalid assembler '$assembly_method'"
                echo "Supported assemblers: spades, megahit, idba"
                exit 1
                ;;
        esac

        "$pipeline_dir/assembly.sh" \
            -a "$assembly_method" \
            -1 "$read_file_1" \
            -2 "$read_file_2" \
            -o "$assembly_folder"
        ;;

    # ----------------------------------------------
    # Gene prediction
    # ----------------------------------------------

    gene-prediction)

        assembled_file="${assembled_file:-${assembly_folder}/poem_contigs.fa}"

        "$pipeline_dir/gene_prediction.sh" \
            -af "$assembled_file" \
            -o "$gene_prediction_folder"
        ;;

    # ----------------------------------------------
    # Function prediction
    # ----------------------------------------------

    function-prediction)

        gene_file="${gene_file:-${gene_prediction_folder}/protein.translations.faa}"

        "$pipeline_dir/function_prediction.sh" \
            -gf "$gene_file" \
            -o "$function_prediction_folder"
        ;;

    # ----------------------------------------------
    # CD-HIT clustering
    # ----------------------------------------------

    cdhit)

        if [ -z "$gene_file" ]; then
            gene_file="${gene_prediction_folder}/protein.translations.faa"
        else
            printf 'Using user-provided protein FASTA: %s\n' "$gene_file"
        fi

        echo "Input: $gene_file"

        "$pipeline_dir/cdhit.sh" \
            -gf "$gene_file" \
            -o "$cdhit_folder"
        ;;
    # ----------------------------------------------
    # Predict operons
    # ----------------------------------------------

    predict-operons)

        gene_coordinates="${gene_coordinates:-${gene_prediction_folder}/gene_coordinates.csv}"
        assembled_file="${assembled_file:-${assembly_folder}/poem_contigs.fa}"

        echo "Starting operon prediction..."
        echo "Gene coordinates: $gene_coordinates"
        echo "Assembled contigs: $assembled_file"

        printf '[1/2] Predicting operons.\n'
        # Step 1: Predict operons
        "$pipeline_dir/predict_operons.sh" \
            -gc "$gene_coordinates" \
            -af "$assembled_file" \
            -o "$operon_prediction_folder" || exit 1

        
        printf '[2/2] Processing predicted operons.\n'

        # Step 2: Process predicted operons
        "$pipeline_dir/process_operons.sh" \
            -fpo "$function_prediction_folder" \
            -opo "$operon_prediction_folder" \
            -cdo "$cdhit_folder" || exit 1

        echo "Operon prediction and processing completed successfully."
        ;;

    # ----------------------------------------------
    # Core operon identification
    # ----------------------------------------------

    core-function)

        operon_prediction_output="${operon_prediction_output:-$operon_prediction_folder}"

        "$pipeline_dir/core_function.sh" \
            -opo "$operon_prediction_output" \
            -o "$core_operon_folder"
        ;;

    # ----------------------------------------------
    # Core operon visualization
    # ----------------------------------------------

    visualize-core-operon)

        core_operon_output="${core_operon_output:-$core_operon_folder}"

        "$pipeline_dir/visualize_core_operon.sh" \
            -coo "$core_operon_output"
        ;;

    esac


# ==================================================
# 13. Check execution status
# ==================================================

status=$?

if [ "$status" -ne 0 ]; then
    echo "Error: Stage '$stage_value' failed (exit code $status)."
    exit "$status"
fi

echo "Stage '$stage_value' completed successfully."