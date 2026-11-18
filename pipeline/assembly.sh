
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
POEM - Metagenomic Assembly

Usage:
  poem -s assembly -1 <R1.fastq> -2 <R2.fastq> -o <poem_results_dir> [-a <assembler>]

Required arguments:
  -1 <file>       Read 1 paired-end FASTQ file
  -2 <file>       Read 2 paired-end FASTQ file
  -o <directory>  Root directory for POEM pipeline results

Optional arguments:
  -a <assembler>  spades (metaSPAdes), megahit (MEGAHIT), or idba (IDBA-UD)
                  Default: spades
  -h, --help      Show this help message

Output:
  POEM creates assembly_output/ inside the results directory and writes
  the assembled contigs to:
    <poem_results_dir>/assembly_output/poem_contigs.fa

  Use the same results directory for gene prediction to pick up these
  contigs automatically.

Examples:
  poem -s assembly -1 R1.fq -2 R2.fq -o results/
  poem -s assembly -1 R1.fq -2 R2.fq -o results/ -a megahit
EOF
}


# ==================================================
# 3. Initialize variables
# ==================================================

assembly_method="spades"

assembly_file_1=""
assembly_file_2=""
output_folder=""


# ==================================================
# 4. Parse command-line arguments
# ==================================================

while [[ $# -gt 0 ]]; do

    case "$1" in

        -a|-1|-2|-o)
            if [ $# -lt 2 ] || [[ "$2" == -* ]]; then
                echo "Error: Missing value for $1"
                exit 1
            fi

            case "$1" in
                -a) assembly_method="$2" ;;
                -1) assembly_file_1="$2" ;;
                -2) assembly_file_2="$2" ;;
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

if [ -z "$assembly_file_1" ] || [ -z "$assembly_file_2" ]; then
    echo "Error: Both paired-end read files are required."
    exit 1
fi

if [ -z "$output_folder" ]; then
    echo "Error: Output directory is required."
    exit 1
fi

if [ ! -f "$assembly_file_1" ]; then
    echo "Error: Read 1 file not found: $assembly_file_1"
    exit 1
fi

if [ ! -f "$assembly_file_2" ]; then
    echo "Error: Read 2 file not found: $assembly_file_2"
    exit 1
fi

case "$assembly_method" in
    spades|megahit|idba)
        ;;
    *)
        echo "Error: Invalid assembler '$assembly_method'"
        echo "Supported assemblers: spades, megahit, idba"
        exit 1
        ;;
esac



# ==================================================
# 6. Prepare output directory
# ==================================================

output_method_folder="${output_folder}/${assembly_method}_output"

# Create main assembly directory if needed
mkdir -p "$output_folder" || exit 1

# Remove previous results for the selected assembler
if [ -d "$output_method_folder" ]; then
    echo "Removing previous $assembly_method assembly results..."
    rm -rf -- "$output_method_folder" || exit 1
fi

echo "Preparing fresh assembly output directory..."

echo ""
echo "POEM - Starting metagenomic assembly"
echo "Assembler: $assembly_method"
echo "Read 1: $assembly_file_1"
echo "Read 2: $assembly_file_2"
echo "Output directory: $output_folder"
echo ""



# ==================================================
# 7. Run assembly
# ==================================================

case "$assembly_method" in

    # ----------------------------------------------
    # MEGAHIT
    # ----------------------------------------------

    megahit)
        printf 'Running MEGAHIT. Log directory: %s\n' "$output_method_folder"

        if ! megahit \
            -1 "$assembly_file_1" \
            -2 "$assembly_file_2" \
            -o "$output_method_folder" >/dev/null 2>&1
        then
            printf 'Error: MEGAHIT failed. Check its log in %s if present.\n' \
                "$output_method_folder" >&2
            exit 1
        fi

        assembled_contigs="${output_method_folder}/final.contigs.fa"
        ;;

    # ----------------------------------------------
    # metaSPAdes
    # ----------------------------------------------

    spades)
        printf '[1/2] Running metaSPAdes. Log directory: %s\n' "$output_method_folder"

        if ! spades.py \
            -1 "$assembly_file_1" \
            -2 "$assembly_file_2" \
            --meta \
            -o "$output_method_folder" >/dev/null
        then
            printf 'Error: SPAdes failed. Check %s/spades.log if present.\n' \
                "$output_method_folder" >&2
            exit 1
        fi

        assembled_contigs="${output_method_folder}/contigs.fasta"
        ;;


    # ----------------------------------------------
    # IDBA-UD
    # ----------------------------------------------

    idba)
        # Keep the converted reads outside IDBA-UD's output directory.
        idba_reads="${output_folder}/idba_read12.fa"

        printf '[1/2] Converting paired-end reads to interleaved FASTA.\n'

        if ! fq2fa --merge \
            "$assembly_file_1" \
            "$assembly_file_2" \
            "$idba_reads" >/dev/null
        then
            printf 'Error: FASTQ conversion failed.\n' >&2
            exit 1
        fi

        printf '[2/2] Running IDBA-UD. Log directory: %s\n' "$output_method_folder"

        if ! idba_ud \
            -r "$idba_reads" \
            --num_threads 8 \
            -o "$output_method_folder" >/dev/null 2>&1
        then
            printf 'Error: IDBA-UD failed. Check its log in %s if present.\n' \
                "$output_method_folder" >&2
            exit 1
        fi

        assembled_contigs="${output_method_folder}/contig.fa"
        ;;

esac


# ==================================================
# 8. Validate and save assembled contigs
# ==================================================

if [ ! -s "$assembled_contigs" ]; then
    echo "Error: Assembly output is missing or empty:"
    echo "  $assembled_contigs"
    exit 1
fi

# Normalize FASTA headers for all assemblers
printf '[2/2] Normalizing contig FASTA headers.\n'

python "$script_path/normalize_fasta_headers.py" \
    -af "$assembled_contigs" \
    -o "$output_folder/poem_contigs.fa" \
    -map "$output_folder/contig_id_mapping.tsv" || {
        echo "Error: FASTA header normalization failed."
        exit 1
    }

# Validate normalized output
if [ ! -s "$output_folder/poem_contigs.fa" ]; then
    echo "Error: Normalized contigs file is missing or empty."
    exit 1
fi

echo ""
echo "Assembly completed successfully."
echo "Final contigs: $output_folder/poem_contigs.fa"