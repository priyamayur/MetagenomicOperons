
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
POEM - Core Operon Visualization

Usage:
  ./poem -s visualize-core-operon -o <poem_results_dir>

Arguments:
  -o <directory>  Root directory for POEM pipeline results
  -h, --help      Show this help message

Input:
  POEM uses the results generated during core-operon identification.

  Required files:
    <poem_results_dir>/core_operon_output/leiden_communities.pkl
    <poem_results_dir>/core_operon_output/path_quasi_hybrid_cores.csv

Visualization:
  Launches an interactive visualization of Leiden communities
  and predicted core operons.

Example:
  ./poem -s visualize-core-operon -o poem_stage_test
EOF
}


# ==================================================
# 3. Initialize variables
# ==================================================

core_operon_path=""


# ==================================================
# 4. Parse command-line arguments
# ==================================================

while [[ $# -gt 0 ]]; do

    case "$1" in

        -coo)
            if [ $# -lt 2 ] || [[ "$2" == -* ]]; then
                echo "Error: Missing value for $1"
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

if [ -z "$core_operon_path" ]; then
    echo "Error: Core-operon directory is required (-coo)."
    exit 1
fi

leiden_graph="${core_operon_path}/leiden_communities.pkl"
core_file="${core_operon_path}/path_quasi_hybrid_cores.csv"
cog_def_file="${cog_path}/cog-20.def.tab"

for file in "$leiden_graph" "$core_file" "$cog_def_file"; do
    if [ ! -s "$file" ]; then
        echo "Error: Required file not found or empty:"
        echo "  $file"
        exit 1
    fi
done


# ==================================================
# 7. Launch visualization
# ==================================================

echo ""
echo "POEM - Starting core-operon visualization"
echo "Core-operon directory: $core_operon_path"
echo ""

python "$script_path/visualize_core_operon.py" \
    -g "$leiden_graph" \
    -cores "$core_file" \
    -cdef "$cog_def_file" \
    -o "$core_operon_path" || {
        echo "Error: Core-operon visualization failed." >&2
        exit 1
    }
