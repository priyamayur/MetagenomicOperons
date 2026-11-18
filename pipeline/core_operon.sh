#! /bin/bash


source config
pipeline_config_path=$pipeline_config
source $pipeline_config_path

echo "visualize"
echo "start"
python ${script_path}/construct_predicted_operons.py -g ${core_operon_path}/leiden_core_graphs_5  -c2g ${operon_output_path}/cogId_to_gene_large.pkl  -o ${core_operon_path}
echo "end"