#! /bin/bash

echo "checking"
conda create -y -n con_env
#eval "$(conda shell.bash hook)"
eval "$(conda shell.bash hook)"
conda activate con_env

conda install -y python=3.12.0
source config

mamba install  -y -c bioconda megahit spades idba prodigal cd-hit blast diamond biopython quast

mamba install -y -c conda-forge pandas scikit-learn tqdm numba
mamba install -y -c conda-forge leidenalg matplotlib igraph networkx dash_cytoscape

#conda install  -y -c bioconda metawrap
#metawrap assembly -h

echo "Download COG database"

#cog_directory="data/cog"
if [ -z "$(ls -A $cog_directory)" ]; then
    echo "Directory is empty. Proceeding with wget..."
    wget -q -c ftp://ftp.ncbi.nih.gov/pub/COG/COG2020/data/*20* -P $cog_path
    gunzip $cog_path/cog-20.fa.gz
else
    echo "Directory is not empty. Skipping wget cog."
fi

conda deactivate

#mamba uninstall bioconda::megahit



