import sys
import argparse
from Bio import SeqIO
import pandas as pd
import pickle

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", required=True, help="clustered proteins")
parser.add_argument("-cf", "--clusteredfile", required=True, help="clustered file")
parser.add_argument("-o", "--output", required=True, default="output", help="output directory name")
# parser.add_argument("-ut", "--upperthreshold", required=False, default="0.80", help="set upper threshold value")


def get_full_gene_name(partial_gene_name, clustered_file_sequences):
  final_gene_full_name = ''
  for gene_name in clustered_file_sequences.keys():
    if partial_gene_name in gene_name:
      final_gene_full_name = gene_name
      break
  return final_gene_full_name
  
def find_clusters(records, clustered_file_sequences):
  cluster_to_genes = {}
  gene_to_clusters ={}
  for record in records:
    # print(record)
    if record.strip().startswith(">"):
      cluster_name = record.strip().split(">")[1]
      
      cluster_to_genes[cluster_name] = []
    else:
      gene_name = record.split(">",1)[1].strip().split("...",1)[0]
      # print(record)
      # print(gene_name)
      full_gene_name = get_full_gene_name(gene_name, clustered_file_sequences)
      # print(gene_name, "-->", full_gene_name)
      cluster_to_genes[cluster_name].append(full_gene_name)
      gene_to_clusters[full_gene_name] = cluster_name
  return cluster_to_genes, gene_to_clusters

def main():
    args = parser.parse_args()
    clustered_proteins_file = args.input
    clustered_file = args.clusteredfile
    output_path = args.output
    
    clustered_file_sequences = SeqIO.to_dict(SeqIO.parse(
        clustered_file,
        'fasta'))
    # print(clustered_file_sequences.keys())
    with open(clustered_proteins_file
            ) as handle:
        clustered_proteins = handle.readlines()
    cluster_to_genes, gene_to_clusters = find_clusters(clustered_proteins, clustered_file_sequences)
    # print(cluster_to_genes)
    # print(gene_to_clusters)
    with open(
            output_path+ '/cluster_to_genes','wb') as file:
        pickle.dump(cluster_to_genes, file)
    with open(
            output_path+ '/gene_to_clusters', 'wb') as file:
        pickle.dump(gene_to_clusters, file)

if __name__ == "__main__":
    main()