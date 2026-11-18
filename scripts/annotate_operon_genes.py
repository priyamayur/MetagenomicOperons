import sys
import argparse
from Bio import SeqIO
import pandas as pd
import pickle
import itertools

parser = argparse.ArgumentParser()
parser.add_argument("-og", "--operon_genes", required=True, help="operon gene pairs")
parser.add_argument("-gc", "--gene_to_cluster", required=True, help="gene_to_cluster")
parser.add_argument("-gcid", "--gene_to_cogId", required=True, help="cluster_to_gene")
parser.add_argument("-o", "--output", required=True, default="output", help="output directory name")
# parser.add_argument("-ut", "--upperthreshold", required=False, default="0.80", help="set upper threshold value")

def get_cog_pairs(filtered_gene_pairs, gene_to_cogId):

  cogIds_to_operon = []
  for i in range(len(filtered_gene_pairs)):

    gene_1 = filtered_gene_pairs[i][0]
    gene_2 = filtered_gene_pairs[i][1]
    operon_id = 1

    # get all the cogs associated with gene 1
    gene1_cogs = []
    if gene_1 in gene_to_cogId.keys():
      gene_1_cogId = gene_to_cogId[gene_1]
      if (gene_1_cogId != gene_1_cogId or gene_1_cogId == ""):
        continue
      if ";" in gene_1_cogId:
        gene1_cogs.extend(gene_1_cogId.split(";"))
      else:
        gene1_cogs.append(gene_1_cogId)

    # get all the cogs associated with gene 2
    gene2_cogs = []
    if gene_2 in gene_to_cogId.keys():
      gene_2_cogId = gene_to_cogId[gene_2]
      if (gene_2_cogId != gene_2_cogId or gene_2_cogId == ""):
        continue
      if ";" in gene_2_cogId:
        gene2_cogs.extend(gene_2_cogId.split(";"))
      else:
        gene2_cogs.append(gene_2_cogId)

    # get all the combinations of gene 1 and gene 2 cogs
    combinations = list(itertools.product(gene1_cogs, gene2_cogs))
    for combination in combinations:
      gene_1_cogId = combination[0]
      gene_2_cogId = combination[1]

      if gene_1_cogId == gene_2_cogId:
        continue
      cogIds_to_operon.append([gene_1_cogId,gene_2_cogId,operon_id])

  cogIds_to_operon_df = pd.DataFrame(cogIds_to_operon, columns=['cog1', 'cog2', 'operon'])

  return cogIds_to_operon_df

def get_filtered_gene_pairs(transitive_gene_pairs, gene_to_clusters):
    gene_cluster_dict = {}
    for gene_pair in transitive_gene_pairs:  # see which is better, only gene-pairs or transitive
        # print(gene_pair)
        gene_cluster_pair = []
        for gene in gene_pair:
            cluster_name = gene_to_clusters[gene]
            gene_cluster_pair.append(cluster_name)
        single_gene_cluster_str = '_'.join(gene_cluster_pair)
        if single_gene_cluster_str in gene_cluster_dict.keys():
            gene_cluster_dict[single_gene_cluster_str].append(gene_pair)
        else:
            gene_cluster_dict[single_gene_cluster_str] = [gene_pair]
    count = 0
    filtered_gene_pairs = []
    for clusters in gene_cluster_dict.keys():
        gene_pair = (gene_cluster_dict[clusters])
        first_gene_pair = list(gene_pair[0])
        first_gene_pair.append(1)
        # print(first_gene_pair)
        filtered_gene_pairs.append(first_gene_pair)
        if len(gene_pair) > 1:
            count += 1
            # print(gene_pair)

    return filtered_gene_pairs

def get_transitive_gene_pairs(gene_pairs_merged_dataset):
  transitive_gene_pairs = []
  operon_class = 1
  for pointer in range(len(gene_pairs_merged_dataset)):
    combinations = list(itertools.combinations(gene_pairs_merged_dataset[pointer], 2))
    combinations_list = [[com[0], com[1]] for com in combinations]

    transitive_gene_pairs.extend(combinations_list)
  # transitive_gene_pairs_df = pd.DataFrame(transitive_gene_pairs, columns=['gene1', 'gene2', 'operon'])
  return transitive_gene_pairs

def merge_gene_pairs_to_operon(gene_names_pairs_to_operon_df, col1 = 'gene1', col2 = 'gene2'):
  gene_names_pairs_to_operon_filter = gene_names_pairs_to_operon_df[gene_names_pairs_to_operon_df['operon'] >= 1]
  dictionary_1 = gene_names_pairs_to_operon_filter.groupby(col1)[col2].apply(list).to_dict()
  dictionary_2 = gene_names_pairs_to_operon_filter.groupby(col2)[col1].apply(list).to_dict()
  gene_names_pairs_to_operon_filter_dataset = gene_names_pairs_to_operon_filter.values

  gene_pairs_to_operon = []
  for i in range(len(gene_names_pairs_to_operon_filter_dataset)):
    gene_1_i = gene_names_pairs_to_operon_filter_dataset[i][0]
    gene_2_i = gene_names_pairs_to_operon_filter_dataset[i][1]
    operon_found = {gene_1_i, gene_2_i}
    operon_checking_set = {gene_1_i, gene_2_i}
    operon_seen = []
    while len(operon_checking_set) > 0:
      # print("start=",operon_checking_set)
      check_gene = operon_checking_set.pop()
      # print("check gene=",check_gene)
      if check_gene in operon_seen:
        continue
      operon_seen.append(check_gene)
      if check_gene in dictionary_1.keys():
        for gene in dictionary_1[check_gene]:
          operon_checking_set.add(gene)
          operon_found.add(gene)
      if check_gene in dictionary_2.keys():
        for gene in dictionary_2[check_gene]:
          operon_checking_set.add(gene)
          operon_found.add(gene)
    gene_pairs_to_operon.append(list(operon_found))
  # Remove duplicates using set and map
  gene_pairs_to_operon = list(map(list, set(map(lambda x: tuple(sorted(x)), gene_pairs_to_operon))))
  return gene_pairs_to_operon


def main():
    #print("Starting annotation of operon gene pairs with COG IDs")
    args = parser.parse_args()
    operon_genes = args.operon_genes
    gene_to_cluster = args.gene_to_cluster
    gene_to_cogId=args.gene_to_cogId
    output_path = args.output

    gene_names_pairs_to_operon_prediction_df = pd.read_csv(
        operon_genes)
    # print("Gene pairs with operon predictions loaded. Sample:")
    # print(gene_names_pairs_to_operon_prediction_df.head())
    gene_names_pairs_to_operon_prediction_df = gene_names_pairs_to_operon_prediction_df[['gene1', 'gene2', 'operon']]


    with open(gene_to_cluster, 'rb') as file:
        gene_to_cluster = pickle.load(file)

    # print(gene_to_cluster.keys())
    with open(gene_to_cogId, 'rb') as file:
        gene_to_cogId = pickle.load(file)

    gene_pairs_merged = merge_gene_pairs_to_operon(gene_names_pairs_to_operon_prediction_df)
    transitive_gene_pairs = get_transitive_gene_pairs(gene_pairs_merged)
    # print(list(gene_to_cluster.keys())[0:5])
    filtered_gene_pairs = get_filtered_gene_pairs(transitive_gene_pairs, gene_to_cluster)
    # filtered_gene_pairs = get_filtered_gene_pairs(gene_pairs_merged, gene_to_cluster)
    cog_pairs_df = get_cog_pairs(filtered_gene_pairs,gene_to_cogId)
    # print("Total COG pairs:", len(cog_pairs_df))
    with open(
            output_path+"/cog_pairs_df",
            'wb') as file:
        pickle.dump(cog_pairs_df, file)
    # print("COG pairs dataframe saved to:", output_path+"/cog_pairs_df")    

if __name__ == "__main__":
    main()