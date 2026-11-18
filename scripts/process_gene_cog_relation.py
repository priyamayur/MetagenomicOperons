import sys
import argparse
from Bio import SeqIO
import pandas as pd
import pickle

parser = argparse.ArgumentParser()
parser.add_argument("-cp", "--cog_prot", required=True, help="cog protein file")
parser.add_argument("-ca", "--cog_align", required=True, help="cog alignment file")
parser.add_argument("-o", "--output", required=True, default="output", help="output directory name")
# parser.add_argument("-ut", "--upperthreshold", required=False, default="0.80", help="set upper threshold value")

def get_gene_to_cog_id(gene_matches, query_to_protId_dict, protId_to_cog_dict):
    queries = list(set(gene_matches['query'].tolist()))
    gene_to_cogId = {}
    count = 0
    for query in queries:
        count += 1
        cogIds = []
        if query in query_to_protId_dict.keys():
            protIds = query_to_protId_dict[query]
            for protId in protIds:
                if protId in protId_to_cog_dict.keys():
                    cogId = protId_to_cog_dict[protId]
                    cogIds.append(cogId)
            cogIds = list(set(cogIds))

            cogId = ';'.join(cogIds)
            gene_to_cogId[query] = cogId
    return gene_to_cogId

def get_cogId_to_gene(gene_matches, query_to_protId_dict, protId_to_cog_dict):
    queries = list(set(gene_matches['query'].tolist()))
    cogId_to_gene = {}
    count = 0
    for query in queries:
        count += 1
        cogIds = []
        if query in query_to_protId_dict.keys():
            protIds = query_to_protId_dict[query]
            for protId in protIds:
                if protId in protId_to_cog_dict.keys():
                    cogId = protId_to_cog_dict[protId]
                    cogIds.append(cogId)
                    if cogId in cogId_to_gene.keys():
                        cogId_to_gene[cogId].append(query)
                        cogId_to_gene[cogId] = list(set(cogId_to_gene[cogId]))
                    else:
                        cogId_to_gene[cogId] = [query]

    return cogId_to_gene

def get_gene_cog_annotations(cog_protein, cog_alignment):
    cog_protein['proteinID'] = [x.replace("_", "") for x in cog_protein['proteinID']]
    cog_protein['proteinID'] = [x.replace(".", "") for x in cog_protein['proteinID']]
    protId_to_cog_dict = cog_protein.set_index('proteinID')['cogID'].to_dict()

    cog_alignment['target'] = [x.replace("_", "") for x in cog_alignment['target']]
    cog_alignment['target'] = [x.replace(".", "") for x in cog_alignment['target']]

    sequence_identity_cutoff = 85
    gene_matches_filtered = cog_alignment[cog_alignment['sequence_identity'] >= sequence_identity_cutoff]
    query_to_protId_dict = gene_matches_filtered.groupby('query')['target'].apply(list).to_dict()
    gene_matches_filtered.head()

    gene_to_cogId = get_gene_to_cog_id(cog_alignment,query_to_protId_dict,protId_to_cog_dict)
    cogId_to_gene = get_cogId_to_gene(cog_alignment,query_to_protId_dict,protId_to_cog_dict)

    return gene_to_cogId, cogId_to_gene

def main():
    args = parser.parse_args()
    cog_prot = args.cog_prot
    cog_align = args.cog_align
    output_path = args.output

    cog_protein = pd.read_csv( cog_prot
        ,
        on_bad_lines='skip', delimiter=',', header=None,
        names=['1', '2', 'proteinID', '4', '5', '6', 'cogID', '8', '9', '10', '11', '12', '13'])

    cog_alignment = pd.read_csv(
        cog_align,
        on_bad_lines='skip', delimiter='\t', header=None,
        names=['query', 'target', 'sequence_identity', '4', '5', '6', 'qstart', 'qend', 'tstart', 'tend', '11', '12'])

    gene_to_cogId, cogId_to_gene = get_gene_cog_annotations(cog_protein, cog_alignment)

    with open(output_path+ '/gene_to_cogId.pkl',
              'wb') as file:
        pickle.dump(gene_to_cogId, file)

    with open(output_path+ '/cogId_to_gene.pkl',
              'wb') as file:
        pickle.dump(cogId_to_gene, file)

if __name__ == "__main__":
    main()