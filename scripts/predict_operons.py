import sys
import argparse
from Bio import SeqIO
import pandas as pd
import pickle
import warnings
from multiprocessing import Pool, cpu_count
from tqdm import tqdm
import numba
import numpy as np

warnings.filterwarnings("ignore")

# -------------------- ARGUMENTS --------------------
parser = argparse.ArgumentParser()
parser.add_argument("-g", "--gene", required=True, help="gene coordinates in csv")
parser.add_argument("-m", "--model", required=True, help="operon prediction model")
# parser.add_argument("-sc", "--scaler", required=True, help="operon prediction model scaler")
parser.add_argument("-af", "--contigs", required=True, help="metagenomic contigs")
parser.add_argument("-o", "--output", required=True, default="output", help="output directory name")

# -------------------- HELPERS --------------------
COMPLEMENT = str.maketrans("atgc", "tacg")

def reverse_complement(seq):
    return seq.translate(COMPLEMENT)[::-1]

# -------------------- NUMBA-ACCELERATED REPEAT RATIOS --------------------
@numba.njit
def compute_repeat_ratios_numba(seq_array, min_forward_len=4, min_reverse_len=4):
    n = len(seq_array)
    forward_total = 0
    reverse_total = 0

    seq_int = np.zeros(n, dtype=np.int8)
    for i in range(n):
        c = seq_array[i]
        if c == 97: seq_int[i] = 0  # a
        elif c == 116: seq_int[i] = 1  # t
        elif c == 103: seq_int[i] = 2  # g
        elif c == 99: seq_int[i] = 3  # c
        else: seq_int[i] = -1

    # Forward repeats
    for i in range(n):
        best_f = 0
        for j in range(i):
            k = 0
            while i + k < n and j + k < n and seq_int[j + k] == seq_int[i + k]:
                k += 1
            if k > best_f:
                best_f = k
        if best_f >= min_forward_len:
            forward_total += best_f

    # Reverse repeats
    for i in range(n):
        best_r = 0
        for j in range(i-1, -1, -1):
            k = 0
            while i + k < n and j - k >= 0:
                s_i = seq_int[i + k]
                s_j = seq_int[j - k]
                match = (s_i == 0 and s_j == 1) or (s_i == 1 and s_j == 0) or (s_i == 2 and s_j == 3) or (s_i == 3 and s_j == 2)
                if not match: break
                k += 1
            if k > best_r:
                best_r = k
        if best_r >= min_reverse_len:
            reverse_total += best_r

    forward_ratio = forward_total / n if n > 0 else 0.0
    reverse_ratio = reverse_total / n if n > 0 else 0.0
    return forward_ratio, reverse_ratio

def compute_repeat_ratios(seq, min_forward_len=4, min_reverse_len=4):
    seq_array = np.frombuffer(seq.encode('ascii'), dtype=np.uint8)
    return compute_repeat_ratios_numba(seq_array, min_forward_len, min_reverse_len)

# -------------------- FEATURE EXTRACTION --------------------
def extract_features_for_gene_pair(
    gene_pair,
    metagenome_contigs,
    min_forward_len=4,
    min_reverse_len=4,
    igd_threshold=300,
    label=None,
    border_window=100
):
    gene1, gene2 = gene_pair

    g1_contig, g1_name, g1_start, g1_end, g1_strand = gene1[1], gene1[2], gene1[3], gene1[4], gene1[5]
    g2_contig, g2_name, g2_start, g2_end, g2_strand = gene2[1], gene2[2], gene2[3], gene2[4], gene2[5]

    contig_seq = metagenome_contigs[g1_contig].seq

    intergenic_seq = str(contig_seq[g1_end:g2_start])
    igd = len(intergenic_seq)
    gene_ratio = min(g1_end - g1_start, g2_end - g2_start) / max(g1_end - g1_start, g2_end - g2_start)

    feats = {
        'igd': igd,
        'gene_ratio': round(gene_ratio, 4),
        'gene1_border_comp_ratio_for': 0,
        'gene1_border_comp_ratio_rev': 0,
        'gene2_border_comp_ratio_for': 0,
        'gene2_border_comp_ratio_rev': 0,
        'compression_ratio_for': 0,
        'compression_ratio_rev': 0,
        'label': label
    }
    igr = str(contig_seq[g1_end:g2_start])
    if len(igr) > 300:
        igr = (igr[:250] )
    igd_for, igd_rev = compute_repeat_ratios(igr, min_forward_len, min_reverse_len)
    if igd == 0:
        igd_for = 1
        igd_rev = 1
    feats["compression_ratio_for"] = igd_for
    feats["compression_ratio_rev"] = igd_rev
    
    g1_border_seq = str(contig_seq[max(0, g1_end - border_window): min(len(contig_seq), g1_end + border_window)])
    g2_border_seq = str(contig_seq[max(0, g2_start - border_window): min(len(contig_seq), g2_start + border_window)])

    g1_for, g1_rev = compute_repeat_ratios(g1_border_seq, min_forward_len, min_reverse_len)
    g2_for, g2_rev = compute_repeat_ratios(g2_border_seq, min_forward_len, min_reverse_len)

    feats.update({
        'gene1_border_comp_ratio_for': round(g1_for, 4),
        'gene1_border_comp_ratio_rev': round(g1_rev, 4),
        'gene2_border_comp_ratio_for': round(g2_for, 4),
        'gene2_border_comp_ratio_rev': round(g2_rev, 4)
    })

    if g1_strand != g2_strand:
        feats['label'] = 0

    return feats

# -------------------- PARALLEL FEATURE EXTRACTION --------------------
def process_gene_pair_wrapper(args):
    return extract_features_for_gene_pair(*args)

def process_all_gene_pairs_parallel(gene_pair_list, metagenome_contigs, min_forward_len=4, min_reverse_len=4, igd_threshold=300, border_window=100):
    n_cores = max(1, cpu_count() - 1)
    pool = Pool(n_cores)
    args_list = [(gp, metagenome_contigs, min_forward_len, min_reverse_len, igd_threshold, None, border_window) for gp in gene_pair_list]

    results = []
    for res in tqdm(pool.imap(process_gene_pair_wrapper, args_list, chunksize=1000), total=len(args_list), desc="Processing gene pairs"):
        results.append(res)

    pool.close()
    pool.join()
    return results

# -------------------- GENE PAIR GENERATION --------------------
def get_gene_pairs(gene_coordinates_df):
    gene_pairs = []
    for contig, group in gene_coordinates_df.groupby('contig_name'):
        group = group.sort_values('start').reset_index(drop=True)
        same_strand = group['strand'].values[:-1] == group['strand'].values[1:]
        pairs = [[group.iloc[i].values, group.iloc[i+1].values] for i in range(len(group)-1) if same_strand[i]]
        gene_pairs.extend(pairs)
    return gene_pairs

# -------------------- ML PREDICTION WITH PROGRESS --------------------
def get_prediction_for_gene_pair_batch_progress(all_encodings, gene_pair_list, ml_model,  batch_size=5000):
    include_cols = [
        'igd', 'gene_ratio',
        'gene1_border_comp_ratio_for','gene2_border_comp_ratio_for',
        'gene1_border_comp_ratio_rev','gene2_border_comp_ratio_rev',
        'compression_ratio_for',
        'compression_ratio_rev'
    ]
    
    n = len(all_encodings)
    gene_names_pairs_to_operon = []

    for start in tqdm(range(0, n, batch_size), desc="Predicting operon labels"):
        end = min(start + batch_size, n)
        X_batch = np.array([[enc[col] for col in include_cols] for enc in all_encodings[start:end]])
        # if scaler is not None:
        #     X_batch = scaler.transform(X_batch)
        labels_batch = ml_model.predict(X_batch)
        for i, label in enumerate(labels_batch):
            gene_names_pairs_to_operon.append([gene_pair_list[start + i][0][2],
                                               gene_pair_list[start + i][1][2],
                                               label])
    return gene_names_pairs_to_operon

# -------------------- MAIN --------------------
def main():
    args = parser.parse_args()
    gene_coord = args.gene
    model = args.model
    # scaler = args.scaler
    contigs = args.contigs
    output_path = args.output

    gene_coordinates_df = pd.read_csv(gene_coord)
    with open(model, 'rb') as file: ml_model = pickle.load(file)
    # with open(scaler, 'rb') as file: scaler_obj = pickle.load(file)
    metagenome_contigs = SeqIO.to_dict(SeqIO.parse(contigs, 'fasta'))
    
    #print(len(metagenome_contigs))
    

    print("Gene coordinates, model, scaler, and contigs loaded successfully.")
    #print(f"Number of gene pairs to process: {len(gene_coordinates_df) - 1}")

    gene_pair_list = get_gene_pairs(gene_coordinates_df)
    #print(f"Generated {len(gene_pair_list)} gene pairs for processing.")

    min_forward_len = 4
    min_reverse_len = 4

    # -------------------- FEATURE EXTRACTION --------------------
    all_encodings = process_all_gene_pairs_parallel(
        gene_pair_list,
        metagenome_contigs,
        min_forward_len,
        min_reverse_len
    )

    # -------------------- BATCH ML PREDICTION WITH PROGRESS --------------------
    gene_names_pairs_to_operon = get_prediction_for_gene_pair_batch_progress(
        all_encodings,
        gene_pair_list,
        ml_model,
        # scaler_obj
    )

    # -------------------- SAVE OUTPUT --------------------
    df_out = pd.DataFrame(gene_names_pairs_to_operon, columns=['gene1', 'gene2', 'operon'])
    df_out.to_csv(f"{output_path}/gene_pair_names_to_operon.csv", index=False)
    print(f"Operon prediction completed. Results saved to {output_path}/gene_pair_names_to_operon.csv")

if __name__ == "__main__":
    main()