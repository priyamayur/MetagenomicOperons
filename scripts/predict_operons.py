"""POEM: predict operon membership of adjacent, same-strand genes.

Workflow:
    1. Load the trained classifier, gene coordinates and assembled contigs.
    2. Identify adjacent genes on the same strand of each contig.
    3. Extract six features for each pair using parallel CPU threads.
    4. Predict operon membership and write results in batches.

Batch processing limits memory usage without changing the trained features.
"""

import sys
import os
from itertools import islice
import argparse
import json
from pathlib import Path
from Bio import SeqIO
import pandas as pd
import pickle
import warnings
from concurrent.futures import ThreadPoolExecutor
from tqdm import tqdm
import numba
import numpy as np

warnings.filterwarnings("ignore")

# 1. Command-line arguments
parser = argparse.ArgumentParser(
    description="POEM: predict operon membership from gene coordinates and contigs."
)
parser.add_argument("-g", "--gene", required=True, help="gene coordinates in csv")
parser.add_argument("--model-manifest", type=Path, required=True,
                    help="training manifest for a named classifier")
parser.add_argument("--model-name", default="monocytogenes",
                    help="classifier key in the training manifest")
parser.add_argument("-af", "--contigs", required=True, help="metagenomic contigs")
parser.add_argument("-o", "--output", required=True, default="output", help="output directory name")

parser.add_argument('--threads', type=int, default=None,
                    help='CPU workers (default: up to 8 allocated CPUs)')
parser.add_argument('--batch-size', type=int, default=5000,
                    help='Number of gene pairs processed and predicted per batch')

# Feature names must match the model training manifest.
BASELINE_FEATURES = ['igd', 'gene_ratio']
COMPRESSION_FEATURES = [
    'gene1_border_comp_ratio_for', 'gene2_border_comp_ratio_for',
    'gene1_border_comp_ratio_rev', 'gene2_border_comp_ratio_rev',
]
FULL_FEATURE_COLUMNS = COMPRESSION_FEATURES + BASELINE_FEATURES


# 2. Repeat detection
# Numba compiles the numerical loops and releases Python's GIL so that
# multiple CPU threads can execute repeat scans concurrently.
# Lowercase ASCII nucleotide codes used by the uint8 sequence arrays.
A, T, C, G = 97, 116, 99, 103


@numba.njit(cache=True, nogil=True)
def find_longest_forward_repeat(seq, i):
    """Return the length of the longest direct repeat starting at i."""
    best_len = 0
    # Compare the sequence starting at i with each earlier position j.
    for j in range(i):
        k = 0
        while i + k < len(seq) and j + k < len(seq) and seq[j + k] == seq[i + k]:
            k += 1
        best_len = max(best_len, k)
    return best_len


@numba.njit(cache=True, nogil=True)
def find_longest_inverted_repeat(seq, i):
    """Return the length of the longest reverse-complement repeat starting at i."""
    best_len = 0
    # Compare rightward bases from i with leftward complementary bases from j.
    for j in range(i - 1, -1, -1):
        k = 0
        while i + k < len(seq) and j - k >= 0:
            left = seq[j - k]
            right = seq[i + k]
            if not (
                (left == A and right == T)
                or (left == T and right == A)
                or (left == C and right == G)
                or (left == G and right == C)
            ):
                break
            k += 1
        best_len = max(best_len, k)
    return best_len


@numba.njit(cache=True, nogil=True)
def lz77_multi_signal_scan(seq, min_forward_len, min_reverse_len):
    """Return the fractions not covered by direct and inverted repeats."""
    n = len(seq)
    if n == 0:
        return 0.0, 0.0

    # Separate boolean masks record bases covered by each type of repeat.
    forward_covered = np.zeros(n, dtype=np.bool_)
    reverse_covered = np.zeros(n, dtype=np.bool_)
    for i in range(n):
        forward_len = find_longest_forward_repeat(seq, i)
        if forward_len >= min_forward_len:
            forward_covered[i:i + forward_len] = True

        reverse_len = find_longest_inverted_repeat(seq, i)
        if reverse_len >= min_reverse_len:
            reverse_covered[i:i + reverse_len] = True

    # The trained compression features are uncovered fractions, not coverage.
    return (n - forward_covered.sum()) / n, (n - reverse_covered.sum()) / n


# 3. Extract features for one gene pair
def extract_features_for_gene_pair(
    gene_pair, metagenome_contigs, min_forward_len=4, min_reverse_len=4, border_window=25,
):
    """Calculate intergenic distance, gene-length ratio and four repeat features."""
    # Input row layout: contig name at index 1, gene ID at 2, coordinates at 3:5.
    gene1, gene2 = gene_pair
    g1_start, g1_end = gene1[3:5]
    g2_start, g2_end = gene2[3:5]
    contig_seq = metagenome_contigs[gene1[1]].seq
    # Two baseline features: intergenic distance and relative gene length.
    gene_lengths = (g1_end - g1_start, g2_end - g2_start)
    features = {
        'igd': len(contig_seq[g1_end:g2_start]),
        'gene_ratio': round(min(gene_lengths) / max(gene_lengths), 4),
    }

    # Extract windows around the end of gene 1 and start of gene 2.
    # Preserve the exact training-time slicing and genomic coordinate order,
    # including for reverse-strand genes and overlapping genes.
    borders = (
        contig_seq[max(g1_start, g1_end - border_window):min(g2_start, g1_end + border_window)],
        contig_seq[max(g1_end, g2_start - border_window):g2_start + border_window],
    )
    # Each border contributes a direct and an inverted repeat feature.
    for index, border in enumerate(borders, 1):
        # Convert lowercase DNA into a numeric array for the Numba functions.
        bases = np.frombuffer(str(border).lower().encode('ascii'), dtype=np.uint8)
        forward, reverse = lz77_multi_signal_scan(
            bases, min_forward_len, min_reverse_len
        )
        features[f'gene{index}_border_comp_ratio_for'] = round(float(forward), 4)
        features[f'gene{index}_border_comp_ratio_rev'] = round(float(reverse), 4)
    return features


# 4. CPU allocation and gene pair selection
def available_cpus():
    """Respect OS CPU affinity and the SLURM allocation, if present."""
    try:
        count = len(os.sched_getaffinity(0))
    except AttributeError:
        count = os.cpu_count() or 1
    allocated = os.environ.get('SLURM_CPUS_PER_TASK')
    if allocated:
        try:
            count = min(count, int(allocated))
        except ValueError:
            pass
    return max(1, count)


def get_gene_pairs(gene_coordinates_df):
    """Yield adjacent, same-strand pairs in the original genomic order."""
    for _, group in gene_coordinates_df.groupby('contig_name'):
        group = group.sort_values('start').reset_index(drop=True)
        # Only adjacent genes on the same strand are eligible.
        same_strand = group['strand'].values[:-1] == group['strand'].values[1:]
        for i in np.flatnonzero(same_strand):
            yield (group.iloc[i].values, group.iloc[i + 1].values)


def count_gene_pairs(gene_coordinates_df):
    """Count eligible pairs for the progress bar without storing them all."""
    total = 0
    for _, group in gene_coordinates_df.groupby('contig_name'):
        strands = group.sort_values('start')['strand'].values
        total += int(np.count_nonzero(strands[:-1] == strands[1:]))
    return total


# 5. Load and validate the trained model
def load_model(args):
    """Load the selected classifier and its training-time feature settings."""
    try:
        manifest = json.loads(args.model_manifest.read_text())
    except FileNotFoundError as exc:
        raise ValueError(
            f"Model manifest is missing: {args.model_manifest}. "
            "Copy model_manifest.json and the trained model files into its directory."
        ) from exc

    models = manifest['models']
    model_name = args.model_name
    if model_name not in models:
        raise ValueError(
            f"Unknown classifier {model_name!r}. Available: {', '.join(sorted(models))}"
        )
    filename = models[model_name]['file']
    if not isinstance(filename, str) or not filename.endswith('.pkl') or Path(filename).name != filename:
        raise ValueError(f"Invalid model filename in the manifest: {filename!r}")
    model_path = args.model_manifest.parent / filename
    # Preserve the feature order expected by the trained classifier.
    feature_columns = manifest['feature_columns']
    if (
        not isinstance(feature_columns, list)
        or len(feature_columns) != len(FULL_FEATURE_COLUMNS)
        or set(feature_columns) != set(FULL_FEATURE_COLUMNS)
    ):
        raise ValueError("Manifest must list all six Full model features exactly once.")
    # Use the same border size and minimum repeat length as model training.
    border_window = manifest['border_size']
    min_repeat_len = manifest['min_repeat_len']
    for name, value in [('border_size', border_window), ('min_repeat_len', min_repeat_len)]:
        if type(value) is not int or value < 1:
            raise ValueError(f"{name} must be a positive integer in the model manifest.")

    from sklearn import __version__ as sklearn_version

    training_version = manifest.get('scikit_learn_version')
    if training_version and training_version != sklearn_version:
        print(
            f"Warning: model trained with scikit-learn {training_version}, "
            f"but POEM is running {sklearn_version}.",
            file=sys.stderr,
        )


    if not model_path.is_file():
        raise ValueError(f"Operon prediction model not found: {model_path}")
    with model_path.open('rb') as handle:
        model = pickle.load(handle)
    expected = getattr(model, 'n_features_in_', len(feature_columns))
    if expected != len(feature_columns):
        raise ValueError(
            f"{model_path.name} expects {expected} features, "
            f"but POEM is configured to supply {len(feature_columns)}."
        )
    print(f"Operon prediction model: {model_path}", flush=True)
    return model, feature_columns, border_window, min_repeat_len


# 6. Batch feature extraction and operon prediction
def process_and_predict(
    gene_coordinates_df, metagenome_contigs, ml_model,
    feature_columns, border_window, min_repeat_len,
    threads, batch_size, output_file,
):
    """Extract features, predict labels and save results one batch at a time."""
    print("Identifying adjacent, same-strand gene pairs...", flush=True)
    total_pairs = count_gene_pairs(gene_coordinates_df)
    print(f"Found {total_pairs:,} eligible gene pairs", flush=True)
    print(f"Using {threads} CPU threads (batch size: {batch_size:,})", flush=True)

    # Warm up Numba before launching threads to avoid concurrent compilation.
    print("Initializing repeat detection...", flush=True)
    example_sequence = np.frombuffer(b'acgtacgtacgt', dtype=np.uint8)
    lz77_multi_signal_scan(example_sequence, min_repeat_len, min_repeat_len)

    def extract_chunk(gene_pairs):
        """Extract features for a small group of pairs on one CPU thread."""
        return [
            extract_features_for_gene_pair(
                pair, metagenome_contigs, min_repeat_len,
                min_repeat_len, border_window,
            )
            for pair in gene_pairs
        ]

    # The generator supplies pairs only when needed, reducing memory usage.
    gene_pair_iterator = iter(get_gene_pairs(gene_coordinates_df))
    temporary_file = output_file.with_suffix('.csv.tmp')
    output_columns = ['gene1', 'gene2', 'operon']

    # Write to a temporary file so incomplete results are not published as final.
    with temporary_file.open('w', newline='') as handle:
        pd.DataFrame(columns=output_columns).to_csv(handle, index=False)

        # Threads share the same contig dictionary; no per-task copying is needed.
        with ThreadPoolExecutor(max_workers=threads) as executor, tqdm(
            total=total_pairs,
            desc='Extracting features and predicting operons',
            unit='pair',
        ) as progress:
            while True:
                # Step 1: Retrieve the next batch of eligible gene pairs.
                batch_pairs = list(islice(gene_pair_iterator, batch_size))
                if not batch_pairs:
                    break

                # Step 2: Split the batch into smaller chunks for load balancing.
                # Aim for about four chunks per thread (rounding up).
                target_chunks = threads * 4
                chunk_size = max(1, (len(batch_pairs) + target_chunks - 1) // target_chunks)
                chunks = (
                    batch_pairs[i:i + chunk_size]
                    for i in range(0, len(batch_pairs), chunk_size)
                )

                # Step 3: Extract the six features in parallel. executor.map
                # preserves input order, so features remain aligned with pairs.
                batch_features = []
                for chunk_features in executor.map(extract_chunk, chunks):
                    batch_features.extend(chunk_features)

                # Step 4: Build the model input in the manifest's feature order.
                X = np.array([
                    [features[column] for column in feature_columns]
                    for features in batch_features
                ])

                # Step 5: Predict whether each pair belongs to the same operon.
                labels = ml_model.predict(X)

                # Step 6: Save gene IDs and predicted labels for this batch.
                rows = [
                    (pair[0][2], pair[1][2], label)
                    for pair, label in zip(batch_pairs, labels)
                ]
                pd.DataFrame(rows, columns=output_columns).to_csv(
                    handle, index=False, header=False,
                )
                handle.flush()

                # Count pairs only after feature extraction, prediction and saving.
                progress.update(len(batch_pairs))

    # Publish the complete CSV only after all batches finish successfully.
    temporary_file.replace(output_file)


# 7. Main workflow
def main():
    args = parser.parse_args()
    if args.threads is not None and args.threads < 1:
        parser.error('--threads must be positive')
    if args.batch_size < 1:
        parser.error('--batch-size must be positive')
    cpus = available_cpus()
    threads = min(args.threads if args.threads is not None else 8, cpus)

    try:
        ml_model, feature_columns, border_window, min_repeat_len = load_model(args)
    except (OSError, KeyError, TypeError, ValueError, pickle.PickleError) as exc:
        parser.exit(1, f'Error loading operon model: {exc}\n')

    # Keep model prediction within the same CPU allocation when supported.
    if hasattr(ml_model, 'n_jobs'):
        ml_model.n_jobs = threads

    # Load input data once and share the contig sequences across CPU threads.
    print("Loading gene coordinates...", flush=True)
    gene_coordinates_df = pd.read_csv(args.gene)
    print(f"Loaded {len(gene_coordinates_df):,} genes", flush=True)

    print("Loading assembled contigs...", flush=True)
    metagenome_contigs = SeqIO.to_dict(SeqIO.parse(args.contigs, 'fasta'))
    print(f"Loaded {len(metagenome_contigs):,} contigs", flush=True)

    # Prepare the output directory and run the batch prediction workflow.
    output_path = Path(args.output)
    output_path.mkdir(parents=True, exist_ok=True)
    output_file = output_path / 'gene_pair_names_to_operon.csv'
    process_and_predict(
        gene_coordinates_df, metagenome_contigs, ml_model,
        feature_columns, border_window, min_repeat_len,
        threads, args.batch_size, output_file,
    )
    print(f'Operon prediction completed. Results saved to {output_file}')


if __name__ == '__main__':
    main()
