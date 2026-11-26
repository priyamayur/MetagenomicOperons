import argparse
import pickle
import pandas as pd
import numpy as np
import networkx as nx
import igraph as ig
import leidenalg
from itertools import combinations, combinations_with_replacement
from scipy.stats import hypergeom
import os

# -------------------------------------------------
# CLI Arguments
# -------------------------------------------------
def parse_args():
    parser = argparse.ArgumentParser(description="Metagenomic Core Operon Discovery")

    parser.add_argument("-cp", "--cog_pairs", required=True,
                        help="Pickle file of COG pairs with 'operon' column")

    parser.add_argument("-o", "--output", required=True,
                        help="Output directory")

    parser.add_argument("-res", "--resolution", type=float, default=1,
                        help="Leiden resolution parameter")

    parser.add_argument("-ms", "--max_size", type=int, default=50,
                        help="Maximum community size before recursive splitting")

    parser.add_argument("-f", "--elbow_fraction", type=float, default=0.7,
                        help="Fraction of elbow threshold")

    parser.add_argument("-k", "--core_length", type=int, default=2,
                        help="Length of core operons")

    return parser.parse_args()


# -------------------------------------------------
# Fisher Significance Scores
# -------------------------------------------------
def compute_fisher_scores(df: pd.DataFrame) -> pd.DataFrame:

    if df.empty:
        return df.assign(sig_weight=pd.Series(dtype=float))

    N = df['operon'].sum()

    cog_totals = pd.concat([
        df[['cog1', 'operon']].rename(columns={'cog1': 'cog'}),
        df[['cog2', 'operon']].rename(columns={'cog2': 'cog'})
    ]).groupby('cog')['operon'].sum()

    def score(row):

        k = row['operon']
        n_a = cog_totals[row['cog1']]
        n_b = cog_totals[row['cog2']]

        p_val = hypergeom.sf(k - 1, N, n_a, n_b)

        return -np.log10(max(p_val, 1e-300))

    df['sig_weight'] = df.apply(score, axis=1)

    return df

# -------------------------------------------------
# Dynamic Elbow Threshold
# -------------------------------------------------
def find_elbow(weights: np.ndarray, fraction: float = 0.7):

    weights = np.sort(weights)[::-1]

    n_points = len(weights)

    if n_points == 0:
        return 0
    if n_points == 1:
        return weights[0] * fraction

    coords = np.vstack((range(n_points), weights)).T

    first = coords[0]
    last = coords[-1]

    line_vec = last - first
    line_vec = line_vec / np.linalg.norm(line_vec)

    vec_from_first = coords - first

    scalar_proj = np.sum(vec_from_first * line_vec, axis=1)
    vec_proj = np.outer(scalar_proj, line_vec)

    dist = np.linalg.norm(vec_from_first - vec_proj, axis=1)

    elbow_idx = np.argmax(dist)

    return weights[elbow_idx] * fraction


def apply_elbow_threshold(df: pd.DataFrame, fraction: float):
    threshold = find_elbow(df['sig_weight'].values, fraction)
    filtered = df[df['sig_weight'] >= threshold].copy()
    print(f"Dynamic elbow threshold: {threshold:.4f} | edges retained: {len(filtered)}/{len(df)}")
    return filtered

# -------------------------------------------------
# Build igraph Graph
# -------------------------------------------------
def build_igraph(df: pd.DataFrame):

    edges = list(df[['cog1', 'cog2']].itertuples(index=False, name=None))

    unique = list(set(sum(edges, ())))

    mapping = {cog: i for i, cog in enumerate(unique)}

    ig_edges = [(mapping[a], mapping[b]) for a, b in edges]

    weights = df['sig_weight'].tolist()

    G = ig.Graph(edges=ig_edges, directed=False)

    G.vs['name'] = unique
    G.es['weight'] = weights

    return G


# -------------------------------------------------
# Recursive Leiden
# -------------------------------------------------
def recursive_leiden(G: ig.Graph, res_param=1.5, max_size=50):

    all_subgraphs = []

    def split_subgraph(subgraph: ig.Graph, current_resolution):

        n = len(subgraph.vs)

        if n <= max_size:

            all_subgraphs.append({
                'vertices': subgraph.vs['name'],
                'edges': [
                    [subgraph.vs[e.source]['name'],
                     subgraph.vs[e.target]['name'],
                     e['weight']]
                    for e in subgraph.es
                ]
            })

            return

        partition = leidenalg.find_partition(
            subgraph,
            leidenalg.CPMVertexPartition,
            weights=subgraph.es['weight'],
            resolution_parameter=current_resolution
        )

        for community in partition:

            sg = subgraph.subgraph(community)

            split_subgraph(sg, current_resolution * 1.5)

    for comp in G.connected_components():

        split_subgraph(G.induced_subgraph(comp), res_param)

    # print("\nCommunity Summary")
    # print("ID   Nodes   Edges")

    # for i, sg in enumerate(all_subgraphs):
    #     print(i, len(sg['vertices']), len(sg['edges']))

    return all_subgraphs


def extract_core_operons_max_size(leiden_results, k=3, density_threshold=0.7, max_core_size=10):
    if k != 2:
        raise ValueError("Supplementary Algorithm 3 requires a two-node seed (-k 2).")
    if max_core_size < 2:
        raise ValueError("Maximum core size must be at least two.")

    rows = []

    for cid, comm in enumerate(leiden_results):
        G = nx.Graph()
        G.add_weighted_edges_from(comm['edges'])
        nodes = list(G.nodes())
        if len(nodes) < k:
            continue

        seeded = set()

        for start in nodes:
            if start in seeded:
                continue

            # 1. Construct a two-node seed using the highest-weight edge.
            neighbors = {n: d['weight'] for n, d in G[start].items()
                         if n != start}
            if not neighbors:
                seeded.add(start)
                continue

            nxt = max(neighbors, key=neighbors.get)
            expanded = {start, nxt}

            # 2. Rank neighboring candidates by total edge weight to the core.
            while len(expanded) < max_core_size:
                candidates = set().union(*(G.neighbors(n) for n in expanded)) - expanded
                if not candidates:
                    break

                ranked = sorted(
                    candidates,
                    key=lambda node: sum(G[node][n]['weight']
                                         for n in expanded if G.has_edge(node, n)),
                    reverse=True,
                )

                added = False
                for node in ranked:
                    test_set = expanded | {node}
                    num_nodes = len(test_set)
                    possible_edges = num_nodes * (num_nodes - 1) / 2
                    actual_edges = sum(G.has_edge(u, v)
                                       for u, v in combinations(test_set, 2))
                    current_density = actual_edges / possible_edges

                    if current_density >= density_threshold:
                        expanded.add(node)
                        added = True
                        break

                if not added:
                    break

            seeded.update(expanded)

            # 3. Final scoring and storage (unchanged).
            final_sub = G.subgraph(expanded)
            total_w = sum(d['weight'] for _, _, d in final_sub.edges(data=True))
            final_density = nx.density(final_sub)
            core_size = final_sub.number_of_nodes()
            possible_pairs = core_size * (core_size - 1) // 2
            # Missing pairs contribute zero. Self-loops are not distinct COG pairs.
            pair_weight = sum(
                data['weight'] for u, v, data in final_sub.edges(data=True) if u != v
            )
            connectivity_score = pair_weight / possible_pairs if possible_pairs else 0.0

            rows.append({
                "community_id": cid,
                "genes": tuple(sorted(expanded)),
                "total_weight": total_w,
                "density": final_density,
                "score": connectivity_score
            })

    columns = ["community_id", "genes", "total_weight", "density", "score"]
    return (
        pd.DataFrame(rows, columns=columns)
        .drop_duplicates("genes")
        .sort_values("score", ascending=False)
    )

def export_network_edges(leiden_results, cores, retained_pairs, output_path):
    """Export intra-community and intra-core edges for Cytoscape."""
    # Use one run-wide scale so the same edge has the same weight in both exports.
    max_weight = retained_pairs["sig_weight"].max()
    scale = max_weight if max_weight > 0 else 1.0
    pair_metrics = {
        (row.cog1, row.cog2): (row.operon, row.sig_weight, row.sig_weight / scale)
        for row in retained_pairs.itertuples(index=False)
    }

    edge_columns = [
        "source", "target", "operon_count", "fisher_weight",
        "fisher_weight_normalized",
    ]
    community_rows = []
    community_edges = []
    for community_id, community in enumerate(leiden_results):
        edges = {}
        for source, target, _ in community["edges"]:
            pair = tuple(sorted((source, target)))
            edge = (*pair, *pair_metrics[pair])
            edges[pair] = edge
            community_rows.append((community_id, *edge))
        community_edges.append(edges)

    core_rows = []
    for core in cores.itertuples(index=False):
        edges = community_edges[core.community_id]
        for pair in combinations_with_replacement(core.genes, 2):
            if pair in edges:
                core_rows.append((core.core_id, core.community_id, *edges[pair]))

    pd.DataFrame(
        community_rows, columns=["community_id", *edge_columns],
    ).to_csv(os.path.join(output_path, "community_edges.csv"), index=False)
    pd.DataFrame(
        core_rows, columns=["core_id", "community_id", *edge_columns],
    ).to_csv(os.path.join(output_path, "core_operon_edges.csv"), index=False)


def main():

    args = parse_args()
    os.makedirs(args.output, exist_ok=True)

    print("Loading COG pairs...")

    with open(args.cog_pairs, 'rb') as f:
        df = pickle.load(f)

    df[['cog1', 'cog2']] = np.sort(df[['cog1', 'cog2']], axis=1)

    df_group = df.groupby(['cog1', 'cog2'], as_index=False)['operon'].sum()

    df_group = compute_fisher_scores(df_group)

    filtered_df = apply_elbow_threshold(df_group, args.elbow_fraction)
    print("Building graph...")

    G = build_igraph(filtered_df)

    print("IGRAPH SUMMARY")
    print("Nodes:", G.vcount())

    print("Edges:", G.ecount())

    print("Running Leiden communities...")

    leiden_results = recursive_leiden(
        G,
        res_param=args.resolution,
        max_size=args.max_size
    )


    print("Extracting core operons...")
    cores = extract_core_operons_max_size(
        leiden_results,
        k=args.core_length,
        density_threshold=0.7,
        max_core_size=6,
    )
    core_lengths = cores["genes"].apply(len)
    # print(f"Core operons: {len(cores)}")
    # print(core_lengths.value_counts().sort_index())

    print("Saving results...")
    with open(os.path.join(args.output, "leiden_communities.pkl"), "wb") as f:
        pickle.dump(leiden_results, f)

    # The visualization stage still reads this file and its original schema.
    
    cores = cores.assign(
        core_id=range(1, len(cores) + 1),
        core_size=cores["genes"].map(len).astype(int),
    )
    cores.to_csv(os.path.join(args.output, "core_operons.csv"), index=False)
    export_network_edges(leiden_results, cores, filtered_df, args.output)

    print("Finished.")


if __name__ == "__main__":
    main()
