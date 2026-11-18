import argparse
import pickle
import pandas as pd
import numpy as np
import networkx as nx
import igraph as ig
import leidenalg
from itertools import combinations
from scipy.stats import hypergeom
from typing import List, Dict, Any
import os
import matplotlib.pyplot as plt
import numpy as np
import networkx as nx

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

# def compute_fisher_scores(df: pd.DataFrame) -> (pd.DataFrame):
#     N = df['operon'].sum()

#     # Calculate total occurrences for each COG
#     cog_totals = pd.concat([
#         df[['cog1', 'operon']].rename(columns={'cog1': 'cog'}),
#         df[['cog2', 'operon']].rename(columns={'cog2': 'cog'})
#     ]).groupby('cog')['operon'].sum()

#     # Identify Promiscuous Genes (Top 20 by total counts)
#     promiscuous_genes = cog_totals.sort_values(ascending=False).head(20)
    
#     print(f"\n{'='*40}")
#     print(f"TOP 20 PROMISCUOUS GENES (High frequency)")
#     print(f"{'COG':<15} | {'Total Pair Counts':<20}")
#     print("-" * 40)
#     for cog, count in promiscuous_genes.items():
#         print(f"{str(cog):<15} | {int(count):<20}")
#     print(f"{'='*40}\n")

#     def score(row):
#         k = row['operon']
#         n_a = cog_totals[row['cog1']]
#         n_b = cog_totals[row['cog2']]
#         p_val = hypergeom.sf(k - 1, N, n_a, n_b)
#         return -np.log10(max(p_val, 1e-300))

#     df['sig_weight'] = df.apply(score, axis=1)

#     return df 


# -------------------------------------------------
# Dynamic Elbow Threshold
# -------------------------------------------------
def find_elbow(weights: np.ndarray, fraction: float = 0.7):

    weights = np.sort(weights)[::-1]

    n_points = len(weights)

    if n_points == 0:
        return 0

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


# def apply_elbow_threshold(df: pd.DataFrame, fraction: float):

#     threshold = find_elbow(df['sig_weight'].values, fraction)

#     filtered = df[df['sig_weight'] >= threshold].copy()

#     print(f"Dynamic elbow threshold: {threshold:.4f}   edges retained: {len(filtered)}/{len(df)}")

#     return filtered

def apply_elbow_threshold(df: pd.DataFrame, fraction: float):
    threshold = find_elbow(df['sig_weight'].values, fraction)
    filtered = df[df['sig_weight'] >= threshold].copy()
    print(f"Dynamic elbow threshold: {threshold:.4f} | edges retained: {len(filtered)}/{len(df)}")
    # Return both so we can plot it
    return filtered, threshold

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
# Recursive Leiden (UNCHANGED)
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

    print("\nCommunity Summary")
    print("ID   Nodes   Edges")

    for i, sg in enumerate(all_subgraphs):
        print(i, len(sg['vertices']), len(sg['edges']))

    return all_subgraphs


def extract_path_quasi_hybrid(leiden_results, k=3, density_threshold=0.7):
    """
    Start with greedy paths of length k, then expand each path into quasi-cliques
    if adding nodes maintains density >= density_threshold.
    """
    rows = []

    for cid, comm in enumerate(leiden_results):
        G = nx.Graph()
        G.add_weighted_edges_from(comm['edges'])
        nodes = list(G.nodes())
        if len(nodes) < k:
            continue

        for start in nodes:
            path = [start]
            current = start

            # Greedy path of length k
            for _ in range(k - 1):
                neighbors = {n: d['weight'] for n, d in G[current].items() if n not in path}
                if not neighbors:
                    break
                nxt = max(neighbors, key=neighbors.get)
                path.append(nxt)
                current = nxt

            if len(path) < k:
                continue

            path_set = set(path)

            # Expand path into quasi-clique
            candidates = set(G.nodes()) - path_set
            expanded = path_set.copy()

            while candidates:
                best_candidate = None
                best_density = 0

                for node in candidates:
                    test_set = expanded | {node}
                    pairs = list(combinations(test_set, 2))
                    edge_weights = [G[a][b]['weight'] if G.has_edge(a, b) else 0 for a, b in pairs]
                    density = sum(1 for w in edge_weights if w > 0) / len(pairs)
                    if density >= density_threshold and density > best_density:
                        best_density = density
                        best_candidate = node

                if best_candidate:
                    expanded.add(best_candidate)
                    candidates.remove(best_candidate)
                else:
                    break

            # Score the hybrid core
            pairs = list(combinations(expanded, 2))
            edge_weights = [G[a][b]['weight'] if G.has_edge(a, b) else 0 for a, b in pairs]
            density = sum(1 for w in edge_weights if w > 0) / len(pairs)
            score = sum(edge_weights) * density

            rows.append({
                "community_id": cid,
                "genes": tuple(sorted(expanded)),
                "total_weight": sum(edge_weights),
                "density": density,
                "score": score
            })

    df = pd.DataFrame(rows)
    return df.drop_duplicates("genes").sort_values("score", ascending=False)

def extract_path_quasi_hybrid_max_size(leiden_results, k=3, density_threshold=0.8, max_core_size=10):
    rows = []

    for cid, comm in enumerate(leiden_results):
        G = nx.Graph()
        G.add_weighted_edges_from(comm['edges'])
        nodes = list(G.nodes())
        if len(nodes) < k:
            continue

        for start in nodes:
            path = [start]
            current = start

            # 1. Build initial greedy path
            for _ in range(k - 1):
                neighbors = {n: d['weight'] for n, d in G[current].items() if n not in path}
                if not neighbors:
                    break
                nxt = max(neighbors, key=neighbors.get)
                path.append(nxt)
                current = nxt

            if len(path) < k:
                continue

            # 2. Expand into quasi-clique: Best-node-first selection
            expanded = set(path)
            candidates = set(G.nodes()) - expanded

            while len(expanded) < max_core_size and candidates:
                best_node = None
                best_new_density = -1

                for node in candidates:
                    test_set = expanded | {node}
                    num_nodes = len(test_set)
                    possible_edges = (num_nodes * (num_nodes - 1)) / 2
                    
                    # Count actual edges in the test set
                    actual_edges = 0
                    for u, v in combinations(test_set, 2):
                        if G.has_edge(u, v):
                            actual_edges += 1
                    
                    current_density = actual_edges / possible_edges

                    # Update best candidate if it meets threshold and is better than others
                    if current_density >= density_threshold:
                        if current_density > best_new_density:
                            best_new_density = current_density
                            best_node = node

                if best_node:
                    expanded.add(best_node)
                    candidates.remove(best_node)
                else:
                    # No more candidates meet the density threshold
                    break

            # 3. Final Scoring and Storage
            final_sub = G.subgraph(expanded)
            total_w = sum(d['weight'] for _, _, d in final_sub.edges(data=True))
            final_density = nx.density(final_sub)
            
            rows.append({
                "community_id": cid,
                "genes": tuple(sorted(expanded)),
                "total_weight": total_w,
                "density": final_density,
                "score": total_w * final_density
            })

    return pd.DataFrame(rows).drop_duplicates("genes").sort_values("score", ascending=False)

def extract_path_quasi_hybrid_weighted(leiden_results, k=3, density_threshold=0.7, max_core_size=10):
    """
    Extract core operons using path + quasi-clique expansion based on weighted density.
    
    Args:
        leiden_results: List of communities from Leiden.
        k: Initial greedy path length.
        density_threshold: Minimum weighted density to add a node.
        max_core_size: Maximum allowed size of a core.
        
    Returns:
        DataFrame with columns: community_id, genes, total_weight, weighted_density, score
    """
    rows = []

    for cid, comm in enumerate(leiden_results):
        G = nx.Graph()
        G.add_weighted_edges_from(comm['edges'])
        nodes = list(G.nodes())
        if len(nodes) < k:
            continue

        max_weight = max([d['weight'] for _, _, d in G.edges(data=True)] or [1])

        for start in nodes:
            path = [start]
            current = start

            # 1. Build initial greedy path
            for _ in range(k - 1):
                neighbors = {n: d['weight'] for n, d in G[current].items() if n not in path}
                if not neighbors:
                    break
                nxt = max(neighbors, key=neighbors.get)
                path.append(nxt)
                current = nxt

            if len(path) < k:
                continue

            # 2. Expand using weighted density
            expanded = set(path)
            candidates = set(G.nodes()) - expanded

            while len(expanded) < max_core_size and candidates:
                best_candidate = None
                best_density = -1

                for node in candidates:
                    test_set = expanded | {node}
                    pairs = list(combinations(test_set, 2))
                    edge_weights = [G[a][b]['weight'] if G.has_edge(a, b) else 0 for a, b in pairs]
                    weighted_density = sum(edge_weights) / (len(pairs) * max_weight)

                    if weighted_density >= density_threshold and weighted_density > best_density:
                        best_density = weighted_density
                        best_candidate = node

                if best_candidate:
                    expanded.add(best_candidate)
                    candidates.remove(best_candidate)
                else:
                    break

            # 3. Compute final metrics
            pairs = list(combinations(expanded, 2))
            edge_weights = [G[a][b]['weight'] if G.has_edge(a, b) else 0 for a, b in pairs]
            total_weight = sum(edge_weights)
            final_density = total_weight / (len(pairs) * max_weight) if pairs else 0
            score = total_weight * final_density

            rows.append({
                "community_id": cid,
                "genes": tuple(sorted(expanded)),
                "total_weight": total_weight,
                "weighted_density": final_density,
                "score": score
            })

    return pd.DataFrame(rows).drop_duplicates("genes").sort_values("score", ascending=False)

# def export_large_communities_to_cytoscape(leiden_results, output_path, min_nodes=20):
#     edge_rows = []
#     node_rows = []
    
#     # Counter for how many communities met the criteria
#     large_comm_count = 0

#     for cid, comm in enumerate(leiden_results):
#         nodes_in_comm = comm['vertices']
#         comm_size = len(nodes_in_comm)
        
#         # Only process communities larger than the threshold
#         if comm_size > min_nodes:
#             large_comm_count += 1
#             print(f"Processing Community {cid} with {comm_size} nodes")
#             # 1. Process Nodes
#             for node in nodes_in_comm:
#                 node_rows.append({
#                     "Node": node,
#                     "Community_ID": cid,
#                     "Community_Size": comm_size
#                 })
                
#             # 2. Process Edges
#             for edge in comm['edges']:
#                 edge_rows.append({
#                     "Source": edge[0],
#                     "Target": edge[1],
#                     "Weight": edge[2],
#                     "Community_ID": cid
#                 })

#     if not node_rows:
#         print(f"No communities found with more than {min_nodes} nodes.")
#         return

#     # Create DataFrames and remove potential duplicates
#     nodes_df = pd.DataFrame(node_rows).drop_duplicates(subset=['Node', 'Community_ID'])
#     edges_df = pd.DataFrame(edge_rows)

#     # Save to CSV
#     nodes_df.to_csv(f"{output_path}/cytoscape_nodes_gt20.csv", index=False)
#     edges_df.to_csv(f"{output_path}/cytoscape_edges_gt20.csv", index=False)
    
#     print(f"Exported {large_comm_count} communities with > {min_nodes} nodes.")
#     print(f"Files saved to {output_path}")
    
    
def export_large_communities_to_cytoscape(leiden_results, output_path, min_nodes=20):
    edge_rows = []
    node_rows = []
    centrality_rows = []
    
    # Create a specific directory for individual community files
    indiv_dir = os.path.join(output_path, "individual_communities")
    os.makedirs(indiv_dir, exist_ok=True)

    print(f"\n{'='*75}")
    print(f"CALCULATING CENTRALITY & SAVING INDIVIDUAL FILES (Size > {min_nodes})")
    print(f"{'='*75}")

    for cid, comm in enumerate(leiden_results):
        nodes_in_comm = comm['vertices']
        if len(nodes_in_comm) > min_nodes:
            
            # 1. Create the subgraph
            G_temp = nx.Graph()
            G_temp.add_weighted_edges_from(comm['edges'])
            
            try:
                # Calculate Betweenness Centrality (The "Gatekeeper" score)
                # normalized=True allows comparison between different community sizes
                betweenness_dict = nx.betweenness_centrality(G_temp, normalized=True, weight='weight')
                degree_dict = dict(G_temp.degree())

                # 2. PRINT RANKINGS TO CONSOLE
                top_central = sorted(betweenness_dict.items(), key=lambda x: x[1], reverse=True)[:5]
                
                print(f"\n[Community ID: {cid} | Size: {len(nodes_in_comm)}]")
                print(f"{'Rank':<6} | {'COG/Gene':<15} | {'Centrality Score':<18}")
                print("-" * 45)
                
                for i, (node, score) in enumerate(top_central, 1):
                    print(f"{i:<6} | {str(node):<15} | {score:.4f}")
                    
                    # Store for the master summary CSV
                    centrality_rows.append({
                        "Community_ID": cid, "Rank": i, "Node": node, 
                        "Betweenness": score, "Degree": degree_dict.get(node)
                    })

                # 3. Prepare and SAVE INDIVIDUAL FILES
                this_comm_nodes = []
                for node in nodes_in_comm:
                    n_data = {
                        "Node": node, 
                        "Betweenness": betweenness_dict.get(node, 0), 
                        "Degree": degree_dict.get(node, 0)
                    }
                    this_comm_nodes.append(n_data)
                    node_rows.append({**n_data, "Community_ID": cid})

                this_comm_edges = []
                for edge in comm['edges']:
                    e_data = {"Source": edge[0], "Target": edge[1], "Weight": edge[2]}
                    this_comm_edges.append(e_data)
                    edge_rows.append({**e_data, "Community_ID": cid})

                pd.DataFrame(this_comm_nodes).to_csv(f"{indiv_dir}/comm_{cid}_nodes.csv", index=False)
                pd.DataFrame(this_comm_edges).to_csv(f"{indiv_dir}/comm_{cid}_edges.csv", index=False)

            except Exception as e:
                print(f"Error processing Community {cid}: {e}")

    # Save Master Summary Files
    if node_rows:
        pd.DataFrame(node_rows).to_csv(f"{output_path}/master_nodes_all.csv", index=False)
        pd.DataFrame(edge_rows).to_csv(f"{output_path}/master_edges_all.csv", index=False)
        pd.DataFrame(centrality_rows).to_csv(f"{output_path}/top_5_central_nodes_summary.csv", index=False)
        print(f"\n{'='*75}")
        print(f"DONE. Individual files in: {indiv_dir}")
    else:
        print("\nNo communities met the size threshold.")
def plot_elbow_threshold(df: pd.DataFrame, threshold: float, output_path: str):
    """
    Visualizes the ranked Fisher scores and the chosen elbow threshold.
    """
    weights = np.sort(df['sig_weight'].values)[::-1]
    
    plt.figure(figsize=(10, 6))
    plt.plot(range(len(weights)), weights, color='#2c3e50', linewidth=2, label='Edge Weights (Ranked)')
    
    # Draw the threshold line
    plt.axhline(y=threshold, color='#e74c3c', linestyle='--', linewidth=2, 
                label=f'Elbow Threshold (W={threshold:.2f})')
    
    # Identify the index where the threshold crosses the data
    # (Finding the first weight that is less than or equal to the threshold)
    intersect_idx = np.where(weights <= threshold)[0][0] if any(weights <= threshold) else 0
    plt.scatter(intersect_idx, threshold, color='#e74c3c', s=100, zorder=5)

    plt.title('Pruning Edges: Edge Weight Distribution & Elbow Detection', fontsize=14)
    plt.xlabel('Edge Rank (Strongest to Weakest Edge Weights)', fontsize=12)
    plt.ylabel('Normalized Edge Weights (-log10 P)', fontsize=12)
    
    # Style the plot
    plt.grid(True, linestyle=':', alpha=0.6)
    plt.legend()
    
    # Save the figure
    plot_file = os.path.join(output_path, "edge_pruning_elbow_plot.png")
    plt.savefig(plot_file, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"Elbow visualization saved to: {plot_file}")
    
def weighted_density_from_comm(comm):
    """
    comm['edges'] is a list of [node1, node2, weight]
    """
    G = nx.Graph()
    # Convert list of lists to tuples for add_weighted_edges_from
    edges = [tuple(e) for e in comm['edges']]
    G.add_weighted_edges_from(edges)
    
    total_weight = sum(d['weight'] for _, _, d in G.edges(data=True))
    n = len(G.nodes())
    max_edges = n * (n - 1) / 2
    return total_weight / max_edges if max_edges > 0 else 0
# -------------------------------------------------
# MAIN PIPELINE
# -------------------------------------------------
def main():

    args = parse_args()

    print("Loading COG pairs...")

    with open(args.cog_pairs, 'rb') as f:
        df = pickle.load(f)

    df[['cog1', 'cog2']] = pd.DataFrame(
        np.sort(df[['cog1', 'cog2']], axis=1)
    )

    df_group = df.groupby(['cog1', 'cog2'], as_index=False)['operon'].sum()

    df_group = compute_fisher_scores(df_group)

    # filtered_df = apply_elbow_threshold(df_group, args.elbow_fraction)
    # print(f"Total unique edges before pruning: {len(df_group)}")
    filtered_df, threshold_val = apply_elbow_threshold(df_group, args.elbow_fraction)
    plot_elbow_threshold(df_group, threshold_val, args.output)
    # print(f"Total unique edges after pruning: {len(filtered_df)}")
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


    densities = []

    for comm in leiden_results:
        G = nx.Graph()
        G.add_weighted_edges_from(comm['edges'])
        densities.append(nx.density(G))

    mean_density = np.mean(densities)
    std_density = np.std(densities)

    
    
    
    
    print("Extracting core operons...")

 
    # -------------------------------
    # PATH + QUASI-CLIQUE HYBRID
    # -------------------------------
    # path_quasi_hybrid_df = extract_path_quasi_hybrid(
    #     leiden_results,
    #     k=args.core_length,
    #     density_threshold=0.7
    # )
    max_core_size = 6
    path_quasi_hybrid_df = extract_path_quasi_hybrid_max_size(
        leiden_results,
        k=args.core_length,
        density_threshold=0.7,
        max_core_size=max_core_size
    )
    # max_core_size = 6
    # path_quasi_hybrid_df = extract_path_quasi_hybrid_weighted(
    # leiden_results,
    # k=args.core_length,
    # density_threshold=0.7,
    # max_core_size=max_core_size
    # )
    # Print distribution of core lengths
    print("\nPath + Quasi-Hybrid core size distribution:")
    print(path_quasi_hybrid_df['genes'].apply(len).value_counts().sort_index())
    core_lengths = path_quasi_hybrid_df['genes'].apply(len)
    print("count of cores", len(core_lengths))
    print("\nPath + Quasi-Hybrid core size distribution:")
    print(core_lengths.value_counts().sort_index())
    
    print("Saving results...")

    with open(f"{args.output}/leiden_communities.pkl", "wb") as f:
        pickle.dump(leiden_results, f)
        

    path_quasi_hybrid_df.to_csv(f"{args.output}/path_quasi_hybrid_cores.csv", index=False)
    
    export_large_communities_to_cytoscape(leiden_results, args.output, min_nodes=20)
    
    print(f"Mean density: {mean_density:.3f} ± {std_density:.3f}")
    
    densities = [weighted_density_from_comm(comm) for comm in leiden_results]

    mean_density = np.mean(densities)
    std_density = np.std(densities)

    print(f"Weighted mean ± std: {mean_density:.3f} ± {std_density:.3f}")
        
    print("Done.")

    print("Finished.")
if __name__ == "__main__":
    main()
    
    