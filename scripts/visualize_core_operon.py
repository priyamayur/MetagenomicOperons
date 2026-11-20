#!/usr/bin/env python3
"""POEM Core Operon Explorer: Leiden communities and predicted cores."""

import argparse
import ast
import json
import pickle
from pathlib import Path

import networkx as nx
import pandas as pd
from dash import Dash, Input, Output, State, dcc, html, ctx
import dash_cytoscape as cyto


# --------------------------------------------------
# 1. Command-line arguments
# --------------------------------------------------

def parse_args():
    parser = argparse.ArgumentParser(
        description="POEM Core Operon Explorer"
    )

    parser.add_argument(
        "-g", "--graph", required=True,
        help="Path to leiden_communities.pkl"
    )

    parser.add_argument(
        "-cores", "--cores", required=True,
        help="Path to core_operons.csv"
    )

    parser.add_argument(
        "-cdef", "--cog_def", required=True,
        help="Path to COG descriptions"
    )

    parser.add_argument(
        "-o", "--output", required=True,
        help="Output directory"
    )

    parser.add_argument("--host", default="127.0.0.1")
    parser.add_argument("--port", type=int, default=8050)

    return parser.parse_args()


# --------------------------------------------------
# 2. Load input data
# --------------------------------------------------

def parse_genes(value):
    if isinstance(value, (tuple, list, set)):
        return [str(v).strip() for v in value]

    parsed = ast.literal_eval(str(value))

    if not isinstance(parsed, (tuple, list, set)):
        raise ValueError(
            f"Expected tuple/list of COG IDs, got {value!r}"
        )

    return [str(v).strip() for v in parsed]


def load_data(args):

    # Load Leiden communities
    with open(args.graph, "rb") as f:
        communities = pickle.load(f)

    # Load final core operons
    cores_df = pd.read_csv(args.cores)

    required = {
        "community_id", "genes", "score", "density"
    }

    missing = required - set(cores_df.columns)

    if missing:
        raise ValueError(
            f"Core CSV missing columns: {sorted(missing)}"
        )

    cores = {}

    for index, row in cores_df.iterrows():

        community_id = int(row["community_id"])

        if not 0 <= community_id < len(communities):
            raise ValueError(
                f"Invalid community ID: {community_id}"
            )

        genes = parse_genes(row["genes"])
        core = {
            "id": str(row.get("core_id", index)),
            "genes": genes,
            "size": int(row.get("core_size", len(genes))),
            "score": float(row["score"]),
            "density": float(row["density"])
        }

        cores.setdefault(community_id, []).append(core)

    # Rank cores by score
    for group in cores.values():
        group.sort(
            key=lambda item: -item["score"]
        )

    # Load COG descriptions
    cog_df = pd.read_csv(
        args.cog_def,
        sep="\t",
        encoding="latin-1",
        header=None,
        on_bad_lines="skip",
        names=[
            "cogID",
            "functionalCategory",
            "cogName",
            "GeneAssociatedCog",
            "FunctionalPathway",
            "pubmedID",
            "pdbID"
        ]
    )

    descriptions = {}

    for _, row in cog_df.iterrows():

        cog = str(row["cogID"]).strip()

        descriptions[cog] = {
            "name": (
                "" if pd.isna(row["cogName"])
                else str(row["cogName"])
            ),
            "category": (
                "" if pd.isna(row["functionalCategory"])
                else str(row["functionalCategory"])
            )
        }

    return communities, cores, descriptions


# --------------------------------------------------
# 3. Prepare and validate Leiden graphs
# --------------------------------------------------

def load_edge_metrics(graph_path):
    """Read edge counts and normalized weights from the community edge table."""
    edge_file = Path(graph_path).with_name("community_edges.csv")
    if not edge_file.exists():
        return {}

    columns = {"community_id", "source", "target", "operon_count", "fisher_weight_normalized"}
    rows = pd.read_csv(edge_file, usecols=lambda column: column in columns)
    missing = columns - set(rows.columns)
    if missing:
        raise ValueError(f"{edge_file} is missing columns: {sorted(missing)}")

    metrics = {}
    for row in rows.itertuples(index=False):
        endpoints = sorted((str(row.source).strip(), str(row.target).strip()))
        metrics[(int(row.community_id), *endpoints)] = {
            "operon_count": int(row.operon_count),
            "fisher_weight_normalized": float(row.fisher_weight_normalized),
        }
    return metrics


def prepare_networks(communities, cores, edge_metrics=None):
    edge_metrics = edge_metrics or {}

    networks = {}

    for cid, community in enumerate(communities):

        G = nx.Graph()

        vertices = [
            str(node).strip()
            for node in community["vertices"]
        ]

        G.add_nodes_from(vertices)

        for source, target, weight in community["edges"]:

            source = str(source).strip()
            target = str(target).strip()

            if source not in G or target not in G:
                raise ValueError(
                    f"Community {cid}: "
                    f"edge endpoint not found in vertices"
                )

            G.add_edge(
                source,
                target,
                weight=float(weight),
                **edge_metrics.get((cid, *sorted((source, target))), {}),
            )

        # Validate community connectivity
        if G.number_of_nodes() > 0:

            if not nx.is_connected(G):

                print(
                    f"WARNING: Community {cid} contains "
                    f"{nx.number_connected_components(G)} "
                    f"connected components.",
                    flush=True
                )

        # Calculate node positions ONCE
        if G.number_of_nodes() == 1:

            positions = {
                next(iter(G.nodes)): (0.0, 0.0)
            }

        elif G.number_of_nodes() > 1:

            positions = nx.spring_layout(
                G,
                seed=42,
                weight=None,
                iterations=200
            )

        else:
            positions = {}

        positions = {
            node: {
                "x": float(coords[0]) * 450,
                "y": float(coords[1]) * 450
            }
            for node, coords in positions.items()
        }

        networks[cid] = {
            "graph": G,
            "positions": positions
        }

        # Validate final cores
        for core in cores.get(cid, []):

            missing = (
                set(core["genes"]) - set(G.nodes)
            )

            if missing:
                raise ValueError(
                    f"Core {core['id']} contains "
                    f"missing COGs: {missing}"
                )

            subgraph = G.subgraph(core["genes"])

            if (
                subgraph.number_of_nodes() > 1
                and not nx.is_connected(subgraph)
            ):
                print(
                    f"WARNING: Core {core['id']} "
                    f"in community {cid} is disconnected.",
                    flush=True
                )

    return networks


# --------------------------------------------------
# 4. Cytoscape stylesheet
# --------------------------------------------------

STYLESHEET = [

    # Regular nodes
    {
        "selector": "node",
        "style": {
            "label": "data(label)",
            "background-color": "#cbd5e1",
            "color": "#1e293b",
            "font-size": 11,
            "width": 42,
            "height": 42,
            "text-valign": "center",
            "text-halign": "center"
        }
    },

    # Nodes belonging to selected core
    {
        "selector": "node.core",
        "style": {
            "background-color": "#f4d58d",
            "border-color": "#a16207",
            "border-width": 3,
            "shape": "hexagon",
            "width": 50,
            "height": 50
        }
    },

    # Regular community edges
    {
        "selector": "edge",
        "style": {
            "line-color": "#64748b",
            "width": 2.5,
            "opacity": 0.75,
            "curve-style": "bezier"
        }
    },

    # Core edges
    {
        "selector": "edge.core",
        "style": {
            "line-color": "#c05621",
            "width": 5,
            "opacity": 1
        }
    }

]


# --------------------------------------------------
# 5. Create Cytoscape graph elements
# --------------------------------------------------

def make_elements(
    community_id,
    network,
    core,
    view,
    descriptions
):

    G = network["graph"]

    selected = (
        set(core["genes"]) if core else set()
    )

    # Full community or final core only
    if view == "core":
        visible = selected
    else:
        visible = set(G.nodes)

    elements = []

    # Unique prefix for each community
    prefix = f"community_{community_id}__"

    # Add nodes
    for cog in sorted(visible):

        annotation = descriptions.get(cog, {})

        elements.append({
            "data": {
                "id": prefix + cog,
                "cog_id": cog,
                "label": cog,
                "description": annotation.get(
                    "name", "Unknown COG"
                ),
                "category": annotation.get(
                    "category", ""
                )
            },

            "position": network["positions"][cog],

            "classes": (
                "core" if cog in selected else ""
            )
        })

    # Add ALL existing edges between visible nodes
    for index, (source, target, attrs) in enumerate(
        G.edges(data=True)
    ):

        if source not in visible:
            continue

        if target not in visible:
            continue

        fisher_weight = float(attrs["weight"])
        raw_count = attrs.get("operon_count")
        label = f"Fisher {fisher_weight:.2f}"
        if raw_count is not None:
            label = f"Raw {raw_count:,} | {label}"

        elements.append({
            "data": {
                "id": f"{prefix}edge_{index}",
                "source": prefix + source,
                "target": prefix + target,
                "weight": fisher_weight,
                "fisher_weight": fisher_weight,
                "operon_count": raw_count,
                "fisher_weight_normalized": attrs.get("fisher_weight_normalized"),
                "edge_label": label,
            },

            "classes": (
                "core"
                if source in selected and target in selected
                else ""
            )
        })

    return elements


# --------------------------------------------------
# 6. Create Dash application
# --------------------------------------------------

def create_app(networks, cores, descriptions):

    app = Dash(__name__)
    all_cores = [core for group in cores.values() for core in group]
    core_sizes = sorted({core["size"] for core in all_cores})
    scores = [core["score"] for core in all_cores]
    min_score = min(scores, default=0.0)
    max_score = max(scores, default=0.0)
    score_span = max_score - min_score

    def matching_cores(cid, sizes, score_range):
        low, high = score_range or (min_score, max_score)
        return [
            core for core in cores.get(cid, [])
            if (not sizes or core["size"] in sizes) and low <= core["score"] <= high
        ]

    def community_choices(sizes, score_range):
        options = []
        for cid, network in networks.items():
            count = len(matching_cores(cid, sizes, score_range))
            if count:
                n_nodes = network["graph"].number_of_nodes()
                options.append({
                    "label": f"Community {cid} ({n_nodes} COGs, {count} cores)",
                    "value": cid,
                })
        return options

    community_options = community_choices([], None)
    first_community = community_options[0]["value"] if community_options else None

    # ----------------------------------------------
    # Layout
    # ----------------------------------------------

    app.layout = html.Div([

        html.H2("POEM — Core Operon Explorer"),

        # Community and core selection
        html.Div([

            html.Div([
                html.Label("Core Operon Size"),
                dcc.Dropdown(
                    id="core-size",
                    options=[{"label": f"{size} COGs", "value": size} for size in core_sizes],
                    value=[],
                    multi=True,
                    placeholder="All core operon sizes",
                ),
            ], style={"flex": 1, "minWidth": "200px"}),

            html.Div([

                html.Label("Leiden Community"),

                dcc.Dropdown(
                    id="community",
                    options=community_options,
                    value=first_community,
                    clearable=False,
                    searchable=True
                )

            ], style={
                "flex": 1,
                "minWidth": "250px"
            }),

            html.Div([

                html.Label("Core Operon"),

                dcc.Dropdown(
                    id="core",
                    clearable=False
                )

            ], style={
                "flex": 1,
                "minWidth": "250px"
            })

        ], style={
            "display": "flex",
            "gap": "18px",
            "flexWrap": "wrap"
        }),

        html.Div([
            html.Label("Core Operon Score Range"),
            html.Div(
                dcc.RangeSlider(
                    id="core-score-range",
                    min=min_score,
                    max=max_score if score_span else min_score + 1.0,
                    step=score_span / 200 if score_span else 1.0,
                    value=[min_score, max_score],
                    marks={min_score: f"{min_score:.4g}", max_score: f"{max_score:.4g}"},
                    allowCross=False,
                    disabled=not score_span,
                ),
                style={} if score_span else {"display": "none"},
            ),
            html.Div(id="score-filter-summary"),
        ], style={"marginTop": "18px"}),

        # Visualization mode
        html.Div([

            html.Label("Visualization Mode"),

            dcc.RadioItems(
                id="mode",
                value="full",
                inline=True,
                options=[
                    {
                        "label": " Full Leiden Community",
                        "value": "full"
                    },
                    {
                        "label": " Core Operon Only",
                        "value": "core"
                    }
                ]
            ),
            dcc.Checklist(
                id="edge-weights",
                options=[{"label": " Show edge weights on graph", "value": "show"}],
                value=[],
            ),

        ], style={"marginTop": "18px"}),

        # Summary
        html.Div(
            id="summary",
            style={
                "margin": "16px 0",
                "fontWeight": "bold"
            }
        ),

        # Cytoscape network
        cyto.Cytoscape(
            id="network",

            elements=[],

            stylesheet=STYLESHEET,

            layout={
                "name": "preset",
                "fit": True,
                "padding": 50
            },

            minZoom=0.1,
            maxZoom=5,

            style={
                "height": "680px",
                "width": "100%",
                "border": "1px solid #cbd5e1"
            }
        ),

        # COG and edge information
        html.Div(
            id="click-details",
            children="Click a node or edge for details.",
            style={
                "padding": "14px",
                "backgroundColor": "#f1f5f9",
                "marginTop": "12px"
            }
        ),

        # Download buttons
        html.Div([

            html.Button(
                "Download JSON",
                id="json-button"
            ),

            dcc.Download(
                id="json-download"
            ),

            html.Button(
                "Download PNG",
                id="png-button",
                style={"marginLeft": "12px"}
            )

        ], style={"marginTop": "16px"})

    ], style={
        "maxWidth": "1250px",
        "margin": "auto",
        "padding": "20px",
        "fontFamily": "Arial"
    })

    # --------------------------------------------------
    # 7. Update core dropdown
    # --------------------------------------------------

    @app.callback(
        Output("community", "options"),
        Output("community", "value"),
        Output("score-filter-summary", "children"),
        Input("core-size", "value"),
        Input("core-score-range", "value"),
        State("community", "value"),
    )
    def update_communities(sizes, score_range, cid):
        options = community_choices(sizes, score_range)
        available = [option["value"] for option in options]
        selected = cid if cid in available else next(iter(available), None)
        if not scores:
            return options, selected, "No predicted cores."

        low, high = score_range or (min_score, max_score)
        matching_count = sum(len(matching_cores(cid, sizes, score_range)) for cid in available)
        summary = (
            f"Core Operon Score Range: {min_score:.4g} to {max_score:.4g} | "
            f"Selected: {low:.4g} to {high:.4g} | "
            f"{matching_count:,} of {len(scores):,} core operons match"
        )
        return options, selected, summary

    @app.callback(
        Output("core", "options"),
        Output("core", "value"),
        Input("community", "value"),
        Input("core-size", "value"),
        Input("core-score-range", "value"),
        State("core", "value"),
    )
    def update_cores(cid, sizes, score_range, core_id):
        matches = matching_cores(cid, sizes, score_range)
        options = [
            {
                "label": f"Core {core['id']} ({core['size']} COGs, score {core['score']:.2f})",
                "value": core["id"],
            }
            for core in matches
        ]
        selected = next((core["id"] for core in matches if core["id"] == core_id), None)
        if selected is None and matches:
            selected = matches[0]["id"]
        return options, selected

    @app.callback(
        Output("network", "stylesheet"),
        Input("edge-weights", "value"),
    )
    def update_edge_labels(selection):
        if "show" not in (selection or []):
            return STYLESHEET
        return [
            *STYLESHEET,
            {
                "selector": "edge",
                "style": {
                    "label": "data(edge_label)",
                    "font-size": 16,
                    "text-background-color": "#ffffff",
                    "text-background-opacity": 0.9,
                    "text-background-padding": "2px",
                },
            },
        ]

    # --------------------------------------------------
    # 8. Update Cytoscape network
    # --------------------------------------------------

    @app.callback(
        Output("network", "elements"),
        Output("summary", "children"),

        Input("community", "value"),
        Input("core", "value"),
        Input("mode", "value"),
        Input("core-size", "value"),
        Input("core-score-range", "value"),
    )
    def update_network(cid, core_id, mode, sizes, score_range):
        matches = matching_cores(cid, sizes, score_range)
        if not matches:
            return [], "No predicted cores match the current selection."

        selected_core = next(
            (core for core in matches if core["id"] == core_id),
            matches[0],
        )

        network = networks[cid]

        elements = make_elements(
            cid,
            network,
            selected_core,
            mode,
            descriptions
        )

        # Validate Cytoscape node references
        node_ids = {
            e["data"]["id"]
            for e in elements
            if "source" not in e["data"]
        }

        missing_endpoints = []

        for e in elements:

            if "source" not in e["data"]:
                continue

            if (
                e["data"]["source"] not in node_ids
                or e["data"]["target"] not in node_ids
            ):
                missing_endpoints.append(
                    e["data"]["id"]
                )

        if missing_endpoints:
            raise RuntimeError(
                f"Community {cid}: edges reference "
                f"missing nodes: {missing_endpoints}"
            )

        # Validate unique Cytoscape IDs
        ids = [
            e["data"]["id"]
            for e in elements
        ]

        if len(ids) != len(set(ids)):
            raise RuntimeError(
                f"Community {cid}: duplicate Cytoscape IDs"
            )

        # Count nodes and edges
        n_nodes = sum(
            "source" not in e["data"]
            for e in elements
        )

        n_edges = len(elements) - n_nodes

        # Full community must retain every saved edge
        if mode == "full":

            expected = network["graph"].number_of_edges()

            if n_edges != expected:
                raise RuntimeError(
                    f"Community {cid}: displayed "
                    f"{n_edges} of {expected} edges"
                )

        summary = (
            f"Community {cid} | "
            f"{n_nodes} COGs | "
            f"{n_edges} edges | "
            f"Core {selected_core['id']}: "
            f"{selected_core['size']} COGs | "
            f"Core operon score {selected_core['score']:.3f}"
        )

        print(summary, flush=True)

        return elements, summary

    # --------------------------------------------------
    # 9. Node and edge click information
    # --------------------------------------------------

    @app.callback(
        Output("click-details", "children"),

        Input("network", "tapNodeData"),
        Input("network", "tapEdgeData")
    )
    def show_details(node, edge):

        triggered = (
            ctx.triggered[0]["prop_id"]
            if ctx.triggered
            else ""
        )

        if triggered.endswith("tapEdgeData") and edge:

            source = edge["source"].split("__", 1)[-1]
            target = edge["target"].split("__", 1)[-1]

            raw_count = edge.get("operon_count")
            normalized = edge.get("fisher_weight_normalized")
            raw_text = f"{raw_count:,}" if raw_count is not None else "Unavailable"
            details = (
                f"Edge: {source} — {target} | Raw operon count: {raw_text} | "
                f"Fisher weight: {edge['fisher_weight']:.4f}"
            )
            if normalized is not None:
                details += f" | Normalized Fisher weight: {normalized:.4f}"
            return details

        if triggered.endswith("tapNodeData") and node:

            return (
                f"COG: {node['cog_id']} | "
                f"{node.get('description', 'Unknown COG')} | "
                f"Category: {node.get('category') or 'N/A'}"
            )

        return "Click a node or edge for details."

    # --------------------------------------------------
    # 10. Download displayed network
    # --------------------------------------------------



    @app.callback(
        Output("json-download", "data"),
        Input("json-button", "n_clicks"),
        State("network", "elements"),
        State("community", "value"),
        State("core", "value"),
        State("mode", "value"),
        prevent_initial_call=True
    )
    def download_json(n_clicks, elements, cid, core_id, mode):

        filename = f"poem_community_{cid}_core_{core_id}_{mode}.json"

        return {
            "content": json.dumps(elements or [], indent=2),
            "filename": filename
        }

    @app.callback(
        Output("network", "generateImage"),
        Input("png-button", "n_clicks"),
        State("community", "value"),
        State("core", "value"),
        State("mode", "value"),
        prevent_initial_call=True
    )
    def download_png(n_clicks, cid, core_id, mode):

        filename = f"poem_community_{cid}_core_{core_id}_{mode}"

        return {
            "type": "png",
            "action": "download",
            "filename": filename
        }

    # Return Dash application
    return app
# --------------------------------------------------
# 11. Main
# --------------------------------------------------

def main():

    args = parse_args()

    Path(args.output).mkdir(
        parents=True,
        exist_ok=True
    )

    communities, cores, descriptions = load_data(args)

    edge_metrics = load_edge_metrics(args.graph)
    networks = prepare_networks(communities, cores, edge_metrics)

    print(
        f"Loaded {len(communities)} Leiden communities "
        f"and {sum(len(v) for v in cores.values())} "
        f"final core operons."
    )

    app = create_app(
        networks,
        cores,
        descriptions
    )

    app.run(
        host=args.host,
        port=args.port,
        debug=False
    )


if __name__ == "__main__":
    main()
