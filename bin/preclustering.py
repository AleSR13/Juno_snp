#!/usr/bin/env python3
import argparse
import os
from pathlib import Path

import networkx as nx
import pandas as pd
import yaml


def read_mash(filepath: str) -> pd.DataFrame:
    df = pd.read_csv(
        filepath,
        sep="\t",
        header=None,
        names=["query", "ref", "dist", "p-value", "matches"],
    )
    df["query_clean"] = df["query"].apply(lambda x: Path(x).stem)
    df["ref_clean"] = df["ref"].apply(lambda x: Path(x).stem)
    df = df[df["query_clean"] != df["ref_clean"]]
    return df


def create_graph(df: pd.DataFrame) -> tuple[nx.Graph, dict]:
    G = nx.from_pandas_edgelist(
        df, source="query_clean", target="ref_clean", edge_attr="dist"
    )
    pos = nx.circular_layout(
        G
    )  # To do: add different layout options e.g. spring_layout
    return G, pos


def split_graph(
    graph: nx.Graph, threshold: float, pos: dict, out_dir: Path
) -> nx.Graph:
    graph_copy = graph.copy()
    edges_above_threshold = [
        [n1, n2] for n1, n2, dist in graph_copy.edges(data="dist") if dist > threshold
    ]
    graph_copy.remove_edges_from(edges_above_threshold)
    return graph_copy


def plot_graphs(
    whole_graph: nx.Graph, splitted_graph: nx.Graph, pos: dict, out_dir: Path
) -> None:
    # print(whole_graph.edges)
    # print(splitted_graph.edges)
    os.makedirs(out_dir, exist_ok=True)
    dist_labels = nx.get_edge_attributes(whole_graph, "dist")
    nx.draw(whole_graph, pos)
    nx.draw_networkx_edge_labels(whole_graph, pos, edge_labels=dist_labels)
    plt.savefig(out_dir.joinpath("pre_filter.png"))
    plt.clf()

    dist_labels = nx.get_edge_attributes(splitted_graph, "dist")
    nx.draw(splitted_graph, pos)
    nx.draw_networkx_edge_labels(splitted_graph, pos, edge_labels=dist_labels)
    plt.savefig(out_dir.joinpath("post_filter.png"))
    plt.clf()


def define_clusters(graph: nx.Graph) -> dict:
    components = nx.connected_components(graph)
    list_subgraphs = [
        subgraph for subgraph in sorted(components, key=len, reverse=True)
    ]
    sample_clusters = {}
    for cluster, sample_set in enumerate(list_subgraphs, start=1):
        for sample in sample_set:
            sample_clusters[sample] = cluster
    return sample_clusters


def write_results(clusters: dict, outpath: Path) -> None:
    with open(outpath, "w") as file:
        yaml.dump(clusters, file, default_flow_style=False)


def main(args) -> None:
    df = read_mash(args.input)
    whole_graph, positions = create_graph(df)
    splitted_graph = split_graph(whole_graph, args.threshold, positions, args.output)
    if args.plot_output is not None:
        plot_graphs(whole_graph, splitted_graph, positions, args.plot_output)
    clusters = define_clusters(splitted_graph)
    write_results(clusters, args.output)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-i",
        "--input",
        type=str,
        metavar="FILE",
        help="Mash dist output file of samples to compare",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        metavar="FILE",
        help="Output path for clusters.yaml",
    )
    parser.add_argument(
        "-p",
        "--plot-output",
        type=Path,
        metavar="DIR",
        help="Output directory for plots",
    )
    parser.add_argument(
        "-t",
        "--threshold",
        type=float,
        metavar="FLOAT",
        help="Threshold of mash distance to group samples by",
        default=0.01,
    )
    args = parser.parse_args()

    if args.plot_output is not None:
        import matplotlib.pyplot as plt

    main(args)
