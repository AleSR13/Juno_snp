#!/usr/bin/env python3
import argparse
import pathlib

import yaml


def main(args):
    mock_clusters = {}
    for assembly_path in args.input:
        assembly_basename = pathlib.PurePosixPath(assembly_path).stem
        mock_clusters[assembly_basename] = 1
    with open(args.output, "w") as file:
        yaml.dump(mock_clusters, file, default_flow_style=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-i",
        "--input",
        type=str,
        metavar="LIST",
        nargs="+",
        help="List of assemblies",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=pathlib.Path,
        metavar="FILE",
        help="Output path for mock_clusters.yaml",
    )
    args = parser.parse_args()

    main(args)
