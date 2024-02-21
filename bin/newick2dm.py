import argparse
import itertools

import pandas as pd
from Bio import Phylo

# Modified from https://gist.github.com/mgalardini/1dbd33a6250241658d35ec279b9b6a79


def newick2dm(newick_file):
    """Converts a newick file to a matrix. Returns a distance matrix (dictionary)"""
    print(f"Reading Newick file {newick_file}...\n")
    tree = Phylo.read(newick_file, "newick")
    dm = {}
    print("Converting newick format to distance matrix...\n")
    for x, y in itertools.combinations(tree.get_terminals(), 2):
        pair_dist = tree.distance(x, y)
        dm[x.name] = dm.get(x.name, {})
        dm[x.name][y.name] = pair_dist
        dm[y.name] = dm.get(y.name, {})
        dm[y.name][x.name] = pair_dist
    for x in tree.get_terminals():
        dm[x.name][x.name] = 0
    return dm


def dm2file(dm, output_file):
    """Converts a distance matrix (as a dictionary) into a csv file"""
    print("Writing distance matrix to a csv file...\n")
    df = pd.DataFrame(dm)
    df.to_csv(output_file, sep="\t")


def main():
    parser = argparse.ArgumentParser(
        description="Code to obtain a distance matrix from a Newick file"
    )
    parser.add_argument(
        "-i",
        "--input",
        help="Newick file to be used as input.",
        metavar="FILE",
        required=True,
    )
    parser.add_argument(
        "-o",
        "--output",
        help="Output file (csv extension).",
        metavar="FILE",
        required=True,
    )
    args = parser.parse_args()
    dm = newick2dm(args.input)
    dm2file(dm, args.output)


if __name__ == "__main__":
    main()
