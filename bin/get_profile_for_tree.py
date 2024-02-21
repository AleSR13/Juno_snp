import argparse
import pathlib

import pandas as pd
from this import d


def make_profile_tab(core_snps_tab: str, output: pathlib.Path = None):
    if output is not None and not output.suffix == ".tab":
        raise ValueError(f"The output file must have tab extension")
    core_snps = pd.read_csv(core_snps_tab, sep="\t", index_col=0)
    samples = [col for col in core_snps.columns if col not in ["REF", "POS"]]
    profile = core_snps.loc[:, samples].transpose()
    profile.columns.name = None
    if output is not None:
        profile.to_csv(output, sep="\t")
    return profile


def main():
    parser = argparse.ArgumentParser(
        description="Making profile for using in grapetree."
    )
    parser.add_argument(
        "-i",
        "--input",
        help="Input core_snps.tab as obtained from snippy-core",
        type=pathlib.Path,
        default=None,
        required=True,
    )
    parser.add_argument(
        "-o",
        "--output",
        help="Output tab file where the profile will be stored.  ",
        type=pathlib.Path,
        default=None,
        required=True,
    )
    args = parser.parse_args()
    make_profile_tab(args.input, output=args.output)


if __name__ == "__main__":
    main()
