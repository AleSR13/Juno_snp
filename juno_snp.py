"""
Juno-SNP pipeline
Authors: Alejandra Hernandez-Segura
Organization: Rijksinstituut voor Volksgezondheid en Milieu (RIVM)
Department: Infektieziekteonderzoek, Diagnostiek en Laboratorium
            Surveillance (IDS), Bacteriologie (BPD)     
Date: 07-02-2022   

Documentation: https://rivm-bioinformatics.github.io/ids_bacteriology_man/juno-snp.html 

"""

import argparse
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable, Union

import yaml

# Dependencies
from juno_library import Pipeline

from version import __description__, __package_name__, __version__


def main() -> None:
    juno_snp = JunoSnp()
    juno_snp.run()


def check_number_within_range(
    minimum: float = 0, maximum: float = 1
) -> Union[Callable[[str], str], argparse.FileType]:
    """
    Creates a function to check whether a numeric value is within a range, inclusive.

    The generated function can be used by the `type` parameter in argparse.ArgumentParser.
    See https://stackoverflow.com/a/53446178.

    Args:
        value: the numeric value to check.
        minimum: minimum of allowed range, inclusive.
        maximum: maximum of allowed range, inclusive.

    Returns:
        A function which takes a single argument and checks this against the range.

    Raises:
        argparse.ArgumentTypeError: if the value is outside the range.
        ValueError: if the value cannot be converted to float.
    """

    def generated_func_check_range(value: str) -> str:
        value_f = float(value)
        if (value_f < minimum) or (value_f > maximum):
            raise argparse.ArgumentTypeError(
                f"Supplied value {value} is not within expected range {minimum} to {maximum}."
            )
        return str(value)

    return generated_func_check_range


@dataclass
class JunoSnp(Pipeline):
    pipeline_name: str = __package_name__
    pipeline_version: str = __version__
    input_type: str = "both"

    def _add_args_to_parser(self) -> None:
        super()._add_args_to_parser()

        self.parser.description = __description__

        self.add_argument(
            "-r",
            "--reference",
            type=Path,
            required=False,
            metavar="FILE",
            help="Relative or absolute path to a reference genome. Can be FASTA or GenBank",
        )
        self.add_argument(
            "-d",
            "--db-dir",
            type=Path,
            metavar="DIR",
            default="/mnt/db/juno/snp",
            help="Relative or absolute path to the database directory. If non is given, /mnt/db/juno/snp"
            " (where the default db resides at the RIVM will be used).",
        )
        self.add_argument(
            "-a",
            "--ani",
            type=check_number_within_range(0, 1),
            metavar="INT",
            default=0.95,
            help="ANI threshold. Passed to referenceseeker",
        )
        self.add_argument(
            "-cd",
            "--conserved-dna",
            type=check_number_within_range(0, 1),
            metavar="INT",
            default=0.69,
            help="Conserved DNA threshold. Passed to referenceseeker",
        )
        self.add_argument(
            "-sw",
            "--sliding-window",
            type=check_number_within_range(100, 1000),
            metavar="INT",
            default=400,
            help="Sliding window - the lower the more accurate but also slower. Passed to referenceseeker",
        )
        self.add_argument(
            "-kl",
            "--kmer-length",
            type=check_number_within_range(1, 32),
            metavar="INT",
            default=21,
            help="K-mer length - longer kmers increase specificity, shorter kmers increase sensitivity. Passed to mash sketch",
        )
        self.add_argument(
            "-ss",
            "--sketch-size",
            type=int,
            metavar="INT",
            default=1000,
            help="Sketch size - larger sketch size better represents the original sequence,"
            " but leads to large files and longer running time. Passed to mash sketch",
        )
        self.add_argument(
            "-mt",
            "--mash-threshold",
            type=check_number_within_range(0, 1),
            metavar="FLOAT",
            default=0.01,
            help="Mash threshold - maximum mash distance to consider genomes similar. Passed to preclustering script.",
        )
        self.add_argument(
            "-t",
            "--tree-algorithm",
            type=str,
            metavar="ALGORITHM",
            default="upgma",
            choices=["upgma", "nj"],
            help="Algorithm to use for making the tree. It can be 'upgma' or 'nj' (neighbor-joining). Default is upgma",
        )
        self.add_argument(
            "--snippy-report",
            action="store_true",
            help="If set, the pipeline will run Snippy with --report option. This requires quite some extra time and storage, especially if there a lot of variants",
        )

    def _parse_args(self) -> argparse.Namespace:
        args = super()._parse_args()

        self.reference: Path = args.reference
        self.db_dir: Path = args.db_dir
        self.ani: float = args.ani
        self.conserved_dna: float = args.conserved_dna
        self.sliding_window: int = args.sliding_window
        self.kmer_length: int = args.kmer_length
        self.sketch_size: int = args.sketch_size
        self.mash_threshold: float = args.mash_threshold
        self.tree_algorithm: str = args.tree_algorithm
        self.snippy_report: bool = args.snippy_report
        self.dryrun: bool = args.dryrun

        return args

    def setup(self) -> None:
        super().setup()
        if self.snakemake_args["use_singularity"]:
            self.snakemake_args["singularity_args"] = " ".join(
                [
                    self.snakemake_args["singularity_args"],
                    f"--bind {self.db_dir}:{self.db_dir}",
                ]  # paths that singularity should be able to read from can be bound by adding to the above list
            )

        with open(
            Path(__file__).parent.joinpath("config/pipeline_parameters.yaml")
        ) as f:
            parameters_dict = yaml.safe_load(f)
        self.snakemake_config.update(parameters_dict)

        if self.snippy_report:
            snippy_report_cmd = "--report"
        else:
            snippy_report_cmd = ""

        self.user_parameters = {
            "input_dir": str(self.input_dir),
            "output_dir": str(self.output_dir),
            "exclusion_file": str(self.exclusion_file),
            "db_dir": str(self.db_dir),
            "reference": str(self.reference),
            "use_singularity": str(self.snakemake_args["use_singularity"]),
            "dryrun": self.dryrun,
            "referenceseeker": {
                "db": str(self.db_dir.joinpath("bacteria-refseq")),
                "ani_threshold": self.ani,
                "conserved_dna_threshold": self.conserved_dna,
                "sliding_window": self.sliding_window,
            },
            "mash": {
                "kmer_length": self.kmer_length,
                "sketch_size": self.sketch_size,
                "threshold": self.mash_threshold,
            },
            "tree": {"algorithm": self.tree_algorithm},
            "snippy": {"report": snippy_report_cmd},
        }


if __name__ == "__main__":
    main()
