"""
Juno-SNP pipeline
Authors: Alejandra Hernandez-Segura
Organization: Rijksinstituut voor Volksgezondheid en Milieu (RIVM)
Department: Infektieziekteonderzoek, Diagnostiek en Laboratorium
            Surveillance (IDS), Bacteriologie (BPD)     
Date: 07-02-2022   

Documentation: https://rivm-bioinformatics.github.io/ids_bacteriology_man/juno-snp.html 

"""

# Dependencies
from base_juno_pipeline import base_juno_pipeline, helper_functions
import argparse
import pathlib
import subprocess
import yaml
from version import __package_name__, __version__, __description__

class JunoSnpRun(
    base_juno_pipeline.PipelineStartup, base_juno_pipeline.RunSnakemake
):
    """Class with the arguments and specifications that are only for the Juno-typing pipeline but inherit from PipelineStartup and RunSnakemake"""
    pipeline_name: str = __package_name__
    pipeline_version: str = __version__

    def __init__(self, 
                input_dir, 
                ref,
                output_dir, 
                db_dir='/mnt/db/juno/snp',
                ani=0.95,
                conserved_dna=0.69,
                sliding_window=400,
                kmer_length=21,
                sketch_size=1000,
                mash_threshold=0.01,
                tree_algorithm='upgma',
                cores=300,
                time_limit=60,
                local=False,
                queue='bio',
                unlock=False,
                rerunincomplete=False,
                dryrun=False,
                run_in_container=True,
                prefix=None,
                **kwargs):
        """Initiating Juno-SNP pipeline"""
        
        # Get proper file paths
        output_dir = pathlib.Path(output_dir).resolve()
        if ref is not None:
            self.ref = pathlib.Path(ref).resolve()
        else:
            self.ref = None
        self.db_dir = pathlib.Path(db_dir).resolve()
        self.ani_threshold=float(ani)
        self.conserved_dna_threshold=float(conserved_dna)
        self.sliding_window=int(sliding_window)
        self.kmer_length=int(kmer_length)
        self.sketch_size=int(sketch_size)
        self.mash_threshold=float(mash_threshold)
        if tree_algorithm not in ['upgma', 'nj']:
            raise ValueError(
                f'The provided tree algorithm {tree_algorithm} is not supported.' \
                ' Please choose upgma or nj'
            )
        self.tree_algorithm=tree_algorithm
        workdir = pathlib.Path(__file__).parent.resolve()
        self.path_to_audit = output_dir.joinpath('audit_trail')
        base_juno_pipeline.PipelineStartup.__init__(
            self,
            input_dir=pathlib.Path(input_dir).resolve(), 
            input_type='both',
            min_num_lines=2
        ) # Min for viable fasta
        base_juno_pipeline.RunSnakemake.__init__(
            self,
            pipeline_name='Juno_SNP',
            pipeline_version='v0.1',
            output_dir=output_dir,
            workdir=workdir,
            cores=cores,
            time_limit=time_limit,
            local=local,
            queue=queue,
            unlock=unlock,
            rerunincomplete=rerunincomplete,
            dryrun=dryrun,
            useconda=not run_in_container,
            conda_prefix=prefix,
            usesingularity=run_in_container,
            singularityargs=f"--bind {self.input_dir}:{self.input_dir} --bind {output_dir}:{output_dir} --bind {self.db_dir}:{self.db_dir}",
            singularity_prefix=prefix,
            restarttimes=1,
            latency_wait=60,
            name_snakemake_report=str(self.path_to_audit.joinpath('juno_snp_report.html')),
            **kwargs
        )

        # Start pipeline
        self.run_juno_snp_pipeline()

    def start_juno_snp_pipeline(self):
        """
        Function to start the pipeline
        """
        self.start_juno_pipeline()
        with open(self.sample_sheet, 'w') as file_:
            yaml.dump(self.sample_dict, file_, default_flow_style=False)
    
    def write_userparameters(self):

        config_params = {
            'input_dir': str(self.input_dir),
            'ref': self.ref,
            'out': str(self.output_dir),
            'db_dir': str(self.db_dir),
            'referenceseeker': {
                'db': str(self.db_dir.joinpath('bacteria-refseq')),
                'ani_threshold': self.ani_threshold,
                'conserved_dna_threshold': self.conserved_dna_threshold,
                'sliding_window': self.sliding_window
            },
            'mash': {
                'kmer_length': self.kmer_length,
                'sketch_size': self.sketch_size,
                'threshold': self.mash_threshold
            },
            'tree': {
                'algorithm' : self.tree_algorithm
            },
            'dryrun': self.dryrun
        }
        
        with open(self.user_parameters, 'w') as file_:
            yaml.dump(config_params, file_, default_flow_style=False)

        return config_params     
        
    def run_juno_snp_pipeline(self):
        self.start_juno_snp_pipeline()
        self.user_params = self.write_userparameters()
        self.get_run_info()
        if not self.dryrun or self.unlock:
            self.path_to_audit.mkdir(parents=True, exist_ok=True)

        self.successful_run = self.run_snakemake()
        assert self.successful_run, f'Please check the log files'
        if not self.dryrun or self.unlock:
            subprocess.run(
                ['find', self.output_dir, '-type', 'f', '-empty', '-exec', 'rm', '{}', ';']
            )
            subprocess.run(
                ['find', self.output_dir, '-type', 'd', '-empty', '-exec', 'rm', '-rf', '{}', ';']
            )
            self.make_snakemake_report()

def check_sliding_window(num):
    num_int = int(num)
    if num_int > 99 and num_int < 1001:
        return num_int
    else:
        raise ValueError(
            f'The provided input value for sliding window {str(num)} is not valid.' \
            ' Please provide a number between 100-1000'
        )

def check_between_zero_and_one(num):
    num_f = float(num)
    if num_f >= 0 and num_f <= 1:
        return num_f
    else:
        raise ValueError(
            f'The provided input value {str(num)} is not valid. '\
            'Please provide a value between 0-1'
        )

def check_kmer_length(num):
    num_int = int(num)
    if num_int >= 1 and num_int <= 32:
        return num_int
    else:
        raise ValueError(
            f'The provided input value for kmer length {str(num)} is not valid.' \
            ' Please provide a number between 1-32'
        )

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description = __description__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        "-i",
        "--input",
        type = pathlib.Path,
        required = True,
        metavar = "DIR",
        help = "Relative or absolute path to the input directory. It must either be the output directory of the Juno-assembly pipeline or it must contain all the raw reads (fastq) and assemblies (fasta) files for all samples to be processed."
    )
    parser.add_argument(
        "-r",
        "--reference",
        type = pathlib.Path,
        required = False,
        metavar = "FILE",
        help = "Relative or absolute path to a reference fasta file."
    )
    parser.add_argument(
        "-o",
        "--output",
        type = pathlib.Path,
        metavar = "DIR",
        default = "output",
        help = "Relative or absolute path to the output directory. If non is given, an 'output' directory will be created in the current directory."
    )
    parser.add_argument(
        "-d",
        "--db-dir",
        type = pathlib.Path,
        metavar = "DIR",
        default = "/mnt/db/juno/snp",
        help = "Relative or absolute path to the database directory. If non is given, /mnt/db/juno/snp"\
            " (where the default db resides at the RIVM will be used)."
    )
    parser.add_argument(
        "-a",
        "--ani",
        type = check_between_zero_and_one,
        metavar = "INT",
        default = 0.95,
        help="ANI threshold. Passed to referenceseeker"
    )
    parser.add_argument(
        "-cd",
        "--conserved-dna",
        type = check_between_zero_and_one,
        metavar = "INT",
        default = 0.69,
        help="Conserved DNA threshold. Passed to referenceseeker"
    )
    parser.add_argument(
        "-sw",
        "--sliding-window",
        type = check_sliding_window,
        metavar = "INT",
        default = 400,
        help="Sliding window - the lower the more accurate but also slower. Passed to referenceseeker"
    )
    parser.add_argument(
        "-kl",
        "--kmer-length",
        type = check_kmer_length,
        metavar = "INT",
        default = 21,
        help="K-mer length - longer kmers increase specificity, shorter kmers increase sensitivity. Passed to mash sketch"
    )
    parser.add_argument(
        "-ss",
        "--sketch-size",
        type = int,
        metavar = "INT",
        default = 1000,
        help="Sketch size - larger sketch size better represents the original sequence," \
            " but leads to large files and longer running time. Passed to mash sketch"
    )
    parser.add_argument(
        "-mt",
        "--mash-threshold",
        type = check_between_zero_and_one,
        metavar = "FLOAT",
        default = 0.01,
        help="Mash threshold - maximum mash distance to consider genomes similar. Passed to preclustering script."
    )
    parser.add_argument(
        "-t",
        "--tree-algorithm",
        type = str,
        metavar = "ALGORITHM",
        default = 'upgma',
        choices=['upgma', 'nj'],
        help="Algorithm to use for making the tree. It can be 'upgma' or 'nj' (neighbor-joining). Default is upgma"
    )
    parser.add_argument(
        "--no-containers",
        action = 'store_false',
        help = "Use conda environments instead of containers."
    )
    parser.add_argument(
        "-p",
        "--prefix",
        type = pathlib.Path,
        metavar="PATH",
        default=None,
        help = "Conda or singularity prefix. Path to the place where you want to store the conda environments or the singularity images."
    )
    parser.add_argument(
        "-c",
        "--cores",
        type = int,
        metavar = "INT",
        default = 300,
        help="Number of cores to use. Default is 300"
    )
    #TODO: Get from ${irods_runsheet_sys__runsheet__lsf_queue} if it exists
    parser.add_argument(
        "-q",
        "--queue",
        type = str,
        metavar = "STR",
        default = 'bio',
        help = 'Name of the queue that the job will be submitted to if working on a cluster.'
    )
    parser.add_argument(
        "-l",
        "--local",
        action='store_true',
        help="Running pipeline locally (instead of in a computer cluster). Default is running it in a cluster."
    )
    parser.add_argument(
        "-w",
        "--time-limit",
        type = int,
        metavar = "INT",
        default = 60,
        help="Time limit per job in minutes (passed as -W argument to bsub). Jobs will be killed if not finished in this time."
    )
    # Snakemake arguments
    parser.add_argument(
        "-u",
        "--unlock",
        action = 'store_true',
        help = "Unlock output directory (passed to snakemake)."
    )
    parser.add_argument(
        "-n",
        "--dryrun",
        action='store_true',
        help="Dry run printing steps to be taken in the pipeline without actually running it (passed to snakemake)."
    )
    parser.add_argument(
        "--rerunincomplete",
        action='store_true',
        help="Re-run jobs if they are marked as incomplete (passed to snakemake)."
    )
    parser.add_argument(
        "--snakemake-args",
        nargs='*',
        default={},
        action=helper_functions.SnakemakeKwargsAction,
        help="Extra arguments to be passed to snakemake API (https://snakemake.readthedocs.io/en/stable/api_reference/snakemake.html)."
    )
    parser.add_argument('--version', action='version', version=f'{__package_name__} {__version__}')
    args = parser.parse_args()
    JunoSnpRun(
        input_dir=args.input, 
        ref=args.reference,
        output_dir=args.output,
        db_dir=args.db_dir,
        ani=args.ani,
        conserved_dna=args.conserved_dna,
        sliding_window=args.sliding_window,
        kmer_length=args.kmer_length,
        sketch_size=args.sketch_size,
        mash_threshold=args.mash_threshold,
        tree_algorithm=args.tree_algorithm,
        cores=args.cores,
        local=args.local,
        time_limit=args.time_limit,
        queue=args.queue,
        unlock=args.unlock,
        rerunincomplete=args.rerunincomplete,
        dryrun=args.dryrun,
        run_in_container=args.no_containers,
        prefix=args.prefix,
        **args.snakemake_args
    )
