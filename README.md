<div align="center">
    <h1>Juno-SNP</h1>
    <br />
    <h2>SNP analysis pipeline.</h2>
    <br />
    <img src="https://via.placeholder.com/150" alt="pipeline logo">
</div>

## Pipeline information

* **Author(s):**            Alejandra Hernández Segura
* **Organization:**         Rijksinstituut voor Volksgezondheid en Milieu (RIVM)
* **Department:**           Infektieziekteonderzoek, Diagnostiek en Laboratorium Surveillance (IDS), Bacteriologie (BPD)
* **Start date:**           07 - 03 - 2022

## About this project

The goal of this pipeline is to perform SNP analysis. For that, it requires:

1. Two ‘.fastq’ files (paired-end sequencing) derived from short-read sequencing. They should be already filtered and trimmed (for instance, with the [Juno-assembly pipeline](https://github.com/RIVM-bioinformatics/juno-assembly)).
2. An assembly from the same sample in the form of a single ‘.fasta’ file.

Importantly, the Juno-SNP pipeline works directly on output generated by the [Juno-Assembly pipeline](https://github.com/RIVM-bioinformatics/juno-assembly).

The Juno-SNP pipeline will then perform the following steps:  

1. If no reference genome is provided, it will look for one in the NCBI database. This is done with the help of the [referenceseeker package](https://github.com/oschwengers/referenceseeker)
2. If step 1 was performed, it will download the appropriate genome from the NCBI database
3. It will perform the SNP calling on every sample that was given using [Snippy](https://github.com/tseemann/snippy)
4. It will calculate the distance matrix and produce a Newick file for the "Neighbor Joining" tree for the given samples using [GrapeTree](https://achtman-lab.github.io/GrapeTree/).

By default, most rules will be run using containers.

## Prerequisities

* **Linux + mamba/conda** A Linux-like environment with at least 'miniconda' installed. 
* Preferentially, **Singularity**. See instructions for running the pipeline if you don't have singularity installed.


## Installation

1. Clone the repository:

```
git clone https://github.com/RIVM-bioinformatics/juno-snp.git
```
Alternatively, you can download it manually as a zip file (you will need to unzip it then).

2. Enter the directory with the pipeline and install the master environment:

```
cd juno-snp
conda env update -f envs/master_env.yaml
```

## Parameters & Usage

### Command for help

* ```-h, --help``` Shows the help of the pipeline

### Required parameters

* ```-i, --input``` Directory with the input (fasta) files. The fasta files should be all in this directory (no subdirectories) and have the extension '.fasta'. 

### Optional parameters

```
  -i DIR, --input DIR   Relative or absolute path to the input directory. It must
                        either be the output directory of the Juno-assembly pipeline
                        or it must contain all the raw reads (fastq) and assemblies
                        (fasta) files for all samples to be processed. (default:
                        None)
  -r FILE, --reference FILE
                        Relative or absolute path to a reference fasta file.
                        (default: None)
  -o DIR, --output DIR  Relative or absolute path to the output directory. If non is
                        given, an 'output' directory will be created in the current
                        directory. (default: output)
  -d DIR, --db-dir DIR  Relative or absolute path to the database directory. If non
                        is given, /mnt/db/juno/snp (where the default db resides at
                        the RIVM will be used). (default: /mnt/db/juno/snp)
  -a INT, --ani INT     ANI threshold. Passed to referenceseeker (default: 0.95)
  -cd INT, --conserved-dna INT
                        Conserved DNA threshold. Passed to referenceseeker (default:
                        0.69)
  -sw INT, --sliding-window INT
                        Sliding window - the lower the more accurate but also slower.
                        Passed to referenceseeker (default: 400)
  -kl INT, --kmer-length INT
                        K-mer length - longer kmers increase specificity, shorter kmers increase sensitivity. Passed to mash sketch (default: 21)
  -ss INT, --sketch-size INT
                        Sketch size - larger sketch size better represents the original sequence, but leads to large files and longer running time. Passed to mash sketch (default: 1000)
  -mt FLOAT, --mash-threshold FLOAT
                        Mash threshold - maximum mash distance to consider genomes similar. Passed to preclustering script. (default: 0.01)
  -t ALGORITHM, --tree-algorithm ALGORITHM
                        Algorithm to use for making the tree. It can be 'upgma' or
                        'nj' (neighbor-joining). Default is upgma (default: upgma)
  --no-containers       Use conda environments instead of containers. (default: True)
  -p PATH, --prefix PATH
                        Conda or singularity prefix. Path to the place where you want
                        to store the conda environments or the singularity images.
                        (default: None)
  -c INT, --cores INT   Number of cores to use. Default is 300 (default: 300)
  -q STR, --queue STR   Name of the queue that the job will be submitted to if
                        working on a cluster. (default: bio)
  -l, --local           Running pipeline locally (instead of in a computer cluster).
                        Default is running it in a cluster. (default: False)
  -w INT, --time-limit INT
                        Time limit per job in minutes (passed as -W argument to
                        bsub). Jobs will be killed if not finished in this time.
                        (default: 60)
  -u, --unlock          Unlock output directory (passed to snakemake). (default:
                        False)
  -n, --dryrun          Dry run printing steps to be taken in the pipeline without
                        actually running it (passed to snakemake). (default: False)
  --rerunincomplete     Re-run jobs if they are marked as incomplete (passed to
                        snakemake). (default: False)
  --snakemake-args [SNAKEMAKE_ARGS ...]
                        Extra arguments to be passed to snakemake API (https://snakem
                        ake.readthedocs.io/en/stable/api_reference/snakemake.html).
                        (default: {})
```

### The base command to run this program. 

If you want the pipeline to choose a reference genome for you, run the command below.

By default, genomes will be clustered using `mash` and genomes that are too different are split into separate clusters.
Mapping will subsequently be performed per cluster. This prevents isolates that are too diverse to be compared.

```
python juno_snp.py -i [path/to/input_directory] 
```

If you want to provide your own reference genome, specify this with `--reference`.
This will force all isolates to be mapped against the provided reference genome.

```
python juno_snp.py -i [path/to/input_directory] --reference [path/to/reference.fasta]
```

### An example on how to run the pipeline.

```
python juno_snp.py -i my_input_files -o my_results --db_dir my_db_dir --local --cores 2
```

## Explanation of the output

* **log:** Log files with output and error files from each Snakemake rule/step that is performed. 
* **audit_trail:** Information about the versions of software and databases used.
* **output per sample:** The pipeline will create one subfolder per each step performed (find_reference, ref_genome_used, snp_analysis, tree).
        
## Issues  

* All default values have been chosen to work with the RIVM Linux environment, therefore, there might not be applicable to other environments (although they should work if the appropriate arguments/parameters are given).
* Any issue can be reported in the [Issues section](https://github.com/RIVM-bioinformatics/juno-snp/issues) of this repository.

## Future ideas for this pipeline

* -

## License
This pipeline is licensed with an AGPL3 license. Detailed information can be found inside the 'LICENSE' file in this repository.

## Contact
* **Contact person:**       Roxanne Wolthuis
* **Email**                 roxanne.wolthuis@rivm.nl
