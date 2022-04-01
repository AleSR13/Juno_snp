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

1. Two ‘.fastq’ files (paired-end sequencing) derived from short-read sequencing. They should be already filtered and trimmed (for instance, with the [Juno-assembly pipeline](https://github.com/RIVM-bioinformatics/Juno_pipeline)).
2. An assembly from the same sample in the form of a single ‘.fasta’ file.

Importantly, the Juno-SNP pipeline works directly on output generated by the [Juno-assembly pipeline](https://github.com/AleSR13/Juno_pipeline).

The Juno-SNP pipeline will then perform the following steps:  

1. If no reference genome is provided, it will look for one in the NCBI database. This is done with the help of the [referenceseeker package](https://github.com/oschwengers/referenceseeker)
2. If step 1 was performed, it will download the appropriate genome from the NCBI database
3. It will perform the SNP calling on every sample that was given using [Snippy](https://github.com/tseemann/snippy)
4. It will calculate the distance matrix and produce a Newick file for the "Neighbor Joining" tree for the given samples using [GrapeTree](https://achtman-lab.github.io/GrapeTree/).

## Prerequisities

* **Linux + conda** A Linux-like environment with at least 'miniconda' installed. 
* Preferentially also **Singularity**. See instructions for running the pipeline if you don't have singularity installed.


## Installation

1. Clone the repository:

```
git clone https://github.com/RIVM-bioinformatics/Juno_snp.git
```
Alternatively, you can download it manually as a zip file (you will need to unzip it then).

2. Enter the directory with the pipeline and install the master environment:

```
cd Juno_snp
conda env update -f envs/master_env.yaml
```

## Parameters & Usage

### Command for help

* ```-h, --help``` Shows the help of the pipeline

### Required parameters

* ```-i, --input``` Directory with the input (fasta) files. The fasta files should be all in this directory (no subdirectories) and have the extension '.fasta'. 

### Optional parameters

* ```-o --output``` Directory (if not existing it will be created) where the output of the pipeline will be collected. The default behavior is to create a folder called 'output' within the pipeline directory. 
* ```-d --db_dir``` Directory (if not existing it will be created) where the databases used by this pipeline will be downloaded or where they are expected to be present. Default is '/mnt/db/juno/typing_db' (internal RIVM path to the databases of the Juno pipelines). It is advisable to provide your own path if you are not working inside the RIVM Linux environment.
* `-a --ani` ANI threshold. Passed to referenceseeker. Default: 0.95
* `-cd INT, --conserved-dna` Conserved DNA threshold. Passed to referenceseeker. Default: 0.69
* `-sw --sliding-window` Sliding window - the lower the more accurate but also slower. Passed to referenceseeker. Default: 400
* `--no-containers` Use conda environments instead of containers. Use this option if you don't have singularity installed
* `-p --prefix` Conda or singularity prefix. Path to the place where you want to store the conda environments or the singularity images.
* `-c --cores` Number of cores to use. Default is 300
* `-q --queue` Name of the queue that the job will be submitted to if working on a cluster.
* `-l --local` Running pipeline locally (instead of in a computer cluster). Default is running it in a cluster.
* `-w --time-limit` Time limit per job in minutes (passed as -W argument to bsub). Jobs will be killed if not finished in this time. Default: 60
* `-u --unlock` Unlock output directory (passed to snakemake).
* `-n --dryrun` Dry run printing steps to be taken in the pipeline without actually running it (passed to snakemake).
* `--rerunincomplete` Re-run jobs if they are marked as incomplete (passed to snakemake).
* `--snakemake-args` Extra arguments to be passed to [snakemake API](https://snakemake.readthedocs.io/en/stable/api_reference/snakemake.html).

### The base command to run this program. 

```
python juno_snp.py -i [path/to/input_directory] 
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
* Any issue can be reported in the [Issues section](https://github.com/RIVM-bioinformatics/Juno-typing/issues) of this repository.

## Future ideas for this pipeline

* -

## License
This pipeline is licensed with an AGPL3 license. Detailed information can be found inside the 'LICENSE' file in this repository.

## Contact
* **Contact person:**       Alejandra Hernández Segura
* **Email**                 alejandra.hernandez.segura@rivm.nl
