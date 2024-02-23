"""
Juno-SNP
Author: Alejandra Hernandez-Segura
Organization: Rijksinstituut voor Volksgezondheid en Milieu (RIVM)
Department: Infektieziekteonderzoek, Diagnostiek en Laboratorium Surveillance (IDS), Bacteriologie (BPD)
Date: 07-03-2022
"""
#################################################################################
##### Import config file, sample_sheet and set output folder names          #####
#################################################################################

from yaml import safe_load
from pathlib import Path
from shutil import copyfile

#################################################################################
#####                   Load samplesheet and config params                  #####
#################################################################################

# Loading sample sheet as dictionary
# ("R1" and "R2" keys for fastq)
sample_sheet = config["sample_sheet"]
SAMPLES = {}
with open(sample_sheet) as sample_sheet_file:
    SAMPLES = safe_load(sample_sheet_file)

GIVEN_REF = config["reference"]

#################################################################################
#####                         Expected output                               #####
#################################################################################
output_dir = Path(config["output_dir"])
log_dir = output_dir.joinpath("log")
db_dir = Path(config["db_dir"])
mash_db = db_dir.joinpath("bacteria-refseq", "db.msh")
referenceseeker_md5 = str(db_dir.joinpath("bacteria-refseq", "downloaded_db.txt"))

if (config["dryrun"] is True) and (GIVEN_REF != "None"):
    ref_genome = GIVEN_REF
else:
    ref_genome = output_dir.joinpath(
        "ref_genomes_used", "cluster_1", "ref_genome.fasta"
    )

# GIVEN_REF is converted to str
if (GIVEN_REF != "None") and (not ref_genome.exists()):
    print(f"Copying reference genome {GIVEN_REF} to {ref_genome}")
    output_dir.mkdir(exist_ok=True, parents=True)
    ref_dir = ref_genome.parent
    ref_dir.mkdir(exist_ok=True, parents=True)
    copyfile(GIVEN_REF, ref_genome)


def get_output_per_cluster(cluster):
    with open(checkpoints.preclustering.get(**cluster).output[0]) as file:
        SAMPLE_CLUSTERS = yaml.safe_load(file)
    CLUSTERS = set([cluster for sample, cluster in SAMPLE_CLUSTERS.items()])
    output_files = expand(
        output_dir.joinpath("tree/cluster_{cluster}/{file}"),
        cluster=CLUSTERS,
        file=["newick_tree.txt", "snp_matrix.csv"],
    )
    output_iqtree = expand(
        output_dir.joinpath("ml_tree", "cluster_{cluster}"), cluster=CLUSTERS
    )
    output_qc = expand(
        output_dir.joinpath("qc", "cluster_{cluster}", "multiqc", "multiqc.html"),
        cluster=CLUSTERS,
    )
    return output_files + output_iqtree + output_qc


#################################################################################
#####                              Processes                                #####
#################################################################################

if GIVEN_REF != "None":

    include: "bin/rules/mock_cluster.smk"

else:

    include: "bin/rules/pre_cluster.smk"
    include: "bin/rules/find_reference.smk"


include: "bin/rules/snp_analysis.smk"
include: "bin/rules/dm_n_viz.smk"
include: "bin/rules/qc.smk"


#################################################################################
#####              Finalize pipeline (error/success)                        #####
#################################################################################


onerror:
    shell(
        """
rm -f tmp*npy
rm -f tmp*_fastme_stat.txt
rm -f tmp*_fastme_tree.nwk
rm -f tmp*dist.list
echo -e "Something went wrong with Juno-SNP pipeline. Please check the logging files in {output_dir}/log/"
    """
    )


#################################################################################
#####                       Specify final output                            #####
#################################################################################


localrules:
    all,


rule all:
    input:
        get_output_per_cluster,
