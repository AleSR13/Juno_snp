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
import pathlib
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

GIVEN_REF=config['ref']

#@################################################################################
#@####                         Expected output                               #####
#@################################################################################
output_dir = pathlib.Path(config["out"])
log_dir = output_dir.joinpath('log')
db_dir = pathlib.Path(config["db_dir"])
mash_db = db_dir.joinpath('bacteria-refseq', 'db.msh')
referenceseeker_md5 = str(db_dir.joinpath('bacteria-refseq', 'downloaded_db.txt'))
#scores_refseq_candidates = output_dir.joinpath('ref_genome_used', 'scores_refseq_candidates.csv')

if config['dryrun'] is True and GIVEN_REF is not None:
    ref_genome = GIVEN_REF
else:
    ref_genome = output_dir.joinpath('ref_genomes_used', 'cluster_{cluster}')

if GIVEN_REF is not None and not ref_genome.exists():
    output_dir.mkdir(exist_ok=True)
    ref_dir = ref_genome.parent
    ref_dir.mkdir(exist_ok=True)
    copyfile(GIVEN_REF, ref_genome)
    ## Force mapping of all samples on GIVEN_REF in this case

def get_vcf_clusters(cluster):
    with open(checkpoints.preclustering.get(**cluster).output[0]) as file:
        SAMPLE_CLUSTERS = yaml.safe_load(file)
    CLUSTERS = set([ cluster for sample, cluster in SAMPLE_CLUSTERS.items() ])
    return expand(output_dir.joinpath('snp_analysis', 'snippy-core', 'cluster_{cluster}', 'core_snps.vcf'), cluster=CLUSTERS)

def get_tree_clusters(cluster):
    with open(checkpoints.preclustering.get(**cluster).output[0]) as file:
        SAMPLE_CLUSTERS = yaml.safe_load(file)
    CLUSTERS = set([ cluster for sample, cluster in SAMPLE_CLUSTERS.items() ])
    return expand(output_dir.joinpath('tree/cluster_{cluster}/newick_tree.txt'), cluster=CLUSTERS)

#@################################################################################
#@####                              Processes                                #####
#@################################################################################

include: "bin/rules/pre_cluster.smk"
include: "bin/rules/find_reference.smk"
include: "bin/rules/snp_analysis.smk"
include: "bin/rules/dm_n_viz.smk"

#@################################################################################
#@####              Finalize pipeline (error/success)                        #####
#@################################################################################

onerror:
    shell("""
rm -f tmp*npy
rm -f tmp*_fastme_stat.txt
rm -f tmp*_fastme_tree.nwk
rm -f tmp*dist.list
echo -e "Something went wrong with Juno-SNP pipeline. Please check the logging files in {output_dir}/log/"
    """)


#################################################################################
#####                       Specify final output                            #####
#################################################################################

localrules:
    all

rule all:
    input:
        #ref_genome,
        # expand(
        #     output_dir.joinpath('snp_analysis', '{sample}', 'snps.tab'), 
        #     sample=SAMPLES
        # ),
        # get_vcf_clusters,
        get_tree_clusters,
        # output_dir.joinpath('snp_analysis', 'core_snps.vcf'),
        # output_dir.joinpath('tree', '{cluster}', 'distance_matrix.csv'),
        # output_dir.joinpath('tree', '{cluster}', 'newick_tree.txt'),
        # output_dir.joinpath('tree', '{cluster}', 'snp_matrix.csv'),
        output_dir.joinpath('preclustering', 'clusters.yaml')

