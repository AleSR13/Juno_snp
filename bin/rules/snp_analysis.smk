import yaml

def get_mapped_per_cluster(cluster):
    with open(checkpoints.preclustering.get(**cluster).output[0]) as file:
        SAMPLE_CLUSTERS = yaml.safe_load(file)
    # print(SAMPLE_CLUSTERS)
    SELECTED_SAMPLES = []
    for sample, sample_cluster in SAMPLE_CLUSTERS.items():
        # print(f"comparing sample_cluster {sample_cluster} to cluster {cluster}")
        if str(sample_cluster) == str(cluster):
            # print("match found")
            SELECTED_SAMPLES.append(sample)
            # print(f"SELECTED_SAMPLES now contains {SELECTED_SAMPLES}")
    # print(expand(output_dir.joinpath('snp_analysis', 'cluster_{cluster}', '{sample}'), sample=SELECTED_SAMPLES, allow_missing=True))
    return expand(output_dir.joinpath('snp_analysis', 'cluster_{cluster}', '{sample}'), sample=SELECTED_SAMPLES, allow_missing=True)

rule snp_analysis:
    input: 
        r1 = lambda wildcards: SAMPLES[wildcards.sample]["R1"],
        r2 = lambda wildcards: SAMPLES[wildcards.sample]["R2"],
        ref = output_dir.joinpath("ref_genomes_used/cluster_{cluster}/ref_genome.fasta")
    output: 
        res = directory(output_dir.joinpath('snp_analysis', 'cluster_{cluster}', '{sample}'))
    message: "Running snippy on sample {wildcards.sample}."
    log:
        log_dir.joinpath('snp_analysis', 'cluster_{cluster}', 'snippy_{sample}.log')
    #container: 'docker://staphb/snippy:4.6.0-SC2'
    conda:
        "../../envs/snippy.yaml"  
    threads: config['threads']['snippy']
    resources: mem_gb=config['mem_gb']['snippy']
    params:
        # out_dir = lambda wildcards: output_dir.joinpath('snp_analysis', wildcards.sample),
        mincov = 10,
        minfrac = 0.9,
        minqual = 100,
        mapqual = 60,
        basequal = 13,
        maxsoft = 'x'
    shell:
        """
snippy --cpus {threads} \
    --outdir {output} \
    --ref {input.ref} \
    --R1 {input.r1} \
    --R2 {input.r2} \
    --report \
    --force &> {log}
        """


checkpoint snp_core:
    input: 
        samples = get_mapped_per_cluster,
        ref = output_dir.joinpath("ref_genomes_used/cluster_{cluster}/ref_genome.fasta")
        # snps = expand(output_dir.joinpath('snp_analysis', '{sample}', 'snps.tab'), sample=SAMPLES),
        # ref = ref_genome
    output: 
        directory(output_dir.joinpath('snp_analysis/snippy-core/cluster_{cluster}'))
        # aln = output_dir.joinpath('snp_analysis', 'core_snps.aln'),
        # full_aln = output_dir.joinpath('snp_analysis', 'core_snps.full.aln'),
        # fa = output_dir.joinpath('snp_analysis', 'core_snps.ref.fa'),
        # tab = output_dir.joinpath('snp_analysis', 'core_snps.tab'),
        # txt = output_dir.joinpath('snp_analysis', 'core_snps.txt'),
        # vcf = output_dir.joinpath('snp_analysis', 'core_snps.vcf')
    message: "Getting SNP core."
    log:
        log_dir.joinpath('snp_analysis', 'snippy_core', 'cluster_{cluster}.log')
    conda:
        "../../envs/snippy.yaml"
    #container: 'docker://staphb/snippy:4.6.0-SC2'
    threads: config['threads']['snippy']
    resources: mem_gb=config['mem_gb']['snippy']
    # params:
    #     prefix = str(output_dir.joinpath('snp_analysis', 'core_snps'))
    shell:
        """
mkdir -p {output}
snippy-core --ref {input.ref} --prefix {output}/core_snps {input.samples} &> {log}
        """
