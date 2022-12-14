rule snp_analysis:
    input: 
        r1 = lambda wildcards: SAMPLES[wildcards.sample]["R1"],
        r2 = lambda wildcards: SAMPLES[wildcards.sample]["R2"]
    output: 
        res = output_dir.joinpath('snp_analysis', '{sample}', 'snps.tab')
    message: "Running reference-free SNP analysis."
    log:
        log_dir.joinpath('snp_analysis', 'snippy_{sample}.log')
    #container: 'docker://biocontainers/parsnp:v1.2dfsg-5-deb_cv1'
    conda:
        "../../envs/snippy.yaml"
    threads: config['threads']['snippy']
    resources: mem_gb=config['mem_gb']['snippy']
    params:
        out_dir = lambda wildcards: output_dir.joinpath('snp_analysis', wildcards.sample),
        mincov = 10,
        minfrac = 0.9,
        minqual = 100,
        mapqual = 60,
        basequal = 13,
        maxsoft = 'x'
    shell:
        """
snippy --cpus {threads} \
    --outdir {params.out_dir} \
    --ref {input.ref} \
    --R1 {input.r1} \
    --R2 {input.r2} \
    --report \
    --force &> {log}
        """


rule snp_core:
    input: 
        snps = expand(output_dir.joinpath('snp_analysis', '{sample}', 'snps.tab'), sample=SAMPLES),
        ref = ref_genome
    output: 
        aln = output_dir.joinpath('snp_analysis', 'core_snps.aln'),
        full_aln = output_dir.joinpath('snp_analysis', 'core_snps.full.aln'),
        fa = output_dir.joinpath('snp_analysis', 'core_snps.ref.fa'),
        tab = output_dir.joinpath('snp_analysis', 'core_snps.tab'),
        txt = output_dir.joinpath('snp_analysis', 'core_snps.txt'),
        vcf = output_dir.joinpath('snp_analysis', 'core_snps.vcf')
    message: "Getting SNP core."
    log:
        log_dir.joinpath('snp_analysis', 'snippy_core.log')
    container: 'docker://biocontainers/parsnp:v1.2dfsg-5-deb_cv1'
    threads: config['threads']['snippy']
    resources: mem_gb=config['mem_gb']['snippy']
    params:
        prefix = str(output_dir.joinpath('snp_analysis', 'core_snps'))
    shell:
        """
input_dirs="{input.snps}"
input_dirs="${{input_dirs//snps.tab}}"
snippy-core --ref {input.ref} --prefix {params.prefix} ${{input_dirs}} &> {log}
        """
 