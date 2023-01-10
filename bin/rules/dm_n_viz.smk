rule make_tree:
    input: output_dir.joinpath('snp_analysis', 'snippy-core', 'cluster_{cluster}')
    output: 
        tree = output_dir.joinpath('tree', 'cluster_{cluster}', 'newick_tree.txt')
    #container: "docker://andersenlab/vcf-kit:20200822175018b7b60d"
    conda:
        "../../envs/vcfkit.yaml"
    message: "Making tree..."
    log:
        log_dir.joinpath('making_tree_cluster_{cluster}.log')
    threads: config['threads']['vcfkit']
    resources: mem_gb=config['mem_gb']['vcfkit']
    params:
        input = lambda wildcards: output_dir.joinpath('snp_analysis', 'snippy-core', f'cluster_{wildcards.cluster}', 'core_snps.vcf'),
        algorithm = config['tree']['algorithm']
    shell:
        """
vk phylo tree {params.algorithm} {params.input} > {output.tree} 2> {log}
        """

rule get_dm:
    input: output_dir.joinpath('tree', 'cluster_{cluster}', 'newick_tree.txt')
    output: output_dir.joinpath('tree', 'cluster_{cluster}', 'distance_matrix.csv')
    message: "Getting distance matrix..."
    log:
        log_dir.joinpath('get_distance_matrix_cluster_{cluster}.log')
    threads: config['threads']['other']
    resources: mem_gb=config['mem_gb']['other']
    params:
        dm = lambda wildcards: output_dir.joinpath('tree', f'cluster_{wildcards.cluster}', 'distance_matrix.tab')
    shell:
        '''
python bin/newick2dm.py -i {input} -o {output}
        '''

rule get_snp_matrix:
    input:
        output_dir.joinpath("snp_analysis", 'snippy-core', 'cluster_{cluster}')
    output:
        snp_matrix = output_dir.joinpath("tree", 'cluster_{cluster}', "snp_matrix.csv")
    conda:
        "../../envs/snp_dists.yaml"
    message: "Making SNP matrix"
    log:
        log_dir.joinpath("snp_matrix_cluster_{cluster}.log")
    threads: config['threads']['other']
    resources: mem_gb=config['mem_gb']['other']
    params:
        input = lambda wildcards: output_dir.joinpath('snp_analysis', 'snippy-core', f'cluster_{wildcards.cluster}', 'core_snps.full.aln'),
    shell:
        '''
snp-dists -c {params.input} 1>{output.snp_matrix} 2>{log}
        '''