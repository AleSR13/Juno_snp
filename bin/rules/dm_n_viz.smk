rule make_tree:
    input: output_dir.joinpath('snp_analysis', 'core_snps.vcf')
    output: 
        tree = output_dir.joinpath('tree', 'newick_tree.txt')
    #container: "docker://andersenlab/vcf-kit:20200822175018b7b60d"
    conda:
        "../../envs/vcfkit.yaml"
    message: "Making tree..."
    log:
        log_dir.joinpath('making_tree.log')
    threads: config['threads']['vcfkit']
    resources: mem_gb=config['mem_gb']['vcfkit']
    params:
        dm = output_dir.joinpath('tree', 'distance_matrix.tab'),
        algorithm = config['tree']['algorithm']
    shell:
        """
vk phylo tree {params.algorithm} {input} > {output.tree} 2> {log}
        """

rule get_dm:
    input: output_dir.joinpath('tree', 'newick_tree.txt')
    output: output_dir.joinpath('tree', 'distance_matrix.csv')
    message: "Getting distance matrix..."
    log:
        log_dir.joinpath('get_distance_matrix.log')
    threads: config['threads']['other']
    resources: mem_gb=config['mem_gb']['other']
    params:
        dm = output_dir.joinpath('tree', 'distance_matrix.tab')
    shell:
        '''
python bin/newick2dm.py -i {input} -o {output}
        '''

rule get_snp_matrix:
    input:
        output_dir.joinpath("snp_analysis", "core_snps.full.aln")
    output:
        snp_matrix = output_dir.joinpath("tree", "snp_matrix.csv")
    conda:
        "../../envs/snp_dists.yaml"
    message: "Making SNP matrix"
    log:
        log_dir.joinpath("snp_matrix.log")
    threads: config['threads']['other']
    resources: mem_gb=config['mem_gb']['other']
    shell:
        '''
        snp-dists -c {input} 1>{output.snp_matrix} 2>{log}
        '''