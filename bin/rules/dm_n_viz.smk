rule make_profile_tab:
    input: output_dir.joinpath('snp_analysis', 'core_snps.tab')
    output: temporary(output_dir.joinpath('tree', 'profile.tab'))
    container: "oras://ghcr.io/alesr13/pandas_ncbidatasets:v0.1"
    message: "Making profile.tab for all samples with core SNPs included..."
    log:
        log_dir.joinpath('grapetree.log')
    threads: config['threads']['other']
    resources: mem_gb=config['mem_gb']['other']
    shell:
        """
python bin/get_profile_for_tree.py -i {input} -o {output} &> {log}
        """

rule make_tree:
    input: output_dir.joinpath('tree', 'profile.tab')
    output: 
        dm = output_dir.joinpath('tree', 'distance_matrix.tab'),
        tree = output_dir.joinpath('tree', 'newick_tree.txt')
    container: "docker://quay.io/biocontainers/grapetree:2.1--pyh3252c3a_0"
    message: "Making tree..."
    log:
        log_dir.joinpath('grapetree.log')
    threads: config['threads']['grapetree']
    resources: mem_gb=config['mem_gb']['grapetree']
    shell:
        """
# Distance matrix (no tree)
grapetree --profile {input} \
    --method "distance" \
    --missing 0 \
    --n_proc {threads} > {output.dm} 2> {log}

# Make tree
grapetree --profile {input} \
    --method "NJ" \
    --matrix "symmetric" \
    --missing 0 \
    --n_proc {threads} > {output.tree} 2> {log}
        """