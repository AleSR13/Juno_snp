rule make_tree:
    input:
        output_dir.joinpath("snp_analysis", "snippy-core", "cluster_{cluster}"),
    output:
        tree=output_dir.joinpath("tree", "cluster_{cluster}", "newick_tree.txt"),
    container:
        "docker://andersenlab/vcf-kit:20200822175018b7b60d"
    conda:
        "../../envs/vcfkit.yaml"
    message:
        "Making tree..."
    log:
        log_dir.joinpath("making_tree_cluster_{cluster}.log"),
    threads: config["threads"]["vcfkit"]
    resources:
        mem_gb=config["mem_gb"]["vcfkit"],
    params:
        input=lambda wildcards: output_dir.joinpath(
            "snp_analysis",
            "snippy-core",
            f"cluster_{wildcards.cluster}",
            "core_snps.vcf",
        ),
        algorithm=config["tree"]["algorithm"],
    shell:
        """
vk phylo tree {params.algorithm} {params.input} > {output.tree} 2> {log}
        """


rule make_ml_tree:
    input:
        output_dir.joinpath("snp_analysis", "snippy-core", "cluster_{cluster}"),
    output:
        directory(output_dir.joinpath("ml_tree", "cluster_{cluster}")),
    container:
        "docker://staphb/iqtree2:2.2.2.6"
    conda:
        "../../envs/iqtree.yaml"
    message:
        "Making ML tree..."
    log:
        log_dir.joinpath("making_ML_tree_cluster_{cluster}.log"),
    threads: config["threads"]["iqtree"]
    resources:
        mem_gb=config["mem_gb"]["iqtree"],
    params:
        prefix="cluster_{cluster}",
    shell:
        """
mkdir -p {output}

NR_SAMPLES=$(grep -c '>' {input}/core_snps.aln)
if [ $NR_SAMPLES -le 2 ]
then
    echo "Not running IQ-tree, does not reach minimal of three samples" > {output}/iqtree_not_started_for_cluster.txt
else
    iqtree \
    -s {input}/core_snps.aln \
    -fconst $(snp-sites -C {input}/core_snps.full.aln) \
    -nt {threads} \
    --prefix {output}/{params.prefix} \
    --seed 1 \
    --mem {resources.mem_gb}G 2>&1>{log}

    if [ ! -f {output}/{params.prefix}.treefile ]
    then
        echo "Treefile is missing, exiting with error now" >>{log}
        exit 1
    fi
fi
        """


rule get_dm:
    input:
        output_dir.joinpath("tree", "cluster_{cluster}", "newick_tree.txt"),
    output:
        output_dir.joinpath("tree", "cluster_{cluster}", "distance_matrix.csv"),
    message:
        "Getting distance matrix..."
    log:
        log_dir.joinpath("get_distance_matrix_cluster_{cluster}.log"),
    threads: config["threads"]["other"]
    resources:
        mem_gb=config["mem_gb"]["other"],
    params:
        dm=lambda wildcards: output_dir.joinpath(
            "tree", f"cluster_{wildcards.cluster}", "distance_matrix.tab"
        ),
    shell:
        """
python bin/newick2dm.py -i {input} -o {output}
        """


rule get_snp_matrix:
    input:
        output_dir.joinpath("snp_analysis", "snippy-core", "cluster_{cluster}"),
    output:
        snp_matrix=output_dir.joinpath("tree", "cluster_{cluster}", "snp_matrix.csv"),
    container:
        "docker://staphb/snp-dists:0.8.2"
    conda:
        "../../envs/snp_dists.yaml"
    message:
        "Making SNP matrix"
    log:
        log_dir.joinpath("snp_matrix_cluster_{cluster}.log"),
    threads: config["threads"]["other"]
    resources:
        mem_gb=config["mem_gb"]["other"],
    params:
        input=lambda wildcards: output_dir.joinpath(
            "snp_analysis",
            "snippy-core",
            f"cluster_{wildcards.cluster}",
            "core_snps.full.aln",
        ),
    shell:
        """
snp-dists -c {params.input} 1>{output.snp_matrix} 2>{log}
        """
