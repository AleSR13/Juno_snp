checkpoint preclustering:
    input:
        [SAMPLES[sample]["assembly"] for sample in SAMPLES],
    output:
        yaml=output_dir.joinpath("preclustering", "mock_clusters.yaml"),
    message:
        "Creating mock cluster to force mapping on single reference genome."
    log:
        log_dir.joinpath("mock_clustering.log"),
    conda:
        "../../envs/preclustering.yaml"
    threads: config["threads"]["other"]
    resources:
        mem_gb=config["mem_gb"]["other"],
    params:
        script=srcdir("../../bin/mock_cluster.py"),
    shell:
        """
python {params.script} --input {input} --output {output} 2>&1>{log}
        """
