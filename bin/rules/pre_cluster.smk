rule sketch_genomes:
    input:
        [SAMPLES[sample]['assembly'] for sample in SAMPLES]
    output:
        output_dir.joinpath('preclustering', 'mash_sketches.msh')
    message: "Sketching genome of using mash."
    log:
        log_dir.joinpath('sketch_genomes.log')
    conda:
        "../../envs/mash.yaml"
    threads: config['threads']['mash']
    resources: mem_gb=config['mem_gb']['mash']
    params:
        kmer_length = config['mash']['kmer_length'],
        sketch_size = config['mash']['sketch_size']
    shell:
        """
mash sketch -p {threads} -k {params.kmer_length} -s {params.sketch_size} -o {output} {input} 2>{log}
        """

rule calculate_mash_distances:
    input:
        output_dir.joinpath('preclustering', 'mash_sketches.msh')
    output:
        output_dir.joinpath('preclustering', 'mash_distances.tsv')
    message: "Calculating distances for all sketches using mash."
    log:
        log_dir.joinpath('calculate_mash_distances.log')
    conda:
        "../../envs/mash.yaml"
    threads: config['threads']['mash']
    resources: mem_gb=config['mem_gb']['mash']
    params:
        kmer_length = config['mash']['kmer_length'],
        sketch_size = config['mash']['sketch_size']
    shell:
        """
mash dist -p {threads} {input} {input} > {output} 2>{log}
        """

checkpoint preclustering:
    input:
        output_dir.joinpath('preclustering', 'mash_distances.tsv')
    output:
        yaml = output_dir.joinpath('preclustering', 'clusters.yaml'),
        plot_dir = directory(output_dir.joinpath('preclustering', 'plots'))
    message: "Defining pre-clusters based on mash distances."
    log:
        log_dir.joinpath('preclustering.log')
    conda:
        "../../envs/preclustering.yaml"
    threads: config['threads']['other']
    resources: mem_gb=config['mem_gb']['other']
    params:
        mash_threshold = config['mash']['threshold'],
        script = srcdir("../../bin/preclustering.py")
    shell:
        """
python {params.script} --input {input} --threshold {params.mash_threshold} \
    --output {output.yaml} --plot-output {output.plot_dir} 2>&1>{log}
        """