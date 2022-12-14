rule download_referenceseeker_db:
    output: 
        md5_db = referenceseeker_md5,
        tar_db = temp(db_dir.joinpath('bacteria-refseq.tar.gz')),
        mash_db = db_dir.joinpath('bacteria-refseq', 'db.msh')
    message: "Downloading the bacterial referenceseeker database."
    log: log_dir.joinpath("download_referenceseeker_db.log")
    threads: config['threads']['other']
    resources: mem_gb=config['threads']['other']
    params:
        url="https://zenodo.org/record/4415843/files/bacteria-refseq.tar.gz"
    shell: 
        """
db_dir="$(dirname {output.tar_db})"
cd "${{db_dir}}"
wget "{params.url}" &> {log}
md5sum_db=$(md5sum bacteria-refseq.tar.gz)
tar -xvzf bacteria-refseq.tar.gz >> {log}
echo "Downloaded database from {params.url}" > {output} 2> {log}
echo "md5sum: ${{md5sum_db}}" >> {output} 2> {log}
        """


rule find_reference:
    input:
        fasta = lambda wildcards: SAMPLES[wildcards.sample]["assembly"],
        md5_db = referenceseeker_md5
    output:
        output_dir.joinpath('find_reference', 'referenceseeker_{sample}.tab')
    message: "Running referenceseeker for sample {wildcards.sample}."
    log:
        log_dir.joinpath('referenceseeker', '{sample}.log')
    conda:
        "../../envs/reference_seeker_env.yaml"
    #container: 'docker://pvstodghill/referenceseeker:1.8.0__2022-01-31'
    threads: config['threads']['referenceseeker']
    resources: mem_gb=config['mem_gb']['referenceseeker']
    params:
        db_referenceseeker = config['referenceseeker']['db'],
        ani_threshold = config['referenceseeker']['ani_threshold'],
        conserved_dna_threshold = config['referenceseeker']['conserved_dna_threshold'],
        sliding_window = config['referenceseeker']['sliding_window']
    shell:
        """
referenceseeker --ani {params.ani_threshold} \
    --conserved-dna {params.conserved_dna_threshold} \
    --sliding-window {params.sliding_window} \
    --verbose \
    --threads {threads} \
    {params.db_referenceseeker} {input.fasta} > {output} 2> {log}
        """


rule get_best_ref:
    input: expand(output_dir.joinpath('find_reference', 'referenceseeker_{sample}.tab'), sample=SAMPLES)
    output: 
        scores = scores_refseq_candidates,
        fasta = ref_genome
    message: "Finding best reference genome for dataset."
    log:
        log_dir.joinpath('find_reference.log')
    conda:
        "../../envs/reference_seeker_env.yaml"    
    #container: 'oras://ghcr.io/alesr13/pandas_ncbidatasets:v0.1'
    threads: config['threads']['other']
    resources: mem_gb=config['mem_gb']['other']
    params:
        output_dir = output_dir.joinpath('ref_genome_used')
    shell:
        """
python3 bin/find_best_ref.py --input-files {input} --output {params.output_dir} &> {log}
        """
