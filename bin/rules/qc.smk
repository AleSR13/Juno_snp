def return_filter_status_per_cluster(cluster):
    with open(checkpoints.preclustering.get(**cluster).output[0]) as file:
        SAMPLE_CLUSTERS = yaml.safe_load(file)
    SELECTED_SAMPLES = []
    for sample, sample_cluster in SAMPLE_CLUSTERS.items():
        if str(sample_cluster) == str(cluster):
            SELECTED_SAMPLES.append(sample)
    return expand(
        output_dir.joinpath(
            "qc",
            "cluster_{cluster}",
            "get_filter_status",
            "{sample}.tsv",
        ),
        sample=SELECTED_SAMPLES,
        allow_missing=True,
    )


def return_multiqc_per_cluster(cluster):
    with open(checkpoints.preclustering.get(**cluster).output[0]) as file:
        SAMPLE_CLUSTERS = yaml.safe_load(file)
    SELECTED_SAMPLES = []
    for sample, sample_cluster in SAMPLE_CLUSTERS.items():
        if str(sample_cluster) == str(cluster):
            SELECTED_SAMPLES.append(sample)
    return (
        expand(
            output_dir.joinpath(
                "qc", "cluster_{cluster}", "insertsize", "{sample}_metrics.txt"
            ),
            sample=SELECTED_SAMPLES,
            allow_missing=True,
        )
        + expand(
            output_dir.joinpath(
                "qc",
                "cluster_{cluster}",
                "CollectAlignmentSummaryMetrics",
                "{sample}.txt",
            ),
            sample=SELECTED_SAMPLES,
            allow_missing=True,
        )
        + expand(
            output_dir.joinpath(
                "qc", "cluster_{cluster}", "CollectWgsMetrics", "{sample}.txt"
            ),
            sample=SELECTED_SAMPLES,
            allow_missing=True,
        )
        + expand(
            output_dir.joinpath(
                "snp_analysis", "cluster_{cluster}", "{sample}", "{sample}.txt"
            ),
            sample=SELECTED_SAMPLES,
            allow_missing=True,
        )
    )


rule picard_CollectInsertSizeMetrics:
    input:
        bam=output_dir.joinpath(
            "snp_analysis", "cluster_{cluster}", "{sample}", "{sample}.bam"
        ),
    output:
        txt=output_dir.joinpath(
            "qc", "cluster_{cluster}", "insertsize", "{sample}_metrics.txt"
        ),
        pdf=output_dir.joinpath(
            "qc", "cluster_{cluster}", "insertsize", "{sample}_report.pdf"
        ),
    message:
        "Calculating insert size for {wildcards.sample}"
    conda:
        "../../envs/gatk_picard.yaml"
    container:
        "docker://broadinstitute/picard:2.27.5"
    log:
        log_dir.joinpath("get_insert_size", "cluster_{cluster}", "{sample}.log"),
    params:
        use_singularity=config["use_singularity"],
    threads: config["threads"]["picard"]
    resources:
        mem_gb=config["mem_gb"]["picard"],
    shell:
        """
if [ {params.use_singularity} == True ]
then
    EXEC=\"java -jar /usr/picard/picard.jar\"
else
    EXEC=picard
fi

$EXEC CollectInsertSizeMetrics \
-I {input.bam} \
-O {output.txt} \
-H {output.pdf} \
--VALIDATION_STRINGENCY LENIENT 2>&1>{log}
        """


rule picard_CollectAlignmentSummaryMetrics:
    input:
        bam=output_dir.joinpath(
            "snp_analysis", "cluster_{cluster}", "{sample}", "{sample}.bam"
        ),
        ref=output_dir.joinpath("ref_genomes_used/cluster_{cluster}/ref_genome.fasta"),
    output:
        txt=output_dir.joinpath(
            "qc",
            "cluster_{cluster}",
            "CollectAlignmentSummaryMetrics",
            "{sample}.txt",
        ),
    conda:
        "../../envs/gatk_picard.yaml"
    container:
        "docker://broadinstitute/picard:2.27.5"
    log:
        log_dir.joinpath(
            "CollectAlignmentSummaryMetrics", "cluster_{cluster}", "{sample}.log"
        ),
    params:
        use_singularity=config["use_singularity"],
    threads: config["threads"]["picard"]
    resources:
        mem_gb=config["mem_gb"]["picard"],
    shell:
        """
if [ {params.use_singularity} == True ]
then
    EXEC=\"java -jar /usr/picard/picard.jar\"
else
    EXEC=picard
fi

$EXEC CollectAlignmentSummaryMetrics \
-I {input.bam} \
-R {input.ref} \
-O {output} \
--VALIDATION_STRINGENCY LENIENT 2>&1>{log}
        """


rule picard_CollectWgsMetrics:
    input:
        bam=output_dir.joinpath(
            "snp_analysis", "cluster_{cluster}", "{sample}", "{sample}.bam"
        ),
        ref=output_dir.joinpath("ref_genomes_used/cluster_{cluster}/ref_genome.fasta"),
    output:
        txt=output_dir.joinpath(
            "qc", "cluster_{cluster}", "CollectWgsMetrics", "{sample}.txt"
        ),
    conda:
        "../../envs/gatk_picard.yaml"
    container:
        "docker://broadinstitute/picard:2.27.5"
    log:
        log_dir.joinpath("CollectWgsMetrics", "cluster_{cluster}", "{sample}.log"),
    params:
        use_singularity=config["use_singularity"],
    threads: config["threads"]["picard"]
    resources:
        mem_gb=config["mem_gb"]["picard"],
    shell:
        """
if [ {params.use_singularity} == True ]
then
    EXEC=\"java -jar /usr/picard/picard.jar\"
else
    EXEC=picard
fi

$EXEC CollectWgsMetrics \
-I {input.bam} \
-R {input.ref} \
-O {output} \
--VALIDATION_STRINGENCY LENIENT 2>&1>{log}
        """


rule multiqc:
    input:
        return_multiqc_per_cluster,
        output_dir.joinpath(
            "snp_analysis", "snippy-core", "cluster_{cluster}", "cluster_{cluster}.txt"
        ),
    output:
        output_dir.joinpath("qc", "cluster_{cluster}", "multiqc", "multiqc.html"),
    message:
        "Making MultiQC report."
    conda:
        "../../envs/multiqc.yaml"
    container:
        "docker://quay.io/biocontainers/multiqc:1.19--pyhdfd78af_0"
    threads: config["threads"]["multiqc"]
    resources:
        mem_gb=config["mem_gb"]["multiqc"],
    params:
        config_file="config/multiqc_config.yaml",
    log:
        log_dir.joinpath("multiqc", "cluster_{cluster}", "multiqc.log"),
    shell:
        """
DIRNAME=$(dirname {output})
multiqc --interactive --force --config {params.config_file} \
    -o $DIRNAME \
    -n multiqc.html {input} &> {log}
        """
