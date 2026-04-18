process SAMTOOLS_STATS {
    tag "$meta.id"
    label 'process_low'

    container 'quay.io/biocontainers/samtools:1.20--h50ea8bc_0'

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.stats")   , emit: stats
    tuple val(meta), path("*.flagstat"), emit: flagstat
    tuple val(meta), path("*.idxstats"), emit: idxstats
    path  "versions.yml"               , emit: versions

    script:
    """
    samtools stats $bam > ${meta.id}.stats
    samtools flagstat $bam > ${meta.id}.flagstat
    samtools idxstats $bam > ${meta.id}.idxstats

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
