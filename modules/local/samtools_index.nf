process SAMTOOLS_INDEX {
    tag "$meta.id"
    label 'process_low'
    container 'quay.io/biocontainers/samtools:1.16.1--h6899075_1'

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path(bam), path("*.bai"), emit: bam_bai
    path "versions.yml"                      , emit: versions

    script:
    """
    samtools index $bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
