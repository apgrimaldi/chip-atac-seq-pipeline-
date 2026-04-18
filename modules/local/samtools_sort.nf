process SAMTOOLS_SORT {
    tag "$meta.id"
    label 'process_medium'
    container 'quay.io/biocontainers/samtools:1.20--h50ea8bc_0'

    input:
    tuple val(meta), path(raw_bam)

    output:
    tuple val(meta), path("*.sorted.bam"), emit: bam
    path "versions.yml"                  , emit: versions

    script:
    def args = task.ext.args ?: '-m 1G' // Usa i parametri del config, altrimenti default 1G
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    samtools sort \\
        $args \\
        -@ $task.cpus \\
        -o ${prefix}.sorted.bam \\
        $raw_bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
