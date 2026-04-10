process TRIMGALORE {
    tag "${meta.id}"
    label 'process_high' 

    // Usiamo un container di Quay.io, molto più stabile
    container 'quay.io/biocontainers/trim-galore:0.6.10--hdfd78af_0'

    publishDir "${params.outdir}/trimgalore", mode: 'copy'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.fq.gz")              , emit: reads
    tuple val(meta), path("*_trimming_report.txt"), emit: log
    path "versions.yml"                           , emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    trim_galore \\
        $args \\
        --paired \\
        --gzip \\
        ${reads[0]} \\
        ${reads[1]}

    cat <<EOF > versions.yml
    "${task.process}":
        trimgalore: \$(echo \$(trim_galore --version 2>&1) | sed 's/^.*version //; s/ .*\$//')
    EOF
    """
}
