process DEEPTOOLS_COMPUTEMATRIX {
    tag "$meta.id"
    label 'process_high'
    // Utilizziamo un container stabile
    container 'quay.io/biocontainers/deeptools:3.5.5--pyhdfd78af_0'

    input:
    tuple val(meta), path(bigwig)
    path  bed

    output:
    tuple val(meta), path("*.mat.gz"), emit: matrix
    path  "versions.yml"              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    
    // Assicuriamoci che numberOfProcessors non superi le risorse allocate
    """
    computeMatrix \\
        $args \\
        --regionsFileName $bed \\
        --scoreFileName $bigwig \\
        --outFileName ${prefix}.computeMatrix.mat.gz \\
        --numberOfProcessors $task.cpus

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        deeptools: \$(computeMatrix --version | sed -e "s/computeMatrix //g")
    END_VERSIONS
    """
}
