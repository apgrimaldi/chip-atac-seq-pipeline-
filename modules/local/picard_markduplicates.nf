process PICARD_MARKDUPLICATES {
    tag "$meta.id"
    label 'process_medium'
    container 'quay.io/biocontainers/picard:2.27.4--hdfd78af_0'
    
    publishDir "${params.outdir}/04_deduplicated", mode: 'copy'

    input:
    tuple val(meta), path(bam)
    path fasta  // Rimosso tuple val(meta2) perché il fasta è un file singolo opzionale
    path fai    // Rimosso tuple val(meta3) perché il fai è un file singolo opzionale

    output:
    tuple val(meta), path("*.removed.bam")       , emit: bam
    tuple val(meta), path("*.{bai,bam.bai}")     , emit: bai, optional: true
    tuple val(meta), path("*.metrics.txt")       , emit: metrics
    path  "versions.yml"                         , emit: versions

    script:
    def args = task.ext.args ?: '--REMOVE_DUPLICATES true --ASSUME_SORTED true --VALIDATION_STRINGENCY LENIENT --CREATE_INDEX true'
    def prefix = task.ext.prefix ?: "${meta.id}"
    
    // Gestione riferimento opzionale
    def reference = fasta ? "--REFERENCE_SEQUENCE ${fasta}" : ""
    
    // Calcolo memoria per Java (80% della memoria assegnata al task)
    def avail_mem = 3072
    if (!task.memory) {
        log.info '[Picard MarkDuplicates] Available memory not known - defaulting to 3GB.'
    } else {
        avail_mem = (task.memory.mega * 0.8).intValue()
    }

    """
    picard \\
        -Xmx${avail_mem}M \\
        MarkDuplicates \\
        $args \\
        --INPUT $bam \\
        --OUTPUT ${prefix}.removed.bam \\
        --METRICS_FILE ${prefix}.metrics.txt \\
        $reference

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        picard: \$(picard MarkDuplicates --version 2>&1 | grep -o 'Version:.*' | cut -f2- -d:)
    END_VERSIONS
    """
}
