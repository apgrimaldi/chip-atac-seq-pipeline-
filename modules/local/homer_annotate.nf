process HOMER_ANNOTATEPEAKS {
    tag "$meta.id"
    label 'process_medium'

    container 'quay.io/biocontainers/homer:4.11--pl526hc9558a2_3'

    input:
    tuple val(meta), path(peak)
    path  fasta
    path  gtf

    output:
    tuple val(meta), path("*.annotatePeaks.txt"), emit: txt
    path  "versions.yml"                        , emit: versions

    script:
    def prefix = "${peak.baseName}"
    def VERSION = '4.11' 
    """
    annotatePeaks.pl \\
        $peak \\
        $fasta \\
        -gtf $gtf \\
        -cpu $task.cpus \\
        > ${prefix}.annotatePeaks.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        homer: $VERSION
    END_VERSIONS
    """
}
