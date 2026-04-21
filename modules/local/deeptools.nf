process DEEPTOOLS {
    tag "$meta.id"
    label 'process_high'
    container 'quay.io/biocontainers/deeptools:3.5.5--pyhdfd78af_0'
    
    publishDir "${params.outdir}/07_advanced_qc", mode: 'copy'

    input:
    tuple val(meta), path(bam), path(bai)
    // Rimuoviamo l'obbligo del file BED qui se non lo usiamo

    output:
    path "*.fingerprint.pdf"    , emit: fingerprint_pdf
    path "*.fingerprint.txt"    , emit: fingerprint_txt 
    path "*.bigWig"             , emit: bw
    path "versions.yml"         , emit: versions
    // Rimuoviamo i profili dagli output

    script:
    def prefix = "${meta.id}"
    """
    # 1. Genera BigWig
    bamCoverage -b $bam -o ${prefix}.bigWig --binSize 10 --normalizeUsing CPM --numberOfProcessors $task.cpus

    # 2. Fingerprint
    plotFingerprint -b $bam --plotFile ${prefix}.fingerprint.pdf --outRawCounts ${prefix}.fingerprint.txt --numberOfProcessors $task.cpus --skipZeros

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        deeptools: \$(deeptools --version | cut -d' ' -f2)
    END_VERSIONS
    """
}
