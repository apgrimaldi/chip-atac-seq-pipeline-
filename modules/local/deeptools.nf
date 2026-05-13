process DEEPTOOLS {
    tag "$meta.id"
    label 'process_high'
    container 'quay.io/biocontainers/deeptools:3.5.5--pyhdfd78af_0'
    
    publishDir "${params.outdir}/06_bigwig", mode: 'copy', pattern: "*.bw"
    publishDir "${params.outdir}/06_bigwig/qc_fingerprint", mode: 'copy', pattern: "*.{pdf,txt}"

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("*.display.bw")           , emit: bw_display
    tuple val(meta), path("*.l6n.bw")               , emit: bw_lanceotron
    tuple val(meta), path("*.plotFingerprint.pdf")  , emit: fingerprint_pdf
    tuple val(meta), path("*.plotFingerprint.raw.txt"), emit: fingerprint_txt 
    tuple val(meta), path("*.plotFingerprint.qcmetrics.txt"), emit: fingerprint_metrics
    path "versions.yml"                             , emit: versions

    script:
    def prefix = "${meta.id}"
    // Gestione estensione reads
    def extend = (meta.single_end && params.fragment_size > 0) ? "--extendReads ${params.fragment_size}" : (params.single_end ? "--extendReads" : "")
    
    """
    # 1. BigWig per Visualizzazione (Normalizzato, BinSize 10bp per fluidità)
    bamCoverage \\
        --bam $bam \\
        --outFileName ${prefix}.display.bw \\
        --binSize 10 \\
        --normalizeUsing CPM \\
        --numberOfProcessors $task.cpus \\
        $extend

    # 2. BigWig per Lanceotron (NON normalizzato, BinSize 1bp - Fondamentale per L6n)
    bamCoverage \\
        --bam $bam \\
        --outFileName ${prefix}.l6n.bw \\
        --binSize 1 \\
        --numberOfProcessors $task.cpus \\
        $extend

    # 3. Fingerprint QC
    plotFingerprint \\
        --bamfiles $bam \\
        --plotFile ${prefix}.plotFingerprint.pdf \\
        --outRawCounts ${prefix}.plotFingerprint.raw.txt \\
        --outQualityMetrics ${prefix}.plotFingerprint.qcmetrics.txt \\
        --numberOfProcessors $task.cpus \\
        $extend \\
        --skipZeros

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        deeptools: \$(plotFingerprint --version | sed -e "s/plotFingerprint //g")
    END_VERSIONS
    """
}
