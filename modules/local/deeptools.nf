process DEEPTOOLS {
    tag "$meta.id"
    label 'process_high'
    container 'quay.io/biocontainers/deeptools:3.5.5--pyhdfd78af_0'
    
    publishDir "${params.outdir}/07_advanced_qc", mode: 'copy'

    input:
    // Questa tupla deve ricevere meta, bam e bai tutti insieme
    tuple val(meta), path(bam), path(bai)
    path genes_bed 

    output:
    path "*.fingerprint.pdf"    , emit: fingerprint_pdf
    path "*.fingerprint.txt"    , emit: fingerprint_txt 
    path "*.bigWig"             , emit: bw
    path "*.profile.pdf"        , emit: profile_pdf
    path "*.profile.data.gz"    , emit: profile_data
    path "versions.yml"         , emit: versions

    script:
    def prefix = "${meta.id}"
    """
    # 1. Genera BigWig
    bamCoverage -b $bam -o ${prefix}.bigWig --binSize 10 --normalizeUsing CPM --numberOfProcessors $task.cpus

    # 2. Fingerprint
    plotFingerprint -b $bam --plotFile ${prefix}.fingerprint.pdf --outRawCounts ${prefix}.fingerprint.txt --numberOfProcessors $task.cpus --skipZeros

    # 3. Compute Matrix
    computeMatrix reference-point \\
        --referencePoint TSS \\
        -b 2000 -a 2000 \\
        -R $genes_bed \\
        -S ${prefix}.bigWig \\
        -o ${prefix}.matrix.gz \\
        --numberOfProcessors $task.cpus

    # 4. Plot Profile
    plotProfile -m ${prefix}.matrix.gz \\
        -out ${prefix}.profile.pdf \\
        --outFileNameData ${prefix}.profile.data.gz \\
        --plotTitle "${prefix} TSS Profile"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        deeptools: \$(deeptools --version | cut -d' ' -f2)
    END_VERSIONS
    """
}
