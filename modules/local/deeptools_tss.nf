process DEEPTOOLS_TSS {
    tag "$meta.id"
    label 'process_high'
    container 'quay.io/biocontainers/deeptools:3.5.5--pyhdfd78af_0'
    
    publishDir "${params.outdir}/07_advanced_qc/tss_profiles", mode: 'copy'

    input:
    tuple val(meta), path(bw)
    path gtf

    output:
    tuple val(meta), path("*.tss.matrix.gz"), emit: matrix
    path "*.tss.png"                       , emit: png
    path "*.tss.tab"                       , emit: table 
    path "versions.yml"                    , emit: versions

    script:
def prefix = "${meta.id}"
"""
# Rimuove righe di intestazione o commenti che rompono deepTools
grep -v '^#' $gtf > clean.gtf

computeMatrix reference-point \\
    --referencePoint TSS \\
    -S $bw \\
    -R clean.gtf \\
    -a 3000 -b 3000 \\
    --skipZeros \\
    -o ${prefix}.tss.matrix.gz \\
    --numberOfProcessors $task.cpus

    # 2. Genera il grafico e i DATI (tab) per MultiQC
    # MultiQC ha bisogno del file .tab per generare il grafico interattivo
    plotProfile \\
        -m ${prefix}.tss.matrix.gz \\
        -out ${prefix}.tss.png \\
        --outFileNameData ${prefix}.tss.tab \\
        --plotTitle "${prefix} TSS Profile"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        deeptools: \$(deeptools --version | cut -d' ' -f2)
    END_VERSIONS
    """
}
