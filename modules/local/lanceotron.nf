process LANCEOTRON {
    tag "${meta.id}"
    label 'process_high'
    container 'quay.io/biocontainers/lanceotron:1.2.7--pyhdfd78af_0'

    publishDir "${params.outdir}/05_peak_calling/lanceotron", mode: 'copy'

    input:
    tuple val(meta), path(bam_ip), path(bw_ip), path(bam_ctrl), path(bw_ctrl)

    output:
    tuple val(meta), path("*_peaks.bed")      , emit: peaks
    tuple val(meta), path("*_counts.txt")     , emit: counts_mqc, optional: true
    path "versions.yml"                        , emit: versions

    script:
    def prefix = "${meta.id}"
    def command = (bam_ctrl && bw_ctrl) ? "callPeaksInput ${bw_ip} -i ${bw_ctrl}" : "callPeaks ${bw_ip}"
    
    """
    lanceotron ${command} \\
        -f . \\
        -t 4 \\
        -w 400

    FOUND_BED=\$(ls *.bed 2>/dev/null | grep -v "filtered" | head -n 1)

    if [ -n "\$FOUND_BED" ]; then
        mv "\$FOUND_BED" ${prefix}_lanceotron_peaks.bed
    else
        touch ${prefix}_lanceotron_peaks.bed
    fi

    echo "Sample Peaks" > ${prefix}.lanceotron_counts.txt
    if [ -s ${prefix}_lanceotron_peaks.bed ]; then
        COUNT=\$(grep -v "^#" ${prefix}_lanceotron_peaks.bed | wc -l)
        echo "${prefix} \$COUNT" >> ${prefix}.lanceotron_counts.txt
    else
        echo "${prefix} 0" >> ${prefix}.lanceotron_counts.txt
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        lanceotron: 1.2.7
    END_VERSIONS
    """
}
