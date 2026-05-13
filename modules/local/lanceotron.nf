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
    path "versions.yml"                       , emit: versions

    script:
    def prefix = "${meta.id}"
    
    // Cambiamo la logica dei comandi in base alla presenza del controllo
    if (bam_ctrl && bw_ctrl) {
        """
        lanceotron callPeaksInput \\
            -f ${bam_ip} \\
            -w ${bw_ip} \\
            -i ${bam_ctrl} \\
            -c ${bw_ctrl} \\
            -o . \\
            -t 0.9 \\
            -s 1000

        mv L_extract_peaks.bed ${prefix}_lanceotron_peaks.bed
        
        echo "Sample Peaks" > ${prefix}.lanceotron_counts.txt
        COUNT=\$(grep -v "^#" ${prefix}_lanceotron_peaks.bed | wc -l)
        echo "${prefix} \$COUNT" >> ${prefix}.lanceotron_counts.txt

        cat <<EOF > versions.yml
        "${task.process}":
            lanceotron: 1.2.7
        EOF
        """
    } else {
        """
        lanceotron callPeaks \\
            -f ${bam_ip} \\
            -w ${bw_ip} \\
            -o . \\
            -t 0.9 \\
            -s 1000

        mv L_extract_peaks.bed ${prefix}_lanceotron_peaks.bed
        
        echo "Sample Peaks" > ${prefix}.lanceotron_counts.txt
        COUNT=\$(grep -v "^#" ${prefix}_lanceotron_peaks.bed | wc -l)
        echo "${prefix} \$COUNT" >> ${prefix}.lanceotron_counts.txt

        cat <<EOF > versions.yml
        "${task.process}":
            lanceotron: 1.2.7
        EOF
        """
    }
}
