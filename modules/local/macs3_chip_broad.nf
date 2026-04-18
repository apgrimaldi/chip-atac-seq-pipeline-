process MACS3_CHIP_BROAD {
    tag "$meta.id"
    label 'process_medium'
    container 'quay.io/biocontainers/macs3:3.0.1--py311h0152c62_3'

    input:
    tuple val(meta), path(ip_bam), path(control_bam)

    output:
    tuple val(meta), path("*.broadPeak") , emit: peaks
    tuple val(meta), path("*.gappedPeak"), emit: gapped_peaks
    path "versions.yml"                  , emit: versions

    script:
    def prefix = "${meta.id}_broad"
    def format = meta.single_end ? 'BAM' : 'BAMPE'
    def m_genome = (params.genome == 'hg38' || params.genome == 'GRCh38') ? 'hs' : params.genome
    """
    macs3 callpeak \\
        -t $ip_bam \\
        -c $control_bam \\
        -f $format \\
        -g $genome \\
        -n $prefix \\
        --broad \\
        --broad-cutoff 0.1

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        macs3: \$(macs3 --version | sed 's/macs3 //g')
    END_VERSIONS
    """
}
