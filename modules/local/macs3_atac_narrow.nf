process MACS3_ATAC_NARROW {
    tag "$meta.id"
    label 'process_medium'
    container 'quay.io/biocontainers/macs3:3.0.1--py311h0152c62_3'

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.narrowPeak"), emit: peaks
    path "versions.yml"                  , emit: versions

    script:
    def prefix = "${meta.id}_atac_narrow"
    def format = meta.single_end ? 'BAM' : 'BAMPE'
    def m_genome = (params.genome == 'hg38' || params.genome == 'GRCh38') ? 'hs' : params.genome
"""
    macs3 callpeak \\
        -t $bam \\
        -f $format \\
        -g $genome \\
        -n $prefix \\
        --nomodel --shift -100 --extsize 200 \\
        --qvalue 0.05

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        macs3: \$(macs3 --version | sed 's/macs3 //g')
    END_VERSIONS
    """
}
